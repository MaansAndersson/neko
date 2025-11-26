!!



module setup
  use neko
  implicit none
contains
  ! Setup rhs
  subroutine set_f(f, dm)
    type(dofmap_t), intent(in) :: dm
    real(kind=rp), intent(inout), dimension(dm%size()) :: f
    integer :: n
    real(kind=rp) :: dx, dy, dz
    real(kind=rp), parameter :: arg = 2d0
    integer :: i, idx(4)

    do i = 1, dm%size()
       idx = nonlinear_index(i, dm%Xh%lx, dm%Xh%ly, dm%Xh%lz)
       dx = dm%x(idx(1), idx(2), idx(3), idx(4)) - 4.0d0
       dy = dm%y(idx(1), idx(2), idx(3), idx(4)) - 4.0d0
       dz = dm%z(idx(1), idx(2), idx(3), idx(4)) - 4.0d0
       f(i) = 500d0*cos(-(dx**arg + dy**arg + dz**arg)/arg)
    end do
  end subroutine set_f

  ! Set Dirichlet conditions
  subroutine set_bc(bc_, msh)
    type(mesh_t), intent(in) :: msh
    type(dirichlet_t), intent(inout) :: bc_
    integer :: i

    do i = 1, msh%nelv
       if (msh%facet_neigh(1, i) .eq. 0) then
         call bc_%mark_facet(1, i)
       end if
       if (msh%facet_neigh(2, i) .eq. 0) then
         call bc_%mark_facet(2, i)
       end if
       if (msh%facet_neigh(3, i) .eq. 0) then
         call bc_%mark_facet(3, i)
       end if
       if (msh%facet_neigh(4, i) .eq. 0) then
         call bc_%mark_facet(4, i)
       end if
       if (msh%facet_neigh(5, i) .eq. 0) then
         call bc_%mark_facet(5, i)
       end if
       if (msh%facet_neigh(6, i) .eq. 0) then
         call bc_%mark_facet(6, i)
       end if
    enddo
  end subroutine set_bc

end module setup


program nekobench
  use neko
  use setup
  use inexact_pc

  implicit none

  character(len=NEKO_FNAME_LEN) :: fname, lxchar
  type(mesh_t) :: msh
  type(file_t) :: nmsh_file, mf
  type(space_t) :: Xh
  type(dofmap_t) :: dm
  type(gs_t) :: gs_h
  type(coef_t) :: coef
  class(ax_t), allocatable :: ax_helm
  type(field_t) :: f1, f2, f3, f4
  real(kind=rp) :: byte, flop, n_tot, t0, t1, time
  real(kind=rp) :: c1, c2, abstol
  class(ksp_t), allocatable :: solver
  type(ksp_monitor_t) :: ksp_mon
  type(dirichlet_t) :: dir_bc
  type(bc_list_t) :: bclst
  integer :: argc, niter, ierr, lx, nelt 
  real(kind=rp) :: tau
  integer :: i, n

  ! Inexact preconditioner
  class(ksp_t), allocatable :: solver_pc
  class(inexact_pc_t), allocatable :: pc
  !class(gs_bcknd_t), allocatable :: bcknd !< Gather-scatter backend
  type(gs_t) :: gs_mpi_straggler

  argc = command_argument_count()

  if ((argc .lt. 4) .or. (argc .gt. 4)) then
     write(*,*) 'Usage: ./solver.bin <neko mesh> <N> <niter> <tau (percentage of messages in Cheby)>'
     write(*,*) 'Use meshes from poisson, e.g. ../poisson/data.512.nmsh'
  endif

  call neko_init

  call get_command_argument(1, fname)
  call get_command_argument(2, lxchar)
  read(lxchar, *) lx
  call get_command_argument(3, lxchar)
  read(lxchar, *) niter
  call get_command_argument(4, lxchar)
  read(lxchar, *) tau

  call nmsh_file%init(fname)
  call nmsh_file%read(msh)

  !Init things
  call Xh%init(GLL, lx, lx, lx)
  call dm%init(msh,Xh)
  call gs_h%init(dm,1,1)
  call gs_mpi_straggler%init(dm,1,5) ! Add percentage
  !call gs_mpi_straggler%comm%gs_set_tau_mpi(tau)
  call f1%init(dm)
  call f2%init(dm)
  call f3%init(dm)
  call f4%init(dm)
  call coef%init(gs_h)

  n_tot = dble(msh%glb_nelv)*dble(niter)*dble(Xh%lxyz) 
  n = dm%size()
  if (pe_rank .eq. 0) then
     write(*,*) 'Straggling Solver '
     write(*,*) 'lx:', lx
     write(*,*) 'N elements tot:', msh%glb_nelv
     write(*,*) 'N ranks:', pe_size
     write(*,*) 'N iterations:', niter
     write(*,*) 'lx^3*N elements tot*N iter:', n_tot 
  end if
  call MPI_Barrier(NEKO_COMM, ierr)


  f1 = 1.0_rp
  f2 = 1.0_rp
  f3 = 1.0_rp

  call device_sync()

  call set_f(f1%x,f1%dof)
  if(NEKO_BCKND_DEVICE .eq. 1) call device_memcpy(f1%x, f1%x_d, n, HOST_TO_DEVICE,sync=.true.) 
  call device_sync()

  call MPI_Barrier(NEKO_COMM, ierr)
  t0 = MPI_Wtime()
  t1 = MPI_Wtime()
  time = t1 - t0

  !example us of cg solver
  !init bcs...
  call dir_bc%init_from_components(coef,real(0.0d0,rp))

  !user specified
  call set_bc(dir_bc, msh)

  call dir_bc%finalize()
  call bclst%init()
  call bclst%append(dir_bc)

  abstol = 1e-4

  call ax_helm_factory(ax_helm, .FALSE.)
  call krylov_solver_factory(solver_pc, dm%size(), 'cheby', niter, abstol)
  ! create a Chebyshev iteration with gs_mpi_straggler
  ! change gs_h such that it is a gs that times out for small ts

  allocate(inexact_pc_t::pc)
  select type (pc => pc)

  type is (inexact_pc_t)
    call pc%init(ax_helm, solver_pc, 100, coef, dm, gs_mpi_straggler, gs_h, bclst)
  end select
  abstol = 1e-2! Something small so we dont converge


  call krylov_solver_factory(solver, dm%size(), 'gmres', niter, abstol, pc)
  f2 = 0.0_rp
  t0 = MPI_Wtime()

  ksp_mon = solver%solve(ax_helm, f2, f1%x, dm%size(), coef, bclst, gs_h, niter)

  t1 = MPI_Wtime()
  time = t1 - t0

  write(*,*) ksp_mon%iter, &
         ksp_mon%res_start, &
         ksp_mon%res_final, &
         time


  if (NEKO_BCKND_DEVICE .eq. 1) &
     call device_memcpy(f2%x, f2%x_d, n, DEVICE_TO_HOST, sync=.true.)
  ! Store the solution
  fname = 'out.fld'
  call mf%init(fname)
  call mf%write(f2)

  call neko_finalize

end program nekobench

