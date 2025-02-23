!> Module for reusing linear solvers as preconditioner.
!!
!! Used with FMGRES or BiCGSTAB.

!> Krylov preconditioner (using Krylov solver)
module inexact_pc
  use neko_config, only : NEKO_BCKND_DEVICE
  use math, only : copy
  use device_math, only : device_copy
  use precon, only : pc_t
  use coefs, only : coef_t
  use num_types, only : c_rp, rp
  use, intrinsic :: iso_c_binding , only : c_ptr, c_int
  use field, only : field_t
  use ax_product, only : ax_t
  use dofmap
  use gather_scatter
  use krylov
  use device, only : device_get_ptr
  use bc_list , only : bc_list_t

  implicit none
  private

  type, public, extends(pc_t) :: inexact_pc_t
    type(gs_t), pointer :: gs_h
    type(dofmap_t), pointer :: dof
    type(coef_t), pointer :: coef
    type(bc_list_t), pointer :: bclst
    class(ax_t), pointer :: Ax
    class(ksp_t), pointer :: M => null() !allocatable?
    type(field_t) :: temp_field
    integer :: inner_iter
    contains
      procedure, pass(this) :: init => inexact_init
      procedure, pass(this) :: solve => inexact_solve
      procedure, pass(this) :: update => inexact_update
      procedure, pass(this) :: free => inexact_free
  end type inexact_pc_t

contains

!> The preconditioner \f$ M z = r \f$ solved with some inexact solver
  subroutine inexact_solve(this, z, r, n)
    class(inexact_pc_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: z
    real(kind=rp), dimension(n), intent(inout) :: r
    type(ksp_monitor_t) :: ksp_mon
    type(c_ptr) :: z_d, r_d

    !? Do I have to sync somewhere?
    this%temp_field = 0.0_rp
    ksp_mon = this%M%solve(this%Ax, this%temp_field, r, n, this%coef, this%bclst, this%gs_h, this%inner_iter)


    write(*,*) ksp_mon%iter, &
         ksp_mon%res_start, &
         ksp_mon%res_final

    if (NEKO_BCKND_DEVICE .eq. 1) then
      z_d = device_get_ptr(z)
      call device_copy(z_d,this%temp_field%x_d,n)
    else
      call copy(z, this%temp_field%x, n)
    endif

    write(*,*) ksp_mon%iter, &
         ksp_mon%res_start, &
         ksp_mon%res_final

  end subroutine inexact_solve

!> Init solver, system, gather-scatter-method (straggler).
  subroutine inexact_init(this, Ax, M, inner_iter, coef, dof, gs_h, bclst)
    class(inexact_pc_t), intent(inout) :: this
    class(ax_t), intent(in), target :: Ax
    class(ksp_t), intent(in), target :: M
    type(gs_t), intent(in), target :: gs_h
    type(dofmap_t), intent(in), target :: dof
    type(coef_t), intent(in), target :: coef
    type(bc_list_t), intent(in), target :: bclst
    integer, intent(in) :: inner_iter

    call this%free()
    call this%temp_field%init(dof)
    this%inner_iter = inner_iter
    this%Ax => Ax
    this%M => M
    this%gs_h => gs_h
    this%dof => dof
    this%coef => coef
    this%bclst => bclst
    call inexact_update(this)

  end subroutine inexact_init

!> Mandatory update routine
!! Should probably update BCs, tolerance etc.
  subroutine inexact_update(this)
    class(inexact_pc_t), intent(inout) :: this
  end subroutine inexact_update

  subroutine inexact_free(this)
    class(inexact_pc_t), intent(inout) ::this
  end subroutine inexact_free

end module inexact_pc
