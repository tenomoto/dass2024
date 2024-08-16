module step_module
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use rk4_module, only: rk4, rk4_tl, rk4_ad
  implicit none

contains

  function step_fom(f, x, nstep, dt, opts) result(x_out)
    interface
      function f(t, x, opts) result(res)
        use, intrinsic :: iso_fortran_env, only: dp=>real64
        real(dp), intent(in) :: t
        real(dp), intent(in) :: x(:), opts(:)
        real(dp) :: res(size(x))
      end function f
    end interface
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: nstep
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: opts(:)
    real(dp) :: x_out(size(x))

    integer :: i
    real(dp) :: t = 0

    x_out = x
    do i = 1, nstep
      x_out = rk4(f, t, x_out, dt, opts)
      t = t + dt
    end do
  end function step_fom

  function step_tlm(f, df, x, dx, nstep, dt, opts) result(dx_out)
    interface
      function f(t, x, opts) result(res)
        use, intrinsic :: iso_fortran_env, only: dp=>real64
        real(dp), intent(in) :: t
        real(dp), intent(in) :: x(:), opts(:)
        real(dp) :: res(size(x))
      end function f
      function df(t, x, dx, opts) result(res)
        use, intrinsic :: iso_fortran_env, only: dp=>real64
        real(dp), intent(in) :: t
        real(dp), intent(in) :: x(:), dx(:), opts(:)
        real(dp) :: res(size(dx))
      end function df
    end interface
    real(dp), intent(in) :: x(:), dx(:)
    integer, intent(in) :: nstep
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: opts(:)
    real(dp) :: dx_out(size(dx))

    real(dp) :: t = 0
    integer :: i

    dx_out = dx
    do i = 1, nstep
      dx_out = rk4_tl(f, df, t, x, dx_out, dt, opts)
      t = t + dt
    end do
  end function step_tlm

  function step_adm(f, fa, x, xa, nstep, dt, opts) result(xa_out)
    interface
      function f(t, x, opts) result(res)
        use, intrinsic :: iso_fortran_env, only: dp=>real64
        real(dp), intent(in) :: t
        real(dp), intent(in) :: x(:), opts(:)
        real(dp) :: res(size(x))
      end function f
      function fa(t, x, xa, opts) result(res)
        use, intrinsic :: iso_fortran_env, only: dp=>real64
        real(dp), intent(in) :: t
        real(dp), intent(in) :: x(:), xa(:), opts(:)
        real(dp) :: res(size(xa))
      end function fa
    end interface
    real(dp), intent(in) :: x(:), xa(:)
    integer, intent(in) :: nstep
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: opts(:)
    real(dp) :: xa_out(size(xa))

    integer :: i
    real(dp) :: t = 0

    t = nstep
    xa_out = xa
    do i = nstep, 1, -1
      xa_out = rk4_ad(f, fa, t, x, xa_out, dt, opts)
      t = t - dt
    end do
  end function step_adm

end module step_module
