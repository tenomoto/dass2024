module l96_module
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use rk4_module, only: rk4, rk4_tl, rk4_ad
  implicit none

contains

  function l96_nl(t, x, opts) result(x_out)
    real(dp), intent(in) :: t
    real(dp), intent(in) :: x(:), opts(:)
    real(dp) :: x_out(size(x))

    x_out = (cshift(x, 1) - cshift(x, -2)) * cshift(x, -1) - x + opts(1)
  end function l96_nl

  function l96_tl(t, x, dx, opts) result(dx_out)
    real(dp), intent(in) :: t
    real(dp), intent(in) :: x(:), dx(:), opts(:)
    real(dp) :: dx_out(size(x))

    dx_out = (cshift(dx, 1) - cshift(dx, -2)) * cshift(x, -1) + &
      (cshift(x, 1) - cshift(x, -2)) * cshift(dx, -1) - dx
  end function l96_tl

  function l96_ad(t, x, dxa, opts) result(xa_out)
    real(dp), intent(in) :: t
    real(dp), intent(in) :: x(:), dxa(:), opts(:)
    real(dp) :: xa_out(size(x))

    xa_out = -cshift(x, 1) * cshift(dxa, 2) + &
      (cshift(x, 2) - cshift(x, -1)) * cshift(dxa, 1) + &
      cshift(x, -2) * cshift(dxa, -1) - dxa
  end function l96_ad

  function l96_fom(x, nstep, dt, F) result(x_out)
    real(dp), intent(in) :: dt, F
    integer, intent(in) :: nstep
    real(dp), intent(in) :: x(:)
    real(dp) :: opts(1)
    real(dp) :: x_out(size(x))

    integer :: i
    real(dp) :: t = 0

    opts(1) = F
    x_out = x
    do i = 1, nstep
      x_out = rk4(l96_nl, t, x_out, dt, opts)
      t = t + dt
    end do
  end function l96_fom

  function l96_tlm(x, dx, nstep, dt, F) result(dx_out)
    real(dp), intent(in) :: dt, F
    integer, intent(in) :: nstep
    real(dp), intent(in) :: x(:), dx(:)
    real(dp) :: opts(1)
    real(dp) :: dx_out(size(dx))

    real(dp) :: t = 0
    integer :: i

    opts(1) = F
    dx_out = dx
    do i = 1, nstep
      dx_out = rk4_tl(l96_nl, l96_tl, t, x, dx_out, dt, opts)
      t = t + dt
    end do
  end function l96_tlm

  function l96_adm(x, xa, nstep, dt, F) result(xa_out)
    real(dp), intent(in) :: dt, F
    integer, intent(in) :: nstep
    real(dp), intent(in) :: x(:), xa(:)
    real(dp) :: opts(1)
    real(dp) :: xa_out(size(xa))

    integer :: i
    real(dp) :: t = 0

    opts(1) = F
    t = nstep
    xa_out = xa
    do i = nstep, 1, -1
      xa_out = rk4_ad(l96_nl, l96_ad, t, x, xa_out, dt, opts)
      t = t - dt
    end do
  end function l96_adm

  subroutine l96_test(a)
    use random_module, only: &
      set_seed => random_set_seed, rnorm => random_normal
    real(dp), intent(in) :: a
    integer, parameter :: &
      ns = 40, nstep_spinup = 1000, seed = 514
    real(dp), parameter :: F = 8.0_dp, dt = 0.05_dp
    real(dp) :: x0(ns), x(ns), x1(ns), dx0(ns), dx(ns), xa(ns)
    real(dp) :: tdxdx, dx0xa, e

    call set_seed(seed)
    x0 = l96_fom(rnorm(ns), nstep_spinup, dt, F)
    x = l96_fom(x0, 1, dt, F)
    
    dx0 = a * rnorm(ns)
    x1 = l96_fom(x0 + dx0, 1, dt, F)
    dx = l96_tlm(x0, dx0, 1, dt, F)
    e = sqrt(sum((x1 - x - dx)**2))
    print *, "TLM:", "a=", a, "l2=", e
    
    xa = l96_adm(x0, dx, 1, dt, F)
    tdxdx = sum(dx * dx)
    dx0xa = sum(dx0 * xa)

    print *, "ADJ: dt^t dx - t(dx0) xa =", tdxdx, "-", dx0xa, "=", tdxdx - dx0xa
  end subroutine l96_test

end module l96_module
