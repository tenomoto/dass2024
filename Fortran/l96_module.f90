module l96_module
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use step_module, only: step_fom, step_tlm, step_adm
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

  subroutine l96_test(a)
    use random_module, only: &
      set_seed => random_set_seed, rnorm => random_normal
    real(dp), intent(in) :: a

    integer, parameter :: &
      ns = 40, nstep_spinup = 1000, seed = 514
    real(dp), parameter :: F = 8.0_dp, dt = 0.05_dp
    real(dp) :: x0(ns), x(ns), x1(ns), dx0(ns), dx(ns), xa(ns)
    real(dp) :: tdxdx, dx0xa, e
    real(dp) :: opts(1) = F

    call set_seed(seed)
    x0 = step_fom(l96_nl, rnorm(ns), nstep_spinup, dt, opts)
    x = step_fom(l96_nl, x0, 1, dt, opts)
    
    dx0 = a * rnorm(ns)
    x1 = step_fom(l96_nl, x0 + dx0, 1, dt, opts)
    dx = step_tlm(l96_nl, l96_tl, x0, dx0, 1, dt, opts)
    e = sqrt(sum((x1 - x - dx)**2))
    print *, "TLM:", "a=", a, "l2=", e
    
    xa = step_adm(l96_nl, l96_ad, x0, dx, 1, dt, opts)
    tdxdx = sum(dx * dx)
    dx0xa = sum(dx0 * xa)

    print *, "ADJ: dt^t dx - t(dx0) xa =", tdxdx, "-", dx0xa, "=", tdxdx - dx0xa

  end subroutine l96_test

end module l96_module
