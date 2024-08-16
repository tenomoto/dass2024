module l63_module
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use step_module, only: step_fom, step_tlm, step_adm
  implicit none

contains

  function l63_nl(t, w, opts) result(dw)
    real(dp), intent(in) :: t 
    real(dp), intent(in) :: w(:), opts(:)
    real(dp) :: dw(size(w))

    real(dp) :: p, r, b

    p = opts(1)
    r = opts(2)
    b = opts(3)
    
    dw(1) =         -p * w(1) + p * w(2)
    dw(2) = (r - w(3)) * w(1) -     w(2)
    dw(3) =              w(1) * w(2) - b * w(3)    

  end function l63_nl

  function l63_tl(t, wb, wt, opts) result(dwt)
    real(dp), intent(in) :: t 
    real(dp), intent(in) :: wb(:), wt(:), opts(:)
    real(dp) :: dwt(size(wb))

    real(dp) :: p, r, b

    p = opts(1)
    r = opts(2)
    b = opts(3)

    dwt(1) =          -p * wt(1) +     p * wt(2)
    dwt(2) = (r - wb(3)) * wt(1) -         wt(2) - wb(1) * wt(3)
    dwt(3) =       wb(2) * wt(1) + wb(1) * wt(2) -     b * wt(3)

  end function l63_tl

  function l63_ad(t, wb, dwa, opts) result(wa)
    real(dp), intent(in) :: t 
    real(dp), intent(in) :: wb(:), dwa(:), opts(:)
    real(dp) :: wa(size(wb))

    real(dp) :: p, r, b

    p = opts(1)
    r = opts(2)
    b = opts(3)

    wa(1) =  -p * dwa(1) + (r - wb(3)) * dwa(2) + wb(2) * dwa(3)
    wa(2) =   p * dwa(1) -               dwa(2) + wb(1) * dwa(3)
    wa(3) =                     -wb(1) * dwa(2) -     b * dwa(3)

  end function l63_ad

  subroutine l63_test(a)
    real(dp), intent(in) :: a

    integer, parameter :: ns = 3
    real(dp), parameter :: &
      p = 10_dp, r = 32_dp, b = 8.0_dp / 3.0_dp, dt = 0.01_dp
    real(dp) :: x0(ns) = [1, 1, 1]
    real(dp) :: x(ns), x1(ns), dx0(ns), dx(ns), xa(ns)
    real(dp) :: e, tdxdx, dx0xa
    real(dp) :: opts(3)

    opts = [p, r, b]
    x = step_fom(l63_nl, x0, 1, dt, opts)
    
    dx0 = a * x0
    x1 = step_fom(l63_nl, x0 + dx0, 1, dt, opts)
    dx = step_tlm(l63_nl, l63_tl, x0, dx0, 1, dt, opts)
    e = sqrt(sum((x1 - x - dx)**2) / ns)
    print *, "TLM:", "a=", a, "l2=", e
    
    xa = step_adm(l63_nl, l63_ad, x0, dx, 1, dt, opts)
    tdxdx = sum(dx * dx) / ns
    dx0xa = sum(dx0 * xa) / ns
    print *, "ADJ: dt^t dx - t(dx0) xa =", tdxdx, "-", dx0xa, "=", tdxdx - dx0xa

  end subroutine l63_test

end module l63_module
