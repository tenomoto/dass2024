module rk4_module
  use iso_fortran_env, only: real64
  implicit none
contains

  function rk4(f, t, y, h) result(y_out)
    interface
      function f(t, y) result(res)
        use iso_fortran_env, only: real64
        real(real64), intent(in) :: t
        real(real64), intent(in) :: y(:)
        real(real64) :: res(size(y))
      end function f
    end interface
    real(real64), intent(in) :: t, h
    real(real64), intent(in) :: y(:)
    real(real64) :: k1(size(y)), k2(size(y)), k3(size(y)), k4(size(y)), y_out(size(y))

    k1 = f(t, y)
    k2 = f(t + 0.5_real64 * h, y + 0.5_real64 * h * k1)
    k3 = f(t + 0.5_real64 * h, y + 0.5_real64 * h * k2)
    k4 = f(t + h, y + h * k3)
    y_out = y + h * (k1 + 2_real64 * k2 + 2_real64 * k3 + k4) / 6.0_real64
  end function rk4

  function rk4_tl(f, df, t, y, dy, h) result(dy_out)
    interface
      function f(t, y) result(res)
        use iso_fortran_env, only: real64
        real(real64), intent(in) :: t
        real(real64), intent(in) :: y(:)
        real(real64) :: res(size(y))
      end function f
      function df(t, y, dy) result(res)
        use iso_fortran_env, only: real64
        real(real64), intent(in) :: t
        real(real64), intent(in) :: y(:), dy(:)
        real(real64) :: res(size(dy))
      end function df
    end interface
    real(real64), intent(in) :: t, h
    real(real64), intent(in) :: y(:), dy(:)
    real(real64) :: k1(size(y)), k2(size(y)), k3(size(y)), &
      dk1(size(dy)), dk2(size(dy)), dk3(size(dy)), dk4(size(dy)), dy_out(size(dy))

    k1 = f(t, y)
    k2 = f(t + 0.5_real64 * h, y + 0.5_real64 * h * k1)
    k3 = f(t + 0.5_real64 * h, y + 0.5_real64 * h * k2)
    dk1 = df(t, y, dy)
    dk2 = df(t + 0.5_real64 * h, y + 0.5_real64 * h * k1, dy + 0.5_real64 * h * dk1)
    dk3 = df(t + 0.5_real64 * h, y + 0.5_real64 * h * k2, dy + 0.5_real64 * h * dk2)
    dk4 = df(t + h, y + h * k3, dy + h * dk3)
    dy_out = dy + h * (dk1 + 2_real64 * dk2 + 2_real64 * dk3 + dk4) / 6.0_real64
  end function rk4_tl

  function rk4_ad(f, fa, t, y, ya, h) result(ya_out)
    interface
      function f(t, y) result(res)
        use iso_fortran_env, only: real64
        real(real64), intent(in) :: t
        real(real64), intent(in) :: y(:)
        real(real64) :: res(size(y))
      end function f
      function fa(t, y, ya) result(res)
        use iso_fortran_env, only: real64
        real(real64), intent(in) :: t
        real(real64), intent(in) :: y(:), ya(:)
        real(real64) :: res(size(ya))
      end function fa
    end interface
    real(real64), intent(in) :: t, h
    real(real64), intent(in) :: y(:), ya(:)
    real(real64) :: k1(size(y)), k2(size(y)), k3(size(y)), &
      y2(size(y)), y3(size(y)), y4(size(y))
    real(real64) :: k1a(size(ya)), k2a(size(ya)), k3a(size(ya)), k4a(size(ya)), &
      dfa(size(ya)), ya_out(size(ya))

    k1 = f(t, y)
    y2 = y + 0.5_real64 * h * k1
    k2 = f(t + 0.5_real64 * h, y2)
    y3 = y + 0.5_real64 * h * k2
    k3 = f(t + 0.5_real64 * h, y3)
    y4 = y + h * k3

    k1a = h * ya / 6.0_real64
    k2a = h * ya / 3.0_real64
    k3a = k2a
    k4a = k1a

    dfa = fa(t + h, y4, k4a)
    ya_out = ya + dfa
    k3a = k3a + h * dfa
    k4a = 0.0_real64
    dfa = fa(t + 0.5_real64 * h, y3, k3a)
    ya_out = ya_out + dfa
    k2a = k2a + 0.5_real64 * h * dfa
    k3a = 0.0_real64
    dfa = fa(t + 0.5_real64 * h, y2, k2a)
    ya_out = ya_out + dfa
    k1a = k1a + 0.5_real64 * h * dfa
    k2a = 0.0_real64
    ya_out = ya_out + fa(t, y, k1a)
    k1a = 0.0_real64
  end function rk4_ad

end module rk4_module

