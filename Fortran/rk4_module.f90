module rk4_module
  use, intrinsic :: iso_fortran_env, only: dp=>real64
  implicit none

contains

  function rk4(f, t, y, h, opts) result(y_out)
    interface
      function f(t, y, opts) result(res)
        use, intrinsic :: iso_fortran_env, only: dp=>real64
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:), opts(:)
        real(dp) :: res(size(y))
      end function f
    end interface
    real(dp), intent(in) :: t, h
    real(dp), intent(in) :: y(:), opts(:)
    real(dp) :: k1(size(y)), k2(size(y)), k3(size(y)), k4(size(y)), y_out(size(y))

    k1 = f(t, y, opts)
    k2 = f(t + 0.5_dp * h, y + 0.5_dp * h * k1, opts)
    k3 = f(t + 0.5_dp * h, y + 0.5_dp * h * k2, opts)
    k4 = f(t + h, y + h * k3, opts)
    y_out = y + h * (k1 + 2_dp * k2 + 2_dp * k3 + k4) / 6.0_dp
  end function rk4

  function rk4_tl(f, df, t, y, dy, h, opts) result(dy_out)
    interface
      function f(t, y, opts) result(res)
        use, intrinsic :: iso_fortran_env, only: dp=>real64
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:), opts(:)
        real(dp) :: res(size(y))
      end function f
      function df(t, y, dy, opts) result(res)
        use, intrinsic :: iso_fortran_env, only: dp=>real64
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:), dy(:), opts(:)
        real(dp) :: res(size(dy))
      end function df
    end interface
    real(dp), intent(in) :: t, h
    real(dp), intent(in) :: y(:), dy(:), opts(:)
    real(dp) :: k1(size(y)), k2(size(y)), k3(size(y)), &
      dk1(size(dy)), dk2(size(dy)), dk3(size(dy)), dk4(size(dy)), dy_out(size(dy))

    k1 = f(t, y, opts)
    k2 = f(t + 0.5_dp * h, y + 0.5_dp * h * k1, opts)
    k3 = f(t + 0.5_dp * h, y + 0.5_dp * h * k2, opts)
    dk1 = df(t, y, dy, opts)
    dk2 = df(t + 0.5_dp * h, y + 0.5_dp * h * k1, dy + 0.5_dp * h * dk1, opts)
    dk3 = df(t + 0.5_dp * h, y + 0.5_dp * h * k2, dy + 0.5_dp * h * dk2, opts)
    dk4 = df(t + h, y + h * k3, dy + h * dk3, opts)
    dy_out = dy + h * (dk1 + 2_dp * dk2 + 2_dp * dk3 + dk4) / 6.0_dp
  end function rk4_tl

  function rk4_ad(f, fa, t, y, ya, h, opts) result(ya_out)
    interface
      function f(t, y, opts) result(res)
        use, intrinsic :: iso_fortran_env, only: dp=>real64
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:), opts(:)
        real(dp) :: res(size(y))
      end function f
      function fa(t, y, ya, opts) result(res)
        use, intrinsic :: iso_fortran_env, only: dp=>real64
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:), ya(:), opts(:)
        real(dp) :: res(size(ya))
      end function fa
    end interface
    real(dp), intent(in) :: t, h
    real(dp), intent(in) :: y(:), ya(:), opts(:)
    real(dp) :: k1(size(y)), k2(size(y)), k3(size(y)), &
      y2(size(y)), y3(size(y)), y4(size(y))
    real(dp) :: k1a(size(ya)), k2a(size(ya)), k3a(size(ya)), k4a(size(ya)), &
      dfa(size(ya)), ya_out(size(ya))

    k1 = f(t, y, opts)
    y2 = y + 0.5_dp * h * k1
    k2 = f(t + 0.5_dp * h, y2, opts)
    y3 = y + 0.5_dp * h * k2
    k3 = f(t + 0.5_dp * h, y3, opts)
    y4 = y + h * k3

    k1a = h * ya / 6.0_dp
    k2a = h * ya / 3.0_dp
    k3a = k2a
    k4a = k1a

    dfa = fa(t + h, y4, k4a, opts)
    ya_out = ya + dfa
    k3a = k3a + h * dfa
    k4a = 0.0_dp
    dfa = fa(t + 0.5_dp * h, y3, k3a, opts)
    ya_out = ya_out + dfa
    k2a = k2a + 0.5_dp * h * dfa
    k3a = 0.0_dp
    dfa = fa(t + 0.5_dp * h, y2, k2a, opts)
    ya_out = ya_out + dfa
    k1a = k1a + 0.5_dp * h * dfa
    k2a = 0.0_dp
    ya_out = ya_out + fa(t, y, k1a, opts)
    k1a = 0.0_dp
  end function rk4_ad

end module rk4_module

