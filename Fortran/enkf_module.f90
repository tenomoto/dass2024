module enkf_module
  use, intrinsic :: iso_fortran_env, dp => real64
  implicit none

contains

  function enkf_analysis(zf, yo, r, loc_inf) result(dz)
    real(dp), intent(in) :: zf(:, :), yo(:), loc_inf(:)
    real(dp), intent(in) :: r
    real(dp), allocatable :: dz(:, :)

    integer :: k, ne!, i
    real(dp), allocatable ::s(:)
    real(dp), allocatable :: zf_anom(:, :)

    k = size(zf, 1)
    ne = size(zf, 2)

    allocate(zf_anom(k, ne), s(k), dz(k, ne))
    zf_anom = zf - spread(sum(zf, 2) / ne, 2, ne)
    s = loc_inf * matmul(zf_anom, zf_anom(k, :)) / (ne - 1)
    dz = spread(s, 2, ne) * spread(yo - zf(k, :), 1, k) / (r + s(k))
    deallocate(zf_anom, s)

  end function enkf_analysis

end module enkf_module
