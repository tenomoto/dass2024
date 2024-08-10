module enkf_module
  use, intrinsic :: iso_fortran_env, dp => real64
  implicit none

contains

  function enkf_analysis(zf, yo, r) result(dz)
    real(dp), dimension(:,:), intent(in) :: zf
    real(dp), dimension(:), intent(in) :: yo
    real(dp), intent(in) :: r
    real(dp), dimension(:,:), allocatable :: dz

    integer :: k, ne!, i
    real(dp), dimension(:), allocatable ::s
    real(dp), dimension(:, :), allocatable :: zf_anom

    k = size(zf, 1)
    ne = size(zf, 2)

    allocate(zf_anom(k, ne), s(k), dz(k, ne))
    zf_anom = zf - spread(sum(zf, 2) / ne, 2, ne)
    s = matmul(zf_anom, zf_anom(k, :)) / (ne - 1)
    dz = spread(s, 2, ne) * spread(yo - zf(k, :), 1, k) / (r + s(k))
    deallocate(zf_anom, s)

  end function enkf_analysis

end module enkf_module
