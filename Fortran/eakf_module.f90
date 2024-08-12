module eakf_module
  use, intrinsic :: iso_fortran_env, dp => real64
  implicit none

contains

  function eakf_analysis(zf, yo, r, loc_inf) result(dz)
    real(dp), intent(in) :: zf(:, :)
    real(dp), intent(in) :: yo, r
    real(dp), intent(in) :: loc_inf(:)

    real(dp), dimension(:,:), allocatable :: dz

    integer :: k, ne, i
    real(dp) :: alpha
    real(dp), dimension(:), allocatable :: zf_mean, s, dz_mean
    real(dp), dimension(:, :), allocatable :: zf_anom

    k = size(zf, 1)
    ne = size(zf, 2)

    allocate(zf_mean(k), zf_anom(k, ne), s(k), &
      dz_mean(k), dz(k, ne))
    zf_mean(:) = sum(zf, 2) / ne
    do i = 1, ne
      zf_anom(:, i) = zf(:, i) - zf_mean(:)
    end do
    s = loc_inf * matmul(zf_anom, zf_anom(k, :)) / (ne - 1)
    alpha = sqrt(r / (r + s(k)))
    dz_mean  = s / (r + s(k)) * (yo - zf_mean(k))
    dz = spread(dz_mean, 2, ne) + (alpha - 1) / s(k) * &
      spread(s, 2, ne) * spread(zf_anom(k, :), 1, k)
    deallocate(zf_mean, zf_anom, s, dz_mean)

  end function eakf_analysis

end module eakf_module
