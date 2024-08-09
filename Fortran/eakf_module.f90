module eakf_module
  use iso_fortran_env
  implicit none

contains

  function eakf_analysis(zf, yo, r) result(dz)
    real(real64), dimension(:,:), intent(in) :: zf
    real(real64), intent(in) :: yo, r

    real(real64), dimension(:,:),  :: dz

    integer(int32) :: k, ne, i
    real(real64) :: yf_mean, alpha
    real(real64), dimension(:), allocatable :: &
      zf_mean, yf_anom, s
    real(real64), dimension(:, :), allocatable :: zf_anom

    k = size(zf, 1)
    ne = size(zf, 2)

    allocate(zf_mean(k), zf_anom(k, ne), &
      yf_anom(ne), s(k))
    zf_mean(:) = sum(zf, 2)
    do i = 1, ne
      zf_anom(:, i) = zf(:, i) - zf_mean(:)
    end do
    yf_anom(:) = zf_anom(k, :)
    yf_mean = zf_mean(k)
    s = matmul(zf_anom, zf_anom(k, :)) / (ne - 1)
    alpha = sqrt(r / (r + s(k)))
    dz = s / (r + s(k)) * (yo - zf_mean(k))
    dz = dz + (alpha - 1) / s(k) * &
      spread(s, 2, size(zf)) * spread(zf_anom(k, :), 1, size(s))
    deallocate(zf_mean, zf_anom, yf_anom, s)

  end function eakf_analysis

end module eakf_module
