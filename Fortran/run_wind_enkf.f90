program run_wind_enkf
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use io_module, only: io_write_binary, io_delete_file
  use random_module, only: random_set_seed, rnorm => random_normal
  use enkf_module, only: enkf_analysis
  use eakf_module, only: eakf_analysis
  implicit none

  integer :: i
  real(dp), allocatable :: yo_p(:), u(:), v(:), zf(:,:), za(:,:)

  character(len=*), parameter :: nml = "wind.nml"
  integer, parameter :: un = 31, k = 3
  integer :: seed, ne
  real(dp) :: sf, uf, vf, yo, so, r
  character(len=4) :: fil
  namelist /random/ seed
  namelist /wind/ ne, sf, uf, vf, yo, so, fil

  open(unit=un, file=nml, status="old")
  read(unit=un, nml=random)
  read(unit=un, nml=wind)
  close(unit=un)
  print random
  print wind

  r = so**2
 allocate(yo_p(ne), u(ne), v(ne), zf(k, ne), za(k, ne))
  u = rnorm(ne, uf, sf)
  zf(1, :) = u - sum(u) / ne + uf
  v = rnorm(ne, vf, sf)
  zf(2, :) = v - sum(v) / ne + vf
  zf(3, :) = calc_speed(u, v)
  if (fil=="enkf") then
    yo_p = rnorm(ne, yo, so)
    yo_p = yo_p - sum(yo_p) / ne + yo
    za = zf + enkf_analysis(zf, yo_p, r)
  else
    za = zf + eakf_analysis(zf, yo, r)
  end if
  call io_delete_file("zf.dat")
  call io_delete_file("za.dat")
  do i = 1, ne 
    call io_write_binary("zf.dat", zf(:, i))
    call io_write_binary("za.dat", za(:, i))
  end do
  deallocate(yo_p, u, v, zf, za)

contains

  function calc_speed(u, v) result(us)
    real(dp), intent(in) :: u(:), v(:)
    real(dp) :: us(size(u))

    us = sqrt(u * u + v * v)

  end function calc_speed


end program run_wind_enkf
