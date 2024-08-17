program run_wind_ens
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use random_module, only: random_set_seed, rnorm => random_normal
  use enkf_module, only: enkf_analysis
  use eakf_module, only: eakf_analysis
  implicit none

  integer :: i
  real(dp), allocatable :: &
    yo_p(:), u(:), v(:), zf(:,:), za(:,:), loc_inf(:)

  character(len=*), parameter :: nml = "wind.nml"
  integer, parameter :: un = 31, k = 3, un_zf=51, un_za=52
  integer :: seed, ne
  real(dp) :: sf, uf, vf, yo, so, r
  character(len=4) :: fil
  namelist /wind/ seed, ne, sf, uf, vf, yo, so, fil

  open(unit=un, file=nml, status="old")
  read(unit=un, nml=wind)
  close(unit=un)
  print wind

  r = so**2
  allocate(yo_p(ne), u(ne), v(ne), zf(k, ne), za(k, ne), loc_inf(3))
  loc_inf = 1.0_dp
  u = rnorm(ne, uf, sf)
  zf(1, :) = u - sum(u) / ne + uf
  v = rnorm(ne, vf, sf)
  zf(2, :) = v - sum(v) / ne + vf
  zf(3, :) = calc_speed(u, v)
  if (fil=="enkf") then
    yo_p = rnorm(ne, yo, so)
    yo_p = yo_p - sum(yo_p) / ne + yo
    za = zf + enkf_analysis(zf, yo_p, r, loc_inf)
  else
    za = zf + eakf_analysis(zf, yo, r, loc_inf)
  end if
  open(unit=un_zf, file="zf.dat", access="stream", form="unformatted", &
        status="replace", action="write")
  open(unit=un_za, file="za.dat", access="stream", form="unformatted", &
        status="replace", action="write")
  do i = 1, ne 
    write(unit=un_zf) zf(:, i)
    write(unit=un_za) za(:, i)
  end do
  deallocate(yo_p, u, v, zf, za)
  close(un_zf)
  close(un_za)

  open(unit=un, file="confplot_wind.R", access="stream", form="formatted", & 
    status="replace", action="write")
  write(un, *) "ne <-", ne
  write(un, *) "k <-", k
  write(un, *) "yo <-", yo
  write(un, *) "so <-", so
  write(un, *) "fil <-'", fil, "'"
  close(un)

contains

  function calc_speed(u, v) result(us)
    real(dp), intent(in) :: u(:), v(:)
    real(dp) :: us(size(u))

    us = sqrt(u * u + v * v)

  end function calc_speed


end program run_wind_ens
