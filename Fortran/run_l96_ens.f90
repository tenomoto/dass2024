program run_l96_ens
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use random_module, only: rnorm => random_normal
  use step_module, only: step_fom
  use l96_module, only: l96_nl
  use enkf_module, only: enkf_analysis
  use eakf_module, only: eakf_analysis
  implicit none

  integer, parameter :: &
    un = 31, un_xt = 51, un_yo = 52, un_xf = 53, seed = 516
  character(len=*), parameter ::  nml = "l96.nml"
  integer :: ns, ne, nt_spinup, nt, i, j, k
  character(len=256) :: xt_fname, yo_fname, xf_fname, fname
  real(dp) :: dt, F, r, e(5), d, c_loc, infl
  real(dp) :: opts(1)
  real(dp), allocatable :: e1(:), e2(:), st(:), &
    xt(:), yo(:), yo_p(:), zf(:, :), dz(:, :), loc_inf(:)
  character(len=4) :: fil

  namelist /l96/ ns, dt, F, nt_spinup, nt, r
  namelist /ens/ ne, fil, c_loc, infl
  namelist /io/ xt_fname, yo_fname, xf_fname

  open(unit=un, file=nml, status="old")
  read(unit=un, nml=l96)
  print l96
  read(unit=un, nml=ens)
  print ens
  read(unit=un, nml=io)
  print io
  close(unit=un)
  opts(1) = F

  allocate(xt(ns), yo(ns), zf(ns + 1, ne), dz(ns + 1, ne), &
    e1(nt), e2(nt), st(nt), loc_inf(ns + 1))
  if (fil == "enkf") then
    allocate(yo_p(ne))
  end if
  open(unit=un_xf, file=xf_fname, access="stream", form="unformatted", &
        status="old", action="read")
  read(unit=un_xf) zf(1:ns, :)
  close(un_xf)
  open(unit=un_xt, file=xt_fname, access="stream", form="unformatted", &
        status="old", action="read")
  open(unit=un_yo, file=yo_fname, access="stream", form="unformatted", &
        status="old", action="read")
  do k = 1, nt
    read(unit=un_xt) xt
    read(unit=un_yo) yo
    dz = 0.0_dp
    do i = 1, ns
      zf(ns + 1, :) = zf(i, :)
      do j = 1, ns
        d = min(abs(j - i), ns - abs(j - i))
        loc_inf(j) = exp(-0.5_dp * (d / c_loc)**2)
      end do
      loc_inf(ns + 1) = infl
      if (fil == "enkf") then
        yo_p = rnorm(ne, yo(i), sqrt(r))
        yo_p = yo_p - sum(yo_p) / ne + yo(i)
        dz = dz + enkf_analysis(zf, yo_p, r, loc_inf)
      else
        dz = dz + eakf_analysis(zf, yo(i), r, loc_inf)
      end if
    end do
    zf = zf + dz
    e = calc_error(zf(1:ns, :), xt)
    e1(k) = e(1)
    e2(k) = e(2)
    st(k) = calc_sd(zf(1:ns, :))
    if (k < nt) then
      do j = 1, ne
        zf(1:ns, j) = step_fom(l96_nl, zf(1:ns, j), 1, dt, opts)
      end do
    end if
  end do
  close(un_xt)
  close(un_yo)

  fname = "l96_" // fil // ".dat"
  open(unit=un, file=fname, access="stream", &
    form="unformatted", status="replace", action="write")
  write(unit=un) e1
  write(unit=un) e2
  write(unit=un) st
  close(un)
  deallocate(xt, yo, zf, e1, e2, st, loc_inf)

  open(unit=un, file="confplot_l96_ens.R", access="stream", &
    form="formatted", status="replace", action="write")
  write(un, *) "nt <-", nt
  write(un, *) "r <-", r
  write(un, *) "fil <-", "'", fil, "'"
  write(un, *) "ne <-", ne
  write(un, *) "infl <-", infl
  write(un, *) "c.loc <-", c_loc
  close(un)

contains

  function calc_error(x, xt) result(e)
    real(dp), intent(in) :: x(:, :), xt(:)
    real(dp) :: e(5)

    real(dp), allocatable :: dx(:, :)
    integer :: ns, ne, j

    ns = size(x, 1)
    ne = size(x, 2)
    allocate(dx(ns, ne))
    
    do j = 1, ne
      dx(:, j) = x(:, j) - xt
    end do
    e(1) = sqrt(sum((sum(dx, 2) / ne)**2) / ns)
    e(2) = sum(sqrt(sum(dx**2, 1) / ns)) / ne
    e(3) = e(1) / e(2)
    e(4) = sqrt(0.5_dp * (ns + 1.0_dp) / ns)
    e(5) = e(3) / e(4)
    deallocate(dx)

  end function calc_error

  function calc_sd(x) result(sd)
    real(dp), intent(in) :: x(:, :)
    real(dp) :: sd

    integer :: ns, ne

    ns = size(x, 1)
    ne = size(x, 2)
    sd = sum(sqrt((sum(x**2, 2) - sum(x, 2)**2 / ne) / (ne - 1.0_dp))) / ns

  end function calc_sd

end program run_l96_ens
