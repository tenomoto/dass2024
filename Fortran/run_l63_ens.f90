program run_l63_ens
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use random_module, only: &
    random_set_seed, rnorm => random_normal
  use l63_module, only: l63_nl
  use step_module, only: step_fom
  use enkf_module, only: enkf_analysis
  use eakf_module, only: eakf_analysis
  implicit none

  character(len=*), parameter :: nml = "l63.nml"
  integer, parameter ::  seed = 514, ns = 3, nopts = 3, un = 31
  real(dp), allocatable :: xt(:, :), yo(:, :), yo_p(:), &
    zf(:, :), dz(:, :), xm(:, :), l2(:)
  integer :: nobs, i, j, k, m
  real(dp) :: loc_inf(1) = 1.0_dp
  character(len=256) :: fname

  character(len=4) :: fil
  integer :: nt, obs_int, ne
  real(dp) :: p, r, b, dt, model_q
  real(dp) ::  xt0(ns), xb0(ns), obs_r(ns), opts(nopts)

  namelist /l63/ p, r, b, dt, model_q, nt, obs_int, xt0, xb0, obs_r
  namelist /ens/ ne, fil

  open(unit=un, file=nml, status="old")
  read(unit=un, nml=l63)
  print l63
  read(unit=un, nml=ens)
  print ens
  close(unit=un)
  opts = [p, r, b]
  nobs = nt / obs_int

  allocate(xt(ns, nt), yo(ns, nobs), yo_p(ne), &
    zf(ns+1, ne), dz(ns+1, ne), xm(ns, nt), l2(nt))

  call random_set_seed(seed)

  xt(:, 1) = xt0
  m = 1
  do k = 2, nt
    xt(:, k) = step_fom(l63_nl, xt(:, k-1), 1, dt, opts)
    if (mod(k, obs_int) == 0) then
      yo(:, m) = xt(:, k) + sqrt(obs_r) * rnorm(ns)
      m = m + 1
    end if
  end do
  do j = 1, ne
    zf(1:ns, j) = xb0 + sqrt(model_q) * rnorm(ns)
  end do

  m = 1
  do k = 1, nt
    if (mod(k, obs_int) == 0) then
      dz = 0.0d0
      do i = 1, ns
        zf(ns+1, :) = zf(i, :)
        if (fil == "enkf") then
          yo_p = rnorm(ns, yo(i, m), obs_r(i))
          yo_p = yo_p - sum(yo_p) / ne + yo(i, m)
          dz = dz + enkf_analysis(zf, yo_p, obs_r(i), loc_inf)
        else
          dz = dz + eakf_analysis(zf, yo(i, m), obs_r(i), loc_inf)
        end if
      end do
      zf = zf + dz
      m = m + 1
    end if
    xm(:, k) = sum(zf(1:ns, :), 2) / ne
    l2(k) = sqrt(sum((xm(:, k) - xt(:, k))**2) / ns)
    print *, k, l2(k)
    if (k < nt) then
      do j = 1, ne
        zf(1:ns, j) = step_fom(l63_nl, zf(1:ns, j), 1, dt, opts)
      end do
    end if
  end do

  fname = "l63_" // fil // ".dat"
  open(unit=un, file=fname, access="stream", &
    form="unformatted", status="replace", action="write")
  write(unit=un) xt
  write(unit=un) yo
  write(unit=un) xm
  close(un)

  open(unit=un, file="confplot_l63_ens.R", access="stream", &
    form="formatted", status="replace", action="write")
  write(un, *) "ns <-", ns 
  write(un, *) "ne <-", ne 
  write(un, *) "nt <-", nt 
  write(un, *) "dt <-", dt 
  write(un, *) "nobs <-", nobs
  write(un, *) "obs.int <-", obs_int
  write(un, *) "fil <-", "'", fil, "'"
  close(un)

  deallocate(xt, yo, yo_p, zf, dz, xm, l2)

end program run_l63_ens
