program run_l63_var
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use random_module, only: &
    random_set_seed, rnorm => random_normal
  use l63_module, only: l63_nl, l63_ad
  use step_module, only: step_fom, step_adm
  implicit none

  character(len=*), parameter :: nml = "l63.nml"
  integer, parameter ::  seed = 514, ns = 3, nopts = 3, un = 31
  real(dp), allocatable :: xt(:, :), xb(:, :), yo(:, :), &
    d(:, :), l2(:), cost(:)
  integer :: nobs, i, j, k, m
  real(dp) :: loc_inf(1) = 1.0_dp
  character(len=256) :: fname

  integer :: nt, obs_int, ni
  real(dp) :: p, r, b, dt, bgd_b, a
  real(dp) ::  xt0(ns), xb0(ns), ad(ns), obs_r(ns), opts(nopts)

  namelist /l63/ p, r, b, dt, nt, obs_int, xt0, xb0, obs_r
  namelist /var/ bgd_b, ni, a

  open(unit=un, file=nml, status="old")
  read(unit=un, nml=l63)
  print l63
  read(unit=un, nml=var)
  print var
  close(unit=un)
  opts = [p, r, b]
  nobs = nt / obs_int

  allocate(xt(ns, nt), xb(ns, nt), yo(ns, nobs), &
    d(ns, nobs), l2(ni), cost(ni))

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

  xb(:, 1) = xb0

  do i = 1, ni
    do j = 2, nt
      xb(:, j) = step_fom(l63_nl, xb(:, j-1), 1, dt, opts)    
    end do
    do m = 1, nobs
      d(:, m) = xb(:, m * obs_int) - yo(:, m)
    end do
    ad = 0.0_dp
    m = nobs
    do j = nt, 1, -1
      ad = step_adm(l63_nl, l63_ad, xb(:, j), ad, 1, dt, opts)
      if (mod(j, obs_int) == 0) then
        ad = ad + d(:, m) / obs_r
        m = m - 1
      end if
    end do
    xb(:, 1) = xb(:, 1) - a * ad
    l2(i) = sqrt(sum((xb(:, 1) - xt0)**2)/ ns)
    cost(i) = calc_cost(xb(:, 1) - xb0, bgd_b, d, obs_r)
    print *, i, l2(i), cost(i)
  end do
  
  fname = "l63_var.dat"
  open(unit=un, file=fname, access="stream", &
    form="unformatted", status="replace", action="write")
  write(unit=un) xt
  write(unit=un) yo
  write(unit=un) xb
  close(un)

  open(unit=un, file="confplot_l63_var.R", access="stream", &
    form="formatted", status="replace", action="write")
  write(un, *) "ns <-", ns 
  write(un, *) "ni <-", ni 
  write(un, *) "nt <-", nt 
  write(un, *) "dt <-", dt 
  write(un, *) "nobs <-", nobs
  write(un, *) "obs.int <-", obs_int
  write(un, *) "a <-", a 
  write(un, *) "l2<-", l2(ni)
  close(un)

  deallocate(xt, xb, yo, d, l2, cost)

contains

  function calc_cost(dx, b, dy, r) result(cost)
    real(dp), intent(in) :: dx(:), dy(:, :)
    real(dp), intent(in) :: b, r(:)
    real(dp) :: cost

    integer :: i
    real(dp) :: Jo

    cost = 0.5_dp * sum(dx * dx) / b
    do i = 1, size(dy, 2)
      cost = cost + 0.5_dp * sum(dy(:, i) * dy(:, i) / r)
    end do

  end function calc_cost

end program run_l63_var
