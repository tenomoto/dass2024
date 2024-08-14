program run_l96_var
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use l96_module, only: l96_fom, l96_adm
  implicit none

  integer, parameter :: &
    un=31, un_xt = 51, un_yo = 52, un_xf = 53
  character(len=*), parameter :: &
    nml = "l96.nml", fname = "l96_var.dat"
  integer :: ns, nt_spinup, nt, nw, nc, ni, i, j, k
  real(dp) :: dt, F, b, r, a, gb, cost, cost_old, gnorm
  real(dp), allocatable ::  xt(:), xt0(:), x0(:), &
    ad(:), yo(:,:), xb(:,:), d(:,:), l2(:)
  character(len=256) :: xt_fname, yo_fname, xf_fname

  namelist /l96/ ns, dt, F, nt_spinup, nt, r
  namelist /var/ b, nw, ni, a, gb
  namelist /io/ xt_fname, yo_fname, xf_fname

  open(unit=un, file=nml, status="old")
  read(unit=un, nml=l96)
  print l96
  read(unit=un, nml=var)
  print var
  read(unit=un, nml=io)
  print io
  close(unit=un)
  nc = nt / nw
  print *, "nc=", nc

  allocate(xt(ns), xt0(ns), x0(ns), &
    ad(ns), yo(ns, nw), xb(ns, nw), d(ns, nw), l2(nc))
  open(unit=un_xt, file=xt_fname, access="stream", &
    form="unformatted", status="old", action="read")
  open(unit=un_yo, file=yo_fname, access="stream", &
    form="unformatted",  status="old", action="read")
  open(unit=un_xf, file=xf_fname, access="stream", &
    form="unformatted", status="old", action="read")
  read(unit=un_xf) x0
  do k = 1, nc
    print "(a, i5)", "cycle=", k
    read(unit=un_xt) xt0
    read(unit=un_yo) yo
    do j = 2, nw
      read(unit=un_xt) xt
    end do 
    xb(:, 1) = x0
    do i = 1, ni
      do j = 2, nw
        xb(:, j) = l96_fom(xb(:, j-1), 1, dt, F)
      end do
      d = xb - yo
      ad = 0.0_dp
      do j = nw, 1, -1
        ad = l96_adm(xb(:, j), ad, 1, dt, F) + d(:, j) / r
      end do
      xb(:, 1) = xb(:, 1) - a * ad
      cost = calc_cost(xb(:, 1) - x0, b, d, r)
      if (i==1) cost_old = cost
      l2(k) = sqrt(sum((xb(:, 1) - xt0)**2) / ns)
      gnorm = sum(ad * ad)
      print "(i5, ' J=', f10.3, ' Jold=', f10.3, ' g=', e10.3, ' l2=', e10.3)", &
        i, cost, cost_old, gnorm, l2(k)
      if (gnorm < gb) exit
      if (cost > cost_old) then               
        xb(:, 1) = xb(:, 1) + a * ad
        l2(k) = sqrt(sum((xb(:, 1) - xt0)**2) / ns)
        exit
      end if
      cost_old = cost  
    end do
    x0 = l96_fom(xb(:, 1), nw + 1, dt, F)
  end do
  close(un_xt)
  close(un_yo)
  close(un_xf)
  open(unit=un, file=fname, access="stream", &
    form="unformatted", status="replace", action="write")
  write(unit=un) l2
  close(un)
  deallocate(xt, xt0, x0, ad, yo, xb, d, l2)

  open(unit=un, file="confplot_l96_var.R", access="stream", &
    form="formatted", status="replace", action="write")
  write(un, *) "r <-", r
  write(un, *) "nw <-", nw
  write(un, *) "nc <-", nc
  write(un, *) "a <-", a
  close(un)

contains

  function calc_cost(dx, b, dy, r) result(cost)
    real(dp), intent(in) :: dx(:), dy(:, :)
    real(dp), intent(in) :: b, r
    real(dp) :: cost

    integer :: i
    real(dp) :: Jo

    cost = 0.5_dp * sum(dx * dx) / b
    do i = 1, size(dy, 2)
      cost = cost + 0.5_dp * sum(dy(:, i) * dy(:, i)) / r
    end do

  end function calc_cost

end program run_l96_var
