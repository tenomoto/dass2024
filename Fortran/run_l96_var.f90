program run_l96
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use io_module, only: io_write_binary, io_read_binary, io_delete_file
  use random_module, only: random_set_seed, rnorm => random_normal
  use l96_module, only: l96_fom, l96_adm
  implicit none

  integer :: seed, ns, nt, nw, nc, ni, i, j, k
  real(dp) :: dt, F, b, r, a, gb, cost, cost_old, gnorm
  real(dp), allocatable ::  xt(:), xt0(:), x0(:), &
    ad(:), yo(:,:), xb(:,:), d(:,:), l2(:)
  logical :: save_state

  character(len=*), parameter :: nml = "l96.nml"
  integer, parameter :: un = 31
  namelist /random/ seed
  namelist /l96/ ns, dt, F, nt
  namelist /var/ b, r, nw, nc
  namelist /opt/ ni, a, gb
  namelist /io/ save_state

  open(unit=un, file=nml, status="old")
  read(unit=un, nml=random)
  print random
  read(unit=un, nml=l96)
  print l96
  read(unit=un, nml=var)
  print var
  read(unit=un, nml=opt)
  print opt
  read(unit=un, nml=io)
  print io
  close(unit=un)

  if (save_state) then
    call io_delete_file("xt.dat")
    call io_delete_file("x0.dat")
    call io_delete_file("xa.dat")
  endif

  allocate(xt(ns), xt0(ns), x0(ns), &
    ad(ns), yo(ns, nw), xb(ns, nw), d(ns, nw), l2(nc))
  call random_set_seed(seed)
  xt = l96_fom(rnorm(ns), nt, dt, F)
  x0 = l96_fom(rnorm(ns), nt, dt, F)
  do k = 1, nc
    print "(a, i5)", "cycle=", k
    xt0 = xt
    if (save_state) then
      call io_write_binary("xt.dat", xt)
      call io_write_binary("x0.dat", x0)
    end if
    do j = 1, nw
      yo(:, j) = xt + sqrt(r) * rnorm(ns)
      xt = l96_fom(xt, 1, dt, F)
    end do 
    if (save_state) then
      call io_write_binary("yo.dat", yo(:, 1))
    end if
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
    if (save_state) then
      call io_write_binary("xa.dat", xb(:, 1))
    end if
    x0 = l96_fom(xb(:, 1), nw + 1, dt, F)
  end do
  call io_delete_file("l2.dat")
  call io_write_binary("l2.dat", l2)
  deallocate(xt, xt0, x0, ad, yo, xb, d, l2)

  open(unit=un, file="config.R", access="stream", form="formatted", & 
    status="replace", action="write")
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

end program run_l96
