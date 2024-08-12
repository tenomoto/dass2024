program prep_l96
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use random_module, only: random_set_seed, rnorm => random_normal
  use l96_module, only: l96_fom
  implicit none

  integer, parameter :: &
    un = 31, un_xt = 51, un_yo = 52, un_xf = 53, seed = 514
  integer :: ns, nt_spinup, nt, ne, i
  character(len=256) :: xt_fname, yo_fname, xf_fname
  real(dp) :: dt, F, r, c_loc, infl
  real(dp), allocatable :: xt(:), xf(:), yo(:)
  character(len=4) :: fil

  character(len=*), parameter :: nml = "l96.nml"
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

  allocate(xt(ns), xf(ns), yo(ns))
  call random_set_seed(seed)
  
  open(unit=un_xt, file=xt_fname, access="stream", &
    form="unformatted", status="replace", action="write")
  open(unit=un_yo, file=yo_fname, access="stream", &
    form="unformatted", status="replace", action="write")
  open(unit=un_xf, file=xf_fname, access="stream", &
    form="unformatted", status="replace", action="write")
  xt = l96_fom(rnorm(ns), nt_spinup, dt, F)
  do i = 1, ne
    xf = l96_fom(rnorm(ns), nt_spinup, dt, F)
    write(unit=un_xf) xf
  end do
  do i = 1, nt
    yo = xt + sqrt(r) * rnorm(ns)
    write(unit=un_xt) xt
    write(unit=un_yo) yo
    if (i < nt) then
      xt = l96_fom(xt, 1, dt, F)
    end if
  end do
  deallocate(xt, xf, yo)
  close(unit=un_xt)
  close(unit=un_yo)
  close(unit=un_xf)

end program prep_l96
