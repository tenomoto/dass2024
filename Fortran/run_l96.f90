program run_l96
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use io_module, only: io_write_binary, io_read_binary, io_delete_file
  use random_module, only: random_set_seed, random_normal
  use l96_module, only: l96_fom
  implicit none

  integer :: ns, seed, nstep, ne, i
  real(dp) :: dt, F
  real(dp), allocatable :: x(:)
  character(len=64) :: buf, fname

  character(len=*), parameter :: nml = "l96.nml"
  integer, parameter :: un = 31
  namelist /random/ seed
  namelist /l96/ ns, dt, F

  if (iargc() < 2) then
    print *, "Usage :: ./run_l96 nstep fname [nens]"
    stop
  end if
  call getarg(1, buf)
  read(unit=buf, fmt=*) nstep
  call getarg(2, fname)
  print *, "nstep=", nstep, " fname=", fname
  if (iargc() > 2) then
    call getarg(3, buf)
    read(unit=buf, fmt=*) ne
  else
    ne = 1
  end if
  print *, "ne=", ne

  open(unit=un, file=nml, status="old")
  read(unit=un, nml=random)
  read(unit=un, nml=l96)
  close(unit=un)
  print random
  print l96

  allocate(x(ns))
  call random_set_seed(seed)
  call io_delete_file(fname)
  do i = 1, ne
    x = random_normal(ns)
    x = l96_fom(x, nstep, dt, F)
    call io_write_binary(fname, x)
  end do
  deallocate(x)

end program run_l96
