use, intrinsic :: iso_fortran_env, only: dp => real64
real(dp) :: x(5) = [1, 2, 3, 4, 5], y(10)
logical e
open(unit=51, file="x.dat", access="stream", form="unformatted", &
  status="replace", action="write")
write(unit=51) x
close(unit=51)
open(unit=51, file="x.dat", access="stream", form="unformatted", &
  status="old", action="write", position="append")
write(unit=51) x(5:1:-1)
close(unit=51)
open(unit=51, file="x.dat", access="stream", form="unformatted", &
  status="old", action="read")
read(unit=51) y
close(unit=51)
print *, y
end
