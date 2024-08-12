module io_module
  use, intrinsic :: iso_fortran_env, only: &
    dp => real64, fss => file_storage_size
  implicit none

  integer, parameter :: u = 51

contains

  subroutine io_write_binary(fname, x)
    character(len=*), intent(in) :: fname
    real(dp), intent(in) :: x(:)
    logical :: e
    
    inquire(file=fname, exist=e)
    if (.not.e) then
      open(unit=u, file=fname, access="stream", form="unformatted", &
        status="new", action="write")
    else
      open(unit=u, file=fname, access="stream", form="unformatted", &
        status="old", action="write", position="append")
    end if
    write(unit=u) x
    close(unit=u)

  end subroutine io_write_binary

  subroutine io_read_binary(fname, x, skip)
    character(len=*), intent(in) :: fname
    real(dp), intent(inout) :: x(:)
    integer, intent(in), optional :: skip

    open(unit=u, file=fname, access="stream", form="unformatted", &
      status="old", action="read")
    if (.not. present(skip)) then
      read(unit=u) x
    else
      read(unit=u, pos=skip*fss+1) x
    end if
    close(unit=u)

  end subroutine io_read_binary

  subroutine io_delete_file(fname)
    character(len=*), intent(in) :: fname

    open(unit=u, file=fname, status="unknown")
    close(unit=u, status="delete")

  end subroutine io_delete_file

  subroutine io_test()

    integer, parameter :: n = 5
    real(dp) :: x(n) = [1, 2, 3, 4, 5], y(n), w(n*2)
    character(len=*), parameter :: fname = "x.dat"

    print *, "x=", x
    call io_delete_file(fname)
    call io_write_binary(fname, x)
    call io_read_binary(fname, y)
    print *, "y=", y
    call io_write_binary(fname, x(5:1:-1))
    call io_read_binary(fname, w)
    print *, "w=", w

  end subroutine io_test

end module io_module

