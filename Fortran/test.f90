program test
  use random_module, only: random_test
!  use l96_module, only: l96_test
!  use io_module, only: io_test

  call random_test()
!  call l96_test(1.0d-5)
!  call io_test()
end program test
