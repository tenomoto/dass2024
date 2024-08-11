program test
  use random_module, only: random_test
  use l96_module, only: l96_test
  use io_module, only: io_test

  print *, "random:"
  call random_test()
  print *, "l96:"
  call l96_test(1.0d-5)
  print *, "io:"
  call io_test()
end program test
