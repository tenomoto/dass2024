program test
  use random_module, only: random_test
  use l63_module, only: l63_test
  use l96_module, only: l96_test
  implicit none

  print *, "random:"
  call random_test()
  print *, "l63:"
  call l63_test(1.0d-5)
  print *, "l96:"
  call l96_test(1.0d-5)

end program test
