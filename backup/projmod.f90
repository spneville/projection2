  module projmod

    implicit none

    save

    integer                                      :: nvec
    real*8, dimension(:,:), allocatable          :: vec,hess,peigvec
    character(len=60), dimension(:), allocatable :: avec
    character(len=120)                           :: axyz,ahess

  end module projmod
