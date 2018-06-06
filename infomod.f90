module infomod

  implicit none

  integer                                      :: nsta,nqc
  real*8                                       :: e0
  real*8, dimension(:,:), allocatable          :: zcoo,xcoo,energy
  character(len=80)                            :: aset,acalc,aprog,&
                                                  asymfile
  character(len=80), dimension(:), allocatable :: aqc
  character(len=16), dimension(:), allocatable :: asym

end module infomod
