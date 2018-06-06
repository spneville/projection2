module globalmod

  save
  
  integer                                     :: natm,ncoo
  real*8, dimension(:), allocatable           :: xcoo0,mass
  real*8, dimension(:,:), allocatable         :: sthess,cartz,zcart
  character(len=2), dimension(:), allocatable :: atlbl
  character(len=120)                          :: ain,adat

contains

!#######################################################################

  subroutine rddat

    implicit none
    
    integer                           :: unit,i,j
    real*8, dimension(:), allocatable :: freq
    logical                           :: lopen
    character(len=120)                :: header
      
!-----------------------------------------------------------------------
! Open the data file
!-----------------------------------------------------------------------
    do i=10,1000
       inquire(unit=i,opened=lopen)
       if (.not.lopen) then
          unit=i
          exit
       endif
    enddo
    
    open(unit,file=adat,form='unformatted',status='old')

!-----------------------------------------------------------------------
! Read the data file
!-----------------------------------------------------------------------
    ! Read the file header
    write(unit) header

    ! Read the system size and allocate arrays
    read(unit) ncoo
      
    natm=ncoo/3

    if (allocated(xcoo0)) deallocate(xcoo0)
    if (allocated(mass)) deallocate(mass)
    if (allocated(atlbl)) deallocate(atlbl)
    if (allocated(sthess)) deallocate(sthess)
    if (allocated(cartz)) deallocate(cartz)
    if (allocated(zcart)) deallocate(zcart)

    allocate(xcoo0(ncoo))
    allocate(mass(ncoo))
    allocate(atlbl(natm))
    allocate(sthess(ncoo,ncoo))
    allocate(cartz(ncoo,ncoo))
    allocate(zcart(ncoo,ncoo))

    ! Reference geometry
    read(unit) xcoo0
    
    ! Masses
    read(unit) mass
    
    ! Atom labels
    read(unit) atlbl
    
    ! Coordinate vectors
    read(unit) zcart

    ! Similarity transformed Hessian
    read(unit) sthess

!------------------------------------------------------------------
! Frequency scale (in eV) and mass scale (in amu) the
! mass-weighted-Cartesian-to-zeta transformation matrix s.t. it
! corresponds to the transformation from non-mass-weighted
! Cartesians to our frequency-scaled coordinates
!
! cartz: transformation FROM non-mass-weighted Cartesians
! zcart: transformation TO non-mass-weighted Cartesians
!------------------------------------------------------------------
    ! 'Frequencies' in eV
    allocate(freq(ncoo))
    do i=1,ncoo
       freq(i)=sqrt(abs(sthess(i,i)))*0.6373641d0
    enddo

    ! Transformation matrices
    do i=1,ncoo
       do j=1,ncoo
          cartz(i,j)=zcart(j,i)*(15.4644d0*sqrt(freq(i)*mass(j)))
          zcart(j,i)=zcart(j,i)/(15.4644d0*sqrt(freq(i)*mass(j)))
       enddo
    enddo

!-----------------------------------------------------------------------
! Close the data file
!-----------------------------------------------------------------------
    close(unit)
      
    return

  end subroutine rddat

!#######################################################################

end module globalmod
