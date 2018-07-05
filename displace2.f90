program displace2

  use globalmod

  implicit none

!-----------------------------------------------------------------------
! Read the input file
!-----------------------------------------------------------------------
  call rdinput

!-----------------------------------------------------------------------
! Read the data file and set up the transformation matrices
!-----------------------------------------------------------------------
  call rddat

!-----------------------------------------------------------------------
! Make the cut
!-----------------------------------------------------------------------
  call calc_cut

contains

!#######################################################################

  subroutine rdinput
      
    use globalmod
    use parsemod
    use dispmod
    
    implicit none

    integer            :: i,k,iin
    character(len=120) :: errmsg

!-----------------------------------------------------------------------
! Set defaults
!-----------------------------------------------------------------------
    ain=''
    adat=''
    astem=''
    cut=-1.0d0
    cut2=-1.0d0
    cutdiag=-1.0d0
    lcut2d=.false.
    lcutdiag=.false.
    
!-----------------------------------------------------------------------
! Read the name of the input file and open this file
!-----------------------------------------------------------------------
    call getarg(1,ain)

    if (ain.eq.'') then
       write(6,'(/,2x,a,/)') 'The name of the input file has not &
            been given'
       STOP
    endif

    k=index(ain,'.inp')
    if (k.eq.0) ain=trim(ain)//'.inp'

    iin=20
    open(iin,file=ain,form='formatted',status='old')

!-----------------------------------------------------------------------
! Read input file
!-----------------------------------------------------------------------
5   continue
    call rdinp(iin)
      
    i=0
    if (keyword(1).ne.'end-input') then
10     continue
       i=i+1
       
       if (keyword(i).eq.'cut') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             read(keyword(i),*) cut(1)
             if (keyword(i+1).eq.',') then
                i=i+2
                read(keyword(i),*) cut(2)
             else
                write(6,'(/,2x,a,/)') 'The no. points has not been &
                     given with the cut keyword'
                STOP
             endif
             if (keyword(i+1).eq.',') then
                i=i+2
                read(keyword(i),*) cut(3)
             else
                write(6,'(/,2x,a,/)') 'The point spacing has not &
                     been given with the cut keyword'
                STOP
             endif
          else
             goto 100
          endif
            
       else if (keyword(i).eq.'cut2d') then
          lcut2d=.true.
          if (keyword(i+1).eq.'=') then
             i=i+2
             read(keyword(i),*) cut(1)
             if (keyword(i+1).eq.',') then
                i=i+2
                read(keyword(i),*) cut(2)
             else
                write(6,'(/,2x,a,/)') 'The 1st no. points has not been &
                     given with the cut2d keyword'
                STOP
             endif
             if (keyword(i+1).eq.',') then
                i=i+2
                read(keyword(i),*) cut(3)
             else
                write(6,'(/,2x,a,/)') 'The 1st point spacing has not &
                     been given with the cut2d keyword'
                STOP
             endif
             if (keyword(i+1).eq.',') then
                i=i+2
                read(keyword(i),*) cut2(1)
             else
                write(6,'(/,2x,a,/)') 'The 2nd mode number has not been &
                     given with the cut2d keyword'
                STOP
             endif
             if (keyword(i+1).eq.',') then
                i=i+2
                read(keyword(i),*) cut2(2)
             else
                write(6,'(/,2x,a,/)') 'The 2nd no. points has not been &
                     given with the cut2d keyword'
                STOP
             endif
             if (keyword(i+1).eq.',') then
                i=i+2
                read(keyword(i),*) cut2(3)
             else
                write(6,'(/,2x,a,/)') 'The 2nd point spacing has not &
                     been given with the cut2d keyword'
                STOP
             endif
          else
             goto 100
          endif
          
       else if (keyword(i).eq.'cutdiag') then
          lcutdiag=.true.
          if (keyword(i+1).eq.'=') then
             i=i+2
             read(keyword(i),*) cutdiag(1)
             if (keyword(i+1).eq.',') then
                i=i+2
                read(keyword(i),*) cutdiag(2)
             else
                write(6,'(/,2x,a,/)') 'The 1st no. point spacing has not &
                     been given with the cutdiag keyword'
                STOP
             endif
             if (keyword(i+1).eq.',') then
                i=i+2
                read(keyword(i),*) cutdiag(3)
             else
                write(6,'(/,2x,a,/)') 'The 2nd mode number has not &
                     been given with the cutdiag keyword'
                STOP
             endif
             if (keyword(i+1).eq.',') then
                i=i+2
                read(keyword(i),*) cutdiag(4)
             else
                write(6,'(/,2x,a,/)') 'The 2nd no. point spacing has not &
                     been given with the cutdiag keyword'
                STOP
             endif
             if (keyword(i+1).eq.',') then
                i=i+2
                read(keyword(i),*) cutdiag(5)
             else
                write(6,'(/,2x,a,/)') 'The no. points has not &
                     been given with the cutdiag keyword'
                STOP
             endif
          else
             goto 100
          endif
          
       else if (keyword(i).eq.'data_file') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             adat=keyword(i)
          else
             goto 100
          endif

       else if (keyword(i).eq.'filestem') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             astem=keyword(i)
          else
             goto 100
          endif

       else
          ! Exit if the keyword is not recognised
          errmsg='Unknown keyword: '//trim(keyword(i))
          write(6,'(/,2x,a,/)') errmsg
          STOP
       endif
       
       ! If there are more keywords to be read on the current line,
       ! then read them, else read the next line
       if (i.lt.inkw) then
          goto 10
       else
          goto 5
       endif
         
       ! Exit if a required argument has not been given with a keyword
100    continue
       errmsg='No argument given with the keyword '//trim(keyword(i))
       write(6,'(/,2x,a,/)') errmsg
       STOP
         
    endif

!-----------------------------------------------------------------------
! Close the input file
!-----------------------------------------------------------------------
    close(iin)

!-----------------------------------------------------------------------
! Check that all required information has been given
!-----------------------------------------------------------------------
    if (cut(1).eq.-1.0d0.and.cutdiag(1).eq.-1.0d0) then
       write(6,'(/,2x,a,/)') 'No cut keywords have been given'
       STOP
    endif

    if (adat.eq.'') then
       write(6,'(/,2x,a,/)') 'The name of the data file has not &
            been given'
       STOP
    endif
      
    if (astem.eq.'') then
       write(6,'(/,2x,a,/)') 'The filestem has not been given'
       STOP
    endif

    return
      
  end subroutine rdinput

!#######################################################################

  subroutine transmat
    
    use globalmod
    
    implicit none
      
    integer                 :: i,j
    real*8, dimension(ncoo) :: freq

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
    do i=1,ncoo
       freq(i)=sqrt(abs(sthess(i,i)))*0.6373641d0
    enddo

    ! Transformation matrices
    do i=1,ncoo
       do j=1,ncoo
          zcart(i,j)=cartz(j,i)*(15.4644d0*sqrt(freq(i)*mass(j)))
          cartz(j,i)=cartz(j,i)/(15.4644d0*sqrt(freq(i)*mass(j)))
       enddo
    enddo

    return

  end subroutine transmat

!#######################################################################

  subroutine calc_cut
      
    use globalmod
    use dispmod

    implicit none

    if (lcut2d) then
       call cut2d
    else if (lcutdiag) then
       call cut2d_diag
    else
       call cut1d
    endif

    return

  end subroutine calc_cut

!#######################################################################

  subroutine cut1d
    
    use globalmod
    use dispmod

    implicit none

    integer                 :: i,j,k,unit1,unit2,ni,nf,iz
    real*8                  :: dz
    real*8, dimension(ncoo) :: z,x
    character(len=120)      :: aout
      
    unit1=20
    unit2=21
      
    z=0.0d0

    ni=-cut(2)
    nf=cut(2)
    iz=cut(1)
    dz=cut(3)

    open(unit1,file='all.xyz',form='formatted',status='unknown')

    do i=ni,nf
         
       ! Calculate the current Cartesian coordinates
       z(iz)=i*dz
       x=xcoo0+matmul(zcart,z)

       ! Write the current file name
       call filename1d(i,iz,aout)

       ! Write the current Cartesian coordinates to file
       open(unit2,file=aout,form='formatted',status='unknown')
       write(unit1,'(i3,/)') natm
       write(unit2,'(i3,/)') natm
       do j=1,natm
          write(unit1,'(a2,3(2x,F10.7))') atlbl(j),(x(k),k=j*3-2,j*3)
          write(unit2,'(a2,3(2x,F10.7))') atlbl(j),(x(k),k=j*3-2,j*3)
       enddo
       close(unit2)

    enddo

    close(unit1)

    return

  end subroutine cut1d

!#######################################################################

  subroutine filename1d(i,iz,aout)

    use globalmod
    use dispmod

    implicit none
      
    integer :: i,iz,k
    character(len=120) :: aout
    
    aout=trim(astem)//'_'
    k=len_trim(astem)+2

    if (iz.lt.10) then
       write(aout(k:k+2),'(a1,i1,a1)') '0',iz,'_'
    else
       write(aout(k:k+2),'(i2,a1)') iz,'_'
    endif
    k=k+3

    if (abs(i).lt.10) then
       write(aout(k:k+1),'(a1,i1)') '0',abs(i)
    else
       write(aout(k:k+1),'(i2)') abs(i)
    endif
    k=k+2
      
    if (i.lt.0) then
       write(aout(k:k+4),'(a5)') 'l.xyz'
    else
       write(aout(k:k+4),'(a5)') 'r.xyz'
    endif
      
    return

  end subroutine filename1d

!#######################################################################

  subroutine cut2d

    use globalmod
    use dispmod

    implicit none
      
    integer                 :: i1,i2,ni1,nf1,iz1,ni2,nf2,iz2,unit,j,k
    real*8                  :: dz1,dz2
    real*8, dimension(ncoo) :: z,x
    character(len=120)      :: aout
    
    unit=20
      
    z=0.0d0

    ni1=-cut(2)
    nf1=cut(2)
    iz1=cut(1)
    dz1=cut(3)

    ni2=-cut2(2)
    nf2=cut2(2)
    iz2=cut2(1)
    dz2=cut2(3)

    do i1=ni1,nf1
       do i2=ni2,nf2
          
          ! Calculate the current Cartesian coordinates
          z(iz1)=i1*dz1
          z(iz2)=i2*dz2
          x=xcoo0+matmul(zcart,z)

          ! Write the current file name
          call filename2d(i1,iz1,i2,iz2,aout)

          ! Write the current Cartesian coordinates to file
          open(unit,file=aout,form='formatted',status='unknown')
          write(unit,'(i3,/)') natm
          do j=1,natm
             write(unit,'(a2,3(2x,F10.7))') atlbl(j),(x(k),k=j*3-2,j*3)
          enddo
          close(unit)

       enddo
    enddo
    
    return

  end subroutine cut2d

!#######################################################################

  subroutine cut2d_diag
    
    use globalmod
    use dispmod

    implicit none

    integer                 :: i,j,k,unit1,unit2,iz1,iz2,nf,i1,i2
    real*8                  :: dz1,dz2
    real*8, dimension(ncoo) :: z,x
    character(len=120)      :: aout
      
    unit1=20
    unit2=21
      
    z=0.0d0
    
    iz1=cutdiag(1)
    dz1=cutdiag(2)    
    iz2=cutdiag(3)
    dz2=cutdiag(4)
    nf=cutdiag(5)
    
    open(unit1,file='all.xyz',form='formatted',status='unknown')

    do i=0,nf

       ! Calculate the current Cartesian coordinates
       z(iz1)=i*dz1
       z(iz2)=i*dz2
       x=xcoo0+matmul(zcart,z)
       
       ! Write the current file name
       i1=int(sign(dble(i),dz1))
       i2=int(sign(dble(i),dz2))
       call filename2d(i1,iz1,i2,iz2,aout)
       
       ! Write the current Cartesian coordinates to file
       open(unit2,file=aout,form='formatted',status='unknown')
       write(unit1,'(i3,/)') natm
       write(unit2,'(i3,/)') natm
       do j=1,natm
          write(unit1,'(a2,3(2x,F10.7))') atlbl(j),(x(k),k=j*3-2,j*3)
          write(unit2,'(a2,3(2x,F10.7))') atlbl(j),(x(k),k=j*3-2,j*3)
       enddo
       close(unit2)

    enddo

    close(unit1)
    
    return
    
  end subroutine cut2d_diag
    
!#######################################################################

  subroutine filename2d(i1,iz1,i2,iz2,aout)

    use globalmod
    use dispmod

    implicit none

    integer            :: i1,iz1,i2,iz2,k
    character(len=120) :: aout
    
    aout=trim(astem)//'_'
    k=len_trim(astem)+2
      
    if (iz1.lt.10) then
       write(aout(k:k+2),'(a1,i1,a1)') '0',iz1,'_'
    else
       write(aout(k:k+2),'(i2,a1)') iz1,'_'
    endif
    k=k+3
    if (abs(i1).lt.10) then
       write(aout(k:k+1),'(a1,i1)') '0',abs(i1)
    else
       write(aout(k:k+1),'(i2)') abs(i1)
    endif
    k=k+2
    if (i1.lt.0) then
       write(aout(k:k+1),'(a2)') 'l_'
    else
       write(aout(k:k+1),'(a2)') 'r_'
    endif
    k=k+2
    
    if (iz2.lt.10) then
       write(aout(k:k+2),'(a1,i1,a1)') '0',iz2,'_'
    else
       write(aout(k:k+2),'(i2,a1)') iz2,'_'
    endif
    k=k+3
    if (abs(i2).lt.10) then
       write(aout(k:k+1),'(a1,i1)') '0',abs(i2)
    else
       write(aout(k:k+1),'(i2)') abs(i2)
    endif
    k=k+2
    if (i2.lt.0) then
       write(aout(k:k+4),'(a5)') 'l.xyz'
    else
       write(aout(k:k+4),'(a5)') 'r.xyz'
    endif
    
    return

  end subroutine filename2d

!#######################################################################

end program displace2
