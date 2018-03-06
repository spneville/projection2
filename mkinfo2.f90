  program mkinfo2

    use globalmod
    use infomod

    implicit none

!-----------------------------------------------------------------------
! Read the input file
!-----------------------------------------------------------------------
    call rdinput

!------------------------------------------------------------------
! Read the data file
!------------------------------------------------------------------
    call rddat

!------------------------------------------------------------------
! Read the symmetry labels
!------------------------------------------------------------------
    ! Temporary
    allocate(asym(ncoo))
    asym='A'

!------------------------------------------------------------------
! Read the set file
!------------------------------------------------------------------
    call rdset

!------------------------------------------------------------------
! Read the QC files
!------------------------------------------------------------------
    call rdqc

!------------------------------------------------------------------
! Write the info file
!------------------------------------------------------------------
    call wrinfo

  contains

!#######################################################################

    subroutine rdinput

      use globalmod
      use infomod
      use parsemod

      implicit none

      integer            :: iin,i
      character(len=250) :: msg

!-----------------------------------------------------------------------
! Set defaults
!-----------------------------------------------------------------------
      adat=''
      aset=''
      acalc=''
      aprog=''
      asymfile=''
      e0=9999d0
      nsta=0
      
!-----------------------------------------------------------------------
! Read the name of the input file
!-----------------------------------------------------------------------
      ain=''
      call getarg(1,ain)
      
      if (ain.eq.'') then
         write(6,'(/,2x,a,/)') 'The name of the input file has not &
              been given'
         STOP
      endif

      if (index(ain,'.inp').eq.0) ain=trim(ain)//'.inp'

!-----------------------------------------------------------------------
! Open the input file
!-----------------------------------------------------------------------
      iin=20
      open(iin,file=ain,form='formatted',status='old')    

!-----------------------------------------------------------------------
! Read the input file
!-----------------------------------------------------------------------
10    continue
      call rdinp(iin)

      i=1
15    continue
      if (keyword(i).ne.'end-input') then

         if (keyword(i).eq.'datafile') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') adat
            else
               write(6,'(/,a,/)') &
                    'No argument given with the datafile keyword'
               STOP
            endif

         else if (keyword(i).eq.'files') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') aset
            else
               write(6,'(/,a,/)') &
                    'No argument given with the files keyword'
               STOP
            endif

         else if (keyword(i).eq.'e0') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) e0
            else
               write(6,'(/,a,/)') &
                    'No argument given with the e0 keyword'
               STOP
            endif

         else if (keyword(i).eq.'calc_type') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') acalc
            else
               write(6,'(/,a,/)') &
                    'No argument given with the calc_type keyword'
               STOP
            endif

         else if (keyword(i).eq.'qc_prog') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') aprog
            else
               write(6,'(/,a,/)') &
                    'No argument given with the qc_prog keyword'
               STOP
            endif
         
         else if (keyword(i).eq.'nstates') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) nsta
            else
               write(6,'(/,a,/)') &
                    'No argument given with the nstates keyword'
               STOP
            endif

         else if (keyword(i).eq.'symmetries') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),'(a)') asymfile
            else
               write(6,'(/,a,/)') &
                    'No argument given with the symmetries keyword'
               STOP
            endif

         else
            write(6,'(/,2(a,2x),/)') 'Unknown keyword:',keyword(i)
            STOP
         endif

         ! Check if more keywords are to be read, else read the next line
         if (i.lt.inkw) then
            i=i+1
            goto 15
         else
            goto 10
         endif
         
      endif

!-----------------------------------------------------------------------
! Close the input file
!-----------------------------------------------------------------------
      close(iin)

!-----------------------------------------------------------------------
! Check that all the required information has been given
!-----------------------------------------------------------------------
      msg=''

      if (adat.eq.'') msg='The name of the data file has not been given'
      
      if (aset.eq.'') msg='The name of the set file has not been given'

      if (acalc.eq.'') msg='The calculation type has not been given'

      if (aprog.eq.'') msg='The QC program name has not been given'

      if (e0.eq.9999d0) msg='The reference energy has not been given'

      if (nsta.eq.0) msg='The number of states has not been given'

      if (msg.ne.'') then
         write(6,'(/,2x,a,/)') trim(msg)
         STOP
      endif

      return

    end subroutine rdinput

!#######################################################################

    subroutine rdset

      use parsemod
      use infomod

      implicit none

      integer :: unit,k

!-----------------------------------------------------------------------
! Open the set file
!-----------------------------------------------------------------------
      unit=30
      open(unit,file=aset,form='formatted',status='old')

!-----------------------------------------------------------------------
! Read the set file
!-----------------------------------------------------------------------
      ! (1) Determine the no. files and allocate arrays
      nqc=0
5     call rdinp(unit)
      if (keyword(1).ne.'end-files') then
         nqc=nqc+1
         goto 5
      endif
      allocate(aqc(nqc))

      ! (2) Read the filenames
      rewind(unit)
      k=0
10    call rdinp(unit)
      if (keyword(1).ne.'end-files') then
         k=k+1
         read(keyword(1),'(a)') aqc(k)
         goto 10
      endif

!-----------------------------------------------------------------------
! Close the set file
!-----------------------------------------------------------------------
      close(unit)

      return

    end subroutine rdset

!#######################################################################

    subroutine rdqc

      use infomod
      use globalmod

      implicit none

      integer                 :: n,unit
      real*8, dimension(ncoo) :: dvec

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(zcoo(nqc,ncoo))
      allocate(xcoo(nqc,ncoo))
      allocate(energy(nqc,ncoo))
      xcoo=0.0d0
      energy=0.0d0

!-----------------------------------------------------------------------
! Read the coordinates and energies for each file
!-----------------------------------------------------------------------
      unit=40

      do n=1,nqc

         ! Open file
         open(unit,file=aqc(n),form='formatted',status='old')

         ! Read the Cartesian coordinates
         call rdcoord(n,unit)

         ! Calculate the coordinates in terms of the projected
         ! normal modes and vectors of choice
         dvec=xcoo(n,:)-xcoo0
         zcoo(n,:)=matmul(cartz,dvec)

         ! Read the state energies
         call rdenergy(n,unit)

         ! Close file
         close(unit)

      enddo

      return

    end subroutine rdqc

!#######################################################################

    subroutine rdcoord(n,unit)

      use infomod

      implicit none

      integer :: n,unit

!------------------------------------------------------------------
! Read the coordinates
!------------------------------------------------------------------
      ! Turbomole
      if (aprog.eq.'dftmrci') then
         call rdcoord_dftmrci(n,unit)
      else
         write(6,'(/,2(2x,a),/)') "File type not supported:",&
              trim(aprog)
         STOP
      endif

      return

    end subroutine rdcoord

!#######################################################################

    subroutine rdcoord_dftmrci(n,unit)

      use globalmod
      use infomod
      use parsemod
      
      implicit none

      integer            :: n,unit,i,j
      character(len=120) :: string

!------------------------------------------------------------------
! Read to the coordinates
!------------------------------------------------------------------
5     read(unit,'(a)') string
      if (index(string,'$coord').eq.0) goto 5

!------------------------------------------------------------------
! Read the Cartesian coordinates
!------------------------------------------------------------------
      do i=1,natm
         call rdinp(unit)
         do j=1,3
            read(keyword(j),*) xcoo(n,i*3-3+j)
         enddo
      enddo

      ! Convert to Angstrom
      xcoo(n,:)=xcoo(n,:)*0.529177249d0

      return

    end subroutine rdcoord_dftmrci

!#######################################################################

    subroutine rdenergy(n,unit)

      use infomod

      implicit none

      integer :: n,unit

      if (aprog.eq.'dftmrci') then
         
         if (acalc.eq.'dftmrci') then
            call rdenergy_dftmrci(n,unit)
         else
            write(6,'(/,2(2x,a),/)') "Calculation type not supported:",&
                 trim(acalc)
         endif

      else
         write(6,'(/,2(2x,a),/)') "File type not supported:",&
              trim(aprog)
         STOP
      endif

      return

    end subroutine rdenergy

!#######################################################################

    subroutine rdenergy_dftmrci(n,unit)

      use infomod

      implicit none

      integer            :: n,unit,nfound,itype
      character(len=120) :: string

      rewind(unit)
               
      nfound=0
5     read(unit,'(a)',end=10) string
      if (index(string,'DFTCI').ne.0) then
         nfound=nfound+1
         if (nfound.gt.nsta) goto 10
         read(string,'(31x,F12.6)') energy(n,nfound)
      endif
      goto 5
         
10    continue
      
      ! Take the state energies relative to the ground state and 
      ! convert to eV
      energy(n,:)=(energy(n,:)-e0)*27.2113845d0
      
      ! Exit if we have not found all state energier
      if (nfound.lt.nsta) then
         write(6,'(/,2x,a,1x,i2,2(1x,a),/)') 'Less than',nsta,&
              'states found in the file',trim(aqc(n))
         STOP
      endif

      return

    end subroutine rdenergy_dftmrci

!#######################################################################

    subroutine wrinfo

      use globalmod
      use infomod

      implicit none

      integer                 :: unit,i,j,k,m,nzero,nblk,j1,j2,s,&
                                 n
      real*8, dimension(ncoo) :: freq
      character(len=120)      :: ainfo
      character(len=8)        :: label

!------------------------------------------------------------------
! Open the .info file
!------------------------------------------------------------------
      ainfo=''
      k=index(ain,'.')
      write(ainfo(1:k),'(a)') ain(1:k)
      write(ainfo(k+1:k+4),'(a)') 'info'

      unit=50
      open(unit,file=ainfo,form='formatted',status='unknown')

!------------------------------------------------------------------
! Version number (dummy number)
!------------------------------------------------------------------
      write(unit,'(a,f8.5)') '# VCHam Version :',10.10050
      write(unit,'(a)') '# All energies and frequencies in eV'
      write(unit,'(a)') ' '  

!------------------------------------------------------------------
! System information
!------------------------------------------------------------------
      ! Temporary: 
      nzero=0

      write(unit,'(a)') '#system_dimensions'
      write(unit,'(a)') '#internal_coordinates'
      write(unit,'(i3)') ncoo
      write(unit,'(a)')'#trivial_modes'
      write(unit,'(i3)') nzero
      write(unit,'(a)') '#num_states'
      write(unit,'(i3)') nsta
      write(unit,'(a)') '#states'
      write(unit,'(20i3)') (i, i=1,nsta)
      write(unit,'(a)') '#end_system_dimensions'

      write(unit,*)
      write(unit,'(a)') &
           '#general_information from the ground state'
      write(unit,'(a)') '#zero_of_energy'
      write(unit,'(F20.10)') e0
      write(unit,'(a)') '#equilibrium_state_geometry'
      do i=1,ncoo,3
         write(unit,'(3F15.10)') (xcoo0(j),j=i,i+2)
      enddo
      write(unit,'(a)') '#masses'
      do i=1,ncoo,3
         write(unit,'(F15.10)') mass(i)
      enddo

!------------------------------------------------------------------
! Normal mode information
!------------------------------------------------------------------
      do i=1,ncoo
         freq(i)=sqrt(abs(sthess(i,i)))*0.6373641d0
      enddo

      write(unit,'(a)') '#vibrational_frequencies'
      do m=1,ncoo-nzero
         write(unit,'(i5,2x,a,F20.10)') m,asym(m),&
              freq(m)*0.6373641d0
      enddo

      write(unit,'(a)') '#normal_modes'
      nblk=((ncoo-1)/5)+1
      do m=1,ncoo-nzero
         do i=1,nblk
            j1=(i-1)*5+1
            j2=min(j1+4,ncoo)
            write(unit,'(5F15.10)') (zcart(j,m),j=j1,j2)
         enddo
         write(unit,*) 
      enddo

      write(unit,'(a)') '#end_general_information'
      write(unit,*) 

!------------------------------------------------------------------
! Coordinates and energies
!------------------------------------------------------------------
      write(unit,'(a)') '#dataset_section'

      nblk=((ncoo-nzero-1)/5)+1
      k=0
      do n=1,nqc
         do s=1,nsta

            k=k+1
            
            call getlabel(k,label)
            write(unit,'(a)') '#dataset_'//label//'    '
            write(unit,'(a)') '#energy'
            write(unit,'(f20.10)') energy(n,s)

            write(unit,'(a)') '#state'
            write(unit,'(i3)') s

            write(unit,'(a)') '#point_in_q'
            do i=1,nblk
               j1=(i-1)*5+1
               j2=min(j1+4,ncoo-nzero)
               write(unit,'(5F15.10)') (zcoo(n,j),j=j1,j2) 
            enddo

            write(unit,'(a)') '#point_in_cartesian_coordinates'
            do i=1,ncoo,3
               write(unit,'(3F16.10)') (xcoo(n,j),j=i,i+2)
            enddo
            
            write(unit,'(a)') '#end_dataset_'//label//'    '

         enddo
      enddo

      return

    end subroutine wrinfo

!##################################################################

    subroutine getlabel(n,string)
        
      implicit none
      
      integer, intent(in)           :: n
      character(len=8), intent(out) :: string
      
      string='     '
      if (n.le.9) then
         write(string(1:1),'(i1)') n
      else if (n.le.99) then
         write(string(1:2),'(i2)') n
      else if (n.le.999) then
         write(string(1:3),'(i3)') n
      else if (n.le.9999) then
         write(string(1:4),'(i4)') n
      else if (n.le.99999) then
         write(string(1:5),'(i5)') n
      else if (n.le.999999) then
         write(string(1:6),'(i6)') n
      else if (n.le.9999999) then
         write(string(1:7),'(i7)') n
      else if (n.le.99999999) then
         write(string(1:8),'(i8)') n
      endif
        
      return

    end subroutine getlabel

!#######################################################################

  end program mkinfo2
