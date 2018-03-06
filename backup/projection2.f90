  program projection2
    
    implicit none

!-----------------------------------------------------------------------
! Read the input file
!-----------------------------------------------------------------------
    call rdinput

!-----------------------------------------------------------------------
! Read the Cartesian coordinates of the reference geometry
!-----------------------------------------------------------------------
    call rdrefgeom

!-----------------------------------------------------------------------
! Read the vectors
!-----------------------------------------------------------------------
    call rdvecs

!-----------------------------------------------------------------------
! Read and mass-weight the Hessian
!-----------------------------------------------------------------------
    call rdhess

!-----------------------------------------------------------------------
! Calculate the eigenvectors of the projected Hessian
!-----------------------------------------------------------------------   
    call projhess

!-----------------------------------------------------------------------
! Calculate and output the Hessian in terms the nuclear coordinates
! of interest
!-----------------------------------------------------------------------
    call trans_hess

!-----------------------------------------------------------------------
! Output the system and coordinate information
!-----------------------------------------------------------------------
    call wrout

    STOP

  contains

!#######################################################################

    subroutine rdinput

      use globalmod
      use projmod
      use parsemod

      implicit none

      integer                           :: iin,k,i
      character(len=120)                :: errmsg
      character(len=60), dimension(100) :: atmp

!-----------------------------------------------------------------------
! Set defaults
!-----------------------------------------------------------------------
      ain=''
      axyz=''
      ahess=''
      nvec=0

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
5     continue
      call rdinp(iin)

      i=0
      if (keyword(1).ne.'end-input') then
10       continue
         i=i+1

         if (keyword(i).eq.'geom_file') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               axyz=keyword(i)
            else
               goto 100
            endif

         else if (keyword(i).eq.'hessian_file') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               ahess=keyword(i)
            else
               goto 100
            endif

         else if (keyword(i).eq.'vec_files') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               nvec=1
               atmp(1)=keyword(i)
15             continue
               if (keyword(i+1).eq.',') then
                  i=i+2
                  nvec=nvec+1
                  atmp(nvec)=keyword(i)
                  goto 15
               endif               
               allocate(avec(nvec))
               avec(1:nvec)=atmp(1:nvec)
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
100      continue
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
        if (axyz.eq.'') then
           errmsg='The xyz file name has not been given'
           write(6,'(/,2x,a,/)') errmsg
        endif

        if (ahess.eq.'') then
           errmsg='The Hessian file name has not been given'
           write(6,'(/,2x,a,/)') errmsg
        endif
        
        if (nvec.eq.0) then
           errmsg='No vectors have been given'
           write(6,'(/,2x,a,/)') errmsg
        endif
        
      return

    end subroutine rdinput

!#######################################################################

    subroutine rdrefgeom

      use globalmod
      use projmod
      use parsemod

      implicit none

      integer :: unit,i,j

!-----------------------------------------------------------------------
! Open the xyz file
!-----------------------------------------------------------------------
      unit=20
      open(unit,file=axyz,form='formatted',status='old')

!-----------------------------------------------------------------------
! Read the no. atoms and allocate arrays
!-----------------------------------------------------------------------
      call rdinp(unit)

      read(keyword(1),*) natm

      ncoo=3*natm

      allocate(xcoo0(ncoo))
      allocate(mass(ncoo))
      allocate(atlbl(natm))

!-----------------------------------------------------------------------
! Read the xyz file
!-----------------------------------------------------------------------
      ! Skip past the comment line
      read(unit,*)

      ! Read the atom labels and Cartesian coordinates
      do i=1,natm

         ! Read the keywords on the current line
         call rdinp(unit)

         ! Atom labels
         atlbl(i)=keyword(1)
         atlbl(i)=to_upper(atlbl(i))

         ! Cartesian coordinates
         do j=1,3
            read(keyword(1+j),*) xcoo0(i*3-3+j)
         enddo

      enddo

!-----------------------------------------------------------------------
! Assign masses
!-----------------------------------------------------------------------
      do i=1,natm
         if (atlbl(i).eq.'N') then
            mass(i*3-2:i*3)=14.0067d0
         else if (atlbl(i).eq.'C') then 
            mass(i*3-2:i*3)=12.0107d0
         else if (atlbl(i).eq.'H') then    
            mass(i*3-2:i*3)=1.00794d0
         else
            write(6,'(/,2(a,x),/)') 'Unknown atom type:',atlbl(i)
            STOP
         endif
      enddo

!-----------------------------------------------------------------------
! Close the xyz file
!-----------------------------------------------------------------------
      close(unit)

      return

    end subroutine rdrefgeom

!#######################################################################

    function to_upper(strIn) result(strOut)
      
      implicit none
      
      character(len=*), intent(in) :: strIn
      character(len=len(strIn)) :: strOut
      integer :: i,j
      
      do i = 1, len(strIn)
         j = iachar(strIn(i:i))
         if (j>= iachar("a") .and. j<=iachar("z") ) then
            strOut(i:i) = achar(iachar(strIn(i:i))-32)
         else
            strOut(i:i) = strIn(i:i)
         end if
     end do
     
     return

   end function to_upper

!#######################################################################

   subroutine rdvecs

     use globalmod
     use projmod
     use parsemod

     implicit none

     integer :: unit,i,j,k
     real*8  :: length

     allocate(vec(nvec,ncoo))

!-----------------------------------------------------------------------
! Read the Cartesian vectors
!-----------------------------------------------------------------------
     unit=20
     do i=1,nvec
        open(unit,file=avec(i),form='formatted',status='old')
        do j=1,natm
           call rdinp(unit)
           do k=1,3
              read(keyword(k),*) vec(i,j*3-3+k)
           enddo
        enddo
        close(unit)
     enddo

!-----------------------------------------------------------------------
! Mass weight and normalise the vectors
!-----------------------------------------------------------------------
     do i=1,nvec
        length=0.0d0
        do j=1,ncoo
           vec(i,j)=vec(i,j)*sqrt(mass(j))
           length=length+vec(i,j)**2
        enddo
        length=sqrt(length)
        vec(i,:)=vec(i,:)/length
     enddo

     return

   end subroutine rdvecs

!#######################################################################

   subroutine rdhess

     use globalmod
     use projmod

     implicit none

     integer :: unit,itype,i,j

     allocate(hess(ncoo,ncoo))

!-----------------------------------------------------------------------
! Open the Hessian file
!-----------------------------------------------------------------------     
     unit=20
     open(unit,file=ahess,form='formatted',status='old')

!-----------------------------------------------------------------------
! Determine the file type!
!-----------------------------------------------------------------------
     call getfiletype(unit,itype)

!-----------------------------------------------------------------------
! Read the Hessian
!-----------------------------------------------------------------------
     rewind(unit)
     if (itype.eq.1) then
        call rdaoforce(unit)
     else if (itype.eq.2) then
        call rdnumforce(unit)
     endif

!-----------------------------------------------------------------------
! Close the Hessian file
!-----------------------------------------------------------------------
     close(unit)

!-----------------------------------------------------------------------
! Mass-weight the Hessian
!-----------------------------------------------------------------------
     do i=1,ncoo
        do j=1,ncoo
           hess(i,j)=hess(i,j)/sqrt(mass(i)*mass(j))
        enddo
     enddo

     return

   end subroutine rdhess

!#######################################################################
!
! itype: 1 <--> Turbomole, aoforce
!        2 <--> Turbomole, numforce
!
!##################################################################

    subroutine getfiletype(unit,itype)

      implicit none
      
      integer           :: unit,itype
      character(len=90) :: string
      
      itype=0

10    continue
      read(unit,'(a)',end=999) string
      
      if (string(1:11).eq.'$nprhessian') then
         itype=1
      else if (string(1:8).eq.'$hessian'.and.&
           string(11:19).ne.'projected') then
         itype=2
      endif

      goto 10

999   continue

      if (itype.eq.0) then
         write(6,'(/,a,/)') 'Hessian file type not recognised'
         STOP
      endif

      return

    end subroutine getfiletype

!#######################################################################

    subroutine rdaoforce(unit)

      use projmod

      implicit none

      integer              :: unit,curr,i,j,k,itmp
      real*8, dimension(5) :: ftmp
      character(len=90)    :: string

!-----------------------------------------------------------------------
! Read to the non-projected Hessian section
!-----------------------------------------------------------------------
10    continue
      read(unit,'(a)') string
      if (string(1:11).ne.'$nprhessian') goto 10

!-----------------------------------------------------------------------
! Read the non-projected Hessian
!-----------------------------------------------------------------------
      curr=1
20    continue
      read(unit,'(a)') string
      if (string(1:15).ne.'$nprvibrational') then
         
         read(string,'(i2,1x,i2,2x,5(F13.10,2x))') k,itmp,&
              (ftmp(i), i=1,5)
         
         if (k.gt.curr) curr=curr+1
         
         j=0
         do i=itmp*5-4,itmp*5
            j=j+1
            hess(i,curr)=ftmp(j)
         enddo
         
         goto 20
         
      endif

      return

    end subroutine rdaoforce

!#######################################################################

    subroutine rdnumforce(unit)

      use globalmod
      use projmod

      implicit none
      
      integer              :: unit,curr,i,j,k,l,itmp,nline,count
      real*8               :: f
      real*8, dimension(5) :: ftmp
      character(len=90)    :: string

!-----------------------------------------------------------------------
! Read to the non-projected Hessian section
!-----------------------------------------------------------------------
      read(unit,*)

!-----------------------------------------------------------------------
! Calculate no. lines per column
!-----------------------------------------------------------------------
      f=ncoo/5.0d0
      nline=ceiling(f)

!-----------------------------------------------------------------------
! Read the non-projected Hessian
!-----------------------------------------------------------------------
      do i=1,ncoo
         do j=1,nline
            
            read(unit,'(7x,5(F13.10,2x))') (ftmp(k), k=1,5)

            count=0
            do l=j*5-4,j*5
               count=count+1
               hess(l,i)=ftmp(count)
            enddo

         enddo
      enddo
      
      return

    end subroutine rdnumforce

!#######################################################################

    subroutine projhess

      use globalmod
      use projmod

      implicit none

      integer                      :: i,j,k,l,e2,error,unit
      real*8, dimension(ncoo,ncoo) :: qproj,phess,tmp2
      real*8, dimension(ncoo)      :: tmp1,lambda
      real*8, dimension(3*ncoo)    :: work
      real*8                       :: length,sum

!-----------------------------------------------------------------------
! Form the projector onto the subspace orthogonal to that spanned by
! the vectors
!-----------------------------------------------------------------------
      qproj=0.0d0
      do i=1,ncoo
         qproj(i,i)=1.0d0
      enddo

      do i=1,nvec
         do j=1,ncoo
            do k=1,ncoo
               qproj(j,k)=qproj(j,k)-vec(i,j)*vec(i,k)
            enddo
         enddo
      enddo

!-----------------------------------------------------------------------
! Project the Hessian onto the subspace spanned by the orthogonal
! complement of the vectors
!-----------------------------------------------------------------------
      phess=0.0d0

      do i=1,ncoo
         do j=1,ncoo
            do k=1,ncoo
               do l=1,ncoo
                  phess(i,j)=phess(i,j)+qproj(i,k)*hess(k,l)*qproj(l,j)
               enddo
            enddo
         enddo
      enddo

!-----------------------------------------------------------------------
! Diagonalise the projected Hessian
!-----------------------------------------------------------------------
      allocate(peigvec(ncoo,ncoo))

      peigvec=phess
      e2=3*ncoo

      call dsyev('V','U',ncoo,peigvec,ncoo,lambda,work,e2,error)

      if (error.ne.0) then
         write(6,'(/,a,2/,a,x,i20,/)') &
              ' Diagonalisation of the projected Hessian failed.'&
              ,' dsyev error:',error
         STOP
      endif

!-----------------------------------------------------------------------
! Determine the eigenvectors of interest, and save these as the first
! 3N-nvec columns of peigvec
!-----------------------------------------------------------------------
      tmp1=0.0d0
      tmp2=0.0d0
      
      k=0
      do i=1,ncoo
         sum=0.0d0
         do j=1,nvec
            sum=sum+abs(dot_product(peigvec(:,i),vec(j,:)))           
         enddo
         if (sum.lt.1d-6) then
            k=k+1
            tmp2(:,k)=peigvec(:,i)
         endif
      enddo

      peigvec=tmp2

!-----------------------------------------------------------------------
! Save the input vectors as the columns 3N-nvec+1,...,3N of peigvec
!-----------------------------------------------------------------------
      k=0
      do i=ncoo-nvec+1,ncoo
         k=k+1
         peigvec(:,i)=vec(k,:)
      enddo

      return

    end subroutine projhess

!#######################################################################

    subroutine trans_hess

      use globalmod
      use projmod

      implicit none

      integer                      :: i,j,k,imax,jmax,unit
      real*8, dimension(ncoo,ncoo) :: tmpmat
      real*8, dimension(ncoo)      :: tmp1,diag
      real*8                       :: max,ftmp,length,sc
      real*8, parameter            :: eh2ev=27.211396d0

!-----------------------------------------------------------------------
! Perform the similarity transformation of the Hessian in terms of
! mass-weighted Cartesians using the transformation matrix peigvec
!-----------------------------------------------------------------------
      allocate(sthess(ncoo,ncoo))

      tmpmat=matmul(hess,peigvec)
      sthess=matmul(transpose(peigvec),tmpmat)

      max=0.0d0
      do i=1,ncoo
         do j=1,ncoo
            if (i.ne.j.and.sthess(i,j)*eh2ev.gt.max) then
               max=sthess(i,j)*eh2ev
               imax=i
               jmax=j
            endif
         enddo
      enddo

!------------------------------------------------------------------
! Output the elements of the similarity transformed Hessian
!------------------------------------------------------------------
      unit=20
      open(unit,file='sthess',form='formatted',status='unknown')

      write(unit,'(a,/,a,2x,2(i3,2x),/)') &
           'All energies in cm-1',&
           'Largest off-diagonal element:',imax,jmax

      do i=1,ncoo
         do j=i,ncoo
            write(unit,'(2x,i3,2x,i3,3x,F13.7)') i,j,&
                 sqrt(abs(sthess(i,j)))*5140.66
         enddo
      enddo

      close(unit)

!------------------------------------------------------------------
! Output the the nuclear coordinate vectors to file
!------------------------------------------------------------------
      unit=20
      open(unit,file='peigvec.xyz',form='formatted',status='unknown')
      
      do i=1,ncoo
         
         ! Convert the ith eigenvector to non-mass-weighted Cartesians
         ! and normalise
         do j=1,ncoo
            tmp1(j)=peigvec(j,i)/sqrt(mass(j))
         enddo
         length=sqrt(dot_product(tmp1,tmp1))
         tmp1=tmp1/length

         write(unit,'(i3)') natm

         ftmp=sqrt(abs(sthess(i,i)))*5140.66
         if (sthess(i,i).lt.0.0d0) then
            write(unit,'(F8.2,a1,x,a)') ftmp,'i','cm-1'
         else
            write(unit,'(F8.2,x,a)') ftmp,'cm-1'
         endif

         do j=1,natm
            write(unit,'(a2,6(2x,F10.7))') atlbl(j),&
                 (xcoo0(k),k=j*3-2,j*3),&
                 (tmp1(k),k=j*3-2,j*3)
         enddo

      enddo

      close(unit)

      return

    end subroutine trans_hess

!#######################################################################

    subroutine wrout

      use globalmod
      use projmod

      implicit none

      integer            :: unit,k
      character(len=120) :: aout

!-----------------------------------------------------------------------
! Open the output file
!-----------------------------------------------------------------------
      k=len_trim(ain)-3
      aout=ain(1:k)//'dat'

      unit=20
      open(unit,file=aout,form='unformatted',status='unknown')

!-----------------------------------------------------------------------
! Output the system information
!-----------------------------------------------------------------------
      ! No. Cartesian dofs
      write(unit) ncoo

      ! Reference geometry
      write(unit) xcoo0

      ! Masses
      write(unit) mass

      ! Atom labels
      write(unit) atlbl

      ! Coordinate vectors
      write(unit) peigvec

      ! Similarity transformed Hessian
      write(unit) sthess

!-----------------------------------------------------------------------
! Close the output file
!-----------------------------------------------------------------------
      close(unit)

      return

    end subroutine wrout

!#######################################################################

  end program projection2
