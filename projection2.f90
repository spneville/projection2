  program projection2
    
    use projmod

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
    if (lproj) call rdvecs

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
      lproj=.true.
      lprojrt=.false.

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

         else if (keyword(i).eq.'noproj') then
            lproj=.false.

         else if (keyword(i).eq.'projrt') then
            lprojrt=.true.

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
        
        if (nvec.eq.0.and.lproj) then
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

     integer                      :: unit,i,j
     real*8                       :: length,fac,dp
     real*8, dimension(nvec,ncoo) :: vec2
     logical, dimension(nvec)     :: l2nd
     logical                      :: lortho

     allocate(vec(nvec,ncoo))

!-----------------------------------------------------------------------
! Read the Cartesian vectors
!-----------------------------------------------------------------------
     unit=20
     l2nd=.false.
     do i=1,nvec
        open(unit,file=avec(i),form='formatted',status='old')

        ! Read the first vector
        call read1vec(vec(i,:),unit)

        ! Check to see if there is another vector in the current file
        ! to be taken as a linear combination with the once just read
        call rdinp(unit)
        if (.not.lend) then
           l2nd(i)=.true.
           ! Determine whether the linear combination is symmetric or
           ! antisymmetric
           if (keyword(1).eq.'plus') then
              fac=1.0d0
           else if (keyword(1).eq.'minus') then
              fac=-1.0d0
           endif
           ! Read the next vector
           call read1vec(vec2(i,:),unit)
        endif

        close(unit)
     enddo

!-----------------------------------------------------------------------
! Mass weight and normalise the vectors
!-----------------------------------------------------------------------
     do i=1,nvec

        if (l2nd(i)) then
           length=dsqrt(dot_product(vec(i,:),vec(i,:)))
           vec(i,:)=vec(i,:)/length
           length=dsqrt(dot_product(vec2(i,:),vec2(i,:)))
           vec2(i,:)=vec2(i,:)/length
           vec(i,:)=(1.0d0/dsqrt(2.0d0))*(vec(i,:)+fac*vec2(i,:))
        endif

        length=0.0d0
        do j=1,ncoo
           vec(i,j)=vec(i,j)*sqrt(mass(j))
           length=length+vec(i,j)**2
        enddo
        length=dsqrt(length)
        vec(i,:)=vec(i,:)/length
        
     enddo

!-----------------------------------------------------------------------
! Check orthogonality
!-----------------------------------------------------------------------
     lortho=.true.
     do i=1,nvec
        do j=i+1,nvec
           dp=dot_product(vec(i,:),vec(j,:))           
           if (abs(dp).gt.1d-5) then
              write(6,'(/,x,2(x,a,x,i2),x,a,/)') 'Vectors',i,'and',j,&
                   'are not orthogonal'
              lortho=.false.
           endif
        enddo
     enddo
     
     if (.not.lortho) STOP

     return

   end subroutine rdvecs

!#######################################################################

   subroutine read1vec(vector,unit)

     use globalmod
     use projmod
     use parsemod

     implicit none

     integer                 :: unit,i,k
     real*8, dimension(ncoo) :: vector,v1,v2
     logical                 :: lvec

     ! Determine whether a vector has been given or else two
     ! geometries from which the vector will be calculated
     call rdinp(unit)
     backspace(unit)
     if (inkw.eq.3) then
        lvec=.true.
     else
        lvec=.false.
     endif

     ! Vector input
     if (lvec) then
        do i=1,natm
           call rdinp(unit)
           do k=1,3
              read(keyword(k),*) vector(i*3-3+k)
           enddo
        enddo
     ! Geometry input
     else
        do i=1,natm
           call rdinp(unit)
           do k=1,3
              read(keyword(k),*) v1(i*3-3+k)
           enddo
           do k=4,6
              read(keyword(k),*) v2(i*3-3+k-3)
           enddo
        enddo
        vector=v2-v1
     endif

     return

   end subroutine read1vec

!#######################################################################

   subroutine rdhess

     use globalmod
     use projmod

     implicit none

     integer                 :: unit,itype,i,j
     real*8, dimension(ncoo) :: xcoo_hess

!-----------------------------------------------------------------------
! Allocate the array that will hold the Hessian
!-----------------------------------------------------------------------
     allocate(hess(ncoo,ncoo))

!-----------------------------------------------------------------------
! Open the Hessian file
!-----------------------------------------------------------------------     
     unit=20
     open(unit,file=ahess,form='formatted',status='old')

!-----------------------------------------------------------------------
! Determine the file type
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
     else if (itype.eq.3) then
        call rdfmstype(unit)
     endif

!-----------------------------------------------------------------------
! Read the geometry used in the calculation of the Hessian
!-----------------------------------------------------------------------
     rewind(unit)
     if (itype.eq.1) then
        call read_hess_geom_aoforce(unit,xcoo_hess)
     else if (itype.eq.2) then
        call read_hess_geom_numforce(unit,xcoo_hess)
     else if (itype.eq.3) then
        call read_hess_geom_fmstype(unit,xcoo_hess)
     else
        write(6,'(/,2x,a,/)') 'You need to write the code to read &
             the geometry used in the calculation of the Hessian...'
        STOP
     endif

!-----------------------------------------------------------------------
! Rotate the Hessian to the same frame as the reference geometry
!-----------------------------------------------------------------------
     call rotate_hess(xcoo_hess)

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
!        3 <--> FMS-type file
!
!##################################################################

    subroutine getfiletype(unit,itype)
      
      use parsemod
      
      implicit none
      
      integer           :: unit,itype
      character(len=90) :: string
      logical           :: lnum
      
      itype=0

!-----------------------------------------------------------------------
! (1) Check to see if the file if od the FMS type...
!-----------------------------------------------------------------------
      call rdinp(unit)
      
      if (inkw.eq.1) then
         lnum=is_numeric(keyword(1))
         if (lnum) itype=3
      endif

      if (itype.ne.0) goto 999

!-----------------------------------------------------------------------
! (2) ...if not, then determine which QC program wrote the file
!-----------------------------------------------------------------------      
      rewind(unit)
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

    function is_numeric(string)
      
      implicit none

      character(len=*), intent(in) :: string
      logical                      :: is_numeric
      real                         :: x
      integer                      :: e

      read(string,*,iostat=e) x

      is_numeric = e == 0

      return

    end function is_numeric
    
!#######################################################################

    subroutine rdaoforce(unit)

      use globalmod
      use projmod

      implicit none

      integer              :: unit,curr,i,j,k,itmp,upper
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

         upper=min(itmp*5,ncoo)         
         do i=itmp*5-4,upper
            j=j+1
            hess(i,curr)=ftmp(j)
         enddo
         
         goto 20
         
      endif

      return

    end subroutine rdaoforce

!#######################################################################

    subroutine read_hess_geom_aoforce(unit,xcoo_hess)

      use globalmod

      implicit none
      
      integer                 :: unit,igeom,i,j,k,n
      real*8, dimension(ncoo) :: xcoo_hess
      character(len=120)      :: string,filename

!-----------------------------------------------------------------------
! Read the name of the coord file
!-----------------------------------------------------------------------
10    read(unit,'(a)') string
      if (index(string,'$coord').eq.0) goto 10

      k=index(string,'file=')+5

      filename=''
      n=0
      do i=k,len_trim(string)
         if (string(i:i).eq.' ') exit
         n=n+1
         filename(n:n)=string(i:i)
      enddo

!-----------------------------------------------------------------------
! Read the Cartesian coordinates from the coord file
!-----------------------------------------------------------------------      
      ! Open the coord file
      igeom=33
      open(igeom,file=filename,form='formatted',status='old')

      ! Skip to the Cartesian coordinates
      read(igeom,*)
      
      ! Read the Cartesian coordinates
      do i=1,natm
         read(igeom,'(3(F20.14,2x))') (xcoo_hess(j), j=i*3-2,i*3)
      enddo

      ! Convert to Angstrom
      xcoo_hess=xcoo_hess*0.529177249d0

      ! Close the coord file
      close(igeom)

      return

    end subroutine read_hess_geom_aoforce

!#######################################################################

    subroutine read_hess_geom_numforce(unit,xcoo_hess)

      use globalmod

      implicit none
      
      integer                 :: unit,igeom,i,j
      real*8, dimension(ncoo) :: xcoo_hess
      character(len=120)      :: string,filename
      logical                 :: lexists

!-----------------------------------------------------------------------
! Make sure that the aoforce.out file is present
!-----------------------------------------------------------------------
      inquire(file='aoforce.out',exist=lexists)

      if (.not.lexists) then
         write(6,'(/,2x,a,/)') 'The aoforce.out file is missing...'
         STOP
      endif

!-----------------------------------------------------------------------
! Read the Cartesian coordinates from aoforce.out
!-----------------------------------------------------------------------
      igeom=33
      open(igeom,file='aoforce.out',form='formatted',status='old')

      ! Read to the geometry section
5     read(igeom,'(a)') string
      if (index(string,'actual cartesian coordinates').eq.0) goto 5
      do i=1,4
         read(igeom,*)
      enddo

      ! Read the Cartesian coordinates
      do i=1,natm
         read(igeom,'(8x,3(4x,F18.14))') (xcoo_hess(j), j=i*3-2,i*3)
      enddo

      ! Convert to Angstrom
      xcoo_hess=xcoo_hess*0.529177249d0

      close(igeom)

      return

    end subroutine read_hess_geom_numforce

!#######################################################################

    subroutine read_hess_geom_fmstype(unit,xcoo_hess)

      use parsemod
      use globalmod

      implicit none
      
      integer                 :: unit,igeom,i,j
      real*8, dimension(ncoo) :: xcoo_hess
      character(len=120)      :: string,filename
      logical                 :: lau

!-----------------------------------------------------------------------
! Open the Geometry.dat file
!-----------------------------------------------------------------------
      igeom=33
      open(igeom,file='Geometry.dat',form='formatted',status='old')

!-----------------------------------------------------------------------
! Determine the units used
!-----------------------------------------------------------------------
      call rdinp(igeom)

      if (keyword(3).eq.'bohr') then
         lau=.true.
      else if (keyword(3).eq.'angstrom') then
         lau=.false.
      endif

!-----------------------------------------------------------------------
! Read the Cartesian coordinates
!-----------------------------------------------------------------------
      read(unit,*)

      do i=1,natm
         call rdinp(unit)
         do j=1,3
            read(keyword(j+1),*) xcoo_hess(i*3-3+j)
         enddo
      enddo

      if (lau) xcoo_hess=xcoo_hess*0.529177249d0

!-----------------------------------------------------------------------
! Close the Geometry.dat file
!-----------------------------------------------------------------------
      close(igeom)

      return

    end subroutine read_hess_geom_fmstype

!#######################################################################

    subroutine rdnumforce(unit)

      use globalmod
      use projmod

      implicit none
      
      integer              :: unit,i,j,k,l,nline,count
      real*8               :: f
      real*8, dimension(5) :: ftmp

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

    subroutine rdfmstype(unit)

      use globalmod
      use projmod

      implicit none

      integer :: unit,i,j

      read(unit,*)

      do i=1,ncoo
         read(unit,*) (hess(i,j), j=1,ncoo)
      enddo

      return

    end subroutine rdfmstype

!#######################################################################

    subroutine rotate_hess(xcoo_hess)

      use globalmod
      use projmod

      implicit none

      integer                      :: i,j,k,l
      real*8, dimension(ncoo)      :: xcoo_hess
      real*8, dimension(ncoo,ncoo) :: rotmat,tmp

!-----------------------------------------------------------------------
! (1) Determine the rotation matrix
!-----------------------------------------------------------------------
      call get_rot_mat(xcoo_hess,rotmat)

!-----------------------------------------------------------------------
! (2) Similarity transform the Hessian
!-----------------------------------------------------------------------
      tmp=0.0d0
      do i=1,ncoo
         do j=1,ncoo
            do k=1,ncoo
               do l=1,ncoo
                  tmp(i,j)=tmp(i,j)+rotmat(k,i)*hess(k,l)*rotmat(l,j)
               enddo
            enddo
         enddo
      enddo

      hess=tmp

      return

    end subroutine rotate_hess

!#######################################################################

    subroutine get_rot_mat(xcoo_hess,rotmat)

      use globalmod

      implicit none

      integer                      :: i,j,k,m,n,ilbl,jlbl
      real*8, dimension(ncoo)      :: xcoo_hess,x1,x2
      real*8, dimension(ncoo,ncoo) :: rotmat
      real*8                       :: totmass,det,d
      real*8, dimension(3)         :: com1,com2,pm
      real*8, dimension(natm,3)    :: x1mat,x2mat
      real*8, dimension(3,3)       :: vut,mat

      integer                      :: lwork,info,dim
      real*8, dimension(3)         :: sigma
      real*8, dimension(3,3)       :: covar,u,vt,orig
      real*8, dimension(15)        :: work

      x1=xcoo0
      x2=xcoo_hess

!-----------------------------------------------------------------------
! Translate the origins of both geometries to the centre of mass
!-----------------------------------------------------------------------
      totmass=0.0d0
      com1=0.0d0
      com2=0.0d0
      do i=1,natm
         totmass=totmass+mass(i*3)
         do j=1,3
            com1(j)=com1(j)+x1(i*3-3+j)*mass(i*3)
            com2(j)=com2(j)+x2(i*3-3+j)*mass(i*3)
         enddo
      enddo

      com1=com1/totmass
      com2=com2/totmass

      do i=1,natm
         do j=1,3
            x1(i*3-3+j)=x1(i*3-3+j)-com1(j)
            x2(i*3-3+j)=x2(i*3-3+j)-com2(j)
         enddo
      enddo
      
      do i=1,natm
         do j=1,3
            x1mat(i,j)=x1(i*3-3+j)
            x2mat(i,j)=x2(i*3-3+j)
         enddo
      enddo

!-----------------------------------------------------------------------
! Calculate the 3X3 covariance matrix
!-----------------------------------------------------------------------
      covar=0.0d0
      do i=1,3
         do j=1,3
            do k=1,natm
               covar(i,j)=covar(i,j)+x1(k*3-3+i)*x2(k*3-3+j)
            enddo
         enddo
      enddo

!-----------------------------------------------------------------------
! Calculate the SVD of the covariance matrix
!-----------------------------------------------------------------------
      orig=covar
      lwork=15
      dim=3

      call dgesvd('A','A',dim,dim,covar,dim,sigma,u,dim,vt,dim,work,&
           lwork,info)

      covar=orig

      if (info.ne.0) then
         write(6,'(/,2x,a,/)') 'SVD of the covariance matrix in &
              subroutine get_rot_mat failed'
         STOP
      endif

!-----------------------------------------------------------------------
! Calculate the rotation matrix
!-----------------------------------------------------------------------
      vut=0.0d0
      do i=1,3
         do j=1,3
            do k=1,3
               vut(i,j)=vut(i,j)+vt(k,i)*u(j,k)
            enddo
         enddo
      enddo

      det=finddet(covar,3)

      if (det.lt.0.0d0) then
         d=-1.0d0
      else
         d=1.0d0
      endif

      pm(1)=1.0d0
      pm(2)=1.0d0
      pm(3)=d

      ! 3x3 rotation matrix
      mat=0.0d0
      do i=1,3
         do j=1,3
            do k=1,3
               mat(i,j)=mat(i,j)+vt(k,i)*pm(k)*u(j,k)
            enddo
         enddo
      enddo

      ! 3Nx3N rotation matrix
      rotmat=0.0d0
      do i=1,natm
         k=i*3-2
         ilbl=0
         do m=k,k+2
            ilbl=ilbl+1
            jlbl=0
            do n=k,k+2
               jlbl=jlbl+1
               rotmat(m,n)=mat(ilbl,jlbl)
            enddo
         enddo
      enddo 

      return

    end subroutine get_rot_mat

!#######################################################################

    function finddet(matrix,n)

      implicit none

      integer*4              :: n,i,j,k,l
      real*8, dimension(n,n) :: matrix
      real*8                 :: m,temp,finddet
      logical(kind=4)        :: detexists = .TRUE.      
      
      l=1

      !Convert to upper triangular form
      do k=1,n-1
         if (matrix(k,k).eq.0) then
            detexists = .false.
            do i=k+1,n
               if (matrix(i,k).ne.0) then
                  do j=1,n
                     temp=matrix(i,j)
                     matrix(i,j)=matrix(k,j)
                     matrix(k,j)=temp
                  enddo
                  detexists=.true.
                  l=-l
                  exit
               endif
            enddo
            if (.not.detexists) then
               finddet=0
               return
            endif
         endif
         do j=k+1,n
            m=matrix(j,k)/matrix(k,k)
            do i=k+1,n
               matrix(j,i)=matrix(j,i)-m*matrix(k,i)
            enddo
         enddo
      enddo
      
      !Calculate determinant by finding product of diagonal elements
      finddet=l
      do i=1,n
         finddet=finddet*matrix(i,i)
      enddo
   
    end function finddet

!#######################################################################

    subroutine projhess

      use globalmod
      use projmod

      implicit none

      integer                      :: i,j,k,l,e2,error
      real*8, dimension(ncoo,ncoo) :: qproj,phess,tmp2
      real*8, dimension(ncoo)      :: tmp1,lambda
      real*8, dimension(3*ncoo)    :: work
      real*8                       :: sum
      real*8, parameter            :: ortho_tol=1d-1

!-----------------------------------------------------------------------
! Form the projector onto the subspace orthogonal to that spanned by
! the vectors
!-----------------------------------------------------------------------
         qproj=0.0d0
         do i=1,ncoo
            qproj(i,i)=1.0d0
         enddo

         if (lproj) then         
            do i=1,nvec
               do j=1,ncoo
                  do k=1,ncoo
                     qproj(j,k)=qproj(j,k)-vec(i,j)*vec(i,k)
                  enddo
               enddo
            enddo
         endif

!-----------------------------------------------------------------------
! Contribution from the projector onto the translational and 
! rotational DOFs
!-----------------------------------------------------------------------
      if (lprojrt) call rtproj(qproj)

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
         if (sum.lt.ortho_tol) then
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

    subroutine rtproj(qproj)
      
      use globalmod
      use projmod

      implicit none

      integer                      :: i,j,k,l
      real*8, dimension(ncoo,ncoo) :: qproj
      real*8, dimension(6,ncoo)    :: rtvec


      real*8, dimension(6,6)       :: smat,invsmat
      real*8, dimension(ncoo,6)    :: bmat
      real*8, dimension(ncoo,ncoo) :: rmat
      logical(kind=4)              :: lcheck

      real*8, dimension(6)         :: work
      integer*4, dimension(6)      :: ipiv
      integer*4                    :: info

!------------------------------------------------------------------
! Initialise arrays
!------------------------------------------------------------------
      rtvec=0.0d0

!------------------------------------------------------------------
! Vectors 1-3: translation along the three Cartesian axes
!------------------------------------------------------------------
      ! Loop over the translational DOFs
      do i=1,3
         ! Construct the vector for the current DOF
         do j=1,natm
            k=j*3-3+i
            rtvec(i,k)=sqrt(mass(j*3))
         enddo
      enddo

!------------------------------------------------------------------
! Vectors 4-6: infinitesimal displacements corresponding to
!              rotation about the three Cartesian axes
!------------------------------------------------------------------
      ! Rotation about the x-axis
      do i=1,natm
         j=i*3-1
         k=i*3
         rtvec(4,j)=sqrt(mass(i*3))*xcoo0(k)
         rtvec(4,k)=-sqrt(mass(i*3))*xcoo0(j)
      enddo

      ! Rotation about the y-axis
      do i=1,ncoo/3
         j=i*3-2
         k=i*3
         rtvec(5,j)=-sqrt(mass(i*3))*xcoo0(k)
         rtvec(5,k)=sqrt(mass(i*3))*xcoo0(j)
      enddo

      ! Rotation about the z-axis
      do i=1,natm
         j=i*3-2
         k=i*3-1
         rtvec(6,j)=sqrt(mass(i*3))*xcoo0(k)
         rtvec(6,k)=-sqrt(mass(i*3))*xcoo0(j)
      enddo

!------------------------------------------------------------------
! Calculate the projector R onto the translational and rotational
! DOFs using R=b*S^-1*b^T, where S=vec^T*vec.
!
! Here, R <-> rmat, S <-> smat, b <-> bmat (matrix of vectors)
!------------------------------------------------------------------
      ! Construct the b-matrix
      bmat=0.0d0
      do i=1,6
         do j=1,ncoo
            bmat(j,i)=rtvec(i,j)
         enddo
      enddo

      ! Calculate the S-matrix
      smat=0.0d0
      do i=1,6
         do j=1,6
            do k=1,ncoo
               smat(i,j)=smat(i,j)+bmat(k,i)*bmat(k,j)
            enddo
         enddo
      enddo

      ! Invert the S-matrix
      invsmat=smat
      call dgetrf(6,6,invsmat,6,ipiv,info)
      if (info.ne.0) then
         write(6,'(/,a,/)') &
              "LU factorisation of the S-matrix failed"
         STOP
      endif
      call dgetri(6,invsmat,6,ipiv,work,6,info)
      if (info.ne.0) then
         write(6,'(/,a,/)') &
              "Diagonalisation of the S-matrix failed"
         STOP
      endif

      ! Calculate the projection matrix R <-> rmat
      rmat=0.0d0
      do i=1,ncoo
         do j=1,ncoo
            do k=1,6
               do l=1,6
                  rmat(i,j)=rmat(i,j)+&
                       bmat(i,k)*invsmat(k,l)*bmat(j,l)
               enddo
            enddo
         enddo
      enddo

      ! Subtract the projection matrix R from the total projector
      qproj=qproj-rmat

      return
      
    end subroutine rtproj

!#######################################################################

    subroutine trans_hess

      use globalmod
      use projmod

      implicit none

      integer                      :: i,j,k,imax,jmax,unit
      real*8, dimension(ncoo,ncoo) :: tmpmat
      real*8, dimension(ncoo)      :: tmp1
      real*8                       :: max,ftmp,length
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
