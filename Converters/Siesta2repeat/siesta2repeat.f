C
C    rho2xsf,  a script to transform 3-dim grid function
C             (i.e. LDOS, RHO, DRHO, etc.) written in SIESTA by subr. iorho
C             into an arbitrarily chosen grid for XCrysden,
C             using linear 4-point interpolation 
C
C   !!! --------------------  IMPORTANT  ---------------------- !!!
C   compile this code with the same compiler switches as Siesta,
C         in what regards using single/double precision,
C       otherwise reading the data from unformatted files 
C              written by iorho.F  might be spoiled.  
C                Don't say you haven't been warned.
C                                !!!
C
C             Written by Andrei Postnikov, Mar 2006   Vers_0.3
C             apostnik@uos.de
C
C             Modified by ekadants to write a gaussian cube
C             electrostatic potential 
C             is scaled by -0.5 to convert into "normal atomic units" (Hartrees) 
C             conform to REPEAT charge sign convention
C
C	      Modified by ekadants Sep.7, 2012 to check for the atoms that 
C             are translationally invariant
C
      program siesta2repeat
      implicit none
      integer ii1,ii2,io1,is1
      parameter (ii1=11,ii2=12,io1=14,is1=13)
      integer mesh0(3),mesh1(3),ip,nspin,is,ii,jj,n1,n2,n3,
     .        iat,nat,nz,iix,iiy,iiz,ind,mn,ibox,nboxun,nbox,
     .        ixmax,iymax,izmax,ixmin,iymin,izmin,ialloc,
     .        ishift,jshift,kshift,n1div,n2div,n3div,i1,i2,i3,indcur
     
      integer limit,maxa,maxo,maxuo,maxnh,maxna,il,ia,nrela(3),idim

      parameter (limit=5)  !  tried translations along each lattice vector
      character inpfil*60,outfil*60,syslab*30,suffix*6,
     .          unitlab*1,labunit*9,labbox*1,owrite*1
      logical unitb,charge,waves,filexist
      double precision b2ang,cc_bohr(3,3),cc_ang(3,3),cc_inv(3,3),
     .                 coort(3),obox(3),rbox(3,3),rinv(3,3),coortii(3),
     .                 cell(3,3),dum,rmaxo,rela,modu,rmesh(3),drela(3)
      double precision tmp1, tmp2, tmp3, smtmp
      double precision coortjj(3)
      parameter (b2ang=0.529177)   !  Bohr to Angstroem
      integer, allocatable :: ityp(:),iz(:),coor_box_lunique(:)
      double precision, allocatable :: mass(:),coor_ang(:,:)
      double precision, allocatable :: coor_box_xyz(:,:)
      character(len=2), allocatable :: label(:)
      real, allocatable :: func(:) !  NB! single precision, as in iorho.F
      real   fintp,fmax,fmin,fau,fmaxau,fminau       !  NB! single precision, as in iorho.F
	  
      real sigma, core_charge
      REAL :: pi

	  
      external test_xv,read_xv,fillbox,inver3,intpl04
      pi = 4 * ATAN(1.)	  
C
C     string manipulation functions in Fortran used below:
C     len_trim(string): returns the length of string 
C                       without trailing blank characters,
C     char(integer)   : returns the character in the specified position
C                       of computer's ASCII table, i.e. char(49)=1
      
	  
	  
	  
	  
	  
C      write (6,701,advance="no")
C  701 format(" Specify  SystemLabel (or 'siesta' if none): ")
C   read in system label
C 
      read (5,*) syslab
      inpfil = syslab(1:len_trim(syslab))//'.XV'
      open (ii1,file=inpfil,form='formatted',status='old',err=801)
	  
	  
	  
	  
      call test_xv(ii1,nat)
      allocate (ityp(nat))
      allocate (iz(nat))
      allocate (mass(nat))
      allocate (label(nat))
      allocate (coor_ang(1:3,1:nat),STAT=ialloc)
      if (ialloc.ne.0) then
         write (6,*) ' Fails to allocate space for ',nat,' atoms.'
         stop
      endif
      call read_xv(ii1,nat,ityp,iz,cc_ang,mass,label,coor_ang)
      call inver3(cc_ang,cc_inv)
      close (ii1)
	  
C --- set up and fill output box:     
C      call makebox(obox,rbox)
C  read in origin makebox output is in Bohr
C      write(6,*) ' cc_ang1x ', cc_ang(1,1)
C      write(6,*) ' cc_ang1y ', cc_ang(2,1)
C      write(6,*) ' cc_ang1z ', cc_ang(3,1)
C      write(6,*) ' cc_ang2x ', cc_ang(1,2)
C      write(6,*) ' cc_ang2y ', cc_ang(2,2)
C      write(6,*) ' cc_ang2z ', cc_ang(3,2)
C      write(6,*) ' cc_ang3x ', cc_ang(1,3)
C      write(6,*) ' cc_ang3y ', cc_ang(2,3)
C      write(6,*) ' cc_ang3z ', cc_ang(3,3)      
      
      read  (5,*) (obox(ii),ii=1,3)
      write(6,*) ' originx ', obox(1)
      write(6,*) ' originy ', obox(2)
      write(6,*) ' originz ', obox(3)
      
      rbox(1,1) = cc_ang(1,1)
      rbox(2,1) = cc_ang(2,1)
      rbox(3,1) = cc_ang(3,1)
      
      rbox(1,2) = cc_ang(1,2)
      rbox(2,2) = cc_ang(2,2)
      rbox(3,2) = cc_ang(3,2)
      
      rbox(1,3) = cc_ang(1,3)
      rbox(2,3) = cc_ang(2,3)
      rbox(3,3) = cc_ang(3,3)
      write(6,*) ' rbox1x ', rbox(1,1)
      write(6,*) ' rbox1y ', rbox(2,1)
      write(6,*) ' rbox1z ', rbox(3,1)
      
      write(6,*) ' rbox2x ', rbox(1,2)
      write(6,*) ' rbox2y ', rbox(2,2)
      write(6,*) ' rbox2z ', rbox(3,2)      
      
      write(6,*) ' rbox3x ', rbox(1,3)
      write(6,*) ' rbox3y ', rbox(2,3)
      write(6,*) ' rbox3z ', rbox(3,3)      

C
C   spanning vectors are copied from  XV file
C 
C
C
C      write (6,705,advance="no") '1st',labunit
C      read  (5,*) (rbox(ii,1),ii=1,3)
C      write (6,705,advance="no") '2nd',labunit
C      read  (5,*) (rbox(ii,2),ii=1,3)
C      write (6,705,advance="no") '3rd',labunit
C      read  (5,*) (rbox(ii,3),ii=1,3)
C	
C --- invert the box vectors; will need it in a minute...
      call inver3(rbox,rinv)
	  
	  
	  
	  
C --- write selected atoms first into a scratch file (is1), for the case
C     there are zero. Then the label 'ATOMS' with no  atoms following
C     will crash XCrySDen.
      open (is1, form='formatted',status='scratch')
      call fillbox(is1,obox,rbox,rinv,cc_ang,nat,coor_ang,nbox)

	  
	  
	  
C --- open output file:
      outfil = syslab(1:len_trim(syslab))//'.cube'
      inquire (file=outfil, exist=filexist)
      if (filexist) then
C        write (6,*) ' File ',outfil(1:len_trim(outfil)),' exists and will be overwritten'
	  write(6,*) 'cube file exists and will be overwritten'
          open (io1,file=outfil,form='formatted',status='REPLACE')
      else
        open (io1,file=outfil,form='formatted',status='NEW')
      endif
	  
	  
	  
	  
C     title	  
      write (io1,*) syslab
      write (io1,*) 'VH'
	  
	  
C     I5,3F12.6   #ATOMS, X-,Y-,Z-COORDINATES OF ORIGIN
C      write (6,*) ' The box contains ',nbox,' atoms.'
C   ekadants: check the atoms for the translational invariance
C   
      rewind (is1) !temporary file with atoms
      allocate (coor_box_xyz(1:3,1:nbox),STAT=ialloc)
      allocate (coor_box_lunique(1:nbox),STAT=ialloc)
      do ibox = 1,nbox
        read  (is1,'(i4,3f20.8)')  iat,     (coort(jj),jj=1,3)
	coor_box_xyz(1, ibox) = coort(1)
	coor_box_xyz(2, ibox) = coort(2)
	coor_box_xyz(3, ibox) = coort(3)
      enddo
   
     
      nboxun = nbox
      do ii=1,nbox
      	coortii(1) = coor_box_xyz(1, ii)
	coortii(2) = coor_box_xyz(2, ii)
	coortii(3) = coor_box_xyz(3, ii)
	coor_box_lunique(ii) = 1
      	do jj=ii+1,nbox
		do i1=-1,1
		do i2=-1,1
		do i3=-1,1
		 coortjj(1)=coor_box_xyz(1,jj)+i1*rbox(1,1)+i2*rbox(1,2)+i3*rbox(1,3)
		 coortjj(2)=coor_box_xyz(2,jj)+i1*rbox(2,1)+i2*rbox(2,2)+i3*rbox(2,3)
		 coortjj(3)=coor_box_xyz(3,jj)+i1*rbox(3,1)+i2*rbox(3,2)+i3*rbox(3,3)
		 tmp1 = (coortjj(1)-coortii(1))*(coortjj(1)-coortii(1))
		 tmp2 = (coortjj(2)-coortii(2))*(coortjj(2)-coortii(2))
		 tmp3 = (coortjj(3)-coortii(3))*(coortjj(3)-coortii(3))
		 smtmp=sqrt(tmp1+tmp2+tmp3+1.0e-15)
		 if(smtmp.le.1.0e-6) then 
		    coor_box_lunique(ii) = 0
	            nboxun = nboxun-1 
	            goto 11     
		 endif
		enddo
		enddo
		enddo
	enddo
   11 continue	
      enddo	
      
      write (io1,'(I5,3F12.6)') nboxun, (obox(ii),ii=1,3)

	  
C	  I5,3F12.6   #GRIDPOINTS, INCREMENT VECTOR
C      write (6,704) 
C  102 write (6,705,advance="no") 
C
C  enter number of grid points  
C
      read (5,*) n1,n2,n3
      if (n1.le.0.or.n2.le.0.or.n3.le.0) then
        write (6,*) ' Numbers must be positive '
        stop
      endif
      write (io1,'(I5,3F12.6)') n1, (rbox(ii,1)/(n1),ii=1,3)	  
      write (io1,'(I5,3F12.6)') n2, (rbox(ii,2)/(n2),ii=1,3)
      write (io1,'(I5,3F12.6)') n3, (rbox(ii,3)/(n3),ii=1,3)

	  

C     I5,4F12.6   ATOM NUMBER, CHARGE, X-,Y-,Z-COORDINATE
      rewind (is1) !temporary file with atoms
      do ibox = 1,nbox
        read  (is1,'(i4,3f20.8)')  iat,     (coort(jj),jj=1,3)
	if(coor_box_lunique(ibox).eq.1) then 
        write (io1,'(i5, 4f12.6)')
     .		iz(iat), iz(iat)*1., (coort(jj),jj=1,3)
         endif
      enddo
      close (is1)
      deallocate(coor_box_xyz)
      deallocate(coor_box_lunique)
	  

C     6E13.5      CUBE DATA (WITH Z INCREMENT MOVING FASTEST, THEN Y AND THEN X)
C --- Look for grid data files to include:
C    Open DATA file
C  103 write (6,706,advance="no")
C      read (5,*) suffix
C
C  
C      suffix = 'VH'
      inpfil = syslab(1:len_trim(syslab))//
     .      '.VH'
      open (ii2,file=inpfil,form='unformatted',status='old',err=806)
      write (6,*) 'Found and opened: ',inpfil(1:len_trim(inpfil))
      read (ii2,err=807) cell
      read (ii2,err=808) mesh0, nspin
      allocate (func(1:mesh0(1)*mesh0(2)*mesh0(3)),STAT=ialloc)
      if (ialloc.ne.0) then
        write (6,*) ' Fails to allocate space for ',
     .                mesh0(1)*mesh0(2)*mesh0(3),' grid points.'
        stop
      endif 
      write (6,*) 'mesh0 = (',mesh0,'),   nspin=',nspin
C      write (io1,"('BEGIN_BLOCK_DATAGRID_',i1,'D')") idim
C      write (io1,*) 'DATA_from:'//inpfil(1:len_trim(inpfil))
      do is=1,nspin
        fmax= -9.999E+10
        fmin=  9.999E+10
C        write (io1,*) 'BEGIN_DATAGRID_'//char(48+idim)//'D_'//
C     .        suffix(1:len_trim(suffix))//':spin_'//char(48+is)
C        if (n1.ne.1) write (io1,'(i6)',advance="no") n1
C        if (n2.ne.1) write (io1,'(i6)',advance="no") n2
C        if (n3.ne.1) write (io1,'(i6)',advance="no") n3
C        write (io1,'()') 
C        write (io1,'(1p,3e15.7)') (obox(ii),ii=1,3)
C        if (n1.ne.1) write (io1,'(1p,3e15.7)') (rbox(ii,1),ii=1,3)
C        if (n2.ne.1) write (io1,'(1p,3e15.7)') (rbox(ii,2),ii=1,3)
C        if (n3.ne.1) write (io1,'(1p,3e15.7)') (rbox(ii,3),ii=1,3)
        
C       read data		
        ind=0                        
        do iiz=1,mesh0(3)
          do iiy=1,mesh0(2)
            read (ii2,err=809,end=810) (func(ind+iix),iix=1,mesh0(1))
            ind = ind + mesh0(1)
          enddo
        enddo

C ---   loop over mesh points
C       avoid division by zero if only 1 point is selected: 
C        n1div=max(n1-1,1)  !  if (n1.eq.1) n1div=1 else n1div=n1-1
C        n2div=max(n2-1,1)
C        n3div=max(n3-1,1)

c       for periodic cells the first point is origin+0*step
c                 the last point is origin+n*step t.e. end-1step i.e. not on the box boundary

        ind = 0
	indcur = 0
        do i1=1,n1
          do i2=1,n2
            do i3=1,n3
              ind = ind+1
	      indcur = indcur + 1
              do ii=1,3
                rmesh(ii) = obox(ii) + 
     +                      rbox(ii,1)*(i1-1)/n1 +
     +                      rbox(ii,2)*(i2-1)/n2 +
     +                      rbox(ii,3)*(i3-1)/n3
              enddo
C ---         rmesh(1..3) are absolute Cartesian coordinates
C             of the mesh point (i1,i2,i3) in Bohr
C             Find its relative coordinates on the unit cell grid:
              do ii=1,3
                rela = 0.0
                do jj=1,3
                  rela = rela + cc_inv(ii,jj)*rmesh(jj)
                enddo
C               take modulo to make sure that it falls within [0,1]:
                modu = modulo( rela*mesh0(ii), dble(mesh0(ii)) )
                nrela(ii) = floor(modu) + 1
                drela(ii) = modu - nrela(ii) + 1
              enddo
C ---         mesh point rmesh(1..3) falls within the grid microcell
C             originating at the grid point nrela(1..3),
C             its relative coordinates within this microcell are drela(1..3)
C             Select neighboring grid points and make the interpolation:
              call intpl04 (func(1),fintp,
     .                      mesh0,nrela,drela)
              if (fintp.gt.fmax) then
                fmax = fintp
                ixmax = i1
                iymax = i2
                izmax = i3
              endif
              if (fintp.lt.fmin) then
                fmin = fintp
                ixmin = i1
                iymin = i2
                izmin = i3
              endif
			  
			  
c --  add Gaussian core charge if no core corrections are available --------------
c             if(suffix.eq.'RHO')then
c              core is abs zeroat about 3/4 of cov_radius
c              gaussian=0 at around 4sigma, ~0.1% at 3sigma
c              sigma=0.25*0.75*1.5
c              do ii=1,nat
c			    if(ii.eq.1) then
c                   core_charge=6.
c			    else
c                  core_charge=0.
c                end if			  
c                fintp=fintp + exp(-((rmesh(1)-coor_ang(1,ii))**2+
c     .                               (rmesh(2)-coor_ang(2,ii))**2+
c     .                               (rmesh(3)-coor_ang(3,ii))**2  )/
c     .                        (2*sigma**2))*
c     .                        (core_charge/(sqrt(2*pi)*sigma)**3) 
c              end do				
c              end if
			  
			  
	      fau = -0.5*fintp	  
              write (io1,'(1p,e13.5)',advance="no")
     $                          fau        ! write without linebreak
              if ((mod(indcur,6).eq.0) .or. (i3.eq.n3)) then 
	           write (io1,'()') ! linebreak after 6 numbers
		   if(i3.eq.n3) indcur = 0
	       endif 	   
            enddo    !  do i3 =
          enddo    !  do i2 =
        enddo    !  do i1 =
        write (io1,'()')  ! linebreak after all numbers
C        write (io1,"(' END_DATAGRID_',i1,'D')") idim
	fmaxau = -0.5*fmin
	fminau = -0.5*fmax
        write (6,206) is,'max',fmaxau,ixmax,iymax,izmax
        write (6,206) is,'min',fminau,ixmin,iymin,izmin
      enddo  ! do is=1,nspin
      deallocate (func)
C      write (io1,"('END_BLOCK_DATAGRID_',i1,'D')") idim
      close (ii2)
      stop                                                                 

	  
	  
	  
	  
	  
	  
	  
	  
	  
  201 format (3f12.6)
  202 format (i4,3f20.8)
  203 format (3i6)
  204 format (3f12.7)
  205 format (1p,6e13.6)
C 205 format (1p,8e10.3)
  206 format (' For is=',i1,': ',a3,'. grid value =',1p,e12.5,
     .        ' at iix,iiy,iiz=',3i4)

  702 format(" The atom section may appear in the XSF file",
     .       " as for periodic structure, that will allow you",
     .       " to replicate the units within the XCrySDen.",
     .       " This only makes sense if your selected box",
     .       " coincides with the true periodic cell. Otherwise",
     .       " the atoms will be written non-periodically, as for",
     .       " molecule. ")
  703 format(" Do you want atom section as for periodic ",
     .       " structure, Y or N ? ")
  704 format (" Now define the grid. If you want it two-dimensional,",/
     .   " give 1 as number of grid points along one spanning vector.") 
  705 format (" Enter number of grid points along three vectors: ")
  706 format (' Add grid property (LDOS, RHO, ...;',
     .        ' or BYE if none): ')

  801 write (6,*) ' Error opening file ',
     .            inpfil(1:len_trim(inpfil)),' as old formatted'
      stop
  806 write (6,*) ' A wild guess! There is no file ',
     .              inpfil(1:len_trim(inpfil)),'; close XSF and quit.'
      close (io1)
      stop
  807 write (6,*) ' Error reading cell vectors'
      stop
  808 write (6,*) ' Error reading n1,n2,n3,nspin '
      stop
  809 write (6,*) ' Error reading function values on the grid'
      stop
  810 write (6,*) ' Unexpected end of data for iiy=',iiy,
     .            ' iiz=',iiz,'  is=',is
      stop
  811 write (6,*) ' Error opening file ',
     .            outfil(1:len_trim(outfil)),' as new unformatted'
      stop
      end
