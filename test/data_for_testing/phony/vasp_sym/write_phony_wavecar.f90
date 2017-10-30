!
!	Program to write a valid WAVECAR file where we control content for testing.
!
!	A. Linscheid, May 24  2017
!	P.s. based on vasp 5.4.1
!                                                            
program write_phony_wavecar
implicit none

integer, parameter :: q = SELECTED_REAL_KIND(10) !VASP convention
integer, parameter :: qs = SELECTED_REAL_KIND(5) !VASP convention
complex(q),allocatable :: 	coeff(:)  !PW coefficients
complex(qs),allocatable :: 	scoeff(:) !PW coefficients in single precision mode
complex(q),allocatable ::       fullFourierGrid(:,:,:)
real(q),allocatable :: 		bnds(:)   !Complex energies from file
real(q), allocatable :: 	occ(:)     !Occupations from file
real(q):: 			a(3,3),b(3,3),k(3), recVol
real(q)::			Rrecl, Rnspin, Rversion, Rnkp, Rnbnd, Rnpw !VASP stores integers as real numers ...
real(q)::			ecut, en
integer::			Irecl, nspin, fversion, nkp, nbnd, npw, i, j
logical::			dprec
integer::			ibnd, ik, isp, ipw, iost, irec, ig, vasp_ngx, vasp_ngy,vasp_ngz
real(q)::			energyWindow, dummy
integer, allocatable::          igfou(:,:)
integer::			igz,igy,igx,igzf,igyf,igxf, nxm,nym,nzm,nb(3),nGm,ng
real(q)::			gc,kg(3),g(3), range_min, range_max
double precision, allocatable ::         ldos(:,:,:)
real(q), allocatable ::         kgrid(:,:)

!Open the file correctly. Get record length, reopen and read header.
Irecl=100
open(unit=12,file='WAVECAR',access='direct',recl=Irecl,iostat=iost,status='new')
if (iost.ne.0) write(*,*) 'Error opening WAVECAR file for reading',iost            
nspin=1
fversion=45200

write(*,*) 'Selected precision = ', SELECTED_REAL_KIND(10)
write(*,*) 'record length =',Irecl,' spins =',nspin,  ' version ',fversion
if ( fversion .eq. 53300 ) then
  write(*,*) "VASP.5.3 WAVECAR created"
  dprec  = .false.
else if( fversion .eq. 53310) then
  write(*,*) "VASP.5.3 double precision WAVECAR created"
  dprec = .true.
elseif ( fversion .eq. 45200 ) then
  write(*,*) "VASP 5.4.1 WAVECAR created"
  dprec = .false.
else
  write(*,*) "double precision WAVECAR created"
  dprec = .true.
endif
write(unit=12,rec=1) real(Irecl,kind=q),real(nspin,kind=q),real(version,kind=q)

open(unit=12,file='WAVECAR',access='direct',recl=Irecl, &
     iostat=iost,status='old')
if (iost.ne.0) &
	write(*,*) 'Error opening file',iost            
read(unit=12,rec=2) Rnkp,Rnbnd,ecut,((a(i,j),i=1,3),j=1,3)
nkp=nint(Rnkp)
nbnd=nint(Rnbnd)
write(*,*) 'Statistics: nk =',nkp,'; nbands =',nbnd, 'PWcut =',ecut
write(*,*) 'Using energy window: [',range_min,',',range_max,']'
write(*,*) 'Lattice basis:'
write(*,*) 'a1 =',(a(j,1),j=1,3)
write(*,*) 'a2 =',(a(j,2),j=1,3)
write(*,*) 'a3 =',(a(j,3),j=1,3)

call crossprod(b(:,1),a(:,2),a(:,3))
call crossprod(b(:,2),a(:,3),a(:,1))
call crossprod(b(:,3),a(:,1),a(:,2))
recVol=DOT_PRODUCT(a(:,1),b(:,1))
write(*,*) "Unit cell volume =",recVol
b(:,:)=8.0d0*atan(1.0d0)/recVol*b(:,:)
write(*,*) 'Reciprocal cell'
write(*,*) 'b1 =',(b(j,1),j=1,3)
write(*,*) 'b2 =',(b(j,2),j=1,3)
write(*,*) 'b3 =',(b(j,3),j=1,3)

!Check the possible G vectors
en = 0.262465831d0
gc = sqrt(ecut*en)
nGm = int(gc) + 2
nzm=nint(2.0*gc/abs(dot_product(a(3,:),b(:,3)))*sqrt(DOT_PRODUCT(a(3,:),a(3,:))))
nGm = max(nGm,nzm)
nym=nint(2.0*gc/abs(dot_product(a(2,:),b(:,2)))*sqrt(DOT_PRODUCT(a(2,:),a(2,:))))
nGm = max(nGm,nym)
nxm=nint(2.0*gc/abs(dot_product(a(1,:),b(:,1)))*sqrt(DOT_PRODUCT(a(1,:),a(1,:))))
nGm = max(nGm,nxm)

nb = 0
do igz = -nGm, nGm
  do igy = -nGm, nGm
    do igx = -nGm, nGm
      g(:) = DBLE(igx) * b(:,1) + DBLE(igy) * b(:,2) + DBLE(igz) * b(:,3)
      if ( DOT_PRODUCT(g,g)/en < ecut) then
        nb(1) = MAX( nb(1), ABS( igx ) )
        nb(2) = MAX( nb(2), ABS( igy ) )
        nb(3) = MAX( nb(3), ABS( igz ) )
      endif
    enddo
  enddo
enddo

nxm = 2*(nb(1)+1)+1
nym = 2*(nb(2)+1)+1
nzm = 2*(nb(3)+1)+1

write(*,*) "Determined grid dimension ", nxm+1, nym+1 , nzm+1
!if ( nxm+1 .gt. vasp_ngx .or. nym+1  .gt. vasp_ngy .or. nzm+1  .gt. vasp_ngz ) then
!  write(*,*) "Error Fourier grid dimensions:",vasp_ngx,vasp_ngy,vasp_ngz," Charge dim smaller than wfct dim"
!  stop
!endif

allocate(fullFourierGrid(vasp_ngx,vasp_ngy,vasp_ngz))
allocate(ldos(vasp_ngx,vasp_ngy,vasp_ngz))
plan = fftw_plan_dft_3d(vasp_ngz,vasp_ngy,vasp_ngx, fullFourierGrid, fullFourierGrid,FFTW_FORWARD,FFTW_ESTIMATE)

ldos = 0

!!Open the kpts file with the weights
allocate(kgrid(4,nkp))
!open(unit=48,file='kpts.dat',iostat=iost,status='old')
!do ik=1,nkp
!  read(48,*) (kgrid(i,ik),i=1,4)
!enddo
!close(48)


!Read the payload
allocate(occ(nbnd))
allocate(bnds(nbnd))
irec=2
do isp=1,nspin
   do ik=1,nkp
      irec=irec+1
      read(unit=12,rec=irec) Rnpw,k,(bnds(ibnd),dummy,occ(ibnd),ibnd=1,nbnd)
!      if ( abs( k(1) - kgrid(1,ik) ) > 1e-6 .or. abs( k(2) - kgrid(2,ik) ) > 1e-6 .or. abs( k(3) - kgrid(3,ik) ) > 1e-6 ) then
!        write(*,*) "Error comparing k grids! "
!        stop
!      endif
      npw=nint(Rnpw)
      allocate(igfou(3,npw))
      ng = 0
      do igz=0,vasp_ngz-1
         igzf=igz
         if (igzf.ge.vasp_ngz/2) igzf=igzf-vasp_ngz
         do igy=0,vasp_ngy-1
            igyf=igy
            if (igyf.ge.vasp_ngy/2) igyf=igyf-vasp_ngy
            do igx=0,vasp_ngx-1
               igxf=igx
               if (igxf.ge.vasp_ngx/2) igxf=igx-vasp_ngx
               kg(:)=(k(1)+igxf)*b(:,1)+ (k(2)+igyf)*b(:,2)+(k(3)+igzf)*b(:,3)
               if (dot_product(kg,kg)/en.lt.ecut) then

if ( ik .eq. 5 ) write (*,*) ng, igx, igy, igz
                  ng=ng+1
                  igfou(1,ng)=igx+1
                  igfou(2,ng)=igy+1
                  igfou(3,ng)=igz+1
               end if
            enddo
         enddo
      enddo
      if ( ng .ne. npw ) then
        write(*,*) "Error matching plane waves. Failed at (ik,ng,npw) : ",ik, ng , npw
        stop
      endif

      allocate(coeff(npw))
      if ( .not. dprec ) then
         allocate(scoeff(npw))
      endif
      do ibnd=1,nbnd
         irec=irec+1
         !we include only states within the energy range
!         if ( ((bnds(ibnd) - range_min) < 0) .or. ((bnds(ibnd) - range_max) > 0) ) cycle
         if ( dprec ) read(unit=12,rec=irec) (coeff(ipw), ipw=1,npw)
         if ( .not. dprec ) then
            read(unit=12,rec=irec) (scoeff(ipw), ipw=1,npw)
            coeff(:) = scoeff(:)
         endif
         fullFourierGrid = 0;
if ( ik .eq. 5 .and. ibnd .eq. 1)then
write(*,*) npw, bnds(ibnd)
write(*,*) k
write(*,*) (Real(coeff(i)),aimag(coeff(i)),i=1,npw)
write(*,*) igfou
endif

         do ig=1,ng
           fullFourierGrid(igfou(1,ig),igfou(2,ig),igfou(3,ig)) = coeff(ig)
         enddo
	 call fftw_execute_dft(plan, fullFourierGrid, fullFourierGrid)
         ldos(:,:,:) = ldos(:,:,:) + kgrid(4,ik)/real(vasp_ngx*vasp_ngy*vasp_ngz)&
                       *CONJG(fullFourierGrid(:,:,:))*fullFourierGrid(:,:,:)
      enddo
      if (allocated(coeff) ) deallocate(coeff)
      if (allocated(scoeff) ) deallocate(scoeff)
      deallocate(igfou)
   enddo
enddo

do igz=1,vasp_ngz
  do igy=1,vasp_ngy
    do igx=1,vasp_ngx
      write (10) ldos(igx,igy,igz)
    enddo
  enddo
enddo

deallocate(occ)
deallocate(bnds)
deallocate(fullFourierGrid)

call fftw_destroy_plan(plan)

close(12)
close(10)
stop
end program read_wavecar

subroutine crossprod(a,b,c)
  implicit none
  real(kind=SELECTED_REAL_KIND(10)), intent(out):: a(3)
  real(kind=SELECTED_REAL_KIND(10)), intent(in) :: b(3),c(3)
  a(1)=b(2)*c(3)-b(3)*c(2)
  a(2)=b(3)*c(1)-b(1)*c(3)
  a(3)=b(1)*c(2)-b(2)*c(1)
  return
end subroutine crossprod
