module global_param
use healpix_types
use healpix_modules

implicit none
character(LEN=200):: datafnamein
integer(I4B) :: nside
real(dp) :: fillholesize,fillislandsize,apowidth ! all in arcminutes

integer(I4B),parameter :: nlheader=80
integer(I4B) :: npixtot
real(dp) :: nullval,fmissval=-1.6375e30
real(dp), allocatable, dimension(:,:) :: dw8
real(dp), dimension(2) :: zbounds=0.d0

character(LEN=200) ::pathout,masknameout,filename,paramfile,prefix,suffix
character(len=80),dimension(1)  :: mapheader(1:nlheader)
character(len=8) :: s1,s2

logical :: anynull,swR2N

! Note these variables are transparents to all the modules. Hence one needs to be very carefully while making edits to this file.

contains

!########################################################################
subroutine pass_param(mynside,mysw_R2N,holesize,islandsize,apo)
implicit none
integer(I4B) :: mynside
logical :: mysw_R2N
character(LEN=200) :: filename
real(dp) :: holesize,islandsize,apo

nside=mynside
npixtot=nside2npix(nside)
fillholesize=(holesize/60.d0)*(pi/180.d0)
fillislandsize=(islandsize/60.d0)*(pi/180.d0)
apowidth=(apo/60.d0)*(pi/180.d0)
apowidth=int(apowidth/sqrt(4.d0*pi/float(npixtot)))*sqrt(4.d0*pi/float(npixtot))
swR2N=mysw_R2N
!print*, "apowidth:",apowidth*180.d0*60.d0/pi,apo

!call write_minimal_header(mapheader,"MAP",nside=nside,ordering="RING",polar=.True.)

end subroutine pass_param
!########################################################################

end module global_param
