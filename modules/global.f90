module global_param
use healpix_types
use healpix_modules

implicit none

integer(I4B),parameter :: nlheader=80
integer(I4B) :: nside,npixtot
real(dp) :: nullval,fmissval=-1.6375e30
real(dp) :: fillholesize,fillislandsize,apowidth ! all in arcminutes
real(dp), allocatable, dimension(:,:) :: dw8
real(dp), dimension(2) :: zbounds=0.d0

character(LEN=200) ::pathout,datafnamein,masknameout,filename,paramfile,prefix,suffix
character(len=80),dimension(1)  :: mapheader(1:nlheader)
character(len=8) :: s1,s2

logical :: anynull,swR2N

! Note these variables are transparents to all the modules. Hence one needs to be very carefully while making edits to this file.

contains

!########################################################################
subroutine read_param(paramfile)

implicit none
character(LEN=200) :: paramfile

open(10,file=paramfile)
! Parameters common to simulations and observed maps.
read(10,*) ! Reading in parameter common to simulations and observed maps.
read(10,*) nside
read(10,*) fillholesize,fillislandsize
read(10,*) apowidth
read(10,*) datafnamein
read(10,*) prefix
read(10,*) suffix
read(10,*) pathout
close(10)

write(s1,"(i8.1)") int(fillholesize) ; s1="_fh"//trim(adjustl(s1))//"a"
write(s2,"(i8.1)") int(apowidth) ; s2="_apo"//trim(adjustl(s2))//"a"
masknameout=trim(adjustl(prefix))//trim(adjustl(s1))//trim(adjustl(s2))//trim(adjustl(suffix))

fillholesize=(fillholesize/60.d0)*(pi/180.d0)
fillislandsize=(fillislandsize/60.d0)*(pi/180.d0)
apowidth=(apowidth/60.d0)*(pi/180.d0)

npixtot=nside2npix(nside)

apowidth=int(apowidth/sqrt(4.d0*pi/float(npixtot)))*sqrt(4.d0*pi/float(npixtot))
print*, apowidth*180.d0*60.d0/pi

call write_minimal_header(mapheader,"MAP",nside=nside,ordering="RING",polar=.True.)

end subroutine read_param
!########################################################################

end module global_param
