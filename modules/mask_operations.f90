module process_mask
use global_param
use healpix_types
use healpix_modules

implicit none

real(dp), allocatable, dimension(:,:) :: apomask,distance,tempmap
integer(i4b), allocatable, dimension(:,:) :: mask

contains

!########################################################################
subroutine allocate_mask_arrays()

allocate(tempmap(0:npixtot-1,1),mask(0:npixtot-1,1),apomask(0:npixtot-1,1),distance(0:npixtot-1,1))

apomask=0.d0
mask=0
distance=0.d0
tempmap=0.d0

end subroutine allocate_mask_arrays

!-------------------------------------------------------------------------

subroutine deallocate_mask_arrays()

deallocate(mask,apomask,distance,tempmap)

end subroutine deallocate_mask_arrays
!########################################################################


!########################################################################
subroutine read_data_mask()
implicit none

character(len=80),dimension(1)  :: mapinheader(1:nlheader)
integer(i4b) :: i

call read_bintab(datafnamein,tempmap,npixtot,1,nullval,anynull,mapinheader)
mask=int(tempmap)

do i=1,nlheader
!   print*, i,mapinheader(i)
   if (index(mapinheader(i),"RING").gt.0) then
      swR2N=.true.
      print*, i,mapinheader(i)
   endif
enddo

if (swR2N) call convert_ring2nest(nside,mask)

print*, "Read in original mask"

end subroutine read_data_mask
!########################################################################

!########################################################################
subroutine fill_holes()

implicit none
integer(i4b) :: numpix2fill

numpix2fill=int(fillholesize**2.d0/((4.d0*pi)/float(npixtot)))
!print*, numpix2fill,"numpix2fill"

if (fillholesize.gt.0.d0) then

!  Remove holes.
   call fill_holes_nest(nside, numpix2fill, mask(:,1), mask(:,1)) 

!  Remove islands.
   numpix2fill=int(fillislandsize**2.d0/((4.d0*pi)/float(npixtot)))
   mask(:,1)=1-mask(:,1)
   call fill_holes_nest(nside, numpix2fill, mask(:,1), mask(:,1))
   mask(:,1)=1-mask(:,1)
!  --------------------------------------------------------------

   tempmap=float(mask) 
   if (swR2N) call convert_nest2ring(nside,tempmap)
   filename=trim(adjustl(pathout))//"mask_withfewerholes.fits"
   call write_data(tempmap)
   print*, "Filled holes in the mask"
endif

end subroutine fill_holes
!########################################################################

!########################################################################
subroutine gen_apo_mask()

implicit none
integer(i4b) :: i
real(dp) :: pixval

call dist2holes_nest(nside, mask(:,1), distance(:,1))

filename=trim(adjustl(pathout))//"distance.fits"
if (swR2N) call convert_nest2ring(nside,distance) 
call write_data(distance)
print*, "Written out the distance map"

! Since the distance map has already been converted to RING format the apodized mask 
! is generated in RING format and hence is not convertd.
do i=0,npixtot-1
   apomask(i,1)=1.d0
   if (distance(i,1).le.apowidth) then
      !pixval=1.d0-cos((distance(i,1)/apowidth)*(pi/2.d0))
      pixval=1.d0-cos((distance(i,1)/apowidth)*(pi/2.d0))**2.d0
      apomask(i,1)=pixval
   endif
enddo

filename=trim(adjustl(pathout))//trim(adjustl(masknameout))
call write_data(apomask)

end subroutine gen_apo_mask
!########################################################################

!########################################################################
subroutine write_data(mapout)

implicit none
real(dp), allocatable, dimension(:,:) :: mapout

call write_bintab(mapout,npixtot,1,mapheader,nlheader,"!"//filename)

end subroutine write_data
!########################################################################

end module process_mask
