module process_mask
use global_param
use healpix_types
use healpix_modules

implicit none

integer(i4b), allocatable, dimension(:,:) :: mask
real(dp), allocatable, dimension(:,:) :: apomask
real(dp), allocatable, dimension(:,:) :: distance,tempmap

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
subroutine edit_data_ordering(mymap)
implicit none
real(dp) :: mymap(0:npixtot-1)

mask(:,1)=mymap(:)
if (swR2N) call convert_ring2nest(nside,mask)

end subroutine edit_data_ordering
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
endif

end subroutine fill_holes
!########################################################################

!########################################################################
subroutine gen_apo_mask(mymap)

implicit none
integer(i4b) :: i
real(dp) :: pixval
real(dp) :: mymap(0:npixtot-1)

call dist2holes_nest(nside, mask(:,1), distance(:,1))

if (swR2N) call convert_nest2ring(nside,distance)

do i=0,npixtot-1
   apomask(i,1)=1.d0
   if (distance(i,1).le.apowidth) then
      !pixval=1.d0-cos((distance(i,1)/apowidth)*(pi/2.d0))
      pixval=1.d0-cos((distance(i,1)/apowidth)*(pi/2.d0))**2.d0
      apomask(i,1)=pixval
   endif
enddo

mymap(:)=apomask(:,1)

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
