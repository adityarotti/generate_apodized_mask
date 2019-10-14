module gen_apodized_mask

use global_param
use process_mask
 

contains

subroutine return_apo_mask(inmask,nside,ring,fillholesize,fillislandsize,apowidth,npixtot,outmask)
logical, intent(in) :: ring
integer*4,intent(in) :: nside,npixtot
real*8,intent(in) :: fillholesize,fillislandsize,apowidth
real*8,intent(in) :: inmask(0:npixtot-1)
real*8,intent(out) :: outmask(0:npixtot-1)

call pass_param(nside,ring,fillholesize,fillislandsize,apowidth)
call allocate_mask_arrays()
call edit_data_ordering(inmask)
call fill_holes()
call gen_apo_mask(outmask)
call deallocate_mask_arrays()
end subroutine return_apo_mask

end module gen_apodized_mask
