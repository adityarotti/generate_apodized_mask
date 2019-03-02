use global_param
use process_mask
 

implicit none

read(*,*) paramfile
call read_param(paramfile)
print*, "Done reading paramfile"

call allocate_mask_arrays()

call read_data_mask()
call fill_holes()
call gen_apo_mask()

call deallocate_mask_arrays()

stop 
end
