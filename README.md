This is the document of the adptive data trnasfer library 1.0 for model coupling.
If you have some problems, please contant Cheng Zhang via zhang-cheng09@mails.tsinghua.edu.cn.

The APIs of the library are written in Fortran because most couplers and models are programmed 
in Fortran, although the source code is mainly written in C++. 

1. data_transfer_register_instance: the API registers one data transfer instance, 
    which takes local communicator of this model, the global rank of the root process of the remote model,
    and the action (0 for send, 1 for recv, and 2 for sendrecv) of this instance as input, and returns the index of this instance.

2. data_transfer_register_decomposition: the API registers one parallel decomposition for one instance, 
    which takes the instance index, the number of grid cells, the number of local cells, 
    and the global index of local cells as input. One data transfer should register one (the instance action is send or recv) 
    or two parallel decompositions (the instance action is sendrecv).

3. data_transfer_register_field: the API registers one coupling field for one instance,
    which takes the instance index, the memory space of one coupling field, the action of one coupling field (true for input field, 
    and false for output field) as input.

4. data_transfer_register_mask: the API registers one array of mask for one instance,
    which takes the instance index, and the mask array as input. Each element of the mask array
    corresponds to one coupling field registered to this data transfer instance.

5. data_transfer_init_instance: the API intializes one data transfer instance,
    which takes the instance index as input.

6. data_transfer_exec_instance: the API executes one data transfer instance,
    which takes the instance index as input.

7. data_transfer_final_instance: the API finalizes one data transfer instance,
    which takes the instance index as input.


Here is an example to implement data transfer with the adaptive data trasfer library

Program X_model

    use data_transfer_library_mod

    implicit none 
    
    integer :: instance_id, decomp_id, end_step
    
    integer :: local_comm, remote_root_global_rank, direction
    
    integer :: num_grid_cells, num_local_cells, local_cells_global_index(:)
    
    logical :: mask_array(:), input
    
    real :: send_buf1(:,:), send_buf2(:,:,:)

    instance_id = data_transfer_register_instance(local_comm, remote_root_global_rank, direction)
    
    call data_transfer_register_decomp(instance_id, num_grid_cells, num_local_cells, local_cells_global_index)
    
    call data_transfer_register_field(instance_id, send_buf1, input)
    
    call data_transfer_register_field(instance_id, send_buf2, input)
    
    call data_transfer_register_mask(instance_id, mask_array)
    
    call data_transfer_init_instance(instance_id)
    
    do i=1, end_step
    
        call data_transfer_exec_instance(instance_id)
    
    end do
    
    call data_transfer_final_instance(instance_id)

End program X_model 
