
!***************************************************************************
!* This is a source file of the Adaptive Data Transfer Library version 1.0
!* This file was initially finished by Cheng Zhang
!* If you have any problem,
!* please contact Cheng Zhang via zhang-cheng09@mails.tsinghua.edu.cn
!***************************************************************************

MODULE data_transfer_library_mod
    IMPLICIT NONE
    public :: data_transfer_register_instance
    public :: data_transfer_register_decomp
    public :: data_transfer_register_mask
    public :: data_transfer_init_instance
    public :: data_transfer_exec_instance
    public :: data_transfer_final_instance
    interface data_transfer_register_field ; module procedure &
        data_transfer_register_double_1D_field, &
        data_transfer_register_float_1D_field, &
        data_transfer_register_integer_1D_field, &
        data_transfer_register_logical_1D_field, &
        data_transfer_register_double_2D_field, &
        data_transfer_register_float_2D_field, &
        data_transfer_register_integer_2D_field, &
        data_transfer_register_logical_2D_field, &
        data_transfer_register_double_3D_field, &
        data_transfer_register_float_3D_field, &
        data_transfer_register_integer_3D_field, &
        data_transfer_register_logical_3D_field, &
        data_transfer_register_double_4D_field, &
        data_transfer_register_float_4D_field, &
        data_transfer_register_integer_4D_field, &
        data_transfer_register_logical_4D_field
    end interface

    integer, parameter :: R8 = selected_real_kind(12)
    integer, parameter :: R4 = selected_real_kind(6)
    integer,parameter :: I8 = selected_int_kind ( 12)
    integer,parameter :: I4 = selected_int_kind ( 6)

    contains

    subroutine data_transfer_register_double_1D_field(instance_id, data_buf, input)
        implicit none
        include "mpif.h"
        integer(I4), intent(in) :: instance_id
        real(R8), dimension(:) :: data_buf
        logical(1), intent(in) :: input
        integer(I4)             :: data_size

        data_size = size(data_buf)
        call register_field(instance_id, data_buf, data_size, "real8"//char(0), input)

    end subroutine data_transfer_register_double_1D_field

    subroutine data_transfer_register_double_2D_field(instance_id, data_buf, input)
        implicit none
        integer(I4), intent(in) :: instance_id
        real(R8), dimension(:,:) :: data_buf
        logical(1), intent(in) :: input
        integer(I4)             :: data_size

        data_size = size(data_buf)
        call register_field(instance_id, data_buf, data_size, "real8"//char(0), input)

    end subroutine data_transfer_register_double_2D_field

    subroutine data_transfer_register_double_3D_field(instance_id, data_buf, input)
        implicit none
        integer(I4), intent(in) :: instance_id
        real(R8), dimension(:,:,:) :: data_buf
        logical(1), intent(in) :: input
        integer(I4)             :: data_size

        data_size = size(data_buf)
        call register_field(instance_id, data_buf, data_size, "real8"//char(0), input)

    end subroutine data_transfer_register_double_3D_field

    subroutine data_transfer_register_double_4D_field(instance_id, data_buf, input)
        implicit none
        integer(I4), intent(in) :: instance_id
        real(R8), dimension(:,:,:,:) :: data_buf
        logical(1), intent(in) :: input
        integer(I4)             :: data_size

        data_size = size(data_buf)
        call register_field(instance_id, data_buf, data_size, "real8"//char(0), input)

    end subroutine data_transfer_register_double_4D_field

    subroutine data_transfer_register_float_1D_field(instance_id, data_buf, input)
        implicit none
        integer(I4), intent(in) :: instance_id
        real(R4), dimension(:) :: data_buf
        logical(1), intent(in) :: input
        integer(I4)             :: data_size

        data_size = size(data_buf)
        call register_field(instance_id, data_buf, data_size, "real4"//char(0), input)

    end subroutine data_transfer_register_float_1D_field

    subroutine data_transfer_register_float_2D_field(instance_id, data_buf, input)
        implicit none
        integer(I4), intent(in) :: instance_id
        real(R4), dimension(:,:) :: data_buf
        logical(1), intent(in) :: input
        integer(I4)             :: data_size

        data_size = size(data_buf)
        call register_field(instance_id, data_buf, data_size, "real4"//char(0), input)

    end subroutine data_transfer_register_float_2D_field

    subroutine data_transfer_register_float_3D_field(instance_id, data_buf, input)
        implicit none
        integer(I4), intent(in) :: instance_id
        real(R4), dimension(:,:,:) :: data_buf
        logical(1), intent(in) :: input
        integer(I4)             :: data_size

        data_size = size(data_buf)
        call register_field(instance_id, data_buf, data_size, "real4"//char(0), input)

    end subroutine data_transfer_register_float_3D_field

    subroutine data_transfer_register_float_4D_field(instance_id, data_buf, input)
        implicit none
        integer(I4), intent(in) :: instance_id
        real(R4), dimension(:,:,:,:) :: data_buf
        logical(1), intent(in) :: input
        integer(I4)             :: data_size

        data_size = size(data_buf)
        call register_field(instance_id, data_buf, data_size, "real4"//char(0), input)

    end subroutine data_transfer_register_float_4D_field

    subroutine data_transfer_register_integer_1D_field(instance_id, data_buf, input)
        implicit none
        integer(I4), intent(in) :: instance_id
        integer(I4), dimension(:) :: data_buf
        logical(1), intent(in) :: input
        integer(I4)             :: data_size

        data_size = size(data_buf)
        call register_field(instance_id, data_buf, data_size, "integer"//char(0), input)

    end subroutine data_transfer_register_integer_1D_field

    subroutine data_transfer_register_integer_2D_field(instance_id, data_buf, input)
        implicit none
        integer(I4), intent(in) :: instance_id
        integer(I4), dimension(:,:) :: data_buf
        logical(1), intent(in) :: input
        integer(I4)             :: data_size

        data_size = size(data_buf)
        call register_field(instance_id, data_buf, data_size, "integer"//char(0), input)

    end subroutine data_transfer_register_integer_2D_field

    subroutine data_transfer_register_integer_3D_field(instance_id, data_buf, input)
        implicit none
        integer(I4), intent(in) :: instance_id
        integer(I4), dimension(:,:,:) :: data_buf
        logical(1), intent(in) :: input
        integer(I4)             :: data_size

        data_size = size(data_buf)
        call register_field(instance_id, data_buf, data_size, "integer"//char(0), input)

    end subroutine data_transfer_register_integer_3D_field

    subroutine data_transfer_register_integer_4D_field(instance_id, data_buf, input)
        implicit none
        integer(I4), intent(in) :: instance_id
        integer(I4), dimension(:,:,:,:) :: data_buf
        logical(1), intent(in) :: input
        integer(I4)             :: data_size

        data_size = size(data_buf)
        call register_field(instance_id, data_buf, data_size, "integer"//char(0), input)

    end subroutine data_transfer_register_integer_4D_field

    subroutine data_transfer_register_logical_1D_field(instance_id, data_buf, input)
        implicit none
        integer(I4), intent(in) :: instance_id
        logical(1), dimension(:) :: data_buf
        logical(1), intent(in) :: input
        integer(i4)             :: data_size

        data_size = size(data_buf)
        call register_field(instance_id, data_buf, data_size, "logical"//char(0), input)

    end subroutine data_transfer_register_logical_1D_field

    subroutine data_transfer_register_logical_2D_field(instance_id, data_buf, input)
        implicit none
        integer(I4), intent(in) :: instance_id
        logical(1), dimension(:,:) :: data_buf
        logical(1), intent(in) :: input
        integer(I4)             :: data_size

        data_size = size(data_buf)
        call register_field(instance_id, data_buf, data_size, "logical"//char(0), input)

    end subroutine data_transfer_register_logical_2D_field

    subroutine data_transfer_register_logical_3D_field(instance_id, data_buf, input)
        implicit none
        integer(I4), intent(in) :: instance_id
        logical(1), dimension(:,:,:) :: data_buf
        logical(1), intent(in) :: input
        integer(I4)             :: data_size

        data_size = size(data_buf)
        call register_field(instance_id, data_buf, data_size, "logical"//char(0), input)

    end subroutine data_transfer_register_logical_3D_field

    subroutine data_transfer_register_logical_4D_field(instance_id, data_buf, input)
        implicit none
        integer(I4), intent(in) :: instance_id
        logical(1), dimension(:,:,:,:) :: data_buf
        logical(1), intent(in) :: input
        integer(I4)             :: data_size

        data_size = size(data_buf)
        call register_field(instance_id, data_buf, data_size, "logical"//char(0), input)

    end subroutine data_transfer_register_logical_4D_field

    integer(I4) function data_transfer_register_instance(local_comm, global_rank_remote_root, action)
        implicit none
        integer(I4), intent(in) :: local_comm
        integer(I4), intent(in) :: global_rank_remote_root
        integer(I4), intent(in) :: action

        call register_instance(local_comm, global_rank_remote_root, action, data_transfer_register_instance)

    end function data_transfer_register_instance

    subroutine data_transfer_register_decomp(instance_id, num_grid_cells, num_local_cells, local_cells_global_index)
        implicit none
        integer(I4), intent(in) :: instance_id
        integer(I4), intent(in) :: num_local_cells
        integer(I4), intent(in) :: num_grid_cells
        integer(I4), dimension(:) :: local_cells_global_index
        
        call register_decomp(instance_id, num_grid_cells, num_local_cells, local_cells_global_index);

    end subroutine data_transfer_register_decomp

    subroutine data_transfer_register_mask(instance_id, mask)
        implicit none
        integer(I4), intent(in) :: instance_id
        logical(1), dimension(:) :: mask

        call register_mask(instance_id, mask)

    end subroutine data_transfer_register_mask

    subroutine data_transfer_init_instance(instance_id)
        implicit none
        integer(I4), intent(in) :: instance_id

        call init_instance(instance_id)
        
    end subroutine data_transfer_init_instance

    subroutine data_transfer_exec_instance(instance_id)
        implicit none
        integer(I4), intent(in) :: instance_id

        call exec_instance(instance_id)
        
    end subroutine data_transfer_exec_instance

    subroutine data_transfer_final_instance(instance_id)
        implicit none
        integer(I4), intent(in) :: instance_id

        call final_instance(instance_id)
        
    end subroutine data_transfer_final_instance

END MODULE data_transfer_library_mod
