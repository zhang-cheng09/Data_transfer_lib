#include "global_data.h"

extern "C" void register_instance_(MPI_Comm * local_comm, int * global_rank_remote_root, int * action, int * instance_id)
{
    if(data_transfer_instance_mgt == NULL) data_transfer_instance_mgt = new Data_transfer_instance_mgt();
    *instance_id = data_transfer_instance_mgt->register_data_transfer_instance(*local_comm, * global_rank_remote_root, * action);
}

extern "C" void register_decomp_(int * instance_id, int * num_grid_cells, int * num_local_cells, int * local_cells_global_index)
{
    data_transfer_instance_mgt->register_decomposition(*instance_id, * num_grid_cells, * num_local_cells, local_cells_global_index);
}

extern "C" void register_field_(int * instance_id, void * data_buf, int * data_size, char * data_type, bool *input)
{
    data_transfer_instance_mgt->register_field(*instance_id, data_buf, * data_size, data_type, * input);
}

extern "C" void register_mask_(int * instance_id, bool * mask)
{
    data_transfer_instance_mgt->register_mask(*instance_id, mask);
}

extern "C" void init_instance_(int * instance_id)
{
    data_transfer_instance_mgt->init_data_transfer_instance(*instance_id);
}

extern "C" void exec_instance_(int * instance_id)
{
    data_transfer_instance_mgt->exec_data_transfer_instance(*instance_id);
}

extern "C" void final_instance_(int * instance_id)
{
    data_transfer_instance_mgt->final_data_transfer_instance(*instance_id);
}
