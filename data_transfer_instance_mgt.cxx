
/***************************************************************************
 * This is a source file of the Adaptive Data Transfer Library version 1.0
 * This file was initially finished by Cheng Zhang
 * If you have any problem,
 * please contact Cheng Zhang via zhang-cheng09@mails.tsinghua.edu.cn
 **************************************************************************/

#include "data_transfer_instance_mgt.h"
#include <cstring>

template <class T> void pack_segment_data(T *mpi_buf, T *field_data_buf, int segment_start, int segment_size, int field_2D_size, int num_lev)
{
    int i, j, offset;
    offset = 0;

    for (i = segment_start; i < segment_size+segment_start; i ++)
        for (j = 0; j < num_lev; j ++)
            mpi_buf[offset++] = field_data_buf[i+j*field_2D_size];
}

template <class T> void unpack_segment_data(T *mpi_buf, T *field_data_buf, int segment_start, int segment_size, int field_2D_size, int num_lev)
{
        int i, j, offset;


        for (i = segment_start, offset = 0; i < segment_size+segment_start; i ++)
            for (j = 0; j < num_lev; j ++)
                field_data_buf[i+j*field_2D_size] = mpi_buf[offset++];
}


bool words_are_the_same(const char * string1, const char * string2)
{
    return strcmp(string1, string2) == 0;
}

int get_data_type_size(const char *data_type)
{
    if (words_are_the_same(data_type, "real8"))
      return sizeof(double);
    else if (words_are_the_same(data_type, "real4"))
      return sizeof(float);
    else if (words_are_the_same(data_type, "long"))
      return sizeof(long);
    else if (words_are_the_same(data_type, "integer"))
      return sizeof(int);
    else if (words_are_the_same(data_type, "logical"))
      return sizeof(bool);
    else if (words_are_the_same(data_type, "char"))
      return sizeof(char);
    else if (words_are_the_same(data_type, "short"))
      return sizeof(short);
    else if (words_are_the_same(data_type, "string"))
      return 1024;

    return 0;
}

Data_transfer_instance::Data_transfer_instance(MPI_Comm local_comm, int remote_root_global_rank, int direction)
{
    this->local_comm = local_comm;
    this->remote_root_global_rank = remote_root_global_rank;
    this->direction = direction;
    mask = NULL;
    num_grid_cells = -1;
    num_local_cells = -1;
    num_remap_grid_cells = -1;
    num_remap_local_cells = -1;
    local_cells_global_index = NULL;
    local_remap_cells_global_index = NULL;
    remote_procs_global_rank = NULL;
    send_cells = NULL;
    recv_cells = NULL;
    num_input_fields = 0;
    num_output_fields = 0;
    routing_info_mgt = NULL;
    butterfly_instance = NULL;
    send_data_buf = NULL;
    recv_data_buf = NULL;
}

Data_transfer_instance::~Data_transfer_instance()
{
    for(int i=0; i<coupling_fields.size(); i++)
        delete coupling_fields[i];
    if(remote_procs_global_rank != NULL) delete [] remote_procs_global_rank;
    if(send_cells != NULL) delete [] send_cells;
    if(recv_cells != NULL) delete [] recv_cells;
    if(send_data_buf != NULL) delete [] send_data_buf;
    if(recv_data_buf != NULL) delete [] recv_data_buf;
    
    if(routing_info_mgt != NULL)
    {
        int tmp_size = remote_model_size;
        if(direction == SENDRECV) tmp_size = 2 * remote_model_size;
        if(num_local_cells > 0)
            for(int i=0; i<tmp_size; i++)
            {
                if(routing_info_mgt[i].num_elements_transferred > 0)
                {
                    delete [] routing_info_mgt[i].local_indx_segment_starts;
                    delete [] routing_info_mgt[i].local_indx_segment_lengths;
                }
            }
        delete [] routing_info_mgt;
    }

    if(butterfly_instance != NULL)
        delete butterfly_instance;
}

bool Data_transfer_instance::register_decomposition(int num_grid_cells, int num_local_cells, int * local_cells_global_index)
{
    if(this->num_local_cells < 0){
        this->num_grid_cells = num_grid_cells;
        this->num_local_cells = (num_local_cells > 0) ? num_local_cells : 0;
        this->local_cells_global_index = local_cells_global_index;
        return true;
    }
    else if(this->num_remap_local_cells < 0){
        this->num_remap_grid_cells = num_grid_cells;
        this->num_remap_local_cells = (num_local_cells > 0) ? num_local_cells : 0;
        this->local_remap_cells_global_index = local_cells_global_index;
        return true;
    }
    else
      return false;
}

bool Data_transfer_instance::register_field(void * data_buf, int buf_size, char * data_type, bool input)
{
    Field_info * field_info = new Field_info;
    field_info->data_buf = data_buf;
    field_info->buf_size = buf_size;
    field_info->data_type = get_data_type_size(data_type);
    //if(input) input_coupling_fields.push_back(field_info);
    //else output_coupling_fields.push_back(field_info);
    field_info->input = input;
    coupling_fields.push_back(field_info);
    if(input) num_input_fields++;
    else num_output_fields++;

    return true;
}

bool Data_transfer_instance::register_mask(bool * mask)
{
    if(this->mask != NULL) return false;
    this->mask = mask;
    if(this->mask != NULL) return true;
    else return false;
}

bool Data_transfer_instance::init()
{
    int local_size, local_rank, global_rank;
    MPI_Request send_req, recv_req;
    MPI_Status status;
    MPI_Comm_size(local_comm, &local_size);
    MPI_Comm_rank(local_comm, &local_rank);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    if(local_rank == 0)
    {
        MPI_Isend(&local_size, 1, MPI_INT, remote_root_global_rank, 1000, MPI_COMM_WORLD, &send_req);
        MPI_Irecv(&remote_model_size, 1, MPI_INT, remote_root_global_rank, 1000, MPI_COMM_WORLD, &recv_req);
        MPI_Wait(&send_req, &status);
        MPI_Wait(&recv_req, &status);
        MPI_Bcast(&remote_model_size, 1, MPI_INT, 0, local_comm);
        int * local_procs_global_rank = new int[local_size];
        remote_procs_global_rank = new int[remote_model_size];
        MPI_Gather(&global_rank, 1, MPI_INT, local_procs_global_rank, 1, MPI_INT, 0, local_comm);
        MPI_Isend(local_procs_global_rank, local_size, MPI_INT, remote_root_global_rank, 1000, MPI_COMM_WORLD, &send_req);
        MPI_Irecv(remote_procs_global_rank, remote_model_size, MPI_INT, remote_root_global_rank, 1000, MPI_COMM_WORLD, &recv_req);
        MPI_Wait(&send_req, &status);
        MPI_Wait(&recv_req, &status);
        MPI_Bcast(remote_procs_global_rank, remote_model_size, MPI_INT, 0, local_comm);
        delete [] local_procs_global_rank;
    }
    else
    {
        MPI_Bcast(&remote_model_size, 1, MPI_INT, 0, local_comm);
        remote_procs_global_rank = new int[remote_model_size];
        MPI_Gather(&global_rank, 1, MPI_INT, NULL, 1, MPI_INT, 0, local_comm);
        MPI_Bcast(remote_procs_global_rank, remote_model_size, MPI_INT, 0, local_comm);
    }

    routing_info_mgt = NULL;
   
    if(direction == SENDRECV) build_2D_self_router();
    else build_2D_remote_router();

    build_data_transfer_info();

    send_data_buf = NULL;
    recv_data_buf = NULL;

    int temp_send_size = 0;
    int temp_recv_size = 0;
    total_send_cells = 0;
    total_recv_cells = 0;

    if(send_cells != NULL)
      for(int i=0; i<remote_model_size; i++)
        total_send_cells += send_cells[i];

    if(recv_cells != NULL)
      for(int i=0; i<remote_model_size; i++)
        total_recv_cells += recv_cells[i];

    for(int i=0; i<coupling_fields.size(); i++)
    {
        if(coupling_fields[i]->input)
            temp_send_size += coupling_fields[i]->buf_size / num_local_cells * coupling_fields[i]->data_type * total_send_cells;
        else if(direction == SENDRECV)
            temp_recv_size += coupling_fields[i]->buf_size / num_remap_local_cells * coupling_fields[i]->data_type * total_recv_cells;
        else
            temp_recv_size += coupling_fields[i]->buf_size / num_local_cells * coupling_fields[i]->data_type * total_recv_cells;
    }

    if(temp_send_size > 0) send_data_buf = new char[temp_send_size];
    if(temp_recv_size > 0) recv_data_buf = new char[temp_recv_size];

    int fields_per_cell = 0;

    if(direction != RECV)
    {
        for(int i=0; i<coupling_fields.size(); i++)
            if(coupling_fields[i]->input)
                fields_per_cell += coupling_fields[i]->buf_size / num_local_cells * coupling_fields[i]->data_type;
    }
    else
    {
        for(int i=0; i<coupling_fields.size(); i++)
            if(!coupling_fields[i]->input)
                fields_per_cell += coupling_fields[i]->buf_size / num_local_cells * coupling_fields[i]->data_type;
    }

    butterfly_instance = new Butterfly(MPI_COMM_WORLD, local_comm, remote_procs_global_rank, send_cells, recv_cells, remote_model_size, direction);

    //if(direction == SEND && (output_coupling_fields.size() > 0 || input_coupling_fields.size() == 0) ) return false;
    //if(direction == RECV && (input_coupling_fields.size() > 0 || output_coupling_fields.size() == 0) ) return false;
    //if(direction == SENDRECV && (input_coupling_fields.size() == 0 || output_coupling_fields.size() == 0) ) return false;
    if(direction == SEND && (num_output_fields > 0 || num_input_fields == 0) ) return false;
    if(direction == RECV && (num_input_fields > 0 || num_output_fields == 0) ) return false;
    if(direction == SENDRECV && (num_input_fields == 0 || num_output_fields == 0) ) return false;
    return true;
}

void Data_transfer_instance::build_data_transfer_info()
{
    if(direction == SENDRECV){
        send_cells = new int[remote_model_size];
        recv_cells = new int[remote_model_size];
    }
    else if(direction == SEND)
        send_cells = new int[remote_model_size];
    else
        recv_cells = new int[remote_model_size];

    if(direction == SENDRECV)
        for(int i=0; i<remote_model_size; i++)
        {
            send_cells[i] = routing_info_mgt[i].num_elements_transferred;
            recv_cells[i] = routing_info_mgt[i+remote_model_size].num_elements_transferred;
        }
    else if(direction == SEND)
        for(int i=0; i<remote_model_size; i++)
        {
            send_cells[i] = routing_info_mgt[i].num_elements_transferred;
        }
    else
        for(int i=0; i<remote_model_size; i++)
        {
            recv_cells[i] = routing_info_mgt[i].num_elements_transferred;
        }
}

void Data_transfer_instance::build_2D_remote_router()
{
    int * num_cells_each_remote_proc = new int[remote_model_size];
    int ** cells_index_each_remote_proc = new int *[remote_model_size];
    routing_info_mgt = new Routing_info_with_one_process[remote_model_size];
    int global_rank, local_rank, local_size;
    MPI_Status status;
    
    MPI_Comm_size(local_comm, &local_size);
    MPI_Comm_rank(local_comm, &local_rank);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    
    int * num_cells_each_local_proc = new int[local_size];
    int * displs = new int[local_size];

    MPI_Allgather(&num_local_cells, 1, MPI_INT, num_cells_each_local_proc, 1, MPI_INT, local_comm);
    if(local_rank == 0){ 
        if(global_rank < remote_procs_global_rank[0])
        {
            MPI_Send(num_cells_each_local_proc, local_size, MPI_INT, remote_procs_global_rank[0], 2000+remote_procs_global_rank[0], MPI_COMM_WORLD);
            MPI_Recv(num_cells_each_remote_proc, remote_model_size, MPI_INT, remote_procs_global_rank[0], 2000+global_rank, MPI_COMM_WORLD, &status);
        }
        else
        {
            MPI_Recv(num_cells_each_remote_proc, remote_model_size, MPI_INT, remote_procs_global_rank[0], 2000+global_rank, MPI_COMM_WORLD, &status);
            MPI_Send(num_cells_each_local_proc, local_size, MPI_INT, remote_procs_global_rank[0], 2000+remote_procs_global_rank[0], MPI_COMM_WORLD);
        }

        displs[0] = 0;
        for(int i=1; i<local_size; i++) displs[i] = displs[i-1] + num_cells_each_local_proc[i-1];
    }
    MPI_Bcast(num_cells_each_remote_proc, remote_model_size, MPI_INT, 0, local_comm);

    int remote_grid_size = 0;
    for(int i=0; i<remote_model_size; i++) remote_grid_size += num_cells_each_remote_proc[i];
    int local_grid_size = 0;
    for(int i=0; i<local_size; i++) local_grid_size += num_cells_each_local_proc[i];
    
    int * cells_index_remote_procs = new int[remote_grid_size];
    int * cells_index_local_procs = new int[local_grid_size];
    
    MPI_Gatherv(local_cells_global_index, num_local_cells, MPI_INT, cells_index_local_procs, num_cells_each_local_proc, displs, MPI_INT, 0, local_comm);
    
    if(local_rank == 0){ 
        if(global_rank < remote_procs_global_rank[0])
        {
            MPI_Send(cells_index_local_procs, local_grid_size, MPI_INT, remote_procs_global_rank[0], 3000+remote_procs_global_rank[0], MPI_COMM_WORLD);
            MPI_Recv(cells_index_remote_procs, remote_grid_size, MPI_INT, remote_procs_global_rank[0], 3000+global_rank, MPI_COMM_WORLD, &status);
        }
        else
        {
            MPI_Recv(cells_index_remote_procs, remote_grid_size, MPI_INT, remote_procs_global_rank[0], 3000+global_rank, MPI_COMM_WORLD, &status);
            MPI_Send(cells_index_local_procs, local_grid_size, MPI_INT, remote_procs_global_rank[0], 3000+remote_procs_global_rank[0], MPI_COMM_WORLD);
        }
    }
    MPI_Bcast(cells_index_remote_procs, remote_grid_size, MPI_INT, 0, local_comm);

    remote_grid_size = 0;
    for(int i=0; i<remote_model_size; i++)
    {
        cells_index_each_remote_proc[i] = cells_index_remote_procs + remote_grid_size;
        remote_grid_size += num_cells_each_remote_proc[i];
    }
    
    if(num_local_cells > 0)
        for(int i=0; i<remote_model_size; i++)
            compute_routing_info_between_decomps(num_cells_each_remote_proc[i], cells_index_each_remote_proc[i], i);

    delete [] cells_index_each_remote_proc;
    delete [] num_cells_each_remote_proc;
    delete [] num_cells_each_local_proc;
    delete [] cells_index_remote_procs;
    delete [] cells_index_local_procs;
}

void Data_transfer_instance::compute_routing_info_between_decomps(int num_local_cells_remote, int *local_cells_global_index_remote, int remote_proc_id)
{
    const int * reference_cell_indx;
    int num_reference_cells;
    
    routing_info_mgt[remote_proc_id].num_elements_transferred = 0;
    routing_info_mgt[remote_proc_id].num_local_indx_segments = 0;

    if(num_local_cells_remote < num_local_cells){
        num_reference_cells = num_local_cells_remote;
        reference_cell_indx = local_cells_global_index_remote;
    }
    else{
        num_reference_cells = num_local_cells;
        reference_cell_indx = local_cells_global_index;
    }

    int * logical_indx_lookup_table_remote = new int[num_grid_cells];
    int * logical_indx_lookup_table_local = new int[num_grid_cells];

    for(int i=0; i<num_grid_cells; i++)
    {
        logical_indx_lookup_table_remote[i] = -1;
        logical_indx_lookup_table_local[i] = -1;
    }
    for(int i=0; i<num_local_cells; i++)
        if(local_cells_global_index[i] >= 0)
            logical_indx_lookup_table_local[local_cells_global_index[i]] = i;
    for(int i=0; i<num_local_cells_remote; i++)
        if(local_cells_global_index_remote[i] >= 0)
            logical_indx_lookup_table_remote[local_cells_global_index_remote[i]] = i;

    int last_local_logical_indx = -100;
    for (int j = 0; j < num_reference_cells; j ++) 
        if (reference_cell_indx[j] >= 0 && logical_indx_lookup_table_local[reference_cell_indx[j]] != -1 && logical_indx_lookup_table_remote[reference_cell_indx[j]] != -1) {
            if (last_local_logical_indx + 1 != logical_indx_lookup_table_local[reference_cell_indx[j]]) 
                routing_info_mgt[remote_proc_id].num_local_indx_segments ++;
            last_local_logical_indx = logical_indx_lookup_table_local[reference_cell_indx[j]];
            routing_info_mgt[remote_proc_id].num_elements_transferred ++;
        }

    last_local_logical_indx = -100;
    if (routing_info_mgt[remote_proc_id].num_elements_transferred > 0) {
        routing_info_mgt[remote_proc_id].local_indx_segment_starts = new int [routing_info_mgt[remote_proc_id].num_local_indx_segments];
        routing_info_mgt[remote_proc_id].local_indx_segment_lengths = new int [routing_info_mgt[remote_proc_id].num_local_indx_segments];
        routing_info_mgt[remote_proc_id].num_local_indx_segments = 0;
        for (int j = 0; j < num_reference_cells; j ++) 
            if (reference_cell_indx[j] >= 0 && logical_indx_lookup_table_local[reference_cell_indx[j]] != -1 && logical_indx_lookup_table_remote[reference_cell_indx[j]] != -1) {
                if (last_local_logical_indx + 1 != logical_indx_lookup_table_local[reference_cell_indx[j]]) {
                    routing_info_mgt[remote_proc_id].local_indx_segment_starts[routing_info_mgt[remote_proc_id].num_local_indx_segments] = logical_indx_lookup_table_local[reference_cell_indx[j]];
                    routing_info_mgt[remote_proc_id].local_indx_segment_lengths[routing_info_mgt[remote_proc_id].num_local_indx_segments] = 1;
                    routing_info_mgt[remote_proc_id].num_local_indx_segments ++;
                }
                else routing_info_mgt[remote_proc_id].local_indx_segment_lengths[routing_info_mgt[remote_proc_id].num_local_indx_segments - 1] ++;
                last_local_logical_indx = logical_indx_lookup_table_local[reference_cell_indx[j]];
            }
    }

    delete [] logical_indx_lookup_table_remote;
    delete [] logical_indx_lookup_table_local;
}

void Data_transfer_instance::compute_remap_routing_info_between_decomps(int num_local_cells_remote, int *local_cells_global_index_remote, int remote_proc_id)
{
    const int * reference_cell_indx;
    int num_reference_cells;
    
    routing_info_mgt[remote_proc_id+remote_model_size].num_elements_transferred = 0;
    routing_info_mgt[remote_proc_id+remote_model_size].num_local_indx_segments = 0;

    if(num_local_cells_remote < num_remap_local_cells){
        num_reference_cells = num_local_cells_remote;
        reference_cell_indx = local_cells_global_index_remote;
    }
    else{
        num_reference_cells = num_remap_local_cells;
        reference_cell_indx = local_remap_cells_global_index;
    }

    int * logical_indx_lookup_table_remote = new int[num_grid_cells];
    int * logical_indx_lookup_table_local = new int[num_grid_cells];

    for(int i=0; i<num_grid_cells; i++)
    {
        logical_indx_lookup_table_remote[i] = -1;
        logical_indx_lookup_table_local[i] = -1;
    }
    for(int i=0; i<num_remap_local_cells; i++)
        if(local_remap_cells_global_index[i] >= 0)
            logical_indx_lookup_table_local[local_remap_cells_global_index[i]] = i;
    for(int i=0; i<num_local_cells_remote; i++)
        if(local_cells_global_index_remote[i] >= 0)
            logical_indx_lookup_table_remote[local_cells_global_index_remote[i]] = i;

    int last_local_logical_indx = -100;
    for (int j = 0; j < num_reference_cells; j ++) 
        if (reference_cell_indx[j] >= 0 && logical_indx_lookup_table_local[reference_cell_indx[j]] != -1 && logical_indx_lookup_table_remote[reference_cell_indx[j]] != -1) {
            if (last_local_logical_indx + 1 != logical_indx_lookup_table_local[reference_cell_indx[j]]) 
                routing_info_mgt[remote_proc_id+remote_model_size].num_local_indx_segments ++;
            last_local_logical_indx = logical_indx_lookup_table_local[reference_cell_indx[j]];
            routing_info_mgt[remote_proc_id+remote_model_size].num_elements_transferred ++;
        }

    last_local_logical_indx = -100;
    if (routing_info_mgt[remote_proc_id+remote_model_size].num_elements_transferred > 0) {
        routing_info_mgt[remote_proc_id+remote_model_size].local_indx_segment_starts = new int [routing_info_mgt[remote_proc_id+remote_model_size].num_local_indx_segments];
        routing_info_mgt[remote_proc_id+remote_model_size].local_indx_segment_lengths = new int [routing_info_mgt[remote_proc_id+remote_model_size].num_local_indx_segments];
        routing_info_mgt[remote_proc_id+remote_model_size].num_local_indx_segments = 0;
        for (int j = 0; j < num_reference_cells; j ++) 
            if (reference_cell_indx[j] >= 0 && logical_indx_lookup_table_local[reference_cell_indx[j]] != -1 && logical_indx_lookup_table_remote[reference_cell_indx[j]] != -1) {
                if (last_local_logical_indx + 1 != logical_indx_lookup_table_local[reference_cell_indx[j]]) {
                    routing_info_mgt[remote_proc_id+remote_model_size].local_indx_segment_starts[routing_info_mgt[remote_proc_id+remote_model_size].num_local_indx_segments] = logical_indx_lookup_table_local[reference_cell_indx[j]];
                    routing_info_mgt[remote_proc_id+remote_model_size].local_indx_segment_lengths[routing_info_mgt[remote_proc_id+remote_model_size].num_local_indx_segments] = 1;
                    routing_info_mgt[remote_proc_id+remote_model_size].num_local_indx_segments ++;
                }
                else routing_info_mgt[remote_proc_id+remote_model_size].local_indx_segment_lengths[routing_info_mgt[remote_proc_id+remote_model_size].num_local_indx_segments - 1] ++;
                last_local_logical_indx = logical_indx_lookup_table_local[reference_cell_indx[j]];
            }
    }

    delete [] logical_indx_lookup_table_remote;
    delete [] logical_indx_lookup_table_local;
}

void Data_transfer_instance::build_2D_self_router()
{
    routing_info_mgt = new Routing_info_with_one_process[2*remote_model_size];
    int global_rank;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    
    int * num_remap_cells_each_remote_proc = new int[remote_model_size];
    int ** remap_cells_index_each_remote_proc = new int *[remote_model_size];
    int * remap_cells_displs = new int[remote_model_size];
    int * remap_cells_index_buffer;
    
    MPI_Allgather(&num_remap_local_cells, 1, MPI_INT, num_remap_cells_each_remote_proc, 1, MPI_INT, local_comm);
    remap_cells_displs[0] = 0;
    for(int i=1; i<remote_model_size; i++)
        remap_cells_displs[i] = remap_cells_displs[i-1] + num_remap_cells_each_remote_proc[i-1];
    remap_cells_index_buffer = new int[remap_cells_displs[remote_model_size-1] + num_remap_cells_each_remote_proc[remote_model_size-1]];
    MPI_Allgatherv(local_remap_cells_global_index, num_remap_local_cells, MPI_INT, remap_cells_index_buffer, num_remap_cells_each_remote_proc, remap_cells_displs, MPI_INT, local_comm);
    for(int i=0; i<remote_model_size; i++)
        remap_cells_index_each_remote_proc[i] = remap_cells_index_buffer + remap_cells_displs[i];

 
    if(num_local_cells > 0)
        for(int i=0; i<remote_model_size; i++)
            compute_routing_info_between_decomps(num_remap_cells_each_remote_proc[i], remap_cells_index_each_remote_proc[i], i);
    
    delete [] remap_cells_displs;
    delete [] remap_cells_index_buffer;
    delete [] remap_cells_index_each_remote_proc;
    delete [] num_remap_cells_each_remote_proc;
    
    int * num_cells_each_remote_proc = new int[remote_model_size];
    int ** cells_index_each_remote_proc = new int *[remote_model_size];
    int * local_cells_displs = new int[remote_model_size];
    int * local_cells_index_buffer;
    
    MPI_Allgather(&num_local_cells, 1, MPI_INT, num_cells_each_remote_proc, 1, MPI_INT, local_comm);
    local_cells_displs[0] = 0;
    for(int i=1; i<remote_model_size; i++)
        local_cells_displs[i] = local_cells_displs[i-1] + num_cells_each_remote_proc[i-1];
    local_cells_index_buffer = new int[local_cells_displs[remote_model_size-1] + num_cells_each_remote_proc[remote_model_size-1]];
    MPI_Allgatherv(local_cells_global_index, num_local_cells, MPI_INT, local_cells_index_buffer, num_cells_each_remote_proc, local_cells_displs, MPI_INT, local_comm);
    for(int i=0; i<remote_model_size; i++)
        cells_index_each_remote_proc[i] = local_cells_index_buffer + local_cells_displs[i];
      
    if(num_remap_local_cells > 0)
        for(int i=0; i<remote_model_size; i++)
            compute_remap_routing_info_between_decomps(num_cells_each_remote_proc[i], cells_index_each_remote_proc[i], i);
    
    delete [] local_cells_displs;
    delete [] local_cells_index_buffer;
    delete [] cells_index_each_remote_proc;
    delete [] num_cells_each_remote_proc;
}

bool Data_transfer_instance::exec()
{
    pack_MD_data();
    
    int fields_per_cell = 0;
    if(direction != RECV)
    {
        for(int i=0; i<coupling_fields.size(); i++)
            if(coupling_fields[i]->input && (mask == NULL || mask[i]))
                fields_per_cell += coupling_fields[i]->buf_size / num_local_cells * coupling_fields[i]->data_type;
    }
    else
    {
        for(int i=0; i<coupling_fields.size(); i++)
            if(!coupling_fields[i]->input && (mask == NULL || mask[i]))
                fields_per_cell += coupling_fields[i]->buf_size / num_local_cells * coupling_fields[i]->data_type;
    }

    butterfly_instance->execute(send_data_buf, recv_data_buf, fields_per_cell);
    unpack_MD_data();
    return true;
}

bool Data_transfer_instance::auto_set_p2p_stage_num()
{
    //pack_MD_data();
    
    int fields_per_cell = 0;
    if(direction != SEND)
    {
        for(int i=0; i<coupling_fields.size(); i++)
            if(coupling_fields[i]->input)
                fields_per_cell += coupling_fields[i]->buf_size / num_local_cells * coupling_fields[i]->data_type;
    }
    else
    {
        for(int i=0; i<coupling_fields.size(); i++)
            if(!coupling_fields[i]->input)
                fields_per_cell += coupling_fields[i]->buf_size / num_local_cells * coupling_fields[i]->data_type;
    }

    butterfly_instance->auto_set_p2p_stage_num(send_data_buf, recv_data_buf, fields_per_cell);
    //unpack_MD_data();
    return true;
}

void Data_transfer_instance::pack_MD_data()
{
    int global_rank;
    
    int offset = 0;
    for(int i=0; i<remote_model_size; i++)
        if(routing_info_mgt[i].num_elements_transferred > 0)
        {
            int num_segments = routing_info_mgt[i].num_local_indx_segments;
            for(int j=0; j<coupling_fields.size(); j++)
            {
                if(coupling_fields[j]->input && (mask == NULL || mask[j]))
                {
                    int num_lev = coupling_fields[j]->buf_size / num_local_cells;
                    
                    for(int k=0; k<num_segments; k++)
                    {
                        switch(coupling_fields[j]->data_type){
                            case 1:
                                pack_segment_data((char*)((char*)send_data_buf+offset), (char*) coupling_fields[j]->data_buf, routing_info_mgt[i].local_indx_segment_starts[k], routing_info_mgt[i].local_indx_segment_lengths[k], num_local_cells, num_lev);
                                break;
                            case 2:
                                pack_segment_data((short*)((char*)send_data_buf+offset), (short*) coupling_fields[j]->data_buf, routing_info_mgt[i].local_indx_segment_starts[k], routing_info_mgt[i].local_indx_segment_lengths[k], num_local_cells, num_lev);
                                break;
                            case 4:
                                pack_segment_data((int*)((char*)send_data_buf+offset), (int*) coupling_fields[j]->data_buf, routing_info_mgt[i].local_indx_segment_starts[k], routing_info_mgt[i].local_indx_segment_lengths[k], num_local_cells, num_lev);
                                break;
                            case 8:
                                pack_segment_data((double*)((char*)send_data_buf+offset), (double*) coupling_fields[j]->data_buf, routing_info_mgt[i].local_indx_segment_starts[k], routing_info_mgt[i].local_indx_segment_lengths[k], num_local_cells, num_lev);
                                break;
                            default:
                                break;
                        }
                        
                        offset += num_lev * routing_info_mgt[i].local_indx_segment_lengths[k] * coupling_fields[j]->data_type;
                    }
                }
            }
        }
}

void Data_transfer_instance::unpack_MD_data()
{
    int offset = 0;
    int tmp_local_cells;
    if(direction == SENDRECV) tmp_local_cells = num_remap_local_cells;
    else tmp_local_cells = num_local_cells;
    int tmp_routing_index = 0;
    if(direction == SENDRECV) tmp_routing_index = remote_model_size;


    for(int i=0; i<remote_model_size; i++)
        if(routing_info_mgt[i+tmp_routing_index].num_elements_transferred > 0)
        {
            int num_segments = routing_info_mgt[i+tmp_routing_index].num_local_indx_segments;
            for(int j=0; j<coupling_fields.size(); j++)
            {
                if(!coupling_fields[j]->input && (mask == NULL || mask[j]))
                {
                    int num_lev = coupling_fields[j]->buf_size / tmp_local_cells;
                    for(int k=0; k<num_segments; k++)
                    {
                        switch(coupling_fields[j]->data_type){
                            case 1:
                                unpack_segment_data((char*)((char*)recv_data_buf+offset), (char*) coupling_fields[j]->data_buf, routing_info_mgt[i+tmp_routing_index].local_indx_segment_starts[k], routing_info_mgt[i+tmp_routing_index].local_indx_segment_lengths[k], tmp_local_cells, num_lev);
                                break;
                            case 2:
                                unpack_segment_data((short*)((char*)recv_data_buf+offset), (short*) coupling_fields[j]->data_buf, routing_info_mgt[i+tmp_routing_index].local_indx_segment_starts[k], routing_info_mgt[i+tmp_routing_index].local_indx_segment_lengths[k], tmp_local_cells, num_lev);
                                break;
                            case 4:
                                unpack_segment_data((int*)((char*)recv_data_buf+offset), (int*) coupling_fields[j]->data_buf, routing_info_mgt[i+tmp_routing_index].local_indx_segment_starts[k], routing_info_mgt[i+tmp_routing_index].local_indx_segment_lengths[k], tmp_local_cells, num_lev);
                                break;
                            case 8:
                                unpack_segment_data((double*)((char*)recv_data_buf+offset), (double*) coupling_fields[j]->data_buf, routing_info_mgt[i+tmp_routing_index].local_indx_segment_starts[k], routing_info_mgt[i+tmp_routing_index].local_indx_segment_lengths[k], tmp_local_cells, num_lev);
                                break;
                            default:
                                break;
                        }
                        offset += num_lev * routing_info_mgt[i+tmp_routing_index].local_indx_segment_lengths[k] * coupling_fields[j]->data_type;
                    }
                }
            }
        }
}

Data_transfer_instance_mgt::Data_transfer_instance_mgt()
{
}

Data_transfer_instance_mgt::~Data_transfer_instance_mgt()
{
    for(int i=0; i<data_transfer_instance_mgt.size(); i++)
        if(data_transfer_instance_mgt[i] != NULL)
            delete data_transfer_instance_mgt[i];
}

int Data_transfer_instance_mgt::register_data_transfer_instance(MPI_Comm local_comm, int remote_root_global_rank, int direction)
{
    Data_transfer_instance * data_transfer_instance = new Data_transfer_instance(local_comm, remote_root_global_rank, direction);
    data_transfer_instance_mgt.push_back(data_transfer_instance);
    return data_transfer_instance_mgt.size()-1;
}

bool Data_transfer_instance_mgt::register_decomposition(int instance_id, int num_grid_cells, int num_local_cells, int * local_cells_global_index)
{
    if(instance_id >= data_transfer_instance_mgt.size()) return false;
    return data_transfer_instance_mgt[instance_id]->register_decomposition(num_grid_cells, num_local_cells, local_cells_global_index);
}

bool Data_transfer_instance_mgt::register_field(int instance_id, void * data_buf, int buf_size, char * data_type, bool input)
{
    if(instance_id >= data_transfer_instance_mgt.size()) return false;
    return data_transfer_instance_mgt[instance_id]->register_field(data_buf, buf_size, data_type, input);
}

bool Data_transfer_instance_mgt::init_data_transfer_instance(int instance_id)
{
    if(instance_id >= data_transfer_instance_mgt.size()) return false;
    else{ 
        data_transfer_instance_mgt[instance_id]->init();
        //data_transfer_instance_mgt[instance_id]->auto_set_p2p_stage_num();
        return true;
    }
}

bool Data_transfer_instance_mgt::exec_data_transfer_instance(int instance_id)
{
    if(instance_id >= data_transfer_instance_mgt.size()) return false;
    return data_transfer_instance_mgt[instance_id]->exec();
}

bool Data_transfer_instance_mgt::auto_set_p2p_stage_num(int instance_id)
{
    if(instance_id >= data_transfer_instance_mgt.size()) return false;
    return data_transfer_instance_mgt[instance_id]->auto_set_p2p_stage_num();
}

bool Data_transfer_instance_mgt::final_data_transfer_instance(int instance_id)
{
    if(instance_id >= data_transfer_instance_mgt.size()) return false;
    delete data_transfer_instance_mgt[instance_id];
    data_transfer_instance_mgt[instance_id] = NULL;
    return true;
}

bool Data_transfer_instance_mgt::register_mask(int instance_id, bool * mask)
{
    if(instance_id >= data_transfer_instance_mgt.size()) return false;
    return data_transfer_instance_mgt[instance_id]->register_mask(mask);
}
