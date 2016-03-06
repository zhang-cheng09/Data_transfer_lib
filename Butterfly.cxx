
/***************************************************************************
 * This is a source file of the Adaptive Data Transfer Library version 1.0
 * This file was initially finished by Cheng Zhang
 * If you have any problem,
 * please contact Cheng Zhang via zhang-cheng09@mails.tsinghua.edu.cn
 **************************************************************************/

#include <mpi.h>
#include "Butterfly.h"
#include <cstring>
#include <cmath>
#include <vector>
#include <climits>
#include <stdio.h>
#include <sys/time.h>
#include <iostream>
using namespace std;

void Butterfly::wtime(double * t)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    *t = tv.tv_sec + 1.0e-6*tv.tv_usec;
}

Butterfly::Butterfly(MPI_Comm global_comm, MPI_Comm local_comm, int * remote_model_global_rank, int * send_cells, int * recv_cells, int remote_model_size, int action)
{
    butterfly_node_num = -1;
    butterfly_grp_node_num = -1;
    butterfly_stage_num = -1;
    last_p2p_stage_num = 0;
    max_cells_cnts = 0;
    pre_fields_per_cell = -1;
    intra_node_gather_init(global_comm, local_comm, remote_model_global_rank, send_cells, recv_cells, remote_model_size, action);
    process_mapping_from_send_to_butterfly_init();
    process_mapping_from_butterfly_to_recv_init();
    butterfly_init();
    butterfly_p2p_init();
    first_p2p_init();
    last_p2p_init();
}

Butterfly::~Butterfly()
{
    if(send_cells != NULL)
      delete [] send_cells;
    if(recv_cells != NULL)
      delete [] recv_cells;

    delete [] send_model_global_rank;
    delete [] recv_model_global_rank;

    delete [] send_node_per_proc_local_rank;
    delete [] recv_node_per_proc_local_rank;
    delete [] total_node_root_proc_global_rank;

    delete [] total_send_cells;
    delete [] total_recv_cells;

    delete [] send_node_cells_cnts;
    delete [] butterfly_node_global_rank;
    delete [] butterfly_node_map_to_send_node;

    delete [] butterfly_node_map_to_recv_node;
    delete [] recv_node_map_to_butterfly_node;
    delete [] recv_proc_node_index;

    if(butterfly_data_buf[0] != NULL)
    {
        delete butterfly_data_buf[0];
        delete butterfly_data_buf[1];
        delete butterfly_data_buf[2];
    }
    
    if(butterfly_node_master)
    {
        delete [] butterfly_node_decompsitions;
        delete [] butterfly_node_decompsitions_cells;
        delete [] butterfly_pair_node_list;
    }

    if(intra_proc_displs_cells != NULL)
      delete [] intra_proc_displs_cells;
    
    if(last_p2p_send_model_proc_lists != NULL)
      delete [] last_p2p_send_model_proc_lists;
    if(last_p2p_send_cells != NULL)
      delete [] last_p2p_send_cells;
    if(last_p2p_recv_cells != NULL)
      delete [] last_p2p_recv_cells;
    if(last_p2p_send_displs != NULL)
      delete [] last_p2p_send_displs;
    if(last_p2p_recv_buf != NULL)
      delete [] last_p2p_recv_buf;

    if(butterfly_p2p_remote_procs != NULL)
      delete [] butterfly_p2p_remote_procs;
    if(butterfly_p2p_send_cells != NULL)
      delete [] butterfly_p2p_send_cells;
    if(butterfly_p2p_recv_cells != NULL)
      delete [] butterfly_p2p_recv_cells;
    if(butterfly_p2p_send_displs != NULL)
      delete [] butterfly_p2p_send_displs;
    if(butterfly_p2p_recv_displs != NULL)
      delete [] butterfly_p2p_recv_displs;
    if(butterfly_p2p_send_mask != NULL)
      delete [] butterfly_p2p_send_mask;
    if(butterfly_p2p_send_mask_cells != NULL)
      delete [] butterfly_p2p_send_mask_cells;

    if(p2p_num_per_butterfly_stage != NULL)
      delete [] p2p_num_per_butterfly_stage;
    if(p2p_depth_per_butterfly_stage != NULL)
      delete [] p2p_depth_per_butterfly_stage;
    
    if(first_p2p_sendcells != NULL) delete [] first_p2p_sendcells;
    if(first_p2p_recvcells != NULL) delete [] first_p2p_recvcells;
    if(first_p2p_senddispls != NULL) delete [] first_p2p_senddispls;
    if(first_p2p_recvdispls != NULL) delete [] first_p2p_recvdispls;
    if(first_p2p_send_mask != NULL) delete [] first_p2p_send_mask;
    if(first_p2p_send_mask_cells != NULL) delete [] first_p2p_send_mask_cells;
    if(first_p2p_send_proc_index != NULL) delete [] first_p2p_send_proc_index;
    if(first_p2p_recv_proc_index != NULL) delete [] first_p2p_recv_proc_index;

    delete [] stage_mask;
}

void Butterfly::reset_p2p_stage_num(int * stage_num)
{
    if(last_p2p_send_model_proc_lists != NULL)
      delete [] last_p2p_send_model_proc_lists;
    if(last_p2p_send_cells != NULL)
      delete [] last_p2p_send_cells;
    if(last_p2p_recv_cells != NULL)
      delete [] last_p2p_recv_cells;
    if(last_p2p_send_displs != NULL)
      delete [] last_p2p_send_displs;
    if(last_p2p_recv_buf != NULL)
      delete [] last_p2p_recv_buf;

    if(butterfly_p2p_remote_procs != NULL)
      delete [] butterfly_p2p_remote_procs;
    if(butterfly_p2p_send_cells != NULL)
      delete [] butterfly_p2p_send_cells;
    if(butterfly_p2p_recv_cells != NULL)
      delete [] butterfly_p2p_recv_cells;
    if(butterfly_p2p_send_displs != NULL)
      delete [] butterfly_p2p_send_displs;
    if(butterfly_p2p_recv_displs != NULL)
      delete [] butterfly_p2p_recv_displs;
    if(butterfly_p2p_send_mask != NULL)
      delete [] butterfly_p2p_send_mask;
    if(butterfly_p2p_send_mask_cells != NULL)
      delete [] butterfly_p2p_send_mask_cells;

    if(p2p_num_per_butterfly_stage != NULL)
      delete [] p2p_num_per_butterfly_stage;
    if(p2p_depth_per_butterfly_stage != NULL)
      delete [] p2p_depth_per_butterfly_stage;
    
    if(first_p2p_sendcells != NULL) delete [] first_p2p_sendcells;
    if(first_p2p_recvcells != NULL) delete [] first_p2p_recvcells;
    if(first_p2p_senddispls != NULL) delete [] first_p2p_senddispls;
    if(first_p2p_recvdispls != NULL) delete [] first_p2p_recvdispls;
    if(first_p2p_send_mask != NULL) delete [] first_p2p_send_mask;
    if(first_p2p_send_mask_cells != NULL) delete [] first_p2p_send_mask_cells;
    if(first_p2p_send_proc_index != NULL) delete [] first_p2p_send_proc_index;
    if(first_p2p_recv_proc_index != NULL) delete [] first_p2p_recv_proc_index;

    butterfly_p2p_init(stage_num);
    first_p2p_init();
    last_p2p_send_to_remote_procs.clear();
    last_p2p_recv_from_remote_procs.clear();
    last_p2p_send_mask.clear();
    last_p2p_send_mask_cells.clear();
    last_p2p_init();
}

void Butterfly::first_p2p_init()
{
    butterfly_data_buf[0] = NULL;
    butterfly_data_buf[1] = NULL;
    butterfly_data_buf[2] = NULL;

    first_p2p_sendcells = NULL;
    first_p2p_recvcells = NULL;
    first_p2p_senddispls = NULL;
    first_p2p_recvdispls = NULL;
    first_p2p_send_proc_num;
    first_p2p_recv_proc_num;
    first_p2p_send_mask = NULL;
    first_p2p_send_mask_cells = NULL;
    first_p2p_send_proc_index = NULL;
    first_p2p_recv_proc_index = NULL;
    
    if(butterfly_p2p_stage_num == 0)
    {
        first_p2p_send_proc_num = recv_model_size;
        first_p2p_recv_proc_num = send_model_size;
        if(is_send)
        {
            first_p2p_sendcells = new int[first_p2p_send_proc_num];
            first_p2p_senddispls = new int[first_p2p_send_proc_num];
            first_p2p_send_proc_index = new int[first_p2p_send_proc_num];
            memcpy(first_p2p_sendcells, total_send_cells+local_rank*recv_model_size, sizeof(int)*recv_model_size);
            memcpy(first_p2p_send_proc_index, recv_model_global_rank, sizeof(int)*recv_model_size);

            first_p2p_senddispls[0] = 0;
            for(int i=0; i<recv_model_size; i++)
              first_p2p_senddispls[i] = first_p2p_senddispls[i-1] + first_p2p_sendcells[i-1];
        }

        if(is_recv)
        {
            first_p2p_recvcells = new int[first_p2p_recv_proc_num];
            first_p2p_recvdispls = new int[first_p2p_recv_proc_num];
            first_p2p_recv_proc_index = new int[first_p2p_recv_proc_num];
            memcpy(first_p2p_recvcells, total_recv_cells+local_rank*send_model_size, sizeof(int)*send_model_size);
            memcpy(first_p2p_recv_proc_index, send_model_global_rank, sizeof(int)*send_model_size);

            first_p2p_recvdispls[0] = 0;
            for(int i=0; i<send_model_size; i++)
              first_p2p_recvdispls[i] = first_p2p_recvdispls[i-1] + first_p2p_recvcells[i-1];
        }
    }
    else
    {
        int tmp_p2p_depth = p2p_depth_per_butterfly_stage[0];
        int * local_data_decompositions = new int[butterfly_node_num];
        int * process_data_decompositions = new int[recv_model_size];
        first_p2p_send_proc_num = tmp_p2p_depth;
        first_p2p_recv_proc_num = 0;
        first_p2p_sendcells = new int[first_p2p_send_proc_num];
        first_p2p_recvcells = NULL;
        first_p2p_senddispls = new int[first_p2p_send_proc_num];
        first_p2p_recvdispls = NULL;
        first_p2p_send_mask = new int[recv_model_size];
        first_p2p_send_mask_cells = new int[recv_model_size];
        first_p2p_send_proc_index = new int[first_p2p_send_proc_num];
        first_p2p_recv_proc_index = NULL;

        memset(first_p2p_sendcells, 0, sizeof(int)*tmp_p2p_depth);
        for(int i=0; i<recv_model_size; i++)
        {
            int tmp_node = recv_node_map_to_butterfly_node[i][0];
            process_data_decompositions[i] = tmp_node;
        }

        if(is_send)
        {
            int tmp_butterfly_node = send_node_map_to_butterfly_node_index;

            for(int i=0; i<tmp_p2p_depth; i++)
            {
                int tmp_remote_node = tmp_butterfly_node / tmp_p2p_depth * tmp_p2p_depth + i;
                for(int j=0; j<butterfly_node_num/tmp_p2p_depth; j++)
                {
                    local_data_decompositions[i*butterfly_node_num/tmp_p2p_depth+j] = tmp_remote_node % tmp_p2p_depth + j * tmp_p2p_depth;
                }
            }

            for(int i=0; i<recv_model_size; i++)
            {
                int tmp_mask = -1;
                for(int j=0; j<butterfly_node_num; j++)
                  if(local_data_decompositions[j] == process_data_decompositions[i])
                  {
                      tmp_mask = j/(butterfly_node_num/tmp_p2p_depth);
                      break;
                  }
                if(tmp_mask != -1){
                    first_p2p_send_mask[i] = tmp_mask;
                    first_p2p_send_mask_cells[i] = total_send_cells[local_rank*recv_model_size+i];
                    first_p2p_sendcells[tmp_mask] += total_send_cells[local_rank*recv_model_size+i];
                }
            }

            for(int i=0; i<tmp_p2p_depth; i++)
            {
                int tmp_remote_node = tmp_butterfly_node / tmp_p2p_depth * tmp_p2p_depth + i;
                first_p2p_send_proc_index[i] = butterfly_node_global_rank[tmp_remote_node];
            }

            first_p2p_senddispls[0] = 0;
            for(int i=1; i<tmp_p2p_depth; i++)
              first_p2p_senddispls[i] = first_p2p_senddispls[i-1] + first_p2p_sendcells[i-1];
            
            int tmp_buf_size = first_p2p_senddispls[tmp_p2p_depth-1] + first_p2p_sendcells[tmp_p2p_depth-1];
            if(max_cells_cnts < tmp_buf_size)
              max_cells_cnts = tmp_buf_size;
        }
      
        if(butterfly_node_master)
        {
            for(int i=0; i<tmp_p2p_depth; i++)
            {
                int tmp_remote_node = butterfly_pair_node_list[i];
                first_p2p_recv_proc_num += butterfly_node_map_to_send_node[tmp_remote_node].size();
            }

            first_p2p_recvcells = new int[first_p2p_recv_proc_num];
            first_p2p_recvdispls = new int[first_p2p_recv_proc_num];
            first_p2p_recv_proc_index = new int[first_p2p_recv_proc_num];
            memset(first_p2p_recvcells, 0, sizeof(int)*first_p2p_recv_proc_num);
            
            for(int i=0; i<butterfly_node_num/tmp_p2p_depth; i++)
            {
                int tmp_data_decomposition = butterfly_node_index % tmp_p2p_depth + i * tmp_p2p_depth;
                for(int j=0; j<recv_model_size; j++)
                {
                    if(tmp_data_decomposition == process_data_decompositions[j])
                    {
                        int tmp_proc_index = 0;
                        for(int m=0; m<tmp_p2p_depth; m++)
                        {
                            int tmp_remote_node = butterfly_pair_node_list[m];
                            for(int n=0; n<butterfly_node_map_to_send_node[tmp_remote_node].size(); n++)
                            {
                                int tmp_send_node = butterfly_node_map_to_send_node[tmp_remote_node][n];
                                int tmp_local_rank = send_node_per_proc_local_rank[tmp_send_node][0];
                                first_p2p_recvcells[tmp_proc_index] += total_send_cells[tmp_local_rank*recv_model_size+j];
                                tmp_proc_index ++;
                            }
                        }
                    }
                }
            }

            int tmp_proc_index = 0;
            for(int i=0; i<tmp_p2p_depth; i++)
            {
                int tmp_remote_node = butterfly_pair_node_list[i];
                for(int j=0; j<butterfly_node_map_to_send_node[tmp_remote_node].size(); j++)
                {
                    int tmp_send_node = butterfly_node_map_to_send_node[tmp_remote_node][j];
                    int tmp_local_rank = send_node_per_proc_local_rank[tmp_send_node][0];
                    int tmp_global_rank = send_model_global_rank[tmp_local_rank];
                    first_p2p_recv_proc_index[tmp_proc_index] = tmp_global_rank;
                    tmp_proc_index ++;
                }
            }

            first_p2p_recvdispls[0] = 0;
            for(int i=1; i<first_p2p_recv_proc_num; i++)
              first_p2p_recvdispls[i] = first_p2p_recvdispls[i-1] + first_p2p_recvcells[i-1];
            
            int tmp_buf_size = first_p2p_recvdispls[first_p2p_recv_proc_num-1] + first_p2p_recvcells[first_p2p_recv_proc_num-1];
            if(max_cells_cnts < tmp_buf_size)
              max_cells_cnts = tmp_buf_size;
        }

        delete [] local_data_decompositions;
        delete [] process_data_decompositions;
    }
}

void Butterfly::first_p2p_execute(char * send_data, char * recv_data)
{
    MPI_Request * send_req, * recv_req;
    MPI_Status status;

    if(butterfly_p2p_stage_num == 0)
    {
        if(is_send){
            send_req = new MPI_Request[recv_model_size];
            for(int i=0; i<recv_model_size; i++)
              if(first_p2p_sendcells[i] > 0)
                MPI_Isend(send_data+first_p2p_senddispls[i]*fields_per_cell, first_p2p_sendcells[i]*fields_per_cell, MPI_CHAR, recv_model_global_rank[i], 1300+global_rank, global_comm, send_req+i);
        }

        if(is_recv){
            recv_req = new MPI_Request[send_model_size];
            for(int i=0; i<send_model_size; i++)
              if(first_p2p_recvcells[i] > 0)
                MPI_Irecv(recv_data+first_p2p_recvdispls[i]*fields_per_cell, first_p2p_recvcells[i]*fields_per_cell, MPI_CHAR, send_model_global_rank[i], 1300+send_model_global_rank[i], global_comm, recv_req+i);
        }

        if(is_send){
            for(int i=0; i<recv_model_size; i++)
              if(first_p2p_sendcells[i] > 0)
                MPI_Wait(send_req+i, &status);

            delete [] send_req;
        }

        if(is_recv){
            for(int i=0; i<send_model_size; i++)
              if(first_p2p_recvcells[i] > 0)
                MPI_Wait(recv_req+i, &status);

            delete [] recv_req;
        }

        return;
    }

    if(max_cells_cnts > 0 && (pre_fields_per_cell != fields_per_cell || butterfly_data_buf[0] == NULL))
    {
        if(butterfly_data_buf[0] != NULL)
        {
            delete [] butterfly_data_buf[0];
            delete [] butterfly_data_buf[1];
            delete [] butterfly_data_buf[2];
        }

        butterfly_data_buf[0] = new char[max_cells_cnts * fields_per_cell];
        butterfly_data_buf[1] = new char[max_cells_cnts * fields_per_cell];
        butterfly_data_buf[2] = new char[max_cells_cnts * fields_per_cell];
    }

    char * intra_node_recvbuf = butterfly_data_buf[0];
    char * intra_node_sendbuf = butterfly_data_buf[1];

    if(is_send)
    {
        send_req = new MPI_Request[first_p2p_send_proc_num];
        int tmp_send_displs[first_p2p_send_proc_num];

        for(int i=0; i<first_p2p_send_proc_num; i++)
            tmp_send_displs[i] = first_p2p_senddispls[i] * fields_per_cell;

        int tmp_displs = 0;
        for(int i=0; i<recv_model_size; i++)
            if(first_p2p_send_mask_cells[i] > 0)
            {
                int tmp_mask = first_p2p_send_mask[i];
                memcpy(intra_node_sendbuf+tmp_send_displs[tmp_mask], send_data+tmp_displs, sizeof(char)*first_p2p_send_mask_cells[i]*fields_per_cell);
                tmp_send_displs[tmp_mask] += first_p2p_send_mask_cells[i]*fields_per_cell;
                tmp_displs += first_p2p_send_mask_cells[i]*fields_per_cell;
            }

        for(int i=0; i<first_p2p_send_proc_num; i++)
          if(first_p2p_sendcells[i] > 0)
            MPI_Isend(intra_node_sendbuf+first_p2p_senddispls[i]*fields_per_cell, first_p2p_sendcells[i]*fields_per_cell, MPI_CHAR, first_p2p_send_proc_index[i], 1300+global_rank, global_comm, send_req+i);
    }

    if(butterfly_node_master)
    {
        recv_req = new MPI_Request[first_p2p_recv_proc_num];
        for(int i=0; i<first_p2p_recv_proc_num; i++)
          if(first_p2p_recvcells[i] > 0)
            MPI_Irecv(intra_node_recvbuf+first_p2p_recvdispls[i]*fields_per_cell, first_p2p_recvcells[i]*fields_per_cell, MPI_CHAR, first_p2p_recv_proc_index[i], 1300+first_p2p_recv_proc_index[i], global_comm, recv_req+i);
    }

    if(is_send)
    {
        for(int i=0; i<first_p2p_send_proc_num; i++)
          if(first_p2p_sendcells[i] > 0)
            MPI_Wait(send_req+i, &status);
        
        delete [] send_req;
    }

    if(butterfly_node_master)
    {
        for(int i=0; i<first_p2p_recv_proc_num; i++)
          if(first_p2p_recvcells[i] > 0)
            MPI_Wait(recv_req+i, &status);
        
        delete [] recv_req;
    }
}

void Butterfly::last_p2p_init()
{
    last_p2p_send_model_proc_lists = NULL;
    last_p2p_send_cells = NULL;
    last_p2p_send_displs = NULL;
    last_p2p_recv_cells = NULL;
    last_p2p_recv_buf = NULL;
    last_p2p_recv_buf_cells = 0;

    if(last_p2p_stage_num == 0) return;

    MPI_Request * send_req, * recv_req;
    MPI_Status status;

    int * local_data_decompositions = NULL;
    int * tmp_send_model_list, * all_send_model_list;
    int * global_rank_intra_comm;
    tmp_send_model_list = new int[send_model_size];
    if(action == SENDRECV) {
        all_send_model_list = new int[send_model_size*send_model_size];
        global_rank_intra_comm = new int[send_model_size];
    }
    else {
        all_send_model_list = new int[send_model_size * (send_model_size + recv_model_size)];
        global_rank_intra_comm = new int[send_model_size + recv_model_size];
    }
      
    if(butterfly_node_master && (butterfly_grp_index == 0))
    {
        int tmp_send_model_size = (int)exp2f(butterfly_stage_num - last_p2p_stage_num);
        for(int i=0; i<send_model_size; i++)
            tmp_send_model_list[i] = MPI_PROC_NULL;

        int tmp_index = 0;
        for(int i=0; i<butterfly_grp_num; i++)
        {
            for(int j=0; j<tmp_send_model_size; j++)
            {
                int tmp_butterfly_node = i * butterfly_grp_node_num + butterfly_pair_node_list[j];
                for(int m=0; m<butterfly_node_map_to_send_node[tmp_butterfly_node].size(); m++)
                {
                    int tmp_send_node = butterfly_node_map_to_send_node[tmp_butterfly_node][m];
                    for(int n=0; n<send_node_per_proc_local_rank[tmp_send_node].size(); n++)
                    {
                        int tmp_send_proc = send_node_per_proc_local_rank[tmp_send_node][n];
                        tmp_send_model_list[tmp_index++] = tmp_send_proc;
                    }
                }
            }
        }
        
        tmp_send_model_size = (int) exp2f(last_p2p_stage_num);
        local_data_decompositions = new int[tmp_send_model_size];
        int stage_size = butterfly_grp_node_num / tmp_send_model_size;
        for(int i=0; i<tmp_send_model_size; i++){
            local_data_decompositions[i] = butterfly_grp_node_index % stage_size + i * stage_size;
            int tmp_butterfly_node = local_data_decompositions[i];
            for(int j=0; j<butterfly_node_map_to_recv_node[tmp_butterfly_node].size(); j++)
            {
                int tmp_recv_node = butterfly_node_map_to_recv_node[tmp_butterfly_node][j];
                for(int k=0; k<recv_node_per_proc_local_rank[tmp_recv_node].size(); k++)
                {
                    int tmp_recv_proc = recv_node_per_proc_local_rank[tmp_recv_node][0];
                    last_p2p_send_to_remote_procs.push_back(recv_model_global_rank[tmp_recv_proc]);
                }
            }
        }

        //send_req = new MPI_Request[last_p2p_send_to_remote_procs.size()];
        //for(int i=0; i<last_p2p_send_to_remote_procs.size(); i++)
        //    MPI_Isend(tmp_send_model_list, send_model_size, MPI_INT, last_p2p_send_to_remote_procs[i], 3200+last_p2p_send_to_remote_procs[i], global_comm, send_req+i);
        
        last_p2p_send_cells = new int[last_p2p_send_to_remote_procs.size()];
        last_p2p_send_displs = new int[last_p2p_send_to_remote_procs.size()];
        int * p2p_recv_proc_local_rank = new int[last_p2p_send_to_remote_procs.size()];
        memset(last_p2p_send_cells, 0, sizeof(int)*last_p2p_send_to_remote_procs.size());
        tmp_index = 0;
        int tmp_displs = 0;
        for(int i=0; i<tmp_send_model_size; i++){
            int tmp_butterfly_node = local_data_decompositions[i];
            for(int j=0; j<butterfly_node_map_to_recv_node[tmp_butterfly_node].size(); j++)
            {
                int tmp_recv_node = butterfly_node_map_to_recv_node[tmp_butterfly_node][j];
                for(int k=0; k<recv_node_per_proc_local_rank[tmp_recv_node].size(); k++)
                {
                    int tmp_recv_proc = recv_node_per_proc_local_rank[tmp_recv_node][0];
                    for(int m=0; m<send_model_size; m++)
                    {
                        int tmp_send_proc = tmp_send_model_list[m];
                        if(tmp_send_proc != MPI_PROC_NULL)
                            last_p2p_send_cells[tmp_index] += total_recv_cells[tmp_recv_proc*send_model_size+tmp_send_proc];
                    }
                    p2p_recv_proc_local_rank[tmp_index] = tmp_recv_proc;
                    last_p2p_send_displs[tmp_index] = tmp_displs;
                    tmp_displs += last_p2p_send_cells[tmp_index];
                    tmp_index ++;
                }
            }
        }
        if(max_cells_cnts < tmp_displs) max_cells_cnts = tmp_displs;

        for(int i=0; i<send_model_size; i++)
        {
            int tmp_send_proc = tmp_send_model_list[i];
            if(tmp_send_proc != MPI_PROC_NULL)
            {
                for(int j=0; j<recv_model_size; j++)
                {
                    bool tmp_exist = false;
                    int k;
                    for(k=0; k<last_p2p_send_to_remote_procs.size(); k++)
                    {
                        if(p2p_recv_proc_local_rank[k] == j)
                        {
                            tmp_exist = true;
                            break;
                        }
                    }

                    if(tmp_exist){
                        last_p2p_send_mask.push_back(k);
                        last_p2p_send_mask_cells.push_back(total_send_cells[tmp_send_proc*recv_model_size+j]);
                    }
                }
            }
        }

        delete [] p2p_recv_proc_local_rank;
    }

    MPI_Allgather(&global_rank, 1, MPI_INT, global_rank_intra_comm, 1, MPI_INT, intra_comm);
    MPI_Allgather(tmp_send_model_list, send_model_size, MPI_INT, all_send_model_list, send_model_size, MPI_INT, intra_comm);

    if(is_recv)
    {
        int tmp_send_model_size = (int)exp2f(last_p2p_stage_num);
        int cur_recv_node = recv_proc_node_index[local_rank];
        for(int m=0; m<recv_node_map_to_butterfly_node[cur_recv_node].size(); m++){
            int cur_butterfly_node = recv_node_map_to_butterfly_node[cur_recv_node][m];
            int stage_size = butterfly_grp_node_num / tmp_send_model_size;
            for(int i=0; i<tmp_send_model_size; i++)
            {
                int tmp_butterfly_node = cur_butterfly_node % stage_size + i * stage_size;
                int tmp_send_proc = butterfly_node_global_rank[tmp_butterfly_node];
                last_p2p_recv_from_remote_procs.push_back(tmp_send_proc);
            }
        }
        last_p2p_send_model_proc_lists = new int[send_model_size * last_p2p_recv_from_remote_procs.size()];

        //recv_req = new MPI_Request[last_p2p_recv_from_remote_procs.size()];
        //for(int i=0; i<last_p2p_recv_from_remote_procs.size(); i++)
        //    MPI_Irecv(last_p2p_send_model_proc_lists+i*send_model_size, send_model_size, MPI_INT, last_p2p_recv_from_remote_procs[i], 3200+global_rank, global_comm, recv_req+i);
        
        for(int i=0; i<last_p2p_recv_from_remote_procs.size(); i++)
        {
            int index = 0;
            while(true)
            {
                if(global_rank_intra_comm[index] == last_p2p_recv_from_remote_procs[i])
                {
                    memcpy(last_p2p_send_model_proc_lists+i*send_model_size, all_send_model_list+index*send_model_size, sizeof(int)*send_model_size);
                    break;
                }
                index ++;
            }
        }
    }

    /*
    if(butterfly_node_master && (butterfly_grp_index == 0))
    {
        for(int i=0; i<last_p2p_send_to_remote_procs.size(); i++)
            MPI_Wait(send_req+i, &status);
        delete [] send_req;
        delete [] tmp_send_model_list;
    }
    */

    if(is_recv)
    {
        last_p2p_recv_cells = new int[last_p2p_recv_from_remote_procs.size()];
        memset(last_p2p_recv_cells, 0, sizeof(int)*last_p2p_recv_from_remote_procs.size());
        /*
        for(int i=0; i<last_p2p_recv_from_remote_procs.size(); i++)
            MPI_Wait(recv_req+i, &status);
        delete [] recv_req;
        */

        for(int i=0; i<last_p2p_recv_from_remote_procs.size(); i++)
        {
            for(int j=0; j<send_model_size; j++)
            {
                int tmp_send_proc = last_p2p_send_model_proc_lists[i*send_model_size+j];
                if(tmp_send_proc != MPI_PROC_NULL)
                    last_p2p_recv_cells[i] += total_recv_cells[local_rank*send_model_size+tmp_send_proc];
                last_p2p_recv_buf_cells += last_p2p_recv_cells[i];
            }
        }
    }

    intra_proc_displs_cells = NULL;

    if(is_recv)
    {
        intra_proc_displs_cells = new int[send_model_size];

        intra_proc_displs_cells[0] = 0;
        for(int i=1; i<send_model_size; i++)
          intra_proc_displs_cells[i] = intra_proc_displs_cells[i-1] + total_recv_cells[local_rank * send_model_size + i-1];
    }

    if (local_data_decompositions != NULL) delete [] local_data_decompositions;
    delete [] tmp_send_model_list;
    delete [] all_send_model_list;
    delete [] global_rank_intra_comm;
}

void Butterfly::last_p2p_execute(char * recv_data)
{
    if(last_p2p_stage_num == 0) return;
    if(butterfly_p2p_stage_num == 0) return;

    butterfly_gather_buf = butterfly_data_buf[cur_buf_index];
    last_p2p_send_buf = butterfly_data_buf[local_buf_index];
    MPI_Request * send_req, * recv_req;
    MPI_Status status;

    if(last_p2p_recv_buf_cells > 0 && (pre_fields_per_cell != fields_per_cell || last_p2p_recv_buf == NULL))
    {
        if(last_p2p_recv_buf != NULL)
        {
            delete [] last_p2p_recv_buf;
            last_p2p_recv_buf = NULL;
        }

        last_p2p_recv_buf = new char[fields_per_cell*last_p2p_recv_buf_cells];
    }

    if(is_recv)
    {
        recv_req = new MPI_Request[last_p2p_recv_from_remote_procs.size()];
        int tmp_displs = 0;
        for(int i=0; i<last_p2p_recv_from_remote_procs.size(); i++)
        {
            if(last_p2p_recv_cells[i] > 0){
                MPI_Irecv(last_p2p_recv_buf+tmp_displs, last_p2p_recv_cells[i]*fields_per_cell, MPI_CHAR, last_p2p_recv_from_remote_procs[i], 3500+last_p2p_recv_from_remote_procs[i], global_comm, recv_req+i);
                tmp_displs += last_p2p_recv_cells[i]*fields_per_cell;
            }
        }
    }

    if(butterfly_node_master && (butterfly_grp_index == 0))
    {
        int * tmp_p2p_send_displs = new int[last_p2p_send_to_remote_procs.size()];
        for(int i=0; i<last_p2p_send_to_remote_procs.size(); i++)
            tmp_p2p_send_displs[i] = last_p2p_send_displs[i] * fields_per_cell;
        int tmp_displs = 0;
        for(int i=0; i<last_p2p_send_mask.size(); i++)
        {
            int tmp_mask = last_p2p_send_mask[i];
            int tmp_size = last_p2p_send_mask_cells[i] * fields_per_cell;
            memcpy(last_p2p_send_buf+tmp_p2p_send_displs[tmp_mask], butterfly_gather_buf+tmp_displs, sizeof(char)*tmp_size);
            tmp_p2p_send_displs[tmp_mask] += tmp_size;
            tmp_displs += tmp_size;
        }
        send_req = new MPI_Request[last_p2p_send_to_remote_procs.size()];
        for(int i=0; i<last_p2p_send_to_remote_procs.size(); i++)
        {
            if(last_p2p_send_cells[i] > 0)
                MPI_Isend(last_p2p_send_buf+last_p2p_send_displs[i]*fields_per_cell, last_p2p_send_cells[i]*fields_per_cell, MPI_CHAR, last_p2p_send_to_remote_procs[i], 3500+global_rank, global_comm, send_req+i);
        }
    }

    if(is_recv)
    {
        for(int i=0; i<last_p2p_recv_from_remote_procs.size(); i++)
            if(last_p2p_recv_cells[i] > 0)
                MPI_Wait(recv_req+i, &status);
        delete [] recv_req;
        
        int intra_proc_displs[send_model_size];
        for(int i=0; i<send_model_size; i++)
          intra_proc_displs[i] = intra_proc_displs_cells[i] * fields_per_cell;

        int tmp_displs = 0;
        for(int i=0; i<last_p2p_recv_from_remote_procs.size(); i++)
        {
            for(int j=0; j<send_model_size; j++)
            {
                int tmp_send_proc = last_p2p_send_model_proc_lists[i*send_model_size+j];
                if(tmp_send_proc != MPI_PROC_NULL)
                {
                    memcpy(recv_data+intra_proc_displs[tmp_send_proc], last_p2p_recv_buf+tmp_displs, sizeof(char)*fields_per_cell*total_recv_cells[local_rank*send_model_size+tmp_send_proc]);
                    tmp_displs += total_recv_cells[local_rank*send_model_size+tmp_send_proc] * fields_per_cell;
                }
            }
        }
    }

    if(butterfly_node_master && (butterfly_grp_index == 0))
    {
        for(int i=0; i<last_p2p_send_to_remote_procs.size(); i++)
            if(last_p2p_send_cells[i] > 0)
                MPI_Wait(send_req+i, &status);
        delete [] send_req;
    }

    pre_fields_per_cell = fields_per_cell;
}

void Butterfly::butterfly_init()
{
    butterfly_stage_num = (int)log2f(butterfly_grp_node_num);
    
    MPI_Comm inter_comm;

    if(action == SENDRECV) intra_comm = local_comm;
    else if(action == SEND){
        MPI_Intercomm_create(local_comm, 0, global_comm, recv_model_global_rank[0], 1000, &inter_comm);
        MPI_Intercomm_merge(inter_comm, true, &intra_comm);
    }
    else{
        MPI_Intercomm_create(local_comm, 0, global_comm, send_model_global_rank[0], 1000, &inter_comm);
        MPI_Intercomm_merge(inter_comm, true, &intra_comm);
    }
    best_timer = -1.0;
    profiling = true;
    current_stage = 0;
    stage_mask = new bool[butterfly_stage_num];
    for(int i = 0; i < butterfly_stage_num; i ++)
        stage_mask[i] = true;

    recv_proc_node_index = new int[recv_model_size];

    MPI_Status status;

    if(is_recv)
    {
        MPI_Gather(&recv_node_index, 1, MPI_INT, recv_proc_node_index, 1, MPI_INT, 0, local_comm);
        if(local_rank == 0)
          if(send_model_global_rank[0] != global_rank)
            MPI_Send(recv_proc_node_index, recv_model_size, MPI_INT, send_model_global_rank[0], 1400+send_model_global_rank[0], global_comm);
    }

    if(is_send)
        if(local_rank == 0)
          if(recv_model_global_rank[0] != global_rank)
            MPI_Recv(recv_proc_node_index, recv_model_size, MPI_INT, recv_model_global_rank[0], 1400+global_rank, global_comm, &status);

    MPI_Bcast(recv_proc_node_index, recv_model_size, MPI_INT, 0, local_comm);
    
    if(!butterfly_node_master) return;

    butterfly_node_decompsitions = new vector<int> [butterfly_grp_node_num];
    butterfly_node_decompsitions_cells = new vector<int> [butterfly_grp_node_num];

    butterfly_grp_index = butterfly_node_index / butterfly_grp_node_num;
    butterfly_grp_node_index = butterfly_node_index % butterfly_grp_node_num;


    butterfly_pair_node_list = new int [butterfly_grp_node_num];
    butterfly_pair_node_list[0] = butterfly_node_index;
    
    int stage_size = 1;
    for(int i=0; i<butterfly_stage_num; i++)
    {
        int pair_node;
        if((butterfly_grp_node_index/stage_size) % 2 == 0)
            pair_node = butterfly_grp_node_index + stage_size;
        else
            pair_node = butterfly_grp_node_index - stage_size;
        compute_butterfly_recv_rank(butterfly_pair_node_list+stage_size, pair_node, stage_size);
        stage_size *= 2;
    }

    for(int i=0; i<butterfly_grp_node_num; i++)
    {
        int tmp_butterfly_node = butterfly_grp_index * butterfly_grp_node_num + i;
        for(int j=0; j<butterfly_node_map_to_send_node[tmp_butterfly_node].size(); j++)
        {
            int tmp_send_node = butterfly_node_map_to_send_node[tmp_butterfly_node][j];
            for(int k=0; k<send_node_per_proc_local_rank[tmp_send_node].size(); k++)
            {
                int tmp_send_proc = send_node_per_proc_local_rank[tmp_send_node][k];
                for(int m=0; m<recv_model_size; m++)
                {
                    int tmp_recv_node = recv_proc_node_index[m];
                    int tmp_butterfly_num = recv_node_map_to_butterfly_node[tmp_recv_node].size();
                    int tmp_total_send_cells = total_send_cells[tmp_send_proc * recv_model_size + m];
                    int tmp_send_cells = tmp_total_send_cells / tmp_butterfly_num;
                    int tmp_send_proc_reminder = tmp_total_send_cells % tmp_butterfly_num;
                    for(int n=0; n<tmp_butterfly_num; n++)
                    {
                        if(n < tmp_send_proc_reminder)
                        {
                            butterfly_node_decompsitions[i].push_back(recv_node_map_to_butterfly_node[tmp_recv_node][n]);
                            butterfly_node_decompsitions_cells[i].push_back(tmp_send_cells+1);
                        }
                        else if(tmp_send_cells > 0)
                        {
                            butterfly_node_decompsitions[i].push_back(recv_node_map_to_butterfly_node[tmp_recv_node][n]);
                            butterfly_node_decompsitions_cells[i].push_back(tmp_send_cells);
                        }
                    }
                }
            }
        }
    }

    butterfly_grp_num = butterfly_node_num / butterfly_grp_node_num;
}

void Butterfly::butterfly_p2p_init(int * stage_num)
{
    butterfly_p2p_remote_procs = NULL;
    butterfly_p2p_send_cells = NULL;
    butterfly_p2p_recv_cells = NULL;
    butterfly_p2p_send_displs = NULL;
    butterfly_p2p_recv_displs = NULL;
    butterfly_p2p_send_mask = NULL;
    butterfly_p2p_send_mask_cells = NULL;
    p2p_num_per_butterfly_stage = NULL;
    p2p_depth_per_butterfly_stage = NULL;

    p2p_num_per_butterfly_stage = new int[butterfly_stage_num];
    p2p_depth_per_butterfly_stage = new int[butterfly_stage_num];

    butterfly_p2p_stage_num = 0;
    p2p_depth = 0;

    if(stage_num != NULL) memcpy(p2p_num_per_butterfly_stage, stage_num, sizeof(int)*butterfly_stage_num);

    int tmp_butterfly_stage = 0;
    for(int i=0; i<butterfly_stage_num; i++)
    {
        if(stage_num == NULL) p2p_num_per_butterfly_stage[i] = 1;
        p2p_depth_per_butterfly_stage[i] = (int) exp2f(p2p_num_per_butterfly_stage[i]);
        if((tmp_butterfly_stage + p2p_num_per_butterfly_stage[i]) >= butterfly_stage_num)
        {
            last_p2p_stage_num = butterfly_stage_num - tmp_butterfly_stage;
            break;
        }
        else
        {
            butterfly_p2p_stage_num ++;
            tmp_butterfly_stage += p2p_num_per_butterfly_stage[i];
            p2p_depth += p2p_depth_per_butterfly_stage[i];
        }
    }

    if(!butterfly_node_master) return;
    if(butterfly_p2p_stage_num == 0) return;

    butterfly_p2p_remote_procs = new int[p2p_depth];
    butterfly_p2p_send_cells = new int[p2p_depth];
    butterfly_p2p_recv_cells = new int[p2p_depth];
    butterfly_p2p_send_displs = new int[p2p_depth];
    butterfly_p2p_recv_displs = new int[p2p_depth];
    for(int i=0; i<p2p_depth; i++){
        butterfly_p2p_remote_procs[i] = MPI_PROC_NULL;
        butterfly_p2p_send_cells[i] =0;
        butterfly_p2p_recv_cells[i] =0;
        butterfly_p2p_send_displs[i] =0;
        butterfly_p2p_recv_displs[i] =0;
    }
    butterfly_p2p_send_mask = new vector<int>[butterfly_p2p_stage_num];
    butterfly_p2p_send_mask_cells = new vector<int>[butterfly_p2p_stage_num];

    int cur_butterfly_stage_num = 0;
    int pre_stage_size = 1;
    int stage_size = (int)exp2f(cur_butterfly_stage_num);
    int tmp_index = 0;
    int * local_data_decompositions = new int[p2p_depth*butterfly_grp_node_num/2];
    int tmp_decomp_index = 0;
    int data_decompositions_num;

    for(int i=0; i<butterfly_p2p_stage_num; i++)
    {
        pre_stage_size = (int) exp2f(cur_butterfly_stage_num);
        cur_butterfly_stage_num += p2p_num_per_butterfly_stage[i];
        stage_size = (int) exp2f(cur_butterfly_stage_num);

        tmp_decomp_index = 0;
        for(int j=0; j<stage_size; j+=pre_stage_size)
        {
            int tmp_butterfly_node = butterfly_pair_node_list[j] - butterfly_grp_index * butterfly_grp_node_num;
            butterfly_p2p_remote_procs[tmp_index] = butterfly_node_global_rank[tmp_butterfly_node];
            data_decompositions_num = butterfly_grp_node_num / stage_size;
            for(int k = 0; k<data_decompositions_num; k++)
            {
                local_data_decompositions[tmp_decomp_index++] = tmp_butterfly_node % stage_size + k * stage_size;
            }

            for(int k=0; k<pre_stage_size; k++)
            {
                int tmp_butterfly_node = butterfly_pair_node_list[j+k] - butterfly_grp_index * butterfly_grp_node_num;
                for(int m=0; m<butterfly_node_decompsitions[tmp_butterfly_node].size(); m++)
                {
                    int tmp_decomposition = butterfly_node_decompsitions[tmp_butterfly_node][m];
                    int tmp_decomposition_cells = butterfly_node_decompsitions_cells[tmp_butterfly_node][m];
                    bool tmp_exist = false;
                    for(int n=0; n<data_decompositions_num; n++)
                    {
                        if(tmp_decomposition == local_data_decompositions[n])
                        {
                            tmp_exist = true;
                            break;
                        }
                    }
                    if(tmp_exist) butterfly_p2p_recv_cells[tmp_index] += tmp_decomposition_cells;
                }
            }
            tmp_index++;
        }

        for(int k=0; k<pre_stage_size; k++)
        {
            int tmp_butterfly_node = butterfly_pair_node_list[k] - butterfly_grp_index * butterfly_grp_node_num;
            for(int m=0; m<butterfly_node_decompsitions[tmp_butterfly_node].size(); m++)
            {
                int tmp_decomposition = butterfly_node_decompsitions[tmp_butterfly_node][m];
                int tmp_decomposition_cells = butterfly_node_decompsitions_cells[tmp_butterfly_node][m];
                int tmp_mask = -1;
                for(int n=0; n<tmp_decomp_index; n++)
                {
                    if(tmp_decomposition == local_data_decompositions[n])
                    {
                        tmp_mask = n/data_decompositions_num;
                        break;
                    }
                }

                if(tmp_mask != -1)
                {
                    butterfly_p2p_send_mask[i].push_back(tmp_mask);
                    butterfly_p2p_send_mask_cells[i].push_back(tmp_decomposition_cells);
                    butterfly_p2p_send_cells[tmp_index-p2p_depth_per_butterfly_stage[i]+tmp_mask] += tmp_decomposition_cells;
                }
            }
        }
    }

    tmp_index = 0;
    for(int i=0; i<butterfly_p2p_stage_num; i++)
    {
        for(int j=1; j<p2p_depth_per_butterfly_stage[i]; j++)
            butterfly_p2p_recv_displs[tmp_index+j] = butterfly_p2p_recv_displs[tmp_index+j-1] + butterfly_p2p_recv_cells[tmp_index+j-1];
        int tmp_size = butterfly_p2p_recv_cells[tmp_index+p2p_depth_per_butterfly_stage[i]-1] + butterfly_p2p_recv_displs[tmp_index+p2p_depth_per_butterfly_stage[i]-1];
        if(max_cells_cnts < tmp_size) max_cells_cnts = tmp_size;
        for(int j=2; j<p2p_depth_per_butterfly_stage[i]; j++)
            butterfly_p2p_send_displs[tmp_index+j] = butterfly_p2p_send_displs[tmp_index+j-1] + butterfly_p2p_send_cells[tmp_index+j-1];
        tmp_size = butterfly_p2p_send_displs[tmp_index+p2p_depth_per_butterfly_stage[i]-1] + butterfly_p2p_send_cells[tmp_index+p2p_depth_per_butterfly_stage[i]-1];
        if(max_cells_cnts < tmp_size) max_cells_cnts = tmp_size;
        tmp_index += p2p_depth_per_butterfly_stage[i];
    }

    delete [] local_data_decompositions;
}

void Butterfly::butterfly_p2p_execute()
{
    cur_buf_index = 0;
    local_buf_index = 1;
    remote_buf_index = 2;
    
    if(!butterfly_node_master) return;
    if(butterfly_p2p_stage_num == 0) return;

    MPI_Status status;
    MPI_Request send_req[p2p_depth], recv_req[p2p_depth];

    int tmp_send_displs[p2p_depth], tmp_recv_displs[p2p_depth];
    for(int i=0; i<p2p_depth; i++)
    {
        tmp_send_displs[i] = butterfly_p2p_send_displs[i] * fields_per_cell;
        tmp_recv_displs[i] = butterfly_p2p_recv_displs[i] * fields_per_cell;
    }

    int tmp_index = p2p_depth_per_butterfly_stage[0];
    for(int i=1; i<butterfly_p2p_stage_num; i++)
    {
        char * cur_buf = butterfly_data_buf[cur_buf_index];
        char * local_buf = butterfly_data_buf[local_buf_index];
        char * remote_buf = butterfly_data_buf[remote_buf_index];

        for(int j=1; j<p2p_depth_per_butterfly_stage[i]; j++)
            if(butterfly_p2p_recv_cells[tmp_index+j] > 0)
                MPI_Irecv(local_buf+tmp_recv_displs[tmp_index+j], butterfly_p2p_recv_cells[tmp_index+j]*fields_per_cell, MPI_CHAR, butterfly_p2p_remote_procs[tmp_index+j], 2000+global_rank, global_comm, recv_req+j);

        int cur_buf_displs = 0;
        int local_buf_displs = 0;
        int remote_buf_displs = 0;
        for(int j=0; j<butterfly_p2p_send_mask[i].size(); j++)
        {
            if(butterfly_p2p_send_mask[i][j] == 0)
            {
                memcpy(local_buf+tmp_send_displs[tmp_index], cur_buf+cur_buf_displs, sizeof(char)*butterfly_p2p_send_mask_cells[i][j]*fields_per_cell);
                tmp_send_displs[tmp_index] += butterfly_p2p_send_mask_cells[i][j]*fields_per_cell;
                cur_buf_displs += butterfly_p2p_send_mask_cells[i][j]*fields_per_cell;
            }
            else
            {
                int tmp_mask = butterfly_p2p_send_mask[i][j];
                memcpy(remote_buf+tmp_send_displs[tmp_index+tmp_mask], cur_buf+cur_buf_displs, sizeof(char)*butterfly_p2p_send_mask_cells[i][j]*fields_per_cell);
                tmp_send_displs[tmp_index+tmp_mask] += butterfly_p2p_send_mask_cells[i][j]*fields_per_cell;
                cur_buf_displs += butterfly_p2p_send_mask_cells[i][j]*fields_per_cell;
            }
        }


        for(int j=1; j<p2p_depth_per_butterfly_stage[i]; j++)
            if(butterfly_p2p_send_cells[tmp_index+j] > 0)
                MPI_Isend(remote_buf+butterfly_p2p_send_displs[tmp_index+j]*fields_per_cell, butterfly_p2p_send_cells[tmp_index+j]*fields_per_cell, MPI_CHAR, butterfly_p2p_remote_procs[tmp_index+j], 2000+butterfly_p2p_remote_procs[tmp_index+j], global_comm, send_req+j);

        for(int j=1; j<p2p_depth_per_butterfly_stage[i]; j++)
            if(butterfly_p2p_send_cells[tmp_index+j] > 0)
                MPI_Wait(send_req+j, &status);
        for(int j=1; j<p2p_depth_per_butterfly_stage[i]; j++)
            if(butterfly_p2p_recv_cells[tmp_index+j] > 0)
                MPI_Wait(recv_req+j, &status);

        int tmp_buf_index = cur_buf_index;
        cur_buf_index = local_buf_index;
        local_buf_index = remote_buf_index;
        remote_buf_index = tmp_buf_index;
        tmp_index += p2p_depth_per_butterfly_stage[i];
    }

}

void Butterfly::process_mapping_from_butterfly_to_recv_init()
{
    MPI_Status status;

    int * tmp_recv_node_data_size = new int[recv_node_num];
    int * tmp_recv_node_index = new int[recv_node_num];
    for(int i=0; i<recv_node_num; i++)
        tmp_recv_node_index[i] = i;
    
    recv_node_master = false;
    if(is_recv)
    {
        int * tmp_recv_node_recv_cells_cnts;
        int tmp_recv_node_recv_cells = 0;
        if(recv_node_per_proc_local_rank[recv_node_index][0] == local_rank)
        {
            recv_node_master = true;
            for(int i=0; i<recv_node_per_proc_local_rank[recv_node_index].size(); i++)
            {
                int tmp_local_rank = recv_node_per_proc_local_rank[recv_node_index][i];
                for(int j=0; j<send_model_size; j++)
                {
                    tmp_recv_node_recv_cells += total_recv_cells[tmp_local_rank*send_model_size+j];
                }
            }
        }
        if(local_rank == 0)
        {
            tmp_recv_node_recv_cells_cnts = new int[recv_model_size];
            MPI_Gather(&tmp_recv_node_recv_cells, 1, MPI_INT, tmp_recv_node_recv_cells_cnts, 1, MPI_INT, 0, local_comm);
            for(int i=0; i<recv_node_num; i++)
              tmp_recv_node_data_size[i] = tmp_recv_node_recv_cells_cnts[recv_node_per_proc_local_rank[i][0]];
            if(global_rank != send_model_global_rank[0])
              MPI_Send(tmp_recv_node_data_size, recv_node_num, MPI_INT, send_model_global_rank[0], 1300+send_model_global_rank[0], global_comm);
            delete tmp_recv_node_recv_cells_cnts;
        }
        else
            MPI_Gather(&tmp_recv_node_recv_cells, 1, MPI_INT, tmp_recv_node_recv_cells_cnts, 1, MPI_INT, 0, local_comm);
    }

    if(is_send)
      if(local_rank == 0)
        if(global_rank != recv_model_global_rank[0])
          MPI_Recv(tmp_recv_node_data_size, recv_node_num, MPI_INT, recv_model_global_rank[0], 1300+global_rank, global_comm, &status);

    MPI_Bcast(tmp_recv_node_data_size, recv_node_num, MPI_INT, 0, local_comm);
    decrease_sort(tmp_recv_node_data_size, tmp_recv_node_index, recv_node_num);

    butterfly_node_map_to_recv_node = new vector<int> [butterfly_grp_node_num];
    recv_node_map_to_butterfly_node = new vector<int> [recv_node_num];

    int padding_loop_num = (int) log2f(recv_node_num);
    if( recv_node_num > (int) exp2f(padding_loop_num) )
      padding_loop_num ++;
    int butterfly_loop_num = (int) log2f(butterfly_grp_node_num);
    if( padding_loop_num < butterfly_loop_num )
      padding_loop_num = butterfly_loop_num;
    int padding_nodes_num = (int) exp2f(padding_loop_num);

    vector<int> * tmp_group = new vector<int> [padding_nodes_num];
    int * tmp_group_size = new int[padding_nodes_num];
    int * tmp_group_index = new int[padding_nodes_num];

    for(int i=0; i<recv_node_num; i++)
    {
        int tmp_node_index = tmp_recv_node_index[i]; 
        tmp_group[i].push_back(tmp_node_index);
        tmp_group_index[i] = i;
        tmp_group_size[i] = tmp_recv_node_data_size[i];
    }

    for(int i=recv_node_num; i<padding_nodes_num; i++)
    {
        tmp_group_index[i] = i;
        tmp_group_size[i] = 0;
    }

    for(int i=padding_loop_num-1; i>=butterfly_loop_num; i--)
    {
        int tmp_group_num = (int) exp2f(i);
        for(int j=0; j<tmp_group_num; j++)
        {
            int first_group_index = tmp_group_index[j];
            int second_group_index = tmp_group_index[2*tmp_group_num-1-j];
            tmp_group_size[j] += tmp_group_size[2*tmp_group_num-1-j];
            tmp_group_size[2*tmp_group_num-1-j] = 0;
            for(int k=0; k<tmp_group[second_group_index].size(); k++)
              tmp_group[first_group_index].push_back(tmp_group[second_group_index][k]);
        }
        decrease_sort(tmp_group_size, tmp_group_index, tmp_group_num);
    }

    vector<int> * tmp_map = new vector<int> [butterfly_grp_node_num];
    int * tmp_data_size  = new int[butterfly_grp_node_num];
    int * tmp_map_index = new int[butterfly_grp_node_num];
 
    for(int i=0; i<butterfly_grp_node_num; i++)
    {
        tmp_map[i].push_back(tmp_group_index[i]);
        tmp_data_size[i] = tmp_group_size[i];
        tmp_map_index[i] = i;
    }

    for(int i=butterfly_loop_num-1; i>=0; i--)
    {
        int tmp_group_num = (int) exp2f(i);
        for(int j=0; j<tmp_group_num; j++)
        {
            int first_map_index = tmp_map_index[j];
            int second_map_index = tmp_map_index[2*tmp_group_num-1-j];
            tmp_data_size[j] += tmp_data_size[2*tmp_group_num-1-j];
            for(int k=0; k<tmp_map[second_map_index].size(); k++)
              tmp_map[first_map_index].push_back(tmp_map[second_map_index][k]);
        }
        decrease_sort(tmp_data_size, tmp_map_index, tmp_group_num);
    }

    int * butterfly_node_index = new int[butterfly_grp_node_num];
    int cur_map_index = tmp_map_index[0];
    int tmp_step_size = butterfly_grp_node_num / 2;

    butterfly_node_index[0] = 0;
    for(int i=0; i<butterfly_loop_num; i++)
    {
        for(int j=0; j<(int) exp2f(i); j++)
          butterfly_node_index[j+(int)exp2f(i)] = butterfly_node_index[j] + tmp_step_size;

        tmp_step_size /= 2;
    }

    for(int i=0; i<butterfly_grp_node_num; i++)
    {
        int cur_group_index = tmp_map[cur_map_index][i];
        int cur_butterfly_node = butterfly_node_index[i];
        for(int j=0; j<tmp_group[cur_group_index].size(); j++)
        {
            int tmp_recv_node = tmp_group[cur_group_index][j];
            butterfly_node_map_to_recv_node[cur_butterfly_node].push_back(tmp_recv_node);
            recv_node_map_to_butterfly_node[tmp_recv_node].push_back(cur_butterfly_node);
        }
    }

    delete [] butterfly_node_index;
    delete [] tmp_map;
    delete [] tmp_map_index;
    delete [] tmp_data_size;
    delete [] tmp_group;
    delete [] tmp_group_size;
    delete [] tmp_group_index;
    delete [] tmp_recv_node_data_size;
    delete [] tmp_recv_node_index;
}

void Butterfly::compute_butterfly_recv_rank(int * buf, int index, int size)
{
    int pair_node;
    int next_pair_node;
    pair_node = index;

    if(size == 1)
        buf[0] = pair_node + butterfly_grp_index * butterfly_grp_node_num;
    else
    {
        size /= 2;
        
        if((pair_node / size) % 2 == 0)
            next_pair_node = pair_node + size;
        else
            next_pair_node = pair_node - size;
        compute_butterfly_recv_rank(buf+size, next_pair_node, size);

        compute_butterfly_recv_rank(buf, pair_node, size);
    }
}

void Butterfly::process_mapping_from_send_to_butterfly_init()
{
    butterfly_data_buf[0] = NULL;
    butterfly_data_buf[1] = NULL;
    butterfly_data_buf[2] = NULL;

    MPI_Status status;

    if(butterfly_node_num == -1)
        butterfly_node_num = (int) exp2f((int)log2f(total_node_num));
    
    butterfly_grp_node_num = butterfly_node_num;

    if(butterfly_node_num > total_node_num)
      MPI_Abort(global_comm, -100);

    send_node_cells_cnts = new int[send_node_num];
    int * tmp_send_node_data_size = new int[send_node_num];
    int * tmp_send_node_index = new int[send_node_num];
    for(int i=0; i<send_node_num; i++)
        tmp_send_node_index[i] = i;

    intra_node_recv_cells_cnts = 0;
    if(is_send)
    {
        int * tmp_intra_node_recv_cells_cnts;
        for(int i=0; i<recv_model_size; i++)
          intra_node_recv_cells_cnts += total_send_cells[local_rank * recv_model_size + i];
          
        if(local_rank == 0)
        {
            tmp_intra_node_recv_cells_cnts = new int[send_model_size];
            MPI_Gather(&intra_node_recv_cells_cnts, 1, MPI_INT, tmp_intra_node_recv_cells_cnts, 1, MPI_INT, 0, local_comm);
            for(int i=0; i<send_node_num; i++)
                tmp_send_node_data_size[i] = tmp_intra_node_recv_cells_cnts[send_node_per_proc_local_rank[i][0]];
            if(global_rank != recv_model_global_rank[0]) 
                MPI_Send(tmp_send_node_data_size, send_node_num, MPI_INT, recv_model_global_rank[0], 1200+recv_model_global_rank[0], global_comm);
            delete tmp_intra_node_recv_cells_cnts;
        }
        else
            MPI_Gather(&intra_node_recv_cells_cnts, 1, MPI_INT, tmp_intra_node_recv_cells_cnts, 1, MPI_INT, 0, local_comm);

    }

    if(is_recv)
        if(local_rank == 0)
            if(global_rank != send_model_global_rank[0]) 
                MPI_Recv(tmp_send_node_data_size, send_node_num, MPI_INT, send_model_global_rank[0], 1200+global_rank, global_comm, &status);

    MPI_Bcast(tmp_send_node_data_size, send_node_num, MPI_INT, 0, local_comm);
    memcpy(send_node_cells_cnts , tmp_send_node_data_size, sizeof(int)*send_node_num);
    decrease_sort(tmp_send_node_data_size, tmp_send_node_index, send_node_num);

    int padding_loop_num = (int) log2f(send_node_num);
    if( send_node_num > (int) exp2f(padding_loop_num) )
      padding_loop_num ++;
    int butterfly_loop_num = (int) log2f(butterfly_node_num);
    if( padding_loop_num < butterfly_loop_num )
      padding_loop_num = butterfly_loop_num;
    int padding_nodes_num = (int) exp2f(padding_loop_num);

    vector<int> *tmp_group = new vector<int> [padding_nodes_num];
    int * tmp_group_size = new int[padding_nodes_num];
    int * tmp_group_index = new int[padding_nodes_num];

    for(int i=0; i<send_node_num; i++)
    {
        int tmp_node_index = tmp_send_node_index[i];
        total_node_root_proc_global_rank[i] = send_model_global_rank[send_node_per_proc_local_rank[tmp_node_index][0]];
        tmp_group[i].push_back(tmp_node_index);
        tmp_group_index[i] = i;
        tmp_group_size[i] = tmp_send_node_data_size[i];
    }

    for(int i=send_node_num; i<padding_nodes_num; i++)
    {
        tmp_group_index[i] = i;
        tmp_group_size[i] = 0;
    }

    for(int i=padding_loop_num-1; i>=butterfly_loop_num; i--)
    {
        int tmp_group_num = (int) exp2f(i);
        for(int j=0; j<tmp_group_num; j++)
        {
            int first_group_index = tmp_group_index[j];
            int second_group_index = tmp_group_index[2*tmp_group_num-1-j];
            tmp_group_size[j] += tmp_group_size[2*tmp_group_num-1-j];
            tmp_group_size[2*tmp_group_num-1-j] = 0;
            for(int k=0; k<tmp_group[second_group_index].size(); k++)
              tmp_group[first_group_index].push_back(tmp_group[second_group_index][k]);
        }
        decrease_sort(tmp_group_size, tmp_group_index, tmp_group_num);
    }

    vector<int> * tmp_map = new vector<int> [butterfly_node_num];
    int * tmp_data_size  = new int[butterfly_node_num];
    int * tmp_map_index = new int[butterfly_node_num];
 
    for(int i=0; i<butterfly_node_num; i++)
    {
        tmp_map[i].push_back(tmp_group_index[i]);
        tmp_data_size[i] = tmp_group_size[i];
        tmp_map_index[i] = i;
    }

    for(int i=butterfly_loop_num-1; i>=0; i--)
    {
        int tmp_group_num = (int) exp2f(i);
        for(int j=0; j<tmp_group_num; j++)
        {
            int first_map_index = tmp_map_index[j];
            int second_map_index = tmp_map_index[2*tmp_group_num-1-j];
            tmp_data_size[j] += tmp_data_size[2*tmp_group_num-1-j];
            for(int k=0; k<tmp_map[second_map_index].size(); k++)
              tmp_map[first_map_index].push_back(tmp_map[second_map_index][k]);
        }
        decrease_sort(tmp_data_size, tmp_map_index, tmp_group_num);
    }
    
    int tmp_recv_node_index = 0;
    int cur_map_index = tmp_map_index[0];
    for(int i=0; i<butterfly_node_num; i++)
    {
        int cur_group_index = tmp_map[cur_map_index][i];
        
        if(tmp_group[cur_group_index].size() == 0)
        {
            total_node_root_proc_global_rank[i] = recv_model_global_rank[recv_node_per_proc_local_rank[tmp_recv_node_index][0]];
            tmp_recv_node_index ++;
        }
        else
        {
            int tmp_local_rank = send_node_per_proc_local_rank[tmp_group[cur_group_index][0]][0];
            total_node_root_proc_global_rank[i] = send_model_global_rank[tmp_local_rank];
        }
    }

    butterfly_node_master = false;
    butterfly_node_index = -1;
    butterfly_node_global_rank = new int [butterfly_node_num];
    butterfly_node_map_to_send_node = new vector<int>[butterfly_node_num];
    butterfly_node_recv_cells_cnts = new vector<int>[butterfly_node_num];

    for(int i=0; i<butterfly_node_num; i++)
    {
        butterfly_node_global_rank[i] = total_node_root_proc_global_rank[i];
                
        int cur_group_index = tmp_map[cur_map_index][i];
        for(int j=0; j<tmp_group[cur_group_index].size(); j++)
        {
            int tmp_local_node = tmp_group[cur_group_index][j];
            int tmp_local_rank = send_node_per_proc_local_rank[tmp_local_node][0];
            int tmp_global_rank = send_model_global_rank[tmp_local_rank];
            if(global_rank == tmp_global_rank){
                send_node_map_to_butterfly_node = butterfly_node_global_rank[i];
                send_node_map_to_butterfly_node_index = i;
            }
            butterfly_node_map_to_send_node[i].push_back(tmp_local_node);
            butterfly_node_recv_cells_cnts[i].push_back(send_node_cells_cnts[tmp_local_node]);
        }
                
        if(global_rank == butterfly_node_global_rank[i])
        {
            butterfly_node_master = true;
            butterfly_node_index = i;
        }

    }

    
    delete [] tmp_map_index;
    delete [] tmp_map;
    delete [] tmp_data_size;
    delete [] tmp_send_node_index;
    delete [] tmp_send_node_data_size;
    delete [] tmp_group;
    delete [] tmp_group_index;
    delete [] tmp_group_size;
}

void Butterfly::intra_node_gather_init(MPI_Comm global_comm, MPI_Comm local_comm, int * remote_model_global_rank, int * send_cells, int * recv_cells, int remote_model_size, int action)
{
    is_send = true;
    is_recv = true;
    if(action == SEND)
      is_recv = false;
    else if(action == RECV)
      is_send = false;

    this->global_comm = global_comm;
    this->local_comm = local_comm;
    this->action = action;
    this->remote_model_size = remote_model_size;

    MPI_Comm_size(local_comm, &local_model_size);
    MPI_Comm_rank(local_comm, &local_rank);
    MPI_Comm_rank(global_comm, &global_rank);

    char filename[MPI_MAX_PROCESSOR_NAME];
    sprintf(filename, "File_%d.txt", global_rank);
    
    this->send_cells = NULL;
    this->recv_cells = NULL;
    if(is_send)
    {
        send_model_size = local_model_size;
        recv_model_size = remote_model_size;
        this->send_cells = new int[recv_model_size];
        memcpy(this->send_cells, send_cells, sizeof(int) * recv_model_size);
    }
    if(is_recv)
    {
        send_model_size = remote_model_size;
        recv_model_size = local_model_size;
        this->recv_cells = new int[send_model_size];
        memcpy(this->recv_cells, recv_cells, sizeof(int) * send_model_size);
    }

    send_model_global_rank = new int[send_model_size]; 
    this->recv_model_global_rank = new int[recv_model_size];
    MPI_Status status;


    if(is_send)
    {
        memcpy(this->recv_model_global_rank, remote_model_global_rank, sizeof(int)*remote_model_size);
        MPI_Allgather(&global_rank, 1, MPI_INT, send_model_global_rank, 1, MPI_INT, local_comm);
    }
    if(is_recv)
    {
        memcpy(this->send_model_global_rank, remote_model_global_rank, sizeof(int)*remote_model_size);
        MPI_Allgather(&global_rank, 1, MPI_INT, recv_model_global_rank, 1, MPI_INT, local_comm);
    }
    
    send_node_num = send_model_size;
    recv_node_num = recv_model_size;
    if(action != SEND && action != RECV)
        total_node_num = send_node_num;
    else
        total_node_num = send_node_num+recv_node_num;

    send_node_index = -1;
    recv_node_index = -1;

    send_node_per_proc_local_rank = new vector<int> [send_node_num];
    recv_node_per_proc_local_rank = new vector<int> [recv_node_num];
    total_node_root_proc_global_rank = new int[total_node_num];

    for(int i=0; i<send_node_num; i++)
      send_node_per_proc_local_rank[i].push_back(i);

    for(int i=0; i<recv_node_num; i++)
      recv_node_per_proc_local_rank[i].push_back(i);

    for(int i=0; i<send_node_num; i++)
      total_node_root_proc_global_rank[i] = send_model_global_rank[i];

    for(int i=send_node_num; i<total_node_num; i++)
      total_node_root_proc_global_rank[i] = recv_model_global_rank[i-send_model_size];

    total_send_cells = new int[send_model_size*recv_model_size];
    total_recv_cells = new int[send_model_size*recv_model_size];

    if(is_send)
    {
        send_node_index = local_rank;
        MPI_Allgather(send_cells, recv_model_size, MPI_INT, total_send_cells, recv_model_size, MPI_INT, local_comm);
        if(local_rank == 0)
        {
            if(recv_model_global_rank[0] != global_rank)
            {
                MPI_Send(total_send_cells, send_model_size*recv_model_size, MPI_INT, recv_model_global_rank[0], 1100+global_rank, global_comm);
                MPI_Recv(total_recv_cells, recv_model_size*send_model_size, MPI_INT, recv_model_global_rank[0], 1100+recv_model_global_rank[0], global_comm, &status);
            }
        }
        MPI_Bcast(total_recv_cells, send_model_size*recv_model_size, MPI_INT, 0, local_comm);
    }
    if(is_recv)
    {
        recv_node_index = local_rank;
        MPI_Allgather(recv_cells, send_model_size, MPI_INT, total_recv_cells, send_model_size, MPI_INT, local_comm);
        if(local_rank == 0)
        {
            if(send_model_global_rank[0] != global_rank)
            {
                MPI_Recv(total_send_cells, send_model_size*recv_model_size, MPI_INT, send_model_global_rank[0], 1100+send_model_global_rank[0], global_comm, &status);
                MPI_Send(total_recv_cells, recv_model_size*send_model_size, MPI_INT, send_model_global_rank[0], 1100+global_rank, global_comm);
            }
        }
        MPI_Bcast(total_send_cells, send_model_size*recv_model_size, MPI_INT, 0, local_comm);
    }


    send_node_master = false;

    if(is_send)
        send_node_master = true;
    
}

void Butterfly::execute(char * send_data, char * recv_data, int fields_per_cell)
{
    this->fields_per_cell = fields_per_cell;
    
    if(profiling)
    {
        double tstart, tend, current_timer, max_timer;
        int best_stage_num[butterfly_stage_num];
        memset(best_stage_num, 0, sizeof(int) * butterfly_stage_num);
        if(current_stage != 0)
        {
            int index = 0;
            stage_mask[current_stage] = false;
            best_stage_num[index] = 1;
            for(int i = 1; i < butterfly_stage_num; i++)
            {
                if(stage_mask[i])
                {
                    index ++;
                    best_stage_num[index] = 1;
                }
                else best_stage_num[index] ++;
            }

            reset_p2p_stage_num(best_stage_num);
        }

        MPI_Barrier(intra_comm);
        wtime(&tstart);
        first_p2p_execute(send_data, recv_data);
        butterfly_p2p_execute();
        last_p2p_execute(recv_data);
        wtime(&tend);
        current_timer = tend - tstart;
        MPI_Allreduce(&current_timer, &max_timer, 1, MPI_DOUBLE, MPI_MAX, intra_comm);
        if(current_stage == 0) best_timer = max_timer;
        else if(max_timer < best_timer) best_timer = max_timer;
        else stage_mask[current_stage] = true;

        current_stage ++;
        if(current_stage >= butterfly_stage_num)
        {
            int index = 0;
            best_stage_num[index] = 1;
            for(int i = 1; i < butterfly_stage_num; i++)
            {
                if(stage_mask[i])
                {
                    index ++;
                    best_stage_num[index] = 1;
                }
                else best_stage_num[index] ++;
            }

            reset_p2p_stage_num(best_stage_num);
            profiling = false;
        }
    }
    else
    {
        first_p2p_execute(send_data, recv_data);
        butterfly_p2p_execute();
        last_p2p_execute(recv_data);
    }
}

void Butterfly:: decrease_sort(int * data, int * index, int size)
{
    for(int i=0; i<size-1; i++)
    {
        for(int j=i+1; j<size; j++)
        {
            if(data[i] < data[j])
            {
                int tmp_data = data[i];
                int tmp_index = index[i];
                data[i] = data[j];
                index[i] = index[j];
                data[j] = tmp_data;
                index[j] = tmp_index;
            }
        }
    }
}

int * Butterfly:: auto_set_p2p_stage_num(char * send_data, char * recv_data, int fields_per_cell)
{
    double tstart, tend, best_timer, cur_timer, temp_timer, max_timer;
    int best_p2p_stage_num[butterfly_stage_num], cur_p2p_stage_num;
    bool best_p2p_stage_mask[butterfly_stage_num];

    profiling = false;
    for(int i=0; i<butterfly_stage_num; i++)
    {
        best_p2p_stage_num[i] = 1;
        best_p2p_stage_mask[i] = true;
    }

    reset_p2p_stage_num(best_p2p_stage_num);
    MPI_Barrier(intra_comm);
    wtime(&tstart);
    execute(send_data, recv_data, fields_per_cell);
    wtime(&tend);
    temp_timer = tend-tstart;
    MPI_Allreduce(&temp_timer, &max_timer, 1, MPI_DOUBLE, MPI_MAX, intra_comm);
    best_timer = max_timer;

    int index;
    for(int i=1; i<butterfly_stage_num; i++)
    {
        memset(best_p2p_stage_num, 0, sizeof(int)*butterfly_stage_num);
        best_p2p_stage_mask[i] = false;
        index = 0;
        best_p2p_stage_num[index] = 1;
        for(int j=1; j<butterfly_stage_num; j++)
        {
            if(best_p2p_stage_mask[j])
            {
                index ++;
                best_p2p_stage_num[index] = 1;
            }
            else
            {
                best_p2p_stage_num[index] ++;
            }
        }

        reset_p2p_stage_num(best_p2p_stage_num);
        MPI_Barrier(intra_comm);
        wtime(&tstart);
        execute(send_data, recv_data, fields_per_cell);
        wtime(&tend);
        temp_timer = tend-tstart;
        MPI_Allreduce(&temp_timer, &max_timer, 1, MPI_DOUBLE, MPI_MAX, intra_comm);
        if(max_timer < best_timer) best_timer = max_timer;
        else best_p2p_stage_mask[i] = true;

    }

    memset(best_p2p_stage_num, 0, sizeof(int)*butterfly_stage_num);
    index = 0;
    best_p2p_stage_num[index] = 1;
    for(int i=1; i<butterfly_stage_num; i++)
    {
        if(best_p2p_stage_mask[i])
        {
            index ++;
            best_p2p_stage_num[index] = 1;
        }
        else
            best_p2p_stage_num[index] ++;
    }

    reset_p2p_stage_num(best_p2p_stage_num);

    return p2p_num_per_butterfly_stage;
}
