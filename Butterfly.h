#ifndef BUTTERFLY_H
#define BUTTERFLY_H
#include "mpi.h"
#include <vector>
#include <iostream>
using namespace std;

#define SEND 0
#define RECV 1
#define SENDRECV 2

class Butterfly
{
private:
    MPI_Comm global_comm;
    MPI_Comm local_comm;
    int * send_cells;
    int * recv_cells;
    int remote_model_size;
    int action;

    int local_model_size;
    int * recv_model_global_rank;
    int * send_model_global_rank;
    int global_rank;
    int local_rank;
    bool is_send;
    bool is_recv;

    int send_model_size;
    int recv_model_size;
    int send_node_index;
    int recv_node_index;
    vector<int> * send_node_per_proc_local_rank;
    vector<int> * recv_node_per_proc_local_rank;
    int * total_node_root_proc_global_rank;
    int send_node_num;
    int recv_node_num;
    int total_node_num;

    int * total_send_cells;
    int * total_recv_cells;

    bool send_node_master;
    bool recv_node_master;

    int fields_per_cell;
    int pre_fields_per_cell;

    int max_cells_cnts;

    int intra_node_recv_cells_cnts;
    int butterfly_node_num;
    int butterfly_grp_node_num;
    int * butterfly_node_global_rank;
    int send_node_map_to_butterfly_node;
    int send_node_send_cells_cnts;
    vector<int> * butterfly_node_map_to_send_node;
    vector<int> * butterfly_node_recv_cells_cnts;
    bool butterfly_node_master;
    int butterfly_node_index;
    int * send_node_cells_cnts;

    vector<int> * butterfly_node_map_to_recv_node;
    vector<int> * recv_node_map_to_butterfly_node;

    char * butterfly_data_buf[3];
    int butterfly_stage_num;
    int butterfly_grp_index;
    int butterfly_grp_node_index;
    int * recv_proc_node_index;
    vector<int> * butterfly_node_decompsitions;
    vector<int> * butterfly_node_decompsitions_cells;
    int * butterfly_pair_node_list;

    int butterfly_grp_num;
    char * butterfly_gather_buf;

    int * intra_proc_displs_cells;

    FILE * file;
    int * last_p2p_send_model_proc_lists;
    int last_p2p_stage_num;
    vector<int> last_p2p_send_to_remote_procs;
    vector<int> last_p2p_recv_from_remote_procs;
    int * last_p2p_send_cells;
    int * last_p2p_send_displs;
    int * last_p2p_recv_cells;
    char * last_p2p_recv_buf;
    char * last_p2p_send_buf;
    int last_p2p_recv_buf_cells;
    vector<int> last_p2p_send_mask;
    vector<int> last_p2p_send_mask_cells;

    int butterfly_p2p_stage_num;
    int p2p_depth;
    int * butterfly_p2p_remote_procs;
    int * butterfly_p2p_send_cells;
    int * butterfly_p2p_recv_cells;
    int * butterfly_p2p_send_displs;
    int * butterfly_p2p_recv_displs;
    vector<int> * butterfly_p2p_send_mask;
    vector<int> * butterfly_p2p_send_mask_cells;

    int * p2p_num_per_butterfly_stage;
    int * p2p_depth_per_butterfly_stage;
    
    int cur_buf_index;
    int local_buf_index;
    int remote_buf_index;

    int send_node_map_to_butterfly_node_index;
    int * first_p2p_sendcells;
    int * first_p2p_recvcells;
    int * first_p2p_senddispls;
    int * first_p2p_recvdispls;
    int first_p2p_send_proc_num;
    int first_p2p_recv_proc_num;
    int * first_p2p_send_mask;
    int * first_p2p_send_mask_cells;
    int * first_p2p_send_proc_index;
    int * first_p2p_recv_proc_index;
    
    void compute_butterfly_recv_rank(int *, int, int);
    void intra_node_gather_init(MPI_Comm, MPI_Comm, int *, int *, int *, int, int);
    void process_mapping_from_send_to_butterfly_init();
    void process_mapping_from_butterfly_to_recv_init();
    void decrease_sort(int *, int *, int);
    void butterfly_init();
    void last_p2p_init();
    void first_p2p_init();
    void wtime(double *);
    void butterfly_p2p_init(int * = NULL);

public:
    void butterfly_p2p_execute();
    void first_p2p_execute(char *, char *);
    void last_p2p_execute(char *);
    void execute(char *, char *, int);
    void reset_p2p_stage_num(int *);
    int * auto_set_p2p_stage_num(char *, char *, int);
    Butterfly(MPI_Comm, MPI_Comm, int *, int *, int *, int, int);
    ~Butterfly();
};
#endif
