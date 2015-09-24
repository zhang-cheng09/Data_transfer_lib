#ifndef DATA_TRANSFER_PACK_H
#define DATA_TRANSFER_PACK_H

#include "mpi.h"
#include <vector>
#include "Butterfly.h"

#define SEND 0
#define RECV 1
#define SENDRECV 2

struct Field_info
{
    void * data_buf;
    int buf_size;
    int data_type;
    bool input;
};

struct Routing_info_with_one_process
{
    int num_elements_transferred;
    int num_local_indx_segments;
    int *local_indx_segment_starts;
    int *local_indx_segment_lengths;
    bool send_or_recv;                           // true is send and false is recv
};

class Data_transfer_instance
{
    private:
        MPI_Comm local_comm;
        int remote_root_global_rank;
        int direction;
        bool * mask;
        int num_grid_cells;
        int num_local_cells; 
        int * local_cells_global_index;
        int num_remap_grid_cells;
        int num_remap_local_cells; 
        int * local_remap_cells_global_index;
        int remote_model_size;
        int * remote_procs_global_rank;
        int * send_cells;
        int * recv_cells;
        std::vector<Field_info *> coupling_fields;
        int num_input_fields;
        int num_output_fields;
        Routing_info_with_one_process * routing_info_mgt;
        Butterfly * butterfly_instance;
        char * send_data_buf;
        char * recv_data_buf;
        int total_send_cells;
        int total_recv_cells;
    public:
        Data_transfer_instance(MPI_Comm, int, int);
        ~Data_transfer_instance();
        bool register_decomposition(int, int, int *);
        bool register_field(void *, int, char *, bool);
        bool register_mask(bool *);
        bool init();
        bool exec();
        bool auto_set_p2p_stage_num();
        int * get_send_cells(){return send_cells;}
        int * get_recv_cells(){return recv_cells;}
        char * get_send_buf(){return send_data_buf;}
        char * get_recv_buf(){return recv_data_buf;}
        void * get_coupling_field(int field_id){return coupling_fields[field_id]->data_buf;}
    private:
        void build_2D_remote_router();
        void build_2D_self_router();
        void compute_routing_info_between_decomps(int, int *, int);
        void compute_remap_routing_info_between_decomps(int, int *, int);
        void build_data_transfer_info();
        void pack_MD_data();
        void unpack_MD_data();
};

class Data_transfer_instance_mgt
{
    private:
        std::vector<Data_transfer_instance *> data_transfer_instance_mgt;  
    public:
        Data_transfer_instance_mgt();
        ~Data_transfer_instance_mgt();
        int register_data_transfer_instance(MPI_Comm, int, int);
        bool register_decomposition(int, int, int, int *);
        bool register_field(int, void *, int, char *, bool);
        bool init_data_transfer_instance(int);
        bool exec_data_transfer_instance(int);
        bool auto_set_p2p_stage_num(int);
        bool final_data_transfer_instance(int);
        bool register_mask(int, bool *);
        int * get_send_cells(int instance_id){return data_transfer_instance_mgt[instance_id]->get_send_cells();}
        int * get_recv_cells(int instance_id){return data_transfer_instance_mgt[instance_id]->get_recv_cells();}
        char * get_send_buf(int instance_id){return data_transfer_instance_mgt[instance_id]->get_send_buf();}
        char * get_recv_buf(int instance_id){return data_transfer_instance_mgt[instance_id]->get_recv_buf();}
        void * get_coupling_field(int instance_id, int field_id){return data_transfer_instance_mgt[instance_id]->get_coupling_field(field_id);}
};
#endif
