#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <string>
#include <sstream>
#include "functions.h"

void distribute_matrix_2d(int m, int n, std::vector<std::pair<std::pair<int, int>, int>> &full_matrix,
                          std::vector<std::pair<std::pair<int, int>, int>> &local_matrix,
                          int root, MPI_Comm comm_2d)
{
    int rank, size,coords_rank[2];
    MPI_Comm_rank(comm_2d,&rank);
    MPI_Comm_size(comm_2d,&size);

    //find the dimension of the 2d topology
    int gridsize=std::sqrt(size);

    int block_size_row = m/gridsize;
    int block_size_column = n/gridsize;


    if (rank==root)
    {
        // sending data from root to other ranks
        for (int i=0;i<size;i++){
            std::vector<std::pair<std::pair<int, int>, int>> send_buffer;
            
            MPI_Cart_coords(comm_2d,i,2,coords_rank); //find the destination coordinates

            //decide which block to send
            int row_start = block_size_row * coords_rank[0];
            int row_end= (coords_rank[0] == gridsize - 1) ? m : (row_start + block_size_row);
            int column_start = block_size_column*coords_rank[1];
            int column_end = (coords_rank[1] == gridsize - 1) ? n : (column_start + block_size_column);

            // 
            for (auto &elem: full_matrix){
                int row=elem.first.first;
                int column=elem.first.second;
                if ( row>=row_start && row < row_end && column>=column_start && column<column_end)
                {
                    send_buffer.push_back(elem);
                }
            }

            //sending data if the node is root, then directly distribute to local_matrix
            if (i==root){
                local_matrix=send_buffer;
            }
            else{
                int buffer_size=send_buffer.size();
                MPI_Send(&buffer_size,1,MPI_INT,i,0,comm_2d);
                MPI_Send(send_buffer.data(),buffer_size*sizeof(send_buffer[0]),MPI_BYTE,i,1,comm_2d);
            }
            
        }
    }
    else
    {
            //receving if not root node
            int recv_size;
            MPI_Recv(&recv_size, 1, MPI_INT, root, 0, comm_2d, MPI_STATUS_IGNORE);
            local_matrix.resize(recv_size);
            MPI_Recv(local_matrix.data(),recv_size*sizeof(local_matrix[0]),MPI_BYTE, root, 1, comm_2d, MPI_STATUS_IGNORE);
    }
}
