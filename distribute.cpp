// #include <algorithm>
// #include <vector>
// #include <iostream>
// #include <fstream>
// #include <cassert>
// #include <cmath>
// #include <string>
// #include <sstream>
// #include "functions.h"
// #include <mpi.h>

// void distribute_matrix_2d(int m, int n, 
//                           std::vector<std::pair<std::pair<int, int>, int>> &full_matrix,
//                           std::vector<std::pair<std::pair<int, int>, int>> &local_matrix,
//                           int root, MPI_Comm comm_2d)
// {
//     int myRank, szProcesses;
//     MPI_Comm_rank(comm_2d, &myRank);
//     MPI_Comm_size(comm_2d, &szProcesses);

//     // Determine grid dimensions given the number of processors (assuming a square grid)
//     int gridDim = static_cast<int>(std::sqrt(szProcesses));
//     int blockRows = m / gridDim;
//     int blockCols = n / gridDim;

//     if (myRank == root) {
//         // Create a vector of buckets corresponding to each processor in the grid.
//         std::vector<std::vector<std::pair<std::pair<int,int>, int>>> procBuckets(szProcesses);

//         // Distribute each entry in full_matrix to the appropriate bucket.
//         for (size_t i = 0; i < full_matrix.size(); i++) {
//             int curRow = full_matrix[i].first.first;
//             int curCol = full_matrix[i].first.second;

//             int coords[2] = {std::min(curRow / blockRows, gridDim - 1), std::min(curCol / blockCols, gridDim - 1)};
//             int targetRank = -1;
//             MPI_Cart_rank(comm_2d, coords, &targetRank);
            
//             procBuckets[targetRank].push_back(full_matrix[i]);
//         }

//         // Send data to each processor.
//         for (int p = 0; p < szProcesses; p++) {
//             if (p == root) {
//                 // Directly copy the bucket for the root processor.
//                 local_matrix = procBuckets[p];
//             } else {
//                 // First, send the number of items in this processor's bucket.
//                 int nItems = procBuckets[p].size();
//                 MPI_Send(&nItems, 1, MPI_INT, p, 0, comm_2d);
//                 // Then send the actual bucket data.
//                 MPI_Send(procBuckets[p].data(), 
//                          nItems * sizeof(procBuckets[p][0]), 
//                          MPI_BYTE, p, 1, comm_2d);
//             }
//         }
//     } else {
//         // Non-root processors receive the block size first and then the bucket.
//         int receivedCount = 0;
//         MPI_Recv(&receivedCount, 1, MPI_INT, root, 0, comm_2d, MPI_STATUS_IGNORE);
//         local_matrix.resize(receivedCount);
//         MPI_Recv(local_matrix.data(), 
//                  receivedCount * sizeof(local_matrix[0]), 
//                  MPI_BYTE, root, 1, comm_2d, MPI_STATUS_IGNORE);
//     }
// }

//---------------------JCCCCCCCCCCCC

// #include <algorithm>
// #include <vector>
// #include <iostream>
// #include <fstream>
// #include <cassert>
// #include <cmath>
// #include <string>
// #include <sstream>
// #include "functions.h"

// void distribute_matrix_2d(int m, int n, std::vector<std::pair<std::pair<int, int>, int>> &full_matrix,
//                           std::vector<std::pair<std::pair<int, int>, int>> &local_matrix,
//                           int root, MPI_Comm comm_2d)
// {
//     int rank, size,coords_rank[2];
//     MPI_Comm_rank(comm_2d,&rank);
//     MPI_Comm_size(comm_2d,&size);

//     //find the dimension of the 2d topology
//     int gridsize=std::sqrt(size);

//     int block_size_row = m/gridsize;
//     int block_size_column = n/gridsize;


//     if (rank==root)
//     {
//         // sending data from root to other ranks
//         for (int i=0;i<size;i++){
//             std::vector<std::pair<std::pair<int, int>, int>> send_buffer;
            
//             MPI_Cart_coords(comm_2d,i,2,coords_rank); //find the destination coordinates

//             //decide which block to send
//             int row_start = block_size_row * coords_rank[0];
//             int row_end= (coords_rank[0] == gridsize - 1) ? m : (row_start + block_size_row);
//             int column_start = block_size_column*coords_rank[1];
//             int column_end = (coords_rank[1] == gridsize - 1) ? n : (column_start + block_size_column);

//             // 
//             for (auto &elem: full_matrix){
//                 int row=elem.first.first;
//                 int column=elem.first.second;
//                 if ( row>=row_start && row < row_end && column>=column_start && column<column_end)
//                 {
//                     send_buffer.push_back(elem);
//                 }
//             }

//             //sending data if the node is root, then directly distribute to local_matrix
//             if (i==root){
//                 local_matrix=send_buffer;
//             }
//             else{
//                 int buffer_size=send_buffer.size();
//                 MPI_Send(&buffer_size,1,MPI_INT,i,0,comm_2d);
//                 MPI_Send(send_buffer.data(),buffer_size*sizeof(send_buffer[0]),MPI_BYTE,i,1,comm_2d);
//             }
            
//         }
//     }
//     else
//     {
//             //receving if not root node
//             int recv_size;
//             MPI_Recv(&recv_size, 1, MPI_INT, root, 0, comm_2d, MPI_STATUS_IGNORE);
//             local_matrix.resize(recv_size);
//             MPI_Recv(local_matrix.data(),recv_size*sizeof(local_matrix[0]),MPI_BYTE, root, 1, comm_2d, MPI_STATUS_IGNORE);
//     }
// }

//========================================FINALLLLLLLLLLLL

#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <string>
#include <sstream>
#include "functions.h"
#include <mpi.h>

void distribute_matrix_2d(int m, int n, 
                          std::vector<std::pair<std::pair<int, int>, int>> &full_matrix,
                          std::vector<std::pair<std::pair<int, int>, int>> &local_matrix,
                          int root, MPI_Comm comm_2d)
{
    int myRank, szProcesses;
    MPI_Comm_rank(comm_2d, &myRank);
    MPI_Comm_size(comm_2d, &szProcesses);

    // Determine grid dimensions given the number of processors (assuming a square grid)
    int gridDim = static_cast<int>(std::sqrt(szProcesses));
    int blockRows = m / gridDim;
    int blockCols = n / gridDim;

    if (myRank == root) {
        // Create a vector of buckets corresponding to each processor in the grid.
        std::vector<std::vector<std::pair<std::pair<int,int>, int>>> procBuckets(szProcesses);

        // Distribute each entry in full_matrix to the appropriate bucket.
        for (size_t i = 0; i < full_matrix.size(); i++) {
            int curRow = full_matrix[i].first.first;
            int curCol = full_matrix[i].first.second;
        
            // Safer handling of edge coordinates
            int coords[2];
            coords[0] = (curRow < (gridDim - 1) * blockRows) ? curRow / blockRows : gridDim - 1;
            coords[1] = (curCol < (gridDim - 1) * blockCols) ? curCol / blockCols : gridDim - 1;
        
            int targetRank = -1;
            MPI_Cart_rank(comm_2d, coords, &targetRank);
        
            procBuckets[targetRank].push_back(full_matrix[i]);
        }

        // Send data to each processor.
        for (int p = 0; p < szProcesses; p++) {
            if (p == root) {
                // Directly copy the bucket for the root processor.
                local_matrix = procBuckets[p];
            } else {
                // First, send the number of items in this processor's bucket.
                int nItems = procBuckets[p].size();
                MPI_Send(&nItems, 1, MPI_INT, p, 0, comm_2d);
                // Then send the actual bucket data.
                MPI_Send(procBuckets[p].data(), 
                         nItems * sizeof(procBuckets[p][0]), 
                         MPI_BYTE, p, 1, comm_2d);
            }
        }
    } else {
        // Non-root processors receive the block size first and then the bucket.
        int receivedCount = 0;
        MPI_Recv(&receivedCount, 1, MPI_INT, root, 0, comm_2d, MPI_STATUS_IGNORE);
        local_matrix.resize(receivedCount);
        MPI_Recv(local_matrix.data(), 
                 receivedCount * sizeof(local_matrix[0]), 
                 MPI_BYTE, root, 1, comm_2d, MPI_STATUS_IGNORE);
    }
}
