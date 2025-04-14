#include <vector>
#include <map>
#include <algorithm>
#include <utility>
#include <iostream>
#include <mpi.h>
#include <cassert>
#include "functions.h"

void spgemm_2d(int m, int p, int n,
               std::vector<std::pair<std::pair<int,int>, int>> &A,
               std::vector<std::pair<std::pair<int,int>, int>> &B,
               std::vector<std::pair<std::pair<int,int>, int>> &C,
               std::function<int(int, int)> plus, std::function<int(int, int)> times,
               MPI_Comm row_comm, MPI_Comm col_comm)
{
    // TODO: Write your code here
    int row_rank, row_size;
    MPI_Comm_rank(row_comm,&row_rank);
    MPI_Comm_size(row_comm,&row_size);
    int column_rank, column_size;
    MPI_Comm_rank(col_comm,&column_rank);
    MPI_Comm_size(col_comm,&column_size);



    std::map<std::pair<int, int>, int> local_C;
    std::vector<std::pair<std::pair<int,int>, int>> A_i;
    std::vector<std::pair<std::pair<int,int>, int>> B_i;
    //grid dimension the same as the row_size
    for (int i=0;i<row_size;i++)
    {   
        
            if (row_rank==i)
            {
                A_i=A;
            }

            int A_count = (row_rank == i) ? A_i.size() : 0; 
            //bcast the size
            MPI_Bcast(&A_count, 1, MPI_INT, i, row_comm);
            
            if (row_rank != i) {
                A_i.resize(A_count);
            }
            MPI_Bcast(A_i.data(), A_count * sizeof(A[0]), MPI_BYTE, i, row_comm);


            
            if (column_rank==i)
            {
                B_i=B;
            }

            int B_count = (column_rank == i) ? B_i.size() : 0; 
            //bcast the size
            MPI_Bcast(&B_count, 1, MPI_INT, i, col_comm);
            
            if (column_rank != i) {
                B_i.resize(B_count);
            }
            MPI_Bcast(B_i.data(), B_count * sizeof(B[0]), MPI_BYTE, i, col_comm);

        std::unordered_map<int, std::vector<std::pair<int, int>>> B_hash;
        for (const auto& b : B_i) {
            int brow = b.first.first;
            int bcol = b.first.second;
            int bval = b.second;
            B_hash[brow].emplace_back(bcol, bval);
        }

        // Multiply A_i with hashed B_i
        for (const auto& a : A_i) {
            int arow = a.first.first;
            int acol = a.first.second;
            int aval = a.second;

            auto it = B_hash.find(acol);
            if (it != B_hash.end()) {
                for (const auto& [bcol, bval] : it->second) {
                    std::pair<int, int> pos(arow, bcol);
                    int prod = times(aval, bval);
                    if (local_C.find(pos) == local_C.end())
                        local_C[pos] = prod;
                    else
                        local_C[pos] = plus(local_C[pos], prod);
                }
            }
        }
    }
    C.clear();
    for (const auto &entry : local_C) {
        C.push_back({entry.first, entry.second});
    }
}
