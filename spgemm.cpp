// #include <vector>
// #include <map>
// #include <unordered_map>
// #include <algorithm>
// #include <utility>
// #include <iostream>
// #include <mpi.h>
// #include <cassert>
// #include <functional>

// struct pair {
//     std::size_t operator()(const std::pair<int, int> &p) const {
//         return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
//     }
// };

// // For each (i,k,a) in A and (k,j,b) in B, compute c = a*b and add to C(i,j) via plus.
// void spgemm(
//     const std::vector<std::pair<std::pair<int,int>, int>> &procA,
//     const std::vector<std::pair<std::pair<int,int>, int>> &procB,
//     std::unordered_map<std::pair<int, int>, int, pair> &resultAccum,
//     std::function<int(int, int)> plus,
//     std::function<int(int, int)> times)
// {
//     // Hash procB by row index (i.e., B(k, j)) so we can quickly match with A(i, k)
//     std::unordered_map<int, std::vector<std::pair<int, int>>> hashedB;
//     for (const auto &entryB : procB) {
//         int bRowIndex = entryB.first.first;
//         int bColIndex = entryB.first.second;
//         int bValue = entryB.second;
//         hashedB[bRowIndex].emplace_back(bColIndex, bValue);
//     }

//     // For every nonzero element in procA, look for matching entries in hashedB
//     for (const auto &entryA : procA) {
//         int aRow = entryA.first.first;
//         int aCol = entryA.first.second;
//         int aVal = entryA.second;

//         auto iter = hashedB.find(aCol);  // match on column of A and row of B
//         if (iter != hashedB.end()) {
//             for (const auto &bPair : iter->second) {
//                 int bCol = bPair.first;
//                 int bVal = bPair.second;
//                 int prod = times(aVal, bVal);  // multiply using semiring 'times'

//                 std::pair<int, int> resultPos = {aRow, bCol};

//                 // Accumulate using semiring 'plus'
//                 auto resultIter = resultAccum.find(resultPos);
//                 if (resultIter == resultAccum.end()) {
//                     resultAccum[resultPos] = prod;
//                 } else {
//                     resultAccum[resultPos] = plus(resultIter->second, prod);
//                 }
//             }
//         }
//     }
// }

// // Serializes a vector of COO entries as bytes (not portable, but matches original design).
// void bc_coo(std::vector<std::pair<std::pair<int,int>, int>> &block, int root, MPI_Comm comm) {
//     int sz = static_cast<int>(block.size());
//     MPI_Bcast(&sz, 1, MPI_INT, root, comm);

//     if (static_cast<int>(block.size()) != sz) {
//         block.resize(sz);
//     }

//     MPI_Bcast(block.data(), sz * sizeof(block[0]), MPI_BYTE, root, comm);
// }

// void spgemm_2d(int m, int p, int n,
//                std::vector<std::pair<std::pair<int,int>, int>> &A,
//                std::vector<std::pair<std::pair<int,int>, int>> &B,
//                std::vector<std::pair<std::pair<int,int>, int>> &C,
//                std::function<int(int, int)> plus, std::function<int(int, int)> times,
//                MPI_Comm row_comm, MPI_Comm col_comm)
// {
//     // Get ranks and sizes
//     int rowRank, colRank;
//     MPI_Comm_rank(row_comm, &rowRank);
//     MPI_Comm_rank(col_comm, &colRank);
    

//     int numRows, numCols;
//     MPI_Comm_size(row_comm, &numRows);
//     MPI_Comm_size(col_comm, &numCols);

//     // Unordered map to accumulate local results, using a key (i, j) and value as the product sum.
//     std::unordered_map<std::pair<int, int>, int, pair> resultAccum;

//     // Temporary containers for the broadcast pieces of A and B.
//     std::vector<std::pair<std::pair<int,int>, int>> procA;
//     std::vector<std::pair<std::pair<int,int>, int>> procB;

//     for (int step = 0; step < numRows; step++) {
//         // The processor in the current row (step) provides its A data.
//         procA = (rowRank == step) ? A : std::vector<std::pair<std::pair<int,int>, int>>();
//         bc_coo(procA, step, row_comm);

//         // The processor in the current column (step) provides its B data.
//         procB = (colRank == step) ? B : std::vector<std::pair<std::pair<int,int>, int>>();
//         bc_coo(procB, step, col_comm);

//         // --- Multiply the broadcast pieces and accumulate results ---
//         spgemm(procA, procB, resultAccum, plus, times);
//     }

//     // Convert the accumulated result into the output COO vector.
//     C.clear();
//     for (const auto &entry : resultAccum) {
//         C.emplace_back(entry.first, entry.second);
//     }
// }

//-------------------JCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC------------

// #include <vector>
// #include <map>
// #include <algorithm>
// #include <utility>
// #include <iostream>
// #include <mpi.h>
// #include <cassert>
// #include "functions.h"

// void spgemm_2d(int m, int p, int n,
//                std::vector<std::pair<std::pair<int,int>, int>> &A,
//                std::vector<std::pair<std::pair<int,int>, int>> &B,
//                std::vector<std::pair<std::pair<int,int>, int>> &C,
//                std::function<int(int, int)> plus, std::function<int(int, int)> times,
//                MPI_Comm row_comm, MPI_Comm col_comm)
// {
//     // TODO: Write your code here
//     int row_rank, row_size;
//     MPI_Comm_rank(row_comm,&row_rank);
//     MPI_Comm_size(row_comm,&row_size);
//     int column_rank, column_size;
//     MPI_Comm_rank(col_comm,&column_rank);
//     MPI_Comm_size(col_comm,&column_size);



//     std::map<std::pair<int, int>, int> local_C;
//     std::vector<std::pair<std::pair<int,int>, int>> A_i;
//     std::vector<std::pair<std::pair<int,int>, int>> B_i;
//     //grid dimension the same as the row_size
//     for (int i=0;i<row_size;i++)
//     {   
        
//             if (row_rank==i)
//             {
//                 A_i=A;
//             }

//             int A_count = (row_rank == i) ? A_i.size() : 0; 
//             //bcast the size
//             MPI_Bcast(&A_count, 1, MPI_INT, i, row_comm);
            
//             if (row_rank != i) {
//                 A_i.resize(A_count);
//             }
//             MPI_Bcast(A_i.data(), A_count * sizeof(A[0]), MPI_BYTE, i, row_comm);


            
//             if (column_rank==i)
//             {
//                 B_i=B;
//             }

//             int B_count = (column_rank == i) ? B_i.size() : 0; 
//             //bcast the size
//             MPI_Bcast(&B_count, 1, MPI_INT, i, col_comm);
            
//             if (column_rank != i) {
//                 B_i.resize(B_count);
//             }
//             MPI_Bcast(B_i.data(), B_count * sizeof(B[0]), MPI_BYTE, i, col_comm);

//         std::unordered_map<int, std::vector<std::pair<int, int>>> B_hash;
//         for (const auto& b : B_i) {
//             int brow = b.first.first;
//             int bcol = b.first.second;
//             int bval = b.second;
//             B_hash[brow].emplace_back(bcol, bval);
//         }

//         // Multiply A_i with hashed B_i
//         for (const auto& a : A_i) {
//             int arow = a.first.first;
//             int acol = a.first.second;
//             int aval = a.second;

//             auto it = B_hash.find(acol);
//             if (it != B_hash.end()) {
//                 for (const auto& [bcol, bval] : it->second) {
//                     std::pair<int, int> pos(arow, bcol);
//                     int prod = times(aval, bval);
//                     if (local_C.find(pos) == local_C.end())
//                         local_C[pos] = prod;
//                     else
//                         local_C[pos] = plus(local_C[pos], prod);
//                 }
//             }
//         }
//     }
//     C.clear();
//     for (const auto &entry : local_C) {
//         C.push_back({entry.first, entry.second});
//     }
// }

//=================================NEWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

#include <vector>
#include <map>
#include <unordered_map>
#include <utility>
#include <iostream>
#include <mpi.h>
#include <cassert>
#include <functional>

using MatrixEntry = std::pair<std::pair<int, int>, int>;
using SparseMatrix = std::vector<MatrixEntry>;

// Helper to broadcast a sparse matrix
void broadcast_matrix(SparseMatrix &matrix, int root, MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);

    int count = matrix.size();
    MPI_Bcast(&count, 1, MPI_INT, root, comm);

    if (rank != root) {
        matrix.resize(count);
    }
    MPI_Bcast(matrix.data(), count * sizeof(MatrixEntry), MPI_BYTE, root, comm);
}

// Hash B matrix rows by key to accelerate matching
std::unordered_map<int, std::vector<std::pair<int, int>>> hash_matrix_by_row(const SparseMatrix &B) {
    std::unordered_map<int, std::vector<std::pair<int, int>>> hash;
    for (const auto &entry : B) {
        int row = entry.first.first;
        int col = entry.first.second;
        int val = entry.second;
        hash[row].emplace_back(col, val);
    }
    return hash;
}

// Perform local multiplication between A and hashed B
void local_spgemm(const SparseMatrix &A,
                  const std::unordered_map<int, std::vector<std::pair<int, int>>> &B_hash,
                  std::map<std::pair<int, int>, int> &local_C,
                  std::function<int(int, int)> plus,
                  std::function<int(int, int)> times) {
    for (const auto &a : A) {
        int arow = a.first.first;
        int acol = a.first.second;
        int aval = a.second;

        auto it = B_hash.find(acol);
        if (it != B_hash.end()) {
            for (const auto &[bcol, bval] : it->second) {
                std::pair<int, int> pos(arow, bcol);
                int prod = times(aval, bval);
                if (local_C.find(pos) == local_C.end()) {
                    local_C[pos] = prod;
                } else {
                    local_C[pos] = plus(local_C[pos], prod);
                }
            }
        }
    }
}

// Main 2D SpGEMM driver
void spgemm_2d(int m, int p, int n,
               SparseMatrix &A,
               SparseMatrix &B,
               SparseMatrix &C,
               std::function<int(int, int)> plus,
               std::function<int(int, int)> times,
               MPI_Comm row_comm, MPI_Comm col_comm)
{
    int row_rank, row_size, col_rank, col_size;
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_size(row_comm, &row_size);
    MPI_Comm_rank(col_comm, &col_rank);
    MPI_Comm_size(col_comm, &col_size);

    std::map<std::pair<int, int>, int> local_C;

    for (int k = 0; k < row_size; ++k) {
        SparseMatrix A_i, B_i;

        if (row_rank == k) A_i = A;
        broadcast_matrix(A_i, k, row_comm);

        if (col_rank == k) B_i = B;
        broadcast_matrix(B_i, k, col_comm);

        auto B_hash = hash_matrix_by_row(B_i);
        local_spgemm(A_i, B_hash, local_C, plus, times);
    }

    C.clear();
    for (const auto &entry : local_C) {
        C.emplace_back(entry.first, entry.second);
    }
}
