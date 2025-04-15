#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include <iostream>
#include <mpi.h>
#include <cassert>
#include <functional>

struct pair {
    std::size_t operator()(const std::pair<int, int> &p) const {
        return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
    }
};

// Serializes a vector of COO entries as bytes (not portable, but matches original design).
void bc_coo(std::vector<std::pair<std::pair<int,int>, int>> &block, int root, MPI_Comm comm) {
    int sz = static_cast<int>(block.size());
    MPI_Bcast(&sz, 1, MPI_INT, root, comm);

    if (static_cast<int>(block.size()) != sz) {
        block.resize(sz);
    }

    MPI_Bcast(block.data(), sz * sizeof(block[0]), MPI_BYTE, root, comm);
}

// For each (i,k,a) in A and (k,j,b) in B, compute c = a*b and add to C(i,j) via plus.
void local_spgemm(
    const std::vector<std::pair<std::pair<int,int>, int>> &procA,
    const std::vector<std::pair<std::pair<int,int>, int>> &procB,
    std::unordered_map<std::pair<int, int>, int, pair> &resultAccum,
    std::function<int(int, int)> plus,
    std::function<int(int, int)> times)
{
    // Hash procB by row index (i.e., B(k, j)) so we can quickly match with A(i, k)
    std::unordered_map<int, std::vector<std::pair<int, int>>> hashedB;
    for (const auto &entryB : procB) {
        int bRowIndex = entryB.first.first;
        int bColIndex = entryB.first.second;
        int bValue = entryB.second;
        hashedB[bRowIndex].emplace_back(bColIndex, bValue);
    }

    // For every nonzero element in procA, look for matching entries in hashedB
    for (const auto &entryA : procA) {
        int aRow = entryA.first.first;
        int aCol = entryA.first.second;
        int aVal = entryA.second;

        auto iter = hashedB.find(aCol);  // match on column of A and row of B
        if (iter != hashedB.end()) {
            for (const auto &bPair : iter->second) {
                int bCol = bPair.first;
                int bVal = bPair.second;
                int prod = times(aVal, bVal);  // multiply using semiring 'times'

                std::pair<int, int> resultPos = {aRow, bCol};

                // Accumulate using semiring 'plus'
                auto resultIter = resultAccum.find(resultPos);
                if (resultIter == resultAccum.end()) {
                    resultAccum[resultPos] = prod;
                } else {
                    resultAccum[resultPos] = plus(resultIter->second, prod);
                }
            }
        }
    }
}

void spgemm_2d(int m, int p, int n,
               std::vector<std::pair<std::pair<int,int>, int>> &A,
               std::vector<std::pair<std::pair<int,int>, int>> &B,
               std::vector<std::pair<std::pair<int,int>, int>> &C,
               std::function<int(int, int)> plus, std::function<int(int, int)> times,
               MPI_Comm row_comm, MPI_Comm col_comm)
{
    // Get ranks and sizes
    int rowRank, colRank;
    MPI_Comm_rank(row_comm, &rowRank);
    MPI_Comm_rank(col_comm, &colRank);
    

    int numRows, numCols;
    MPI_Comm_size(row_comm, &numRows);
    MPI_Comm_size(col_comm, &numCols);

    // Unordered map to accumulate local results, using a key (i, j) and value as the product sum.
    std::unordered_map<std::pair<int, int>, int, pair> resultAccum;

    // Temporary containers for the broadcast pieces of A and B.
    std::vector<std::pair<std::pair<int,int>, int>> procA;
    std::vector<std::pair<std::pair<int,int>, int>> procB;

    for (int step = 0; step < numRows; step++) {
        // The processor in the current row (step) provides its A data.
        procA = (rowRank == step) ? A : std::vector<std::pair<std::pair<int,int>, int>>();
        bc_coo(procA, step, row_comm);

        // The processor in the current column (step) provides its B data.
        procB = (colRank == step) ? B : std::vector<std::pair<std::pair<int,int>, int>>();
        bc_coo(procB, step, col_comm);

        // --- Multiply the broadcast pieces and accumulate results ---
        local_spgemm(procA, procB, resultAccum, plus, times);
    }

    // Convert the accumulated result into the output COO vector.
    C.clear();
    for (const auto &entry : resultAccum) {
        C.emplace_back(entry.first, entry.second);
    }
}
