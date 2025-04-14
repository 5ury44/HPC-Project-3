#include <vector>
#include <map>
#include <unordered_map>
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
    // Get the ranks and sizes from the communicator objects.
    int rowRank, numRows;
    MPI_Comm_rank(row_comm, &rowRank);
    MPI_Comm_size(row_comm, &numRows);

    int colRank, numCols;
    MPI_Comm_rank(col_comm, &colRank);
    MPI_Comm_size(col_comm, &numCols);

    // Map to accumulate local results, using a key (i, j) for position and value as the accumulated product.
    std::map<std::pair<int, int>, int> resultAccum;

    // Temporary containers for the broadcast pieces of A and B.
    std::vector<std::pair<std::pair<int,int>, int>> procA;
    std::vector<std::pair<std::pair<int,int>, int>> procB;

    // Each processor loop: for each step in the grid, broadcast the required portions.
    for (int step = 0; step < numRows; step++) {
        // --- For matrix A ---
        // The processor in the current row (step) provides its A data.
        if (rowRank == step) {
            procA = A;  // this process owns the current block of A
        }
        // Determine how many elements to expect.
        int numAElements = (rowRank == step) ? static_cast<int>(procA.size()) : 0;
        MPI_Bcast(&numAElements, 1, MPI_INT, step, row_comm);
        // Allocate space if not the owner.
        if (rowRank != step) {
            procA.resize(numAElements);
        }
        // Broadcast the actual data from the row communicator.
        MPI_Bcast(procA.data(), numAElements * sizeof(procA[0]), MPI_BYTE, step, row_comm);

        // --- For matrix B ---
        // Similarly, the processor in the current column (step) will share its B block.
        if (colRank == step) {
            procB = B;
        }
        int numBElements = (colRank == step) ? static_cast<int>(procB.size()) : 0;
        MPI_Bcast(&numBElements, 1, MPI_INT, step, col_comm);
        if (colRank != step) {
            procB.resize(numBElements);
        }
        MPI_Bcast(procB.data(), numBElements * sizeof(procB[0]), MPI_BYTE, step, col_comm);

        // --- Multiply the broadcast pieces ---
        // Create a hash map for the current procB using its "row index" as the key.
        std::unordered_map<int, std::vector<std::pair<int, int>>> hashedB;
        for (const auto &entryB : procB) {
            int bRowIndex = entryB.first.first;
            int bColIndex = entryB.first.second;
            int bValue = entryB.second;
            hashedB[bRowIndex].push_back({bColIndex, bValue});
        }

        // For every nonzero element in procA, look for matching entries in hashedB.
        for (const auto &entryA : procA) {
            int aRow = entryA.first.first;
            int aCol = entryA.first.second;
            int aVal = entryA.second;

            // Find if there is a matching row in B that corresponds to a's column.
            auto iter = hashedB.find(aCol);
            if (iter != hashedB.end()) {
                // Multiply with each matching element from B.
                for (const auto &pairB : iter->second) {
                    int bCol = pairB.first;
                    int bVal = pairB.second;
                    int prod = times(aVal, bVal);
                    std::pair<int, int> posKey = {aRow, bCol};
                    // If an entry already exists, combine them using the semiring plus.
                    if (resultAccum.find(posKey) == resultAccum.end()) {
                        resultAccum[posKey] = prod;
                    } else {
                        resultAccum[posKey] = plus(resultAccum[posKey], prod);
                    }
                }
            }
        }
    }

    // Finally, convert the accumulated result into the COO formatted output vector.
    C.clear();
    for (const auto &entry : resultAccum) {
        C.push_back({entry.first, entry.second});
    }
}
