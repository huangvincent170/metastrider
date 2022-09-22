#include "matrix.hpp"
#include "metastrider.hpp"

#include <utility>
#include <queue>
#include <vector>
#include <iostream>
#include <iterator>

int main() {
    using idxT = int;
    using numT = int;
    Matrix<idxT, numT> a({
        {1, 0, 2, 0},
        {0, 0, 3, 0},
        {4, 0, 5, 6},
        {0, 0, 0, 0}
    });
    Matrix<idxT, numT> b({
        {6, 0, 2, 0},
        {0, 0, 3, 0},
        {4, 0, 2, 1},
        {0, 0, 0, 0}
    });

//     std::cout << a << std::endl;
    auto nonzero_col_a_iter = a.nonzero_cols.begin();
    auto nonzero_row_b_iter = b.nonzero_rows.begin();
    

    // minheap
    using keyT = std::pair<idxT, idxT>;
    using pqT = std::pair<keyT, numT>;
    std::priority_queue<pqT, std::vector<pqT>, std::greater<pqT>> pq;

    Metastrider<keyT, numT> metastrider(
        [=](numT num1, numT num2) {
            return num1 + num2;
        }
    );

    while (nonzero_col_a_iter != a.nonzero_cols.end() && nonzero_row_b_iter != b.nonzero_rows.end()) {
        idxT a_col = *nonzero_col_a_iter;
        idxT b_row = *nonzero_row_b_iter;
        if (a_col == b_row) {
            for (idxT a_row = 0; a_row < a.rows; a_row++) {
                for (idxT b_col = 0; b_col < b.cols; b_col++) {
                    numT val = a.mat[a_row][a_col] * b.mat[b_row][b_col];
//                     std::cout << "row: " << a_row << " col: " << b_col << " val: " << val << std::endl;
                    pq.push(std::make_pair(std::make_pair(a_row, b_col), val));

                    if (pq.size() >= metastrider.K) {
                        while (pq.size() > 0) {
                            metastrider.enq_reduce(pq.top());
                            pq.pop();
                        }
                        metastrider.add_vec();
                    }
                }
            }

            metastrider.global_reduce();
            // metastrider.lookup(); ?
            std::advance(nonzero_col_a_iter, 1);
            std::advance(nonzero_row_b_iter, 1);
        } else if (a_col < b_row) {
            std::advance(nonzero_col_a_iter, 1);
        } else { // a_col > b_row
            std::advance(nonzero_row_b_iter, 1);
        }
    }
    
    if (pq.size() > 0) {
        while (pq.size() > 0) {
            metastrider.enq_reduce(pq.top());
            pq.pop();
        }
        metastrider.add_vec();
        metastrider.global_reduce();
    }
    
    Matrix<idxT, numT> c(metastrider.data, std::make_pair(a.rows, b.cols));
    std::cout << c << std::endl;

    return 0;
}