#include "matrix.hpp"
#include "metastrider.hpp"

#include <utility>
#include <queue>
#include <vector>
#include <iostream>
#include <iterator>
#include <optional>
#include <ostream>

template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1, T2>& pair) {
    os << "(" << pair.first << ", " << pair.second << ")";
    return os;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: ./spgemm [Matrix market file 1] [Matrix market file 2]" << std::endl;
        return EXIT_FAILURE;
    }

    using numT = double;
    using idxT = size_t;
    using keyT = std::pair<idxT, idxT>;

    Matrix<numT> a;
    Matrix<numT> b;

    //TODO multiple file input

    a.read_matrix_market_file(argv[1]);
    b.read_matrix_market_file(argv[1]);

    std::cout << a << std::endl;
    std::cout << "\n------------------------------------\n";
    std::cout << b << std::endl;
    std::cout << "\n------------------------------------\n";
    for (idxT i = 0; i < b.nonzero_rows.size(); i++) {
        std::cout << b.nonzero_rows.at(i) << " ";
    }
    std::cout << std::endl;
    std::cout << "\n------------------------------------\n";
    for (idxT i = 0; i < a.nonzero_cols.size(); i++) {
        std::cout << a.nonzero_cols.at(i) << " ";
    }
    std::cout << std::endl;
    std::cout << "\n------------------------------------\n";
    auto nonzero_col_a_iter = a.nonzero_cols.begin();
    auto nonzero_row_b_iter = b.nonzero_rows.begin();
    

    // minheap
    
    // using pqT = std::pair<keyT, numT>;
    std::priority_queue<keyT, std::vector<keyT>, std::greater<keyT>> pq;
    std::map<keyT, numT> pq_map;

    Metastrider<keyT, numT> metastrider;

    auto reduce_func = [=](numT num1, numT num2) {
        return num1 + num2;
    };

    metastrider.init(1, reduce_func);

    while (nonzero_col_a_iter != a.nonzero_cols.end() && nonzero_row_b_iter != b.nonzero_rows.end()) {
        idxT a_col = *nonzero_col_a_iter;
        idxT b_row = *nonzero_row_b_iter;
        if (a_col == b_row) {
            for (idxT a_row = 0; a_row < a.rows; a_row++) {
                for (idxT b_col = 0; b_col < b.cols; b_col++) {
                    {
                        numT val = a.mat[a_row][a_col] * b.mat[b_row][b_col];
                        if (val == 0) {
                            continue;
                        }
                        // std::cout << "row: " << a_row << " col: " << b_col << " val: " << val << std::endl;
                        
                        keyT key = std::make_pair(a_row, b_col);
                        auto key_find = pq_map.find(key);
                        if (key_find != pq_map.end()) {
                            val = reduce_func(val, key_find->second);
                        } else {
                            pq.push(key);
                        }
                        
                        pq_map[key] = val;
                    }

                    if (pq.size() >= metastrider.K) {
                        while (pq.size() > 0) {
                            auto key_find = pq_map.find(pq.top());
                            assert(key_find != pq_map.end());
                            auto pair = std::make_pair(key_find->first, key_find->second);
                            std::cout << "insert " << pair << std::endl;
                            metastrider.enq_reduce(pair);
                            pq.pop();
                        }
                        assert(pq.size() == 0);
                        pq_map.clear();
                        metastrider.add_vec();
                        metastrider.print_data("post add_vec");
                        metastrider.validate();
                        
                    }
                }
            }
    
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
            auto key_find = pq_map.find(pq.top());
            assert(key_find != pq_map.end());
            auto pair = std::make_pair(key_find->first, key_find->second);
            std::cout << "insert " << pair << std::endl;
            metastrider.enq_reduce(pair);
            pq.pop();
        }
        metastrider.add_vec();
        metastrider.print_data("post final add_vec");
        metastrider.validate();
    }

    metastrider.global_reduce();
    metastrider.print_data("post global_reduce");
    metastrider.validate(true); // strict validation

    std::map<keyT, numT> mat;
    for (std::optional<std::vector<std::pair<keyT, numT>>> data_node : metastrider.data) {
        if (data_node.has_value()) {
            for (std::pair<keyT, numT> data_pair : data_node.value()) {
                assert(mat.find(data_pair.first) == mat.end());
                mat.emplace(data_pair);
            }
        }
    }
    
    Matrix<numT> c(mat, std::make_pair(a.rows, b.cols));
    std::cout << c << std::endl;

    return 0;
}