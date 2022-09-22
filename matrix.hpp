#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <algorithm>
#include <vector>
#include <map>
#include <iostream>
#include <set>
#include <iostream>

template<typename IndexT, typename NumericT>
class Matrix {
    public:
    std::vector<std::map<IndexT, NumericT>> mat;
    IndexT rows;
    IndexT cols;
    std::vector<IndexT> nonzero_rows;
    std::vector<IndexT> nonzero_cols;

    Matrix(std::vector<std::vector<NumericT>> const& vec_mat) {
        rows = vec_mat.size();
        cols = vec_mat.at(0).size();

        for (IndexT i = 0; i < rows; i++) {
            std::map<IndexT, NumericT> cur_row;
            for (IndexT j = 0; j < cols; j++) {
                cur_row.insert(std::make_pair(j, vec_mat.at(i).at(j)));
            }
            if (cur_row.size() > 0) {
                nonzero_rows.push_back(i);
            }
            mat.push_back(cur_row);
        }

        for (IndexT j = 0; j < cols; j++) {
            for (IndexT i = 0; i < rows; i++) {
                if (mat.at(i).find(j) != mat.at(i).end()) {
                    nonzero_cols.push_back(j);
                    break;
                }
            }
        }
    };
    
    Matrix(std::map<std::pair<IndexT, IndexT>, NumericT> const& map_mat, std::pair<IndexT, IndexT> mat_size) {
        rows = mat_size.first;
        cols = mat_size.second;
        mat.resize(rows);
        std::set<IndexT> nonzero_row_set;
        std::set<IndexT> nonzero_col_set;
        
        for (const auto& key_value : map_mat) {
            IndexT row = key_value.first.first;
            IndexT col = key_value.first.second;
            NumericT val = key_value.second;
            
            mat[row][col] = val;
            nonzero_row_set.insert(row);
            nonzero_col_set.insert(col);
        }
        
        
        nonzero_rows.resize(nonzero_row_set.size());
        nonzero_cols.resize(nonzero_col_set.size());
        std::copy(nonzero_row_set.begin(), nonzero_row_set.end(), nonzero_rows.begin());
        std::copy(nonzero_col_set.begin(), nonzero_col_set.end(), nonzero_cols.begin());
        std::sort(nonzero_rows.begin(), nonzero_rows.end());
        std::sort(nonzero_cols.begin(), nonzero_cols.end());

    };
    // Matrix(const char *file);

    
};

template<typename IndexT, typename NumericT>
    std::ostream& operator<< (std::ostream &out, Matrix<IndexT, NumericT> const& matrix) {
        for (auto row : matrix.mat) {
            for (IndexT col = 0; col < matrix.cols; col++) {
                if (row.find(col) != row.end()) {
                    out << row[col] << "\t";
                } else {
                    out << "0\t";
                }
            }
            out << std::endl;
        }
        return out;
    };

#endif
