#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <algorithm>
#include <vector>
#include <map>
#include <iostream>
#include <set>
#include <string>
#include <fstream>
#include <sstream>
#include <cctype>

template<typename numT>
class Matrix {
    public:
    using idxT = size_t;

    std::vector<std::map<idxT, numT>> mat;
    idxT rows;
    idxT cols;
    std::vector<idxT> nonzero_rows;
    std::vector<idxT> nonzero_cols;

    Matrix() {};

    // todo copy constructor

    Matrix(std::vector<std::vector<numT>> const& vec_mat) {
        rows = vec_mat.size();
        cols = vec_mat.at(0).size();

        for (idxT i = 0; i < rows; i++) {
        std::map<idxT, numT> cur_row;
        for (idxT j = 0; j < cols; j++) {
            cur_row.insert(std::make_pair(j, vec_mat.at(i).at(j)));
        }
        if (cur_row.size() > 0) {
            nonzero_rows.push_back(i);
        }
        mat.push_back(cur_row);
    }

      for (idxT j = 0; j < cols; j++) {
          for (idxT i = 0; i < rows; i++) {
              if (mat.at(i).find(j) != mat.at(i).end()) {
                  nonzero_cols.push_back(j);
                  break;
              }
          }
      }
  };
    
    Matrix(std::map<std::pair<idxT, idxT>, numT> const& map_mat, std::pair<idxT, idxT> mat_size) {
        rows = mat_size.first;
        cols = mat_size.second;
        mat.resize(rows);
        std::set<idxT> nonzero_row_set;
        std::set<idxT> nonzero_col_set;
        
        for (const auto& key_value : map_mat) {
            idxT row = key_value.first.first;
            idxT col = key_value.first.second;
            numT val = key_value.second;
            
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

    // source: https://github.com/karlrupp/spgemm-mkl-benchmark/blob/master/matrix_market.hpp
    /** Trims whitespace from the beginning and the end of the string. */
    inline void my_trim(char * buffer, long max_size) {
        //trim at beginning of string
        long start = 0;
        for (long i=0; i<max_size; ++i)
        {
            if (buffer[i] == ' ')
            ++start;
            else
            break;
        }

        //trim at end of string
        long stop = start;
        for (long i=stop; i<max_size; ++i)
        {
            if (buffer[i] == 0)   //end of string
            break;

            if (buffer[i] != ' ')
            stop = i;
        }

        for (long i=0; i<=stop - start; ++i)
        {
            buffer[i] = buffer[start + i];
        }

        if (buffer[0] != ' ')
            buffer[stop - start + 1] = 0; //terminate string
        else
            buffer[0] = 0;
    }

    // source: https://github.com/karlrupp/spgemm-mkl-benchmark/blob/master/matrix_market.hpp
    /** Converts a string to lowercase */
    inline std::string my_tolower(std::string & s) {
        std::transform(s.begin(), s.end(), s.begin(), static_cast < int(*)(int) > (std::tolower));
        return s;
    }

    // source: https://github.com/karlrupp/spgemm-mkl-benchmark/blob/master/matrix_market.hpp
    /** Reads a matrix from file. Returns the number of lines read, or 0 in case of an error. */
    long read_matrix_market_file(const char * file) {
        char buffer[2049];
        std::ifstream reader(file);
        std::string token;
        long linenum = 0;
        bool symmetric = false;
        bool dense_format = false;
        bool is_header = true;
        bool pattern_matrix = false;
        idxT cur_row = 0;
        idxT cur_col = 0;
        long valid_entries = 0;
        long nnz = 0;

        if (!reader) {
            std::cerr << "Matrix Market Reader: Cannot open file " << file << std::endl;
            return EXIT_FAILURE;
        }

        while (reader.good()) {
            // get a non-empty line
            do {
                reader.getline(buffer, 2048);
                ++linenum;
                my_trim(buffer, 2048);
            } while (reader.good() && buffer[0] == 0);

            if (buffer[0] == '%') {
                if (buffer[1] == '%') {
                    //parse header:
                    std::stringstream line(std::string(buffer + 2));
                    line >> token;
                    if (my_tolower(token) != "matrixmarket") {
                        std::cerr << "Error in file " << file << " at line " << linenum << " in file " << file << ": Expected 'MatrixMarket', got '" << token << "'" << std::endl;
                        return 0;
                    }

                    line >> token;
                    if (my_tolower(token) != "matrix") {
                        std::cerr << "Error in file " << file << " at line " << linenum << " in file " << file << ": Expected 'matrix', got '" << token << "'" << std::endl;
                        return 0;
                    }

                    line >> token;
                    if (my_tolower(token) != "coordinate") {
                        if (my_tolower(token) == "array") {
                            dense_format = true;
                            std::cerr << "Error in file " << file << " at line " << linenum << " in file " << file << ": 'array' type is not supported yet!" << std::endl;
                            return 0;
                        } else {
                            std::cerr << "Error in file " << file << " at line " << linenum << " in file " << file << ": Expected 'array' or 'coordinate', got '" << token << "'" << std::endl;
                            return 0;
                        }
                    }

                    line >> token;
                    if (my_tolower(token) == "pattern") {
                        pattern_matrix = true;
                    } else if (my_tolower(token) == "complex") {
                        //is_complex = true;
                    } else if (my_tolower(token) != "real" && my_tolower(token) != "integer") {
                        std::cerr << "Error in file " << file << ": The MatrixMarket reader supports only real valued floating point arithmetic or pattern type matrices." << std::endl;
                        return 0;
                    }

                    line >> token;
                    if (my_tolower(token) == "general") {

                    } else if (my_tolower(token) == "symmetric") {
                        symmetric = true;
                    } else {
                        std::cerr << "Error in file " << file << ": The MatrixMarket reader supports only general or symmetric matrices." << std::endl;
                        return 0;
                    }
                }
            } else {
                std::stringstream line(std::stringstream::in | std::stringstream::out);
                line << std::string(buffer);

                if (is_header) {
                    //read header line
                    if (line.good()) {
                        line >> rows;
                    } else {
                        std::cerr << "Error in file " << file << ": Could not get matrix dimensions (rows) in line " << linenum << std::endl;
                        return 0;
                    }

                    if (line.good()) {
                        line >> cols;
                    } else {
                        std::cerr << "Error in file " << file << ": Could not get matrix dimensions (columns) in line " << linenum << std::endl;
                        return 0;
                    }

                    if (!dense_format) {
                        if (line.good()) {
                            line >> nnz;
                        } else {
                            std::cerr << "Error in file " << file << ": Could not get matrix dimensions (columns) in line " << linenum << std::endl;
                            return 0;
                        }
                    }

                    if (rows > 0 && cols > 0)
                        mat.resize(rows);

                    is_header = false;
                } else {
                    //read data
                    if (dense_format) {
                        numT value;
                        line >> value;
                        mat[cur_row][cur_col] = value;

                        if (++cur_row == mat.size()) {
                            //next column
                            ++cur_col;
                            cur_row = 0;
                        }
                    } else { //sparse format
                        idxT row;
                        idxT col;
                        numT value = numT(1);

                        //parse data:
                        if (line.good()) {
                            line >> row;
                        } else {
                            std::cerr << "Error in file " << file << ": Parse error for matrix row entry in line " << linenum << std::endl;
                            return 0;
                        }

                        if (line.good()) {
                            line >> col;
                        } else {
                            std::cerr << "Error in file " << file << ": Parse error for matrix col entry in line " << linenum << std::endl;
                            return 0;
                        }

                        //take index base 1 into account:
                        row -= 1;
                        col -= 1;

                        if (!pattern_matrix) { // value for pattern matrix is implicitly 1, so we only need to read data for 'normal' matrices
                            if (line.good()) {
                                line >> value;
                            } else {
                                std::cerr << "Error in file " << file << ": Parse error for matrix entry in line " << linenum << std::endl;
                                return 0;
                            }
                        }

                        if (row >= mat.size() || row < 0) {
                            std::cerr << "Error in file " << file << " at line " << linenum << ": Row index out of bounds: " << row << " (matrix dim: " << mat.size() << " x " << mat.size() << ")" << std::endl;
                            return 0;
                        }

                        if (col >= mat.size() || col < 0) {
                            std::cerr << "Error in file " << file << " at line " << linenum << ": Column index out of bounds: " << col << " (matrix dim: " << mat.size() << " x " << mat.size() << ")" << std::endl;
                            return 0;
                        }

                        mat[row][col] = value;
                        if (symmetric) {
                            mat[col][row] = value;
                        }

                        if (++valid_entries == nnz) {
                            break;
                        }
                    } //else dense_format
                }
            }
        }

        reader.close();

        std::set<idxT> nonzero_col_set;
        for (idxT i = 0; i < rows; i++) {
            if (mat[i].size() > 0) {
                nonzero_rows.push_back(i);
            }
            for (idxT j = 0; j < cols; j++) {
                if (nonzero_col_set.find(j) == nonzero_col_set.end() && mat[i].find(j) != mat[i].end()) {
                    nonzero_col_set.insert(j);
                }
            }
        }
        nonzero_cols.resize(nonzero_col_set.size());
        std::copy(nonzero_col_set.begin(), nonzero_col_set.end(), nonzero_cols.begin());
        std::sort(nonzero_rows.begin(), nonzero_rows.end());
        std::sort(nonzero_cols.begin(), nonzero_cols.end());

        return linenum;
    };
};

template<typename numT>
std::ostream& operator<< (std::ostream &out, Matrix<numT> const& matrix) {
    for (auto row : matrix.mat) {
        for (size_t col = 0; col < matrix.cols; col++) {
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
