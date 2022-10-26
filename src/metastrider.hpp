#include <vector>
#include <utility>
#include <map>
#include <functional>
#include <iostream>
#include <optional>
#include <cassert>

// size of rows in bytes
// assume 2KB row
#define ROW_SIZE 2048

#define AVL_IMPL

#ifdef AVL_IMPL
    #define get_left_child_idx(idx) ((2*idx)+1)
    #define get_right_child_idx(idx) ((2*idx)+2)
#endif

template<typename keyT, typename valT>
class Metastrider {

    using pairT = std::pair<keyT, valT>;
    using idxT = size_t;

    public:
    idxT K;
    std::vector<pairT> q;
    std::function<valT(valT, valT)> reduce;

#ifdef MAP_IMPL
    std::map<keyT, valT> data;
#endif

#ifdef AVL_IMPL
    class MetadataNode {
        public:
        keyT pivot;
        keyT min;
        keyT max;
        keyT maxLeft;
        keyT minRight;

        MetadataNode(const std::vector<pairT> &data) {
            min = data.at(0).first;
            max = data.at(data.size()-1).first;
            pivot = data.at(data.size()/2).first;
        }
    };

    std::vector<std::optional<std::vector<pairT>>> data;
    std::vector<std::optional<MetadataNode>> metadata;
#endif


    /**
     * Init metastrider
     * 
     * @param n     number of cores, currently unimplmented and ignored
     * @param func  The reduce function to be applied
     * @return uint K, the number of key-value pairs in a row
     */
    uint init(
        uint n,
        std::function<valT(valT, valT)> reduce_func
    ) {
        (void)n;
        K = ROW_SIZE / (sizeof(keyT) + sizeof(valT));
        q.reserve(K);
        reduce = reduce_func;
        return K;
    }

    /**
     * Enqueue a key-value pair to be reduced by Metastrider.
     * Key-value pairs should be enqueued in sorted order,
     * and at most K key-value pairs should be enqueued,
     * after which add_vec() must be called.
     * 
     * @param key_value 
     */
    void enq_reduce(pairT key_value) {
        q.push_back(key_value);
    }

    /**
     * Reduces and de-duplicates a node's data with input data.
     * 
     * @param node 
     * @param insert_data 
     * @return std::tuple<std::vector<pairT>, std::vector<pairT>, boolean> ret_val
     *     ret_val.get<0> the leftover data (size at most K) 
     *     ret_val.get<1> if the leftover data should be sent to the right node
     *                    (0 = send to left node, 1 = send to right node) 
     */
    std::optional<std::pair<std::vector<pairT>, bool>> reduce_and_dedup(
        MetadataNode &metadata_node,
        std::vector<pairT> &data_node,
        std::vector<pairT> &insert_data
    ) {
        if (insert_data.size() == 0) {
            return {};
        }

        auto node_data_iter = data_node.begin();
        auto insert_data_iter = insert_data.begin();
        std::vector<pairT> merged_data;
        merged_data.reserve(2 * K);

        // merge the existing key-value pairs with the inserting key-value pairs

        // number of elements in merged_data which are less than pivot
        idxT LT = 0;
        // number of elements in merged_data which are greater than or equal to pivot
        idxT GE = 0;
        while (node_data_iter != data_node.end() || insert_data_iter != insert_data.end()) {
            pairT next_data;
            // gets the next key-value pair in sorted order
            // will take from existing node data if it is smaller
            // will take from inserting data is it is smaller
            // will create new reduced value if existing key and inserting key are the same
            if (insert_data_iter == insert_data.end()) {
                next_data = *node_data_iter;
                std::advance(node_data_iter, 1);
            } else if (node_data_iter == data_node.end()) {
                next_data = *insert_data_iter;
                std::advance(insert_data_iter, 1);
            } else if ((*insert_data_iter).first == (*node_data_iter).first) {
                next_data = std::make_pair(
                    (*insert_data_iter).first,
                    reduce((*insert_data_iter).second, (*node_data_iter).second)
                );
                std::advance(node_data_iter, 1);
                std::advance(insert_data_iter, 1);
            } else if ((*insert_data_iter).first < (*node_data_iter).first) {
                next_data = *insert_data_iter;
                std::advance(insert_data_iter, 1);
            } else if((*node_data_iter).first < (*insert_data_iter).first) {
                next_data = *node_data_iter;
                std::advance(node_data_iter, 1);
            } else {
                // should never reach
                assert(false);
            }

            if (next_data.first < metadata_node.pivot) {
                LT += 1;
            } else {
                GE += 1;
            }
            merged_data.push_back(next_data);
        }

        // assert(old_pivot_idx.has_value());

        

        
        std::optional<std::pair<std::vector<pairT>, bool>> ret = {};
        // if there is no leftover data, then add it to the node
        if (merged_data.size() <= K) {
            // todo move instead of copy?
            data_node = merged_data;
        } else {
            // there exists leftover data which cannot fit in current node, should be passed down
            data_node.clear();
            data_node.reserve(K);
            

            if (LT > merged_data.size() / 2) {
                // send first min(K, elems less than pivot) to left child, write remaining to current node
                idxT num_elems = std::min(K, LT);
                data_node.resize(merged_data.size() - num_elems);
                auto mid = merged_data.begin() + num_elems;
                std::vector<pairT> child_insert_data(merged_data.begin(), mid);
                std::copy(mid, merged_data.end(), data_node.begin());
                ret = std::make_pair(child_insert_data, 0);
            } else {
                // send last min(K, elems greater or equal to pivot) to right child, write remaining to current node
                idxT num_elems = std::min(K, GE);
                data_node.resize(merged_data.size() - num_elems);
                auto mid = merged_data.begin() + (merged_data.size() - num_elems);
                std::vector<pairT> child_insert_data(mid, merged_data.end());
                std::copy(merged_data.begin(), mid, data_node.begin());
                ret = std::make_pair(child_insert_data, 1);
            }
        }

        // update metadata
        metadata_node.min = data_node.at(0).first;
        metadata_node.max = data_node.at(data_node.size()-1).first;
        // metadata_node.pivot = data_node.at(data_node.size()/2).first;

        return ret;
    } 

#ifdef MAP_IMPL
    void add_vec() {
        for (auto key_value : q) {
            if (data.find(key_value.first) == data.end()) {
                data.emplace(key_value.first, key_value.second);
            } else {
//                 std::cout << "key: " << key_value.first.first << " " << key_value.first.second << " old val: " << data.at(key_value.first) << " item: " << key_value.second << " new val: " << reduce(data.at(key_value.first), key_value.second) << std::endl;
                data[key_value.first] = reduce(data.at(key_value.first), key_value.second);
            }
        }
        q.clear();
        q.reserve(K);
    }
#endif

    void add_vec() {
        if (q.size() <= 0) {
            return;
        }

        add_vec(q, 0);
        q.clear();
        q.reserve(K);
    }

    void add_vec(std::vector<pairT> vec, idxT idx) {
        std::optional<std::pair<std::vector<pairT>, bool>> next_insert
            = std::make_pair(vec, false);
        while (next_insert.has_value()) {
            if (idx >= metadata.size()) {
                data.resize(idx+1);
                metadata.resize(idx+1);
            }
            
            if (!metadata.at(idx).has_value()) {
                // new node
                assert(!data.at(idx).has_value());
                data.at(idx) = next_insert.value().first; //todo replace with move?
                metadata.at(idx) = MetadataNode(next_insert.value().first);
                next_insert = {};
                return;
            } else {
                assert(data.at(idx).has_value());
                next_insert = reduce_and_dedup(
                    metadata.at(idx).value(),
                    data.at(idx).value(),
                    next_insert.value().first
                );
            }

            if (!next_insert.has_value()) {
                return;
            }

            idx = next_insert.value().second ? get_right_child_idx(idx) : get_left_child_idx(idx);  
        }
    }


    void global_reduce(idxT node_idx = 0) {
        // if node does not exist, nothing to do.
        if (node_idx >= metadata.size() || !metadata.at(node_idx).has_value()) {
            return;
        }

        #define has_left_child(idx) (get_left_child_idx(idx) < metadata.size() && metadata.at(get_left_child_idx(idx)).has_value())
        #define has_right_child(idx) (get_right_child_idx(idx) < metadata.size() && metadata.at(get_right_child_idx(idx)).has_value())

        // if node's children do not exist, then it is a leaf node. nothing to do
        if (!has_left_child(node_idx) && !has_right_child(node_idx)) {
            return;
        }

        global_reduce(2 * node_idx + 1); // reduce left subtree
        global_reduce(2 * node_idx + 2); // reduce right subtree

        // find leaf node w/ max key in left subtree
        if (has_left_child(node_idx)) {
            idxT max_child_idx = get_left_child_idx(node_idx);
            while (has_right_child(max_child_idx)) {
                max_child_idx = get_right_child_idx(max_child_idx);
            }

            // add vec with max node at parent
            std::vector<pairT> max_data_node(
                data.at(max_child_idx).value().begin(),
                data.at(max_child_idx).value().end()
            );
            data.at(max_child_idx) = {};
            metadata.at(max_child_idx) = {};
            add_vec(max_data_node, node_idx);
        }
        

        // find leaf node w/ min key in right subtree
        if (has_right_child(node_idx)) {
            idxT min_child_idx = get_right_child_idx(node_idx);
            while (has_left_child(min_child_idx)) {
                min_child_idx = get_left_child_idx(min_child_idx);
            }

            // add vec with min node at parent
            std::vector<pairT> min_data_node(
                data.at(min_child_idx).value().begin(),
                data.at(min_child_idx).value().end()
            );
            data.at(min_child_idx) = {};
            metadata.at(min_child_idx) = {};
            add_vec(min_data_node, node_idx);
        }
        

        std::string state_str("post global reduce idx ");
        state_str.append(std::to_string(node_idx));
        print_data(state_str.c_str(), false);
        validate(node_idx, {}, {}, true);
    }

    void validate(bool strict) {
        validate(0, {}, {}, strict);
    }

    // if not strict, we should have all keys in left subtree are < pivot of current node
    // and all keys in right are > pivot of current node
    // if strict, we should have all keys in left subtree are < smallest key of current node
    // and all keys in right are >= largest key of current node
    void validate(
        idxT idx = 0,
        std::optional<keyT> min_key_valid = {},
        std::optional<keyT> max_key_valid = {},
        bool strict = false
    ) {
        if (idx >= data.size() || !data.at(idx).has_value()) {
            assert(idx >= metadata.size() || !metadata.at(idx).has_value());
            return;
        }

        assert(metadata.at(idx).has_value());
        assert(
            !min_key_valid.has_value() || !max_key_valid.has_value()
            || min_key_valid.value() < max_key_valid.value()
        );

        std::vector<pairT> data_node = data.at(idx).value();
        MetadataNode metadata_node = metadata.at(idx).value();

        keyT min_key = data_node.at(0).first;
        keyT max_key = data_node.at(data_node.size() - 1).first;

        assert(min_key == metadata_node.min);
        assert(max_key == metadata_node.max);
        // assert(data_node.at(data_node.size()/2).first == metadata_node.pivot);

        for (idxT i = 0; i < data_node.size(); i++) {
            if (i > 0) {
                assert(min_key < data_node.at(i).first);
                // std::cout << data_node.at(i).first.first << " " << data_node.at(i).first.second
                //     << ", " << data_node.at(i-1).first.first << " "
                //         << data_node.at(i-1).first.second << std::endl;
                assert(data_node.at(i-1).first < data_node.at(i).first);
            }
            if (i < data_node.size() - 1) {
                assert(max_key > data_node.at(i).first);
            }
        }

        assert(!min_key_valid.has_value() || min_key_valid.value() <= min_key);
        assert(!max_key_valid.has_value() || max_key < max_key_valid.value());

        if (strict) {
            // everything in left subtree must be less than current node's smallest key
            validate(get_left_child_idx(idx), min_key_valid, min_key, strict);
            // everything in right subtree must be greater than current node's largest key
            validate(get_right_child_idx(idx), max_key, max_key_valid, strict);
        } else {
            // everything in left subtree must be less than current node's pivot
            validate(get_left_child_idx(idx), min_key_valid, metadata_node.pivot, strict);
            // everything in right subtree must be greater than current node's pivot
            validate(get_right_child_idx(idx), metadata_node.pivot, max_key_valid, strict);  
        }
        
    }

    void print_data(const char *state = "", bool print_data = true) {
        std::cout << "---------------------\nstate ";
        std::cout << state << "\n---------------------\n";
        for (idxT i = 0; i < data.size(); i++) {
            if (data.at(i).has_value()) {
                std::cout << "node " << i << std::endl;
                std::cout << "min, pivot, max: (" << metadata.at(i).value().min << ", " << metadata.at(i).value().pivot << ", " << metadata.at(i).value().max << ")\n";
                if (print_data) {
                    for (auto pair : data.at(i).value()) {
                        std::cout << pair << ", ";
                    }
                    std::cout << std::endl;
                }
                std::cout << std::endl;
            }
        }
    }
};