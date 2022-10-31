#include <vector>
#include <utility>
#include <map>
#include <functional>
#include <iostream>
#include <optional>
#include <cassert>
#include <queue>

// size of rows in bytes
// assume 2KB row
#define ROW_SIZE 2048

#define AVL_IMPL

#ifdef AVL_IMPL
    template<typename keyT, typename valT>
    class MetadataNode {
        using pairT = std::pair<keyT, valT>;

        public:
        keyT pivot;
        std::optional<keyT> min; // min elem in current node
        std::optional<keyT> max; // max elem in current node
        // keyT maxLeft;
        // keyT minRight;

        MetadataNode(const std::vector<pairT> &data) {
            assert(data.size() > 0);
            min = data.at(0).first;
            max = data.at(data.size()-1).first;
            pivot = data.at(data.size()/2).first;
        }
    };

    template<typename keyT, typename valT>
    std::ostream& operator<<(std::ostream& os, const MetadataNode<keyT, valT> &node) {
        os << "[";
        if (node.min.has_value()) {
            os << node.min.value();
        } else {
            os << "{}";
        }
        os << ", " << node.pivot << ", ";
        if (node.max.has_value()) {
            os << node.max.value();
        } else {
            os << "{}";
        }
        os << "]";
        return os;
    }
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
    std::vector<std::optional<std::vector<pairT>>> data;
    std::vector<std::optional<MetadataNode<keyT, valT>>> metadata;

    #define get_left_child_idx(idx) ((2*idx)+1)
    #define get_right_child_idx(idx) ((2*idx)+2)
    #define has_left_child(idx) (get_left_child_idx(idx) < metadata.size() && metadata.at(get_left_child_idx(idx)).has_value())
    #define has_right_child(idx) (get_right_child_idx(idx) < metadata.size() && metadata.at(get_right_child_idx(idx)).has_value())
    #define is_defined_node(idx) (idx < metadata.size() && metadata.at(idx).has_value())

    void update_metadata(idxT idx) {
        if (!is_defined_node(idx)) {
            return;
        }
        if (data.at(idx).value().size() > 0) {
            metadata.at(idx).value().min = data.at(idx).value().at(0).first;
            metadata.at(idx).value().max = data.at(idx).value().at(data.at(idx).value().size() - 1).first;
        } else {
            metadata.at(idx).value().min = {};
            metadata.at(idx).value().max = {};
        }
    }
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
        MetadataNode<keyT, valT> &metadata_node,
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
                auto mid = merged_data.begin() + num_elems;
                std::vector<pairT> child_insert_data(merged_data.begin(), mid);
                data_node.insert(data_node.begin(), mid, merged_data.end());
                ret = std::make_pair(child_insert_data, 0);
            } else {
                // send last min(K, elems greater or equal to pivot) to right child, write remaining to current node
                idxT num_elems = std::min(K, GE);
                auto mid = merged_data.begin() + (merged_data.size() - num_elems);
                std::vector<pairT> child_insert_data(mid, merged_data.end());
                data_node.insert(data_node.begin(), merged_data.begin(), mid);
                ret = std::make_pair(child_insert_data, 1);
            }
        }

        // update metadata
        metadata_node.min = data_node.at(0).first;
        metadata_node.max = data_node.at(data_node.size()-1).first;

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
        if (vec.size() == 0) {
            return;
        }

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
                metadata.at(idx) = MetadataNode<keyT, valT>(next_insert.value().first);
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

    
    // TODO make more efficient
    template<class Compare>
    void get_extrema_keys(std::vector<pairT> &vec, idxT node_idx, idxT n, bool min) {
        std::priority_queue<keyT, std::vector<keyT>, Compare> pq;

        // maps key to vector of (kv-pair, node_idx)
        std::map<keyT, std::vector<std::pair<pairT, idxT>>> key_map;

        // TODO make more efficient
        for (pairT kvpair : data.at(node_idx).value()) {
            keyT k = kvpair.first;
            // if we are getting the minimum from right subtree,
            // then we also want to include current node's >= pivot keys
            // and visa versa for max from left subtree.
            if ((!min && k >= metadata.at(node_idx).value().pivot)
                || (min && k < metadata.at(node_idx).value().pivot)) {
                continue;
            }
            auto k_map_find = key_map.find(k);
            assert(k_map_find == key_map.end());
            pq.push(kvpair.first);
            key_map[k].push_back(std::make_pair(kvpair, node_idx));
        }

        while (pq.size() > n) {
            key_map.erase(pq.top());
            pq.pop();
        }

        std::queue<idxT> q;
        // we want min of right subtree or max of left subtree
        if (min) {
            q.push(get_right_child_idx(node_idx));
        } else {
            q.push(get_left_child_idx(node_idx));
        }
        
        while (q.size() > 0) {
            idxT cur_idx = q.front();
            q.pop();
            if (is_defined_node(cur_idx)) {
                for (pairT kvpair : data.at(cur_idx).value()) {
                    keyT k = kvpair.first;
                    auto k_map_find = key_map.find(k);
                    if (k_map_find == key_map.end()) {
                        pq.push(kvpair.first);
                    }
                    key_map[k].push_back(std::make_pair(kvpair, cur_idx));
                }
                while (pq.size() > n) {
                    key_map.erase(pq.top());
                    pq.pop();
                }
                if (has_left_child(cur_idx)) {
                    q.push(get_left_child_idx(cur_idx));
                }
                if (has_right_child(cur_idx)) {
                    q.push(get_right_child_idx(cur_idx));
                }
            }
        }
        for (idxT i = 0; i < n; i++) {
            if (pq.size() == 0) {
                break;
            }
            // extrema key
            keyT ext_key = pq.top();
            auto find_ext_key = key_map.find(ext_key);
            assert(find_ext_key != key_map.end());
            std::optional<pairT> reduced_kvpair;
            for (std::pair<pairT, idxT> kvpair_node_idx_pair: key_map[ext_key]) {
                pairT kvpair = kvpair_node_idx_pair.first;
                idxT rem_node_idx = kvpair_node_idx_pair.second;

                if (reduced_kvpair.has_value()) {
                    reduced_kvpair.value().second = reduce(reduced_kvpair.value().second, kvpair.second);
                } else {
                    reduced_kvpair = kvpair;
                }

                std::vector<pairT> &data_node = data.at(rem_node_idx).value();
                auto rem_node_find = std::find(data_node.begin(), data_node.end(), kvpair);
                assert(rem_node_find != data_node.end());
                data_node.erase(rem_node_find);
                update_metadata(rem_node_idx);
            }
            assert(reduced_kvpair.has_value());
            vec.push_back(reduced_kvpair.value());
            pq.pop();
        }
        std::sort(vec.begin(), vec.end());
    }

    void get_max_keys(std::vector<pairT> &vec, idxT node_idx, idxT n) {
        get_extrema_keys<std::greater<keyT>>(vec, node_idx, n, false);
    }

    void get_min_keys(std::vector<pairT> &vec, idxT node_idx, idxT n) {
        get_extrema_keys<std::less<keyT>>(vec, node_idx, n, true);
    }

    void global_reduce(idxT node_idx = 0) {
        // if node does not exist, nothing to do.
        if (node_idx >= metadata.size() || !metadata.at(node_idx).has_value()) {
            return;
        }

        // if node's children do not exist, then it is a leaf node. nothing to do
        if (!has_left_child(node_idx) && !has_right_child(node_idx)) {
            return;
        }

        std::vector<pairT> &data_node = data.at(node_idx).value();
        MetadataNode<keyT, valT> &metadata_node = metadata.at(node_idx).value();

        // TODO include LT and GE in metadata?
        idxT LT;
        idxT GE;
        
        if (has_left_child(node_idx)) {
            // replace curr node's keys which are less than pivot with left subtree's
            // greatest values.
            // This way, all of the left subtree's keys are still less than the current node's
            // keys.

            LT = 0;
            GE = 0;
            for (pairT pair : data_node) {
                if (pair.first < metadata_node.pivot) {
                    LT += 1;
                } else {
                    GE += 1;
                }
            }

            std::vector<pairT> max_pairs;
            get_max_keys(max_pairs, node_idx, std::max(LT, K - GE));

            LT = 0;
            GE = 0;
            for (pairT pair : data_node) {
                if (pair.first < metadata_node.pivot) {
                    LT += 1;
                } else {
                    GE += 1;
                }
            }
            
            

            std::vector<pairT> cur_LT_pairs(data_node.begin(), data_node.begin() + LT);
            data_node.erase(data_node.begin(), data_node.begin() + LT);
            data_node.insert(data_node.begin(), max_pairs.begin(), max_pairs.end());
            update_metadata(node_idx);
            add_vec(cur_LT_pairs, get_left_child_idx(node_idx));
            // print_data("post left subtree max vec");
        }
        

        // find leaf node w/ min key in right subtree
        if (has_right_child(node_idx)) {
            LT = 0;
            GE = 0;
            for (pairT pair : data_node) {
                if (pair.first < metadata_node.pivot) {
                    LT += 1;
                } else {
                    GE += 1;
                }
            }

            std::vector<pairT> min_pairs;
            get_min_keys(min_pairs, node_idx, std::max(GE, K - LT));

            LT = 0;
            GE = 0;
            for (pairT pair : data_node) {
                if (pair.first < metadata_node.pivot) {
                    LT += 1;
                } else {
                    GE += 1;
                }
            }
            std::vector<pairT> cur_GE_pairs(data_node.begin() + (data_node.size() - GE), data_node.end());
            data_node.erase(data_node.begin() + (data_node.size() - GE), data_node.end());
            data_node.insert(data_node.end(), min_pairs.begin(), min_pairs.end());
            update_metadata(node_idx);
            add_vec(cur_GE_pairs, get_right_child_idx(node_idx));
        }

        // std::string state_str("post global reduce idx ");
        // state_str.append(std::to_string(node_idx));
        // print_data(state_str.c_str(), true, node_idx);
        validate_node(node_idx);

        global_reduce(2 * node_idx + 1); // reduce left subtree
        global_reduce(2 * node_idx + 2); // reduce right subtree
    }

    void validate_node(idxT idx) {
        std::optional<keyT> min_k = {};
        std::optional<keyT> max_k = {};
        validate_node(idx, min_k, max_k);
    }

    void validate_node(idxT idx, std::optional<keyT> &min_key, std::optional<keyT> &max_key) {
        if (idx >= data.size() || !data.at(idx).has_value()) {
            assert(idx >= metadata.size() || !metadata.at(idx).has_value());
            return;
        }

        assert(metadata.at(idx).has_value());

        const std::vector<pairT> &data_node = data.at(idx).value();
        const MetadataNode<keyT, valT> &metadata_node = metadata.at(idx).value();

        if (data_node.size() != 0) {
            min_key = data_node.at(0).first;
            max_key = data_node.at(data_node.size() - 1).first;

            for (idxT i = 0; i < data_node.size(); i++) {
                if (i > 0) {
                    assert(min_key.value() < data_node.at(i).first);
                    assert(data_node.at(i-1).first < data_node.at(i).first);
                }
                if (i < data_node.size() - 1) {
                    assert(max_key.value() > data_node.at(i).first);
                }
            }
        }
        assert(min_key == metadata_node.min);
        assert(max_key == metadata_node.max);
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

        // assert(metadata.at(idx).has_value());
        assert(
            !min_key_valid.has_value() || !max_key_valid.has_value()
            || min_key_valid.value() < max_key_valid.value()
        );

        // std::vector<pairT> data_node = data.at(idx).value();
        // const MetadataNode<keyT, valT> &metadata_node = metadata.at(idx).value();

        std::optional<keyT> min_key;
        std::optional<keyT> max_key;

        validate_node(idx, min_key, max_key);
        // if (data_node.size() != 0) {
        //     min_key = data_node.at(0).first;
        //     max_key = data_node.at(data_node.size() - 1).first;

        //     for (idxT i = 0; i < data_node.size(); i++) {
        //         if (i > 0) {
        //             assert(min_key.value() < data_node.at(i).first);
        //             assert(data_node.at(i-1).first < data_node.at(i).first);
        //         }
        //         if (i < data_node.size() - 1) {
        //             assert(max_key.value() > data_node.at(i).first);
        //         }
        //     }
        // }
        // assert(min_key == metadata_node.min);
        // assert(max_key == metadata_node.max);
        

        assert(!min_key_valid.has_value() || !min_key.has_value() || min_key_valid.value() <= min_key.value());
        assert(!max_key_valid.has_value() || !max_key.has_value() || max_key.value() < max_key_valid.value());

        if (strict) {
            // everything in left subtree must be less than current node's smallest key
            validate(
                get_left_child_idx(idx),
                min_key_valid, min_key.has_value() ? min_key.value() : max_key_valid,
                strict
            );
            // everything in right subtree must be greater than current node's largest key
            validate(
                get_right_child_idx(idx),
                max_key.has_value() ? max_key.value() : min_key_valid,
                max_key_valid,
                strict
            );
        } else {
            // everything in left subtree must be less than current node's pivot
            validate(get_left_child_idx(idx), min_key_valid, metadata.at(idx).value().pivot, strict);
            // everything in right subtree must be greater than current node's pivot
            validate(get_right_child_idx(idx), metadata.at(idx).value().pivot, max_key_valid, strict);  
        }
        
    }

    void print_data(const char *state = "", bool print_data = true, std::optional<idxT> max_idx = {}) {
        std::cout << "---------------------\nstate ";
        std::cout << state << "\n---------------------\n";
        for (idxT i = 0; i < data.size(); i++) {
            if (max_idx.has_value() && i > max_idx) {
                return;
            }
            if (data.at(i).has_value()) {
                std::cout << "node " << i << std::endl;
                std::cout << "meta: " << metadata.at(i).value() << std::endl;
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