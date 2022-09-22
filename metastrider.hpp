#include <vector>
#include <utility>
#include <map>
#include <functional>
#include <iostream>

template<typename keyT, typename valT>
class Metastrider {

    using pairT = std::pair<keyT, valT>;

    public:
    uint K;
    std::vector<pairT> q;
    std::map<keyT, valT> data;
    std::function<valT(valT, valT)> reduce_func;

    // TODO how do we configure according to available resouces and decide K?
    // MetaStrider configures its K, MLP-strategy (Section 5), reserves corresponding segment(s) of memory
    // TODO figure out what sz means
    // sz = row size? key + value size?
    // n = number of trees default 1 for now
    // f = function
//     Metastrider(uint sz, uint n, int f) {
//         // K = maximum key value pairs that can fit in one row
//         // K = row size / (key size + value size)
//         K = 12;
//     }

    Metastrider(std::function<valT(valT, valT)> func) {
        K = 12;
        q.reserve(K);
        reduce_func = func;
    }

    void enq_reduce(pairT key_value) {
        q.push_back(key_value);
    }

    void add_vec() {
        for (auto key_value : q) {
            if (data.find(key_value.first) == data.end()) {
                data.emplace(key_value.first, key_value.second);
            } else {
//                 std::cout << "key: " << key_value.first.first << " " << key_value.first.second << " old val: " << data.at(key_value.first) << " item: " << key_value.second << " new val: " << reduce_func(data.at(key_value.first), key_value.second) << std::endl;
                data[key_value.first] = reduce_func(data.at(key_value.first), key_value.second);
            }
        }
        q.clear();
        q.reserve(K);
    }

    void global_reduce() {

    }
};