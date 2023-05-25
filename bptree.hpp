#ifndef BPTREE_HPP_BPTREE2_HPP
#define BPTREE_HPP_BPTREE2_HPP
#include <fstream>
#include <iostream>
//#include <queue>
#include "vector.hpp"
template <class Key, class T, int M = 100, int L = 100>
class BPTree {
private:
    std::string name;
    const std::string FILE[5] = { "info", "tree", "data", "vector_t", "vector_d" };
    std::fstream file[5];

    int INT = sizeof(int), KEY = sizeof(Key), _T = sizeof(T);
    int Size = 0, last_t = 0, last_d = 0, root = 0;
    int VALUE = sizeof(std::pair<Key, T>);

    /******************************************************** /
                   Structure of node and data

    info: (size: 4 int)
     __________      ________      ________      __________
    |          |    |        |    |        |    |          |
    |   Size   | -> | last_t | -> | last_d | -> |   root   |
    |__________|    |________|    |________|    |__________|

    node of tree: (size: 4 int and 1 Key in total)
     _________      __________      _________      _________      ___________
    |         |    |          |    |         |    |         |    |           |
    |  count  | -> |  parent  | -> |   key   | -> | pointer | -> | next node | -> ...
    |_________|    |__________|    |_________|    |_________|    |___________|

    in node of tree, each piece of tree_node needs a head node, with an empty key

    data_node: (size: 4 int and 1 Value in total)
     _________      __________      _________      __________      __________
    |         |    |          |    |         |    |          |    |          |
    |  count  | -> |  parent  | -> |  value  | -> | next key | -> | next val | -> ...
    |_________|    |__________|    |_________|    |__________|    |__________|

    in data, value can be of the same key, thus, "next val" is used to 
        represent the next data with the same key.
    "count" is only used in the first piece of data.

    vector_t or vector_d:
     _________      _________      _________
    |         |    |         |    |         |
    |  count  | -> |  value  | -> |  value  | -> ...
    |_________|    |_________|    |_________|

    this "file vector" is used to record the deleted position for using.

    initial state of root: (size: 4 int and 1 Key in total)
     _________      __________      _________      _________      ___________
    |  count  |    |  parent  |    |   key   |    | pointer |    | next node |
    |    0    | -> |    -1    | -> |   NULL  | -> |    -1   | -> |     -1    | -> ...
    |_________|    |__________|    |_________|    |_________|    |___________|

    / ********************************************************/

    int NODE_DATA = 4 * INT + VALUE;
    int NODE_TREE = 4 * INT + KEY;

    using value_type = std::pair<Key, T>;
    using fvector = std::fstream;

    enum node_type { tree_node, data_node };
    enum insert_type { to_key, to_value };

    Key temp_key;
    value_type temp_value;
    int zero = 0, nega_one = -1, one = 1;

    // open the existed files
    void openFile() {
        for (int i = 0; i < 5; ++i) {
            file[i].open(FILE[i], std::ios::binary 
                | std::ios::in | std::ios::out);
        }
    }

    // close already opened files
    void closeFile() {
        for (int i = 0; i < 5; ++i) {
            file[i].close();
        }
    }

    // clear content of all the files
    void clearFile() {
        for (int i = 0; i < 5; ++i) {
            file[i].open(FILE[i], std::ios::out | std::ios::trunc);
        }
    }

    // to judge whether a file has existed
    bool exist(std::string name) {
        std::fstream test(name, std::ios::binary 
            | std::ios::in | std::ios::out);
        bool temp = test.is_open();
        test.close();
        return temp;
    }

    // to update Size, last_t, last_d, root in file "info"
    void updateinfo() {
        file[0].seekg(std::ios::beg);
        file[0].write((char*)&Size, INT);
        file[0].write((char*)&last_t, INT);
        file[0].write((char*)&last_d, INT);
        file[0].write((char*)&root, INT);
    }

    // to get a deleted position in a specific file for using
    int getInsertPos(fvector& f, node_type type) {
        int count, pos;
        f.seekg(std::ios::beg);
        f.read((char*)&count, INT);
        if (!count) {
            switch (type) {
            case tree_node:
                last_t += NODE_TREE;
                return last_t - NODE_TREE;
            case data_node:
                last_d += NODE_DATA;
                return last_d - NODE_DATA;
            }
        }

        f.seekg(count * INT, std::ios::beg);
        f.read((char*)&pos, INT);
        --count; f.seekp(std::ios::beg);
        f.write((char*)&count, INT);
        return pos;
    }

    // to add a deleted position for insert later to a specific file
    void addInsertPos(fvector& f, int pos) {
        int count;
        f.seekg(std::ios::beg);
        f.read((char*)&count, INT);
        
        ++count;
        f.seekp(count * INT, std::ios::beg);
        f.write((char*)&pos, INT);
        f.seekp(std::ios::beg);
        f.write((char*)&count, INT);
    }

    // to get the ancestor tree_node of a specific tree_node
    // "pos" shouldn't be -1 initially
    int getTreeNodeAncestor(std::fstream& f, int pos) {
        int count;
        f.seekg(pos, std::ios::beg);
        f.read((char*)&count, INT);

        if (count > 0) return pos;
        else return -count - 1;
    }

    // to update parent-son relationship between prev and cur among data nodes
    void updateDataNode(std::fstream& f, int prev, int cur, insert_type type = to_key) {
        int next, num_of_int = 0;

        switch (type) {
        case to_key: num_of_int = 2; break;
        case to_value: num_of_int = 3; break;
        }

        f.seekg(prev + num_of_int * INT + VALUE, std::ios::beg);
        f.read((char*)&next, INT);
        f.seekp(prev + num_of_int * INT + VALUE, std::ios::beg);
        f.write((char*)&cur, INT);

        f.seekp(cur + INT, std::ios::beg);
        f.write((char*)&prev, INT);
        f.seekp(cur + num_of_int * INT + VALUE, std::ios::beg);
        f.write((char*)&next, INT);

        if (next != -1) {
            f.seekp(next + INT, std::ios::beg);
            f.write((char*)&cur, INT);
        }
    }

    // to update parent-son relationship between prev and cur among tree nodes
    void updateTreeNode(std::fstream& f, int prev, int cur) {
        int next;

        f.seekg(prev + 3 * INT + KEY, std::ios::beg);
        f.read((char*)&next, INT);
        f.seekp(prev + 3 * INT + KEY, std::ios::beg);
        f.write((char*)&cur, INT);

        f.seekp(cur + INT, std::ios::beg);
        f.write((char*)&prev, INT);
        f.seekp(cur + 3 * INT + KEY, std::ios::beg);
        f.write((char*)&next, INT);

        if (next != -1) {
            f.seekp(next + INT, std::ios::beg);
            f.write((char*)&cur, INT);
        }
    }

    // insert a data_node to file "data", all parameters included
    void insertDataNode(std::fstream& f, int pos, const value_type& val, 
        int parent = -1, int next_key = -1, int next_val = -1, int count = -1) {
        f.seekp(pos, std::ios::beg);
        f.write((char*)&count, INT);
        f.write((char*)&parent, INT);
        f.write((char*)&val, VALUE);
        f.write((char*)&next_key, INT);
        f.write((char*)&next_val, INT);
    }

    // insert a tree_node to file "tree", all parameters included
    void insertTreeNode(std::fstream& f, int pos, const Key& key, int pointer,
        int parent = -1, int next_node = -1, int count = -1) {
        f.seekg(pos, std::ios::beg);
        f.write((char*)&count, INT);
        f.write((char*)&parent, INT);
        f.write((char*)&key, KEY);
        f.write((char*)&pointer, INT);
        f.write((char*)&next_node, INT);
    }

    // use the pos acquired from another operation to insert val to the data file
    void insertToData(int pos, const value_type& val) {
        int count; file[2].seekg(pos, std::ios::beg);
        file[2].read((char*)&count, INT);

        int next_pos, old_pos = pos;
        file[2].seekg(pos + 2 * INT + VALUE, std::ios::beg);
        file[2].read((char*)&next_pos, INT);

        // search sequently until finding the one with a larger key
        value_type temp;
        while (next_pos != -1) {
            file[2].seekg(next_pos + 2 * INT, std::ios::beg);
            file[2].read((char*)&temp, VALUE);

            // judge whether there's a value with the same key
            if (temp.first == val.first) {
                int first_pos = next_pos;

                while (true) {
                    if (temp.second < val.second) {
                        old_pos = next_pos;
                        file[2].seekg(INT, std::ios::cur);
                        file[2].read((char*)&next_pos, INT);

                        if (next_pos == -1) {
                            int posForInsert = getInsertPos(file[4], data_node);
                            insertDataNode(file[2], posForInsert, val);
                            updateDataNode(file[2], old_pos, posForInsert, to_value);
                            checkDataNode(file[2], pos);
                            return;
                        }

                        file[2].seekg(next_pos + 2 * INT, std::ios::beg);
                        file[2].read((char*)&temp, VALUE);
                    }
                    else if (temp.second == val.second) {
                        return;
                    }
                    else {
                        if (next_pos == first_pos) {
                            int posForInsert = getInsertPos(file[4], data_node), next;
                            file[2].seekg(next_pos + 2 * INT + VALUE, std::ios::beg);
                            file[2].read((char*)&next, INT);

                            insertDataNode(file[2], posForInsert, val, 
                                old_pos, next, next_pos);

                            file[2].seekp(old_pos + 2 * INT + VALUE, std::ios::beg);
                            file[2].write((char*)&posForInsert, INT);
                            if (next != -1) {
                                file[2].seekp(next + INT, std::ios::beg);
                                file[2].write((char*)&posForInsert, INT);
                            }
                            return;
                        }
                        else {
                            int posForInsert = getInsertPos(file[4], data_node);
                            insertDataNode(file[2], posForInsert, val);
                            updateDataNode(file[2], old_pos, posForInsert, to_value);
                            return;
                        }
                    }
                }
            }

            // already find the right place
            if (temp.first > val.first) {
                int posForInsert = getInsertPos(file[4], data_node);
                insertDataNode(file[2], posForInsert, val);
                updateDataNode(file[2], old_pos, posForInsert);

                ++Size; ++count;
                file[2].seekp(pos, std::ios::beg);
                file[2].write((char*)&count, INT);
                checkDataNode(file[2], pos);
                return;
            }

            old_pos = next_pos;
            file[2].read((char*)&next_pos, INT);
        }

        // when next_pos is equal to -1
        int posForInsert = getInsertPos(file[4], data_node);
        file[2].seekp(old_pos + 2 * INT + VALUE, std::ios::beg);
        file[2].write((char*)&posForInsert, INT);

        int test;
        file[2].seekg(old_pos + 2 * INT + VALUE, std::ios::beg);
        file[2].read((char*)&test, INT);

        insertDataNode(file[2], posForInsert, val, old_pos);

        ++Size; ++count; file[2].seekp(pos, std::ios::beg);
        file[2].write((char*)&count, INT);

        checkDataNode(file[2], pos);
    }

    // to remove the data same as val
    void removeData(int pos, const value_type& val) {
        int count;
        file[2].seekg(pos, std::ios::beg);
        file[2].read((char*)&count, INT);

        int cur, old_pos = pos;
        file[2].seekg(pos + 2 * INT + VALUE, std::ios::beg);
        file[2].read((char*)&cur, INT);

        value_type temp;
        while (cur != -1) {
            file[2].seekg(cur + 2 * INT, std::ios::beg);
            file[2].read((char*)&temp, VALUE);

            if (temp.first == val.first && temp.second == val.second) {
                int next_val, next_key;
                file[2].seekg(cur + 2 * INT + VALUE, std::ios::beg);
                file[2].read((char*)&next_key, INT);
                file[2].read((char*)&next_val, INT);
                addInsertPos(file[4], cur);

                if (next_val == -1) {
                    file[2].seekp(old_pos + 2 * INT + VALUE, std::ios::beg);
                    file[2].write((char*)&next_key, INT);

                    if (next_key != -1) {
                        file[2].seekp(next_key + INT, std::ios::beg);
                        file[2].write((char*)&old_pos, INT);
                    }

                    --count; --Size;
                    file[2].seekp(pos, std::ios::beg);
                    file[2].write((char*)&count, INT);

                    checkDataNode(file[2], pos);
                    return;
                }
                else {
                    file[2].seekp(old_pos + 2 * INT + VALUE, std::ios::beg);
                    file[2].write((char*)&next_key, INT);
                    updateDataNode(file[2], old_pos, next_val, to_key);
                    return;
                }
            }

            if (temp.first == val.first) {
                old_pos = cur;
                file[2].seekg(cur + 3 * INT + VALUE, std::ios::beg);
                file[2].read((char*)&cur, INT);

                while (cur != -1) {
                    file[2].seekg(cur + 2 * INT, std::ios::beg);
                    file[2].read((char*)&temp, VALUE);

                    if (temp.second == val.second) {
                        addInsertPos(file[4], cur);

                        int next_val;
                        file[2].seekg(INT, std::ios::cur);
                        file[2].read((char*)&next_val, INT);
                        file[2].seekp(old_pos + 3 * INT + VALUE, std::ios::beg);
                        file[2].write((char*)&next_val, INT);

                        if (next_val != -1) {
                            file[2].seekp(next_val + INT, std::ios::beg);
                            file[2].write((char*)&old_pos, INT);
                        }
                        return;
                    }
                    if (temp.second > val.second) return;

                    old_pos = cur;
                    file[2].seekg(cur + 3 * INT + VALUE, std::ios::beg);
                    file[2].read((char*)&cur, INT);
                }
                return;
            }

            if (temp.first > val.first) return;

            old_pos = cur;
            file[2].seekg(cur + 2 * INT + VALUE, std::ios::beg);
            file[2].read((char*)&cur, INT);
        }
    }

    // to find all the value with the samw key, return a sorted vector
    sjtu::vector<T> findValue(int pos, const Key& key) {
        sjtu::vector<T> vec; int cur;

        file[2].seekg(pos + 2 * INT + VALUE, std::ios::beg);
        file[2].read((char*)&cur, INT);

        value_type temp;
        while (cur != -1) {
            file[2].seekg(cur + 2 * INT, std::ios::beg);
            file[2].read((char*)&temp, VALUE);
            file[2].read((char*)&cur, INT);

            if (temp.first == key) {
                vec.push_back(temp.second);
                file[2].read((char*)&cur, INT);

                while (cur != -1) {
                    file[2].seekg(cur + 2 * INT, std::ios::beg);
                    file[2].read((char*)&temp, VALUE);
                    vec.push_back(temp.second);

                    file[2].seekg(INT, std::ios::cur);
                    file[2].read((char*)&cur, INT);
                }

                //vec.sort();
                return vec;
            }

            if (temp.first > key) return vec;
        }
        return vec;
    }

    // to find the right place to perform an operation
    int findInTree(int pos, const value_type& val) {
        // the pointers of the last tree_nodes pointing at data are set to be negative
        if (pos < 0) {
            return -pos - 1;
        }

        int next_pos, old_pos = pos;
        file[1].seekg(pos + 3 * INT + KEY, std::ios::beg);
        file[1].read((char*)&next_pos, INT);

        while (next_pos != -1) {
            Key cur;
            file[1].seekg(next_pos + 2 * INT, std::ios::beg);
            file[1].read((char*)&cur, KEY);
            if (cur > val.first) {
                break;
            }

            old_pos = next_pos;
            file[1].seekg(INT, std::ios::cur);
            file[1].read((char*)&next_pos, INT);
        }

        int next_find_pos;
        file[1].seekg(old_pos + 2 * INT + KEY, std::ios::beg);
        file[1].read((char*)&next_find_pos, INT);
        return findInTree(next_find_pos, val);
    }
    
    // combine the adjacent data nodes
    void checkDataNode(std::fstream& f, int pos) {
        int count, parent; f.seekg(pos, std::ios::beg);
        f.read((char*)&count, INT);
        f.read((char*)&parent, INT);

        int splitPos = (L + 1) >> 1, curPos, oldPos = pos;
        f.seekg(pos + 2 * INT + VALUE, std::ios::beg);
        f.read((char*)&curPos, INT);

        if (count > L) {
            for (int i = 1; i <= splitPos; ++i) {
                oldPos = curPos;
                f.seekg(curPos + 2 * INT + VALUE, std::ios::beg);
                f.read((char*)&curPos, INT);
            }

            // split the node into two, setting their "count"s
            f.seekp(pos, std::ios::beg);
            f.write((char*)&splitPos, INT);
            f.seekp(oldPos + 2 * INT + VALUE, std::ios::beg);
            f.write((char*)&nega_one, INT);

            // create head_node for the second node
            int head = getInsertPos(file[4], data_node);
            insertDataNode(f, head, value_type(Key(), T()), 
                -1, curPos, -1, count - splitPos);
            f.seekp(curPos + INT, std::ios::beg);
            f.write((char*)&head, INT);

            // get the key to insert to the parent node
            value_type value;
            f.seekg(curPos + 2 * INT, std::ios::beg);
            f.read((char*)&value, VALUE);

            // insert operation
            int posForInserted = getInsertPos(file[3], tree_node);
            int ancestor = getTreeNodeAncestor(file[1], parent);
            insertTreeNode(file[1], posForInserted, value.first, 
                -head - 1, -1, -1, -ancestor - 1);
            updateTreeNode(file[1], parent, posForInserted);

            f.seekp(head + INT, std::ios::beg);
            f.write((char*)&posForInserted, INT);

            int cnt;
            file[1].seekg(ancestor, std::ios::beg);
            file[1].read((char*)&cnt, INT);
            ++cnt; file[1].seekp(ancestor, std::ios::beg);
            file[1].write((char*)&cnt, INT);
            
            if (cnt > M) {
                checkTreeNode(file[1], ancestor);
            }
        }
        else if (count < splitPos) {
            int left_p, right_p, left_d, right_d;
            file[1].seekg(parent + INT, std::ios::beg);
            file[1].read((char*)&left_p, INT);
            file[1].seekg(INT + KEY, std::ios::cur);
            file[1].read((char*)&right_p, INT);

            if (right_p != -1) {
                // find the end node among the current data_nodes
                int end, p = pos;
                while (p != -1) {
                    end = p;
                    f.seekg(p + 2 * INT + VALUE, std::ios::beg);
                    f.read((char*)&p, INT);
                }

                file[1].seekg(right_p + 2 * INT + KEY, std::ios::beg);
                file[1].read((char*)&right_d, INT);
                right_d = -right_d - 1;

                int right_count;
                f.seekg(right_d, std::ios::beg);
                f.read((char*)&right_count, INT);

                if (right_count > splitPos) {
                    int first, second;
                    f.seekg(right_d + 2 * INT + VALUE, std::ios::beg);
                    f.read((char*)&first, INT);

                    value_type temp;
                    f.seekg(first + 2 * INT + VALUE, std::ios::beg);
                    f.read((char*)&second, INT);
                    f.seekg(second + 2 * INT, std::ios::beg);
                    f.read((char*)&temp, VALUE);
                    Key second_key(temp.first);

                    --right_count;
                    f.seekp(right_d, std::ios::beg);
                    f.write((char*)&right_count, INT);
                    f.seekp(INT + VALUE, std::ios::cur);
                    f.write((char*)&second, INT);

                    ++count;
                    f.seekp(pos, std::ios::beg);
                    f.write((char*)&count, INT);

                    f.seekp(second + INT, std::ios::beg);
                    f.write((char*)&right_d, INT);

                    f.seekp(end + 2 * INT + VALUE, std::ios::beg);
                    f.write((char*)&first, INT);
                    f.seekp(first + INT, std::ios::beg);
                    f.write((char*)&end, INT);
                    f.seekp(VALUE, std::ios::cur);
                    f.write((char*)&nega_one, INT);

                    file[1].seekp(right_p + 2 * INT, std::ios::beg);
                    file[1].write((char*)&second_key, KEY);
                    return;
                }
                else {
                    int next_node, parent_count;
                    file[1].seekg(right_p + 3 * INT + KEY, std::ios::beg);
                    file[1].read((char*)&next_node, INT);
                    addInsertPos(file[3], right_p);

                    file[1].seekp(parent + 3 * INT + KEY, std::ios::beg);
                    file[1].write((char*)&next_node, INT);

                    if (next_node != -1) {
                        file[1].seekp(next_node + INT, std::ios::beg);
                        file[1].write((char*)&parent, INT);
                    }

                    int ancestor = getTreeNodeAncestor(file[1], parent);
                    file[1].seekg(ancestor, std::ios::beg);
                    file[1].read((char*)&parent_count, INT);
                    --parent_count; file[1].seekp(ancestor, std::ios::beg);
                    file[1].write((char*)&parent_count, INT);

                    int first;
                    f.seekg(right_d + 2 * INT + VALUE, std::ios::beg);
                    f.read((char*)&first, INT);
                    addInsertPos(file[4], right_d);

                    f.seekp(end + 2 * INT + VALUE, std::ios::beg);
                    f.write((char*)&first, INT);
                    f.seekp(first + INT, std::ios::beg);
                    f.write((char*)&end, INT);

                    count += right_count;
                    f.seekp(pos, std::ios::beg);
                    f.write((char*)&count, INT);

                    if (parent_count < splitPos) {
                        checkTreeNode(file[1], ancestor);
                    }
                    return;
                }
            }
            else if (left_p != -1) {
                file[1].seekg(left_p + 2 * INT + KEY, std::ios::beg);
                file[1].read((char*)&left_d, INT);
                left_d = -left_d - 1;

                // find the end node among the current data_nodes
                int p = left_d, end;
                while (p != -1) {
                    end = p;
                    f.seekg(p + 2 * INT + VALUE, std::ios::beg);
                    f.read((char*)&p, INT);
                }

                int left_count;
                f.seekg(left_d, std::ios::beg);
                f.read((char*)&left_count, INT);

                if (left_count > splitPos) {
                    f.seekg(end + INT, std::ios::beg);
                    f.read((char*)&p, INT);
                    f.seekp(p + 2 * INT + VALUE, std::ios::beg);
                    f.write((char*)&nega_one, INT);

                    --left_count;
                    f.seekp(left_d, std::ios::beg);
                    f.write((char*)&left_count, INT);

                    ++count;
                    f.seekp(pos, std::ios::beg);
                    f.write((char*)&count, INT);

                    updateDataNode(file[2], pos, end, to_key);

                    value_type temp;
                    f.seekg(end + 2 * INT, std::ios::beg);
                    f.read((char*)&temp, VALUE);
                    Key new_key(temp.first);

                    file[1].seekp(parent + 2 * INT, std::ios::beg);
                    file[1].write((char*)&new_key, KEY);
                    return;
                }
                else {
                    addInsertPos(file[3], parent);
                    file[1].seekp(left_p + 3 * INT + KEY, std::ios::beg);
                    file[1].write((char*)&nega_one, INT);

                    int ancestor = getTreeNodeAncestor(file[1], left_p), parent_count;
                    file[1].seekg(ancestor, std::ios::beg);
                    file[1].read((char*)&parent_count, INT);
                    --parent_count; file[1].seekp(ancestor, std::ios::beg);
                    file[1].write((char*)&parent_count, INT);

                    int first;
                    f.seekg(pos + 2 * INT + VALUE, std::ios::beg);
                    f.read((char*)&first, INT);
                    addInsertPos(file[4], pos);

                    f.seekp(end + 2 * INT + VALUE, std::ios::beg);
                    f.write((char*)&first, INT);
                    f.seekp(first + INT, std::ios::beg);
                    f.write((char*)&end, INT);

                    count += left_count;
                    f.seekp(left_d, std::ios::beg);
                    f.write((char*)&count, INT);

                    if (parent_count < splitPos) {
                        checkTreeNode(file[1], ancestor);
                        return;
                    }
                }
            }
        }
    }

    // combine the adjacent nodes if needed
    void checkTreeNode(std::fstream& f, int pos) {
        int count, parent;
        f.seekg(pos, std::ios::beg);
        f.read((char*)&count, INT);
        f.read((char*)&parent, INT);

        int splitPos = (M + 1) >> 1, curPos = pos, oldPos = pos;
        int count_second = count - splitPos;

        if (count > M) {
            for (int i = 1; i <= splitPos; ++i) {
                oldPos = curPos;
                f.seekg(curPos + 3 * INT + KEY, std::ios::beg);
                f.read((char*)&curPos, INT);
            }
            // split the node into two and set their "count"s
            f.seekp(pos, std::ios::beg);
            f.write((char*)&splitPos, INT);
            f.seekp(oldPos + 3 * INT + KEY, std::ios::beg);
            f.write((char*)&nega_one, INT);

            f.seekp(curPos, std::ios::beg);
            f.write((char*)&count_second, INT);

            // revise parents for later use
            int p = curPos, pa = -curPos - 1;
            f.seekg(p + 3 * INT + KEY, std::ios::beg);
            f.read((char*)&p, INT);
            while (p != -1) {
                f.seekp(p, std::ios::beg);
                f.write((char*)&pa, INT);
                f.seekg(p + 3 * INT + KEY, std::ios::beg);
                f.read((char*)&p, INT);
            }

            // get the key to be inserted to the parent node
            Key key;
            f.seekg(curPos + 2 * INT, std::ios::beg);
            f.read((char*)&key, KEY);

            // if the node to be split is the root
            if (parent == -1) {
                parent = root = getInsertPos(file[3], tree_node);
                insertTreeNode(file[1], parent, Key(), pos, -1, -1, 1);
                f.seekg(pos + INT, std::ios::beg);
                f.write((char*)&parent, INT);
            }

            // insert operation
            int posForInserted = getInsertPos(file[3], tree_node);
            int ancestor = getTreeNodeAncestor(file[1], parent);
            insertTreeNode(file[1], posForInserted, key, 
                curPos, -1, -1, -ancestor - 1);
            updateTreeNode(file[1], parent, posForInserted);

            f.seekp(curPos + INT, std::ios::beg);
            f.write((char*)&posForInserted, INT);

            int cnt;
            f.seekg(ancestor, std::ios::beg);
            f.read((char*)&cnt, INT);
            ++cnt; f.seekp(ancestor, std::ios::beg);
            f.write((char*)&cnt, INT);

            if (cnt > M) {
                checkTreeNode(f, ancestor);
            }
        }
        else if (count < splitPos && pos != root) {
            int left_p, right_p, left_d, right_d;
            f.seekg(parent + INT, std::ios::beg);
            f.read((char*)&left_p, INT);
            f.seekg(INT + KEY, std::ios::cur);
            f.read((char*)&right_p, INT);

            if (right_p != -1) {
                // find the end node among the current data_nodes
                int end, p = pos;
                while (p != -1) {
                    end = p;
                    f.seekg(p + 3 * INT + KEY, std::ios::beg);
                    f.read((char*)&p, INT);
                }

                f.seekg(right_p + 2 * INT + KEY, std::ios::beg);
                f.read((char*)&right_d, INT);

                int right_count;
                f.seekg(right_d, std::ios::beg);
                f.read((char*)&right_count, INT);

                if (right_count > splitPos) {
                    int first;
                    f.seekg(right_d + 3 * INT + KEY, std::ios::beg);
                    f.read((char*)&first, INT);

                    --right_count;
                    f.seekp(first, std::ios::beg);
                    f.write((char*)&right_count, INT);
                    f.write((char*)&right_p, INT);

                    ++count;
                    f.seekp(pos, std::ios::beg);
                    f.write((char*)&count, INT);

                    Key temp;
                    f.seekg(first + 2 * INT, std::ios::beg);
                    f.read((char*)&temp, KEY);
                    f.seekp(right_p + 2 * INT, std::ios::beg);
                    f.write((char*)&temp, KEY);
                    f.write((char*)&first, INT);

                    int pa = -pos - 1;
                    f.seekp(end + 3 * INT + KEY, std::ios::beg);
                    f.write((char*)&right_d, INT);
                    f.seekp(right_d, std::ios::beg);
                    f.write((char*)&pa, INT);
                    f.write((char*)&end, INT);
                    f.seekp(KEY + INT, std::ios::cur);
                    f.write((char*)&nega_one, INT);

                    int _pa = -first - 1, p;
                    f.seekg(first + 3 * INT + KEY, std::ios::beg);
                    f.read((char*)&p, INT);
                    while (p != -1) {
                        f.seekp(p, std::ios::beg);
                        f.write((char*)&_pa, INT);
                        f.seekg(p + 3 * INT + KEY, std::ios::beg);
                        f.read((char*)&p, INT);
                    }
                    return;
                }
                else {
                    int next_node, parent_count;
                    f.seekg(right_p + 3 * INT + KEY, std::ios::beg);
                    f.read((char*)&next_node, INT);
                    addInsertPos(file[3], right_p);

                    f.seekp(parent + 3 * INT + KEY, std::ios::beg);
                    f.write((char*)&next_node, INT);

                    if (next_node != -1) {
                        f.seekp(next_node + INT, std::ios::beg);
                        f.write((char*)&parent, INT);
                    }

                    int ancestor = getTreeNodeAncestor(file[1], parent);
                    f.seekg(ancestor, std::ios::beg);
                    f.read((char*)&parent_count, INT);
                    --parent_count; f.seekp(ancestor, std::ios::beg);
                    f.write((char*)&parent_count, INT);

                    // modify their ancestors, which is saved in "count"
                    int p = right_d, pa = -pos - 1;
                    while (p != -1) {
                        f.seekp(p, std::ios::beg);
                        f.write((char*)&pa, INT);
                        f.seekg(p + 3 * INT + KEY, std::ios::beg);
                        f.read((char*)&p, INT);
                    }
                    f.seekp(right_d + INT, std::ios::beg);
                    f.write((char*)&end, INT);

                    f.seekp(end + 3 * INT + KEY, std::ios::beg);
                    f.write((char*)&right_d, INT);

                    count += right_count;
                    f.seekp(pos, std::ios::beg);
                    f.write((char*)&count, INT);

                    if (parent_count < splitPos) {
                        checkTreeNode(file[1], ancestor);
                    }
                    return;
                }
            }
            else if (left_p != -1) {
                f.seekg(left_p + 2 * INT + KEY, std::ios::beg);
                f.read((char*)&left_d, INT);

                // find the end node among the current data_nodes
                int p = left_d, end;
                while (p != -1) {
                    end = p;
                    f.seekg(p + 3 * INT + KEY, std::ios::beg);
                    f.read((char*)&p, INT);
                }

                int left_count;
                f.seekg(left_d, std::ios::beg);
                f.read((char*)&left_count, INT);

                if (left_count > splitPos) {
                    f.seekg(end + INT, std::ios::beg);
                    f.read((char*)&p, INT);
                    f.seekp(p + 3 * INT + KEY, std::ios::beg);
                    f.write((char*)&nega_one, INT);
                    --left_count; f.seekp(left_d, std::ios::beg);
                    f.write((char*)&left_count, INT);

                    ++count;
                    f.seekp(end, std::ios::beg);
                    f.write((char*)&count, INT);
                    f.write((char*)&parent, INT);
                    f.seekp(KEY + INT, std::ios::cur);
                    f.write((char*)&pos, INT);

                    Key temp;
                    f.seekg(end + 2 * INT, std::ios::beg);
                    f.read((char*)&temp, KEY);

                    f.seekg(parent + 2 * INT, std::ios::beg);
                    f.write((char*)&temp, KEY);
                    f.write((char*)&end, INT);

                    f.seekp(pos + INT, std::ios::beg);
                    f.write((char*)&end, INT);

                    int p = pos, pa = -end - 1;
                    while (p != -1) {
                        f.seekp(p, std::ios::beg);
                        f.write((char*)&pa, INT);
                        f.seekg(p + 3 * INT + KEY, std::ios::beg);
                        f.read((char*)&p, INT);
                    }
                    return;
                }
                else {
                    addInsertPos(file[3], parent);
                    f.seekp(left_p + 3 * INT + KEY, std::ios::beg);
                    f.write((char*)&nega_one, INT);

                    int ancestor = getTreeNodeAncestor(file[1], left_p), parent_count;
                    f.seekg(ancestor, std::ios::beg);
                    f.read((char*)&parent_count, INT);
                    --parent_count; f.seekp(ancestor, std::ios::beg);
                    f.write((char*)&parent_count, INT);

                    f.seekp(end + 3 * INT + KEY, std::ios::beg);
                    f.write((char*)&pos, INT);
                    f.seekp(pos + INT, std::ios::beg);
                    f.write((char*)&end, INT);

                    int p = pos, pa = -left_d - 1;
                    while (p != -1) {
                        f.seekp(p, std::ios::beg);
                        f.write((char*)&pa, INT);
                        f.seekg(p + 3 * INT + KEY, std::ios::beg);
                        f.read((char*)&p, INT);
                    }

                    count += left_count;
                    f.seekp(left_d, std::ios::beg);
                    f.write((char*)&count, INT);

                    if (parent_count < splitPos) {
                        checkTreeNode(file[1], ancestor);
                        return;
                    }
                }
            }
            else {
                addInsertPos(file[3], parent);

                int parent_parent;
                f.seekg(parent + INT, std::ios::beg);
                f.read((char*)&parent_parent, INT);
                f.seekp(pos + INT, std::ios::beg);
                f.write((char*)&parent_parent, INT);

                if (parent_parent >= 0) {
                    f.seekp(parent_parent + 2 * INT + KEY, std::ios::beg);
                    f.write((char*)&pos, INT);
                }
                else {
                    root = pos;
                }
            }
        }
    }
    
public:
    explicit BPTree(const std::string& name) : name(name) {
        file[1].open(FILE[1], std::ios::binary |
            std::ios::in | std::ios::out | std::ios::app);
        file[2].open(FILE[2], std::ios::binary |
            std::ios::in | std::ios::out | std::ios::app);

        if (exist(FILE[0])) {
            file[0].open(FILE[0], std::ios::binary |
                std::ios::in | std::ios::out | std::ios::app);
            file[0].seekg(std::ios::beg);
            file[0].read((char*)&Size, INT);
            file[0].read((char*)&last_t, INT);
            file[0].read((char*)&last_d, INT);
            file[0].read((char*)&root, INT);
        }
        else {
            last_t += NODE_TREE;
            last_d += NODE_DATA;

            file[0].open(FILE[0], std::ios::binary | 
                std::ios::in | std::ios::out | std::ios::app);
            file[0].seekp(std::ios::beg);
            file[0].write((char*)&Size, INT);   // Size
            file[0].write((char*)&last_t, INT); // last_t
            file[0].write((char*)&last_d, INT); // last_d
            file[0].write((char*)&root, INT);   // root
            
            file[1].seekp(std::ios::beg);
            file[1].write((char*)&one, INT);      // count
            file[1].write((char*)&nega_one, INT);  // parent
            file[1].write((char*)&temp_key, KEY);  // key(NULL)
            file[1].write((char*)&nega_one, INT);  // pointer
            file[1].write((char*)&nega_one, INT);  // next_node

            file[2].seekp(std::ios::beg);
            file[2].write((char*)&zero, INT);  // count
            file[2].write((char*)&zero, INT);  // parent
            file[2].write((char*)&temp_value, VALUE);
            file[2].write((char*)&nega_one, INT);  // next_key
            file[2].write((char*)&nega_one, INT);  // mext_value
        }

        bool VECTOR_T = exist(FILE[3]), VECTOR_D = exist(FILE[4]);
        file[3].open(FILE[3], std::ios::binary |
            std::ios::in | std::ios::out | std::ios::app);
        file[4].open(FILE[4], std::ios::binary |
            std::ios::in | std::ios::out | std::ios::app);
        if (!VECTOR_T) file[3].write((char*)&zero, INT);
        if (!VECTOR_D) file[4].write((char*)&zero, INT);

        closeFile();
        openFile();
    }

    ~BPTree() {
        file[0].close();
        file[1].close();
        file[2].close();
    }

    int size() {
        return Size;
    }

    void insert(const std::pair<Key, T>& val) {
        int pos = findInTree(root, val);
        insertToData(pos, val);
        updateinfo();
    }

    sjtu::vector<T> Find(const Key& key) {
        int pos = findInTree(root, value_type(key, 0));
        return findValue(pos, key);
    }

    void remove(const std::pair<Key, T>& val) {
        int pos = findInTree(root, val);
        removeData(pos, val);
        updateinfo();
    }

    // reset the state of BPT
    void clear() {
        closeFile();
        clearFile();
        closeFile();

        Size = last_t = last_d = root = 0;
        last_t += NODE_TREE;
        last_d += NODE_DATA;

        openFile();
        file[0].seekp(std::ios::beg);
        file[0].write((char*)&Size, INT);   // Size
        file[0].write((char*)&last_t, INT); // last_t
        file[0].write((char*)&last_d, INT); // last_d
        file[0].write((char*)&root, INT);   // root

        file[1].seekp(std::ios::beg);
        file[1].write((char*)&one, INT);      // count
        file[1].write((char*)&nega_one, INT);  // parent
        file[1].write((char*)&temp_key, KEY);  // key(NULL)
        file[1].write((char*)&nega_one, INT);  // pointer
        file[1].write((char*)&nega_one, INT);  // next_node

        file[2].seekp(std::ios::beg);
        file[2].write((char*)&zero, INT);  // count
        file[2].write((char*)&zero, INT);  // parent
        file[2].write((char*)&temp_value, VALUE);
        file[2].write((char*)&nega_one, INT);  // next_key
        file[2].write((char*)&nega_one, INT);  // mext_value

        file[3].write((char*)&zero, INT);
        file[4].write((char*)&zero, INT);
    }

    // TODO: modify the second value of "val" with new_val
    void modify(const std::pair<Key, T>& val, T new_val) {
        
    }

    // TODO: simplify the "find" operation
    // return the first inserted value of a specific key
    std::pair<bool, T> find(const Key& key) {
        sjtu::vector<T> temp = Find(key);
        if (temp.empty()) return std::pair<bool, T>(false, T());
        else return std::pair<bool, T>(true, temp.front());
    }

    /*void print() {
        std::cout << std::endl << Size << std::endl;

        std::queue<int> que;
        que.push(root);

        while (que.front() >= 0) {
            Key temp;

            int cur = que.front(), count;
            que.pop();

            file[1].seekg(cur, std::ios::beg);
            file[1].read((char*)&count, INT);
            std::cout << count << std::endl;

            int next_pos = cur, pointer;
            do {
                file[1].seekg(next_pos + 2 * INT, std::ios::beg);
                file[1].read((char*)&temp, KEY);
                file[1].read((char*)&pointer, INT);
                file[1].read((char*)&next_pos, INT);
                que.push(pointer);
                std::cout << temp << ' ';
            } while (next_pos != -1);
            std::cout << std::endl;
        }
        std::cout << std::endl;

        while (!que.empty()) {
            value_type temp;

            int cur = -que.front() - 1, count;
            que.pop();

            file[2].seekg(cur, std::ios::beg);
            file[2].read((char*)&count, INT);
            file[2].seekg(cur + 2 * INT + VALUE, std::ios::beg);
            file[2].read((char*)&cur, INT);
            std::cout << count << std::endl;

            if (count) {
                do {
                    file[2].seekg(cur + 2 * INT, std::ios::beg);
                    file[2].read((char*)&temp, VALUE);
                    file[2].read((char*)&cur, INT);
                    std::cout << temp.first << ' ';
                } while (cur != -1);
            }
            std::cout << std::endl;
        }
    }*/
};

#endif //BPTREE_HPP_BPTREE2_HPP
