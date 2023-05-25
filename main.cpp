#include <iostream>
#include <cstdio>
#include <cstring>
#include "vector.hpp"
#include "bptree.hpp"

struct String {
    char index[65];

    String(const String& other) {
        for (int i = 0; i < 65; i++) index[i] = other.index[i];
    }

    String() = default;

    friend bool operator>(const String& lhs, const String& rhs) {
        return std::string(lhs.index) > std::string(rhs.index);
    }

    friend bool operator>=(const String& lhs, const String& rhs) {
        return std::string(lhs.index) >= std::string(rhs.index);
    }

    friend bool operator<(const String& lhs, const String& rhs) {
        return std::string(lhs.index) < std::string(rhs.index);
    }

    friend bool operator<=(const String& lhs, const String& rhs) {
        return std::string(lhs.index) <= std::string(rhs.index);
    }

    friend bool operator==(const String& lhs, const String& rhs) {
        return std::string(lhs.index) == std::string(rhs.index);
    }

    friend bool operator!=(const String& lhs, const String& rhs) {
        return std::string(lhs.index) != std::string(rhs.index);
    }

    friend std::ostream& operator<<(std::ostream& os, const String& obj) {
        os << obj.index;
        return os;
    }
};

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0); std::cout.tie(0);

    BPTree<String, int> bpTree("test");
    std::pair<String, int> val;
    int cnt;
    char cmd[10];
    std::cin >> cnt;

    for (int i = 1; i <= cnt; i++) {
        std::cin >> cmd;
        if (cmd[0] == 'i') {
            std::cin >> val.first.index >> val.second;
            bpTree.insert(val);
        }
        else if (cmd[0] == 'f') {
            std::cin >> val.first.index;
            sjtu::vector<int> ans = bpTree.Find(val.first);
            if (!ans.empty()) {
                for (int i = 0; i < ans.size() - 1; i++) std::cout << ans[i] << ' ';
                std::cout << ans[ans.size() - 1] << std::endl;
            }
            else puts("null");
        }
        else if (cmd[0] == 'd') {
            std::cin >> val.first.index >> val.second;
            bpTree.remove(val);
        }
        //bpTree.print();
    }
    return 0;
}
