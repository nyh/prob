#undef TRACE

#ifdef TRACE
#include <iostream>
#include <list>
#include <vector>

namespace std {

template <typename T>
static inline
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    bool first = true;
    os << "{";
    for (auto&& elem : v) {
        if (!first) {
            os << ", ";
        } else {
            first = false;
        }
        os << elem;
    }
    os << "}";
    return os;
}

template <typename T>
static inline
std::ostream& operator<<(std::ostream& os, const std::list<T>& v) {
    bool first = true;
    os << "{";
    for (auto&& elem : v) {
        if (!first) {
            os << ", ";
        } else {
            first = false;
        }
        os << elem;
    }
    os << "}";
    return os;
}

template <typename T1, typename T2>
static inline
std::ostream& operator<<(std::ostream& os, const std::pair<T1, T2>& v) {
    os << "{" << v.first << ", " << v.second << "}";
    return os;
}
}
#endif
