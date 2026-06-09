#pragma once

#include "connectivity.h"

inline void ConnectivityList::add(const Connectivity& conn) {
    connections_.push_back(conn);
}

inline Int ConnectivityList::count() const {
    return static_cast<Int>(connections_.size());
}

inline const Connectivity& ConnectivityList::operator[](Int i) const {
    return connections_[static_cast<std::size_t>(i)];
}

inline Connectivity& ConnectivityList::operator[](Int i) {
    return connections_[static_cast<std::size_t>(i)];
}
