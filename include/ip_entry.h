#ifndef CDFCI_IP_ENTRY_H
#define CDFCI_IP_ENTRY_H

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

namespace ipentry {
template<class T, size_t N> inline T c(const std::array<T,N>& v) { return v[0]; }
template<class T, size_t N> inline T d(const std::array<T,N>&) { return T(0); }
template<class T, size_t N> inline T b(const std::array<T,N>& v) { return v[N/2]; }
template<class T, size_t N> inline bool is_c(const std::array<T,N>& v) { return v[0] != T(0); }
template<class T, size_t N> inline bool has_hii(const std::array<T,N>&) { return false; }
template<class T, size_t N> inline T hii(const std::array<T,N>&) { return std::numeric_limits<T>::quiet_NaN(); }
template<class T, size_t N> inline void set_hii(std::array<T,N>&, T, bool) {}
template<class T, size_t N> inline void set_owner(std::array<T,N>&, bool) {}

inline uint64_t bits(double x) { uint64_t u; std::memcpy(&u,&x,8); return u; }
inline double real(uint64_t u) { double x; std::memcpy(&x,&u,8); return x; }
inline bool has_hii(const std::array<double,3>& v) { return (bits(v[2]) & uint64_t(2)) != 0; }
inline bool is_c(const std::array<double,3>& v) {
    return has_hii(v) ? (bits(v[2]) & 1u) != 0 : v[0] != 0.0;
}
inline double c(const std::array<double,3>& v) { return is_c(v) ? v[0] : 0.0; }
inline double d(const std::array<double,3>& v) { return is_c(v) ? 0.0 : v[0]; }
inline double b(const std::array<double,3>& v) { return v[1]; }
inline double hii(const std::array<double,3>& v) { return real(bits(v[2]) & ~uint64_t(3)); }
inline void set_hii(std::array<double,3>& v, double h, bool owner) {
    v[2]=real((bits(h)&~uint64_t(3))|uint64_t(2)|uint64_t(owner));
}
inline void set_owner(std::array<double,3>& v, bool owner) {
    if(has_hii(v)) set_hii(v,hii(v),owner);
}
}
#endif
