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
template<class T, size_t N> inline bool olsen_pending(const std::array<T,N>&) { return false; }
template<class T, size_t N> inline void mark_olsen_pending(std::array<T,N>&) {}
template<class T, size_t N> inline void clear_olsen_pending(std::array<T,N>&) {}

inline uint64_t bits(double x) { uint64_t u; std::memcpy(&u,&x,8); return u; }
inline double real(uint64_t u) { double x; std::memcpy(&x,&u,8); return x; }
// The two existing metadata bits encode: 00 uncached, 10 external,
// 11 internal/current, 01 internal/queued for the next Olsen report.
// Both internal states keep the same decoded H_ii and c/d ownership.
inline bool has_hii(const std::array<double,3>& v) { return (bits(v[2]) & uint64_t(3)) != 0; }
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
    if(has_hii(v)) {
        const uint64_t u = bits(v[2]);
        const uint64_t flags = u & uint64_t(3);
        // A selected internal coordinate may already be queued. Preserve its
        // marker so subsequent b updates do not enqueue the same row again.
        const uint64_t new_flags = owner ? ((flags & 1u) ? flags : 3u) : 2u;
        v[2] = real((u & ~uint64_t(3)) | new_flags);
    }
}
inline bool olsen_pending(const std::array<double,3>& v) {
    return (bits(v[2]) & uint64_t(3)) == 1u;
}
inline void mark_olsen_pending(std::array<double,3>& v) {
    if ((bits(v[2]) & uint64_t(3)) == 3u)
        v[2] = real(bits(v[2]) & ~uint64_t(2));
}
inline void clear_olsen_pending(std::array<double,3>& v) {
    if (olsen_pending(v)) v[2] = real(bits(v[2]) | uint64_t(2));
}
}
#endif
