//   ______  _______   _______   ______  __
//  /      ||       \ |   ____| /      ||  |
// |  ,----'|  .--.  ||  |__   |  ,----'|  |
// |  |     |  |  |  ||   __|  |  |     |  |
// |  `----.|  '--'  ||  |     |  `----.|  |
//  \______||_______/ |__|      \______||__|
//
// Coordinate Descent Full Configuration Interaction (CDFCI) package in C++17
// https://github.com/quan-tum/CDFCI
//
// Copyright (c) 2019-2026, CDFCI Developers and Contributors
// All rights reserved.
//
// This source code is licensed under the BSD 3-Clause License found in the
// LICENSE file in the root directory of this source tree.

#ifndef CDFCI_WAVEFUNCTION_H
#define CDFCI_WAVEFUNCTION_H 1

#include "config.h"
#include "vecmath.h"
#include "container.h"
#include "determinant.h"
#include "hamiltonian.h"
#include "ip_entry.h"
#include <cmath>
#include <type_traits>
#include <unordered_map>
#include <vector>

#define DEFINE_ACCESSORS(NAME, MEMBER)                                 \
    data_matrix_type get_##NAME##_double() const {                     \
         return to_double_precision(MEMBER); }                         \
    quad_data_matrix_type get_##NAME() const { return MEMBER; }        \
    void set_##NAME(quad_data_matrix_type val) { MEMBER = val; }       \
    void set_##NAME(data_matrix_type val) {                            \
        MEMBER = to_quad_precision(val); }                             \
    void update_##NAME(quad_data_matrix_type inc) {                    \
         add_assign(MEMBER, inc); }                                    \
    void update_##NAME(data_matrix_type inc) {                         \
         add_assign(MEMBER, to_quad_precision(inc)); }

#define IDX(i, j) ((i) + (j) * (NSTATES))
#define X(i) ((i))
#define Z(i) ((i) + (NSTATES))

using namespace vecmath;

template <typename Container, int NSTATES = 1>
class WaveFunctionBase
{
public:
    using key_type = typename Container::key_type;
    using mapped_type = typename Container::mapped_type; // {x, z}
    using hasher = typename Container::hasher;
    using key_equal = typename Container::key_equal;
    using iterator = typename Container::iterator;
    using const_iterator = typename Container::const_iterator;

    using data_type = std::conditional_t<(NSTATES == 1), NumericalType, std::array<NumericalType, NSTATES>>;
    using data_matrix_type = std::conditional_t<(NSTATES == 1), NumericalType, std::array<NumericalType, NSTATES * NSTATES>>;
    using quad_data_type = std::conditional_t<(NSTATES == 1), QUAD_PRECISION, std::array<QUAD_PRECISION, NSTATES>>;
    using quad_data_matrix_type = std::conditional_t<(NSTATES == 1), QUAD_PRECISION, std::array<QUAD_PRECISION, NSTATES * NSTATES>>;

    static constexpr int NSTATES_val = NSTATES; // number of states stored
    int num_states = 1; // number of states desired

protected:
    /* Data */
    // Store the wavefunction {i, {x_i, z_i}} where z approximates Hx
    Container data_;
    // {x'x, z'z}
    std::array<quad_data_matrix_type, 2> norm_square_ = {};
    // x'z
    quad_data_matrix_type dot_product_ = {};
    // scale factor
    ScaleFactorType scale_factor_ = static_cast<ScaleFactorType>(1.0);
    // energy estimate
    data_type energy_ = {};

    // number of nonzeros
    std::array<size_t, 2> size_ = {};

    size_t capacity_ = 0;

public:
    /* Constructors */
    WaveFunctionBase(size_t capacity = 16)
    {
        data_.reserve(capacity);
        capacity_ = capacity;
    }

    /* Destructor */
    virtual ~WaveFunctionBase() {}

    virtual void clear()
    {
        data_.clear();
        norm_square_ = {};
        dot_product_ = {};
        scale_factor_ = static_cast<ScaleFactorType>(1.0);
        size_ = {};
    }

    void reserve_capacity(size_t capacity = 16)
    {
        //TODO: return a flag indicating whether reserve is successful.
        data_.reserve(capacity);
        capacity_ = capacity;
    }

    /* Iterators */
    // WARNING:
    //   Iterator may be invalidated by other threads if
    // running in parallel.
    const_iterator begin() const { return data_.begin(); }

    const_iterator end() const { return data_.end(); }

    iterator begin() { return data_.begin(); }

    iterator end() { return data_.end(); }

    /* Other member functions */
    size_t size_x() const { return size_[0]; }
    void update_size_x(size_t n) { size_[0] += n; }

    size_t size_z() const { return size_[1]; }
    void update_size_z(size_t n) { size_[1] += n; }

    size_t size() const { return data_.size(); }

    DEFINE_ACCESSORS(xx, norm_square_[0])
    DEFINE_ACCESSORS(zz, norm_square_[1])
    DEFINE_ACCESSORS(xz, dot_product_)

    ScaleFactorType get_scale() const { return scale_factor_; }
    void set_scale(ScaleFactorType s) { scale_factor_ = s; }
    void set_scale(NumericalType s) { scale_factor_ = (ScaleFactorType)s; }

    size_t capacity() const { return capacity_; }

    void compute_variational_energy()
    {
        if constexpr (NSTATES == 1) {
            energy_ = dot_product_ / norm_square_[0];
        }
        else {
            Eigen::MatrixXd A(num_states, num_states);
            Eigen::MatrixXd B(num_states, num_states);

            for (int i = 0; i < num_states; ++i)
            for (int j = 0; j < num_states; ++j)
            {
                A(i, j) = dot_product_[IDX(i, j)];
                B(i, j) = norm_square_[0][IDX(i, j)];
            }

            Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(A, B);

            if (solver.info() != Eigen::Success) {
                throw std::runtime_error("General Eigensolver failed");
            }

            for (int i = 0; i < num_states; ++i)
                energy_[i] = solver.eigenvalues()[i];
        }
    }

    data_type get_variational_energy() const {  return energy_; }

    data_type get_x(const key_type &det) const
    {
        data_type x = {};
        auto val = data_.find_val(det);
        if (val)
            if constexpr (NSTATES == 1) x = ipentry::c(*val);
            else x = slice_first_half(*val);
        return x;
    }

    data_type get_z(const key_type &det) const
    {
        data_type z = {};
        auto val = data_.find_val(det);
        if (val)
            z = slice_second_half(*val);
        return z;
    }

    data_type get_grad(const data_type &x, const data_type &z) const
    {
        data_matrix_type xx = get_xx_double();
        if constexpr (NSTATES == 1) {
            return (x * xx + z);
        } else {
            data_type g = {};
            for (int i = 0; i < NSTATES; ++i) {
                g[i] = z[i];
                for (int k = 0; k < NSTATES; ++k)
                    g[i] += x[k] * xx[IDX(k, i)];
            }
            return g;
        }
    }

    /* Helper function */
    quad_data_matrix_type to_quad_precision(const data_matrix_type &input) const {
        if constexpr (NSTATES == 1) {
            return static_cast<QUAD_PRECISION>(input);
        } else {
            quad_data_matrix_type result = {};
            for (size_t i = 0; i < NSTATES * NSTATES; ++i) {
                result[i] = static_cast<QUAD_PRECISION>(input[i]);
            }
            return result;
        }
    }

    data_matrix_type to_double_precision(const quad_data_matrix_type &input) const {
        if constexpr (NSTATES == 1) {
            return static_cast<NumericalType>(input);
        } else {
            data_matrix_type result = {};
            for (size_t i = 0; i < NSTATES * NSTATES; ++i) {
                result[i] = static_cast<NumericalType>(input[i]);
            }
            return result;
        }
    }

    /* Debug function */
    void print_xx_xz(data_matrix_type &xx, data_matrix_type &xz) {
        if constexpr (NSTATES == 1) {
            std::cout << "xx = " << static_cast<NumericalType>(xx)
                      << std::endl;
            std::cout << "xz = " << static_cast<NumericalType>(xz)
                      << std::endl;
        } else {
            std::cout << "xx = " << std::endl;
            for (size_t i = 0; i < NSTATES; i++) {
                for (size_t j = 0; j < NSTATES; j++) {
                    std::cout << static_cast<NumericalType>(xx[IDX(i, j)]) << " ";
                }
                std::cout << std::endl;
            }
            std::cout << "xz = " << std::endl;
            for (size_t i = 0; i < NSTATES; i++) {
                for (size_t j = 0; j < NSTATES; j++) {
                    std::cout << static_cast<NumericalType>(xz[IDX(i, j)]) << " ";
                }
                std::cout << std::endl;
            }
        }
    }

    void print_xx_xz() {
        auto xx = to_double_precision(norm_square_[0]);
        auto xz = to_double_precision(dot_product_);
        print_xx_xz(xx, xz);
    }

    void calculate_xx_xz() {
        data_matrix_type xx = {};
        data_matrix_type xz = {};
        for (auto it = begin(); it != end(); ++it) {
            if constexpr (NSTATES == 1) {
                const auto x = ipentry::c(it->second);
                xx += x * x;
                xz += x * ipentry::b(it->second);
            } else {
                for (size_t i = 0; i < NSTATES; i++)
                for (size_t j = 0; j < NSTATES; j++)
                {
                    xx[IDX(i, j)] += it->second[X(i)] * it->second[X(j)];
                    xz[IDX(i, j)] += it->second[X(i)] * it->second[Z(j)];
                }
            }
        }
        print_xx_xz(xx, xz);
    }
};

template <typename Container, int NSTATES = 1>
class WaveFunctionFragment : public WaveFunctionBase<Container, NSTATES>
{
public:
    double ip_pt_delta = 0;
    size_t ip_diagonal_evaluations = 0;
    using Base = WaveFunctionBase<Container, NSTATES>;
    using typename Base::key_type;
    using typename Base::mapped_type;
    using typename Base::hasher;
    using typename Base::key_equal;
    using typename Base::iterator;
    using typename Base::const_iterator;
    using value_type = typename Container::value_type;

protected:
    using Base::data_;
    using Base::norm_square_;
    using Base::dot_product_;
    using Base::scale_factor_;
    using Base::size_;
    using Base::capacity_;

public:
    WaveFunctionFragment(size_t capacity_ = 16)
        : Base(capacity_){}

    void clear()
    {
        Base::clear();
        ip_pt_delta = 0;
        ip_diagonal_evaluations = 0;
    }

    void push_back(const value_type &val) { data_.push_back(val); } // const reference
    void push_back(value_type &&val) { data_.push_back(std::move(val)); } // rvalue reference

    void push_back(const key_type &det, const mapped_type &val) { data_.emplace_back(det, val); }
    void push_back(key_type &&det, mapped_type &&val) { data_.emplace_back(std::move(det), std::move(val)); }

    void append(WaveFunctionFragment<Container, NSTATES> &sub_xz)
    {
        data_.append(sub_xz);
        ip_pt_delta += sub_xz.ip_pt_delta;
        ip_diagonal_evaluations += sub_xz.ip_diagonal_evaluations;
    }
};

template <typename Container, int NSTATES = 1, bool COMPACT_IP = false>
class WaveFunction : public WaveFunctionBase<Container, NSTATES>
{
public:
    static constexpr bool uses_compact_ip = COMPACT_IP;
    static constexpr bool stores_full_hii =
        std::tuple_size<typename Container::mapped_type>::value == 3;
    static constexpr bool supports_ip_correction =
        COMPACT_IP || stores_full_hii;
    using Base = WaveFunctionBase<Container, NSTATES>;
    using typename Base::key_type;
    using typename Base::mapped_type;
    using typename Base::hasher;
    using typename Base::key_equal;
    using typename Base::iterator;
    using typename Base::const_iterator;
    using ContainerVectorType =
        ContainerVector<key_type, mapped_type, hasher, key_equal>;
    using wff_type = WaveFunctionFragment<ContainerVectorType, NSTATES>;
    struct InternalCorrectionRow {
        NumericalType c, b, hii;
    };

    using typename Base::data_type;
    using typename Base::data_matrix_type;
    using typename Base::quad_data_type;
    using typename Base::quad_data_matrix_type;

protected:
    /* Data */
    using Base::data_;
    using Base::norm_square_;
    using Base::dot_product_;
    using Base::scale_factor_;
    using Base::size_;
    using Base::capacity_;

    NumericalType max_load_factor_;

    bool ip_enabled_ = false;
    bool ip_internal_cache_enabled_ = false;
    QUAD_PRECISION ip_pt_energy_ = 0;
    size_t ip_diagonal_evaluations_ = 0;
    std::vector<InternalCorrectionRow> ip_internal_rows_;
    std::unordered_map<key_type, size_t, hasher, key_equal> ip_internal_index_;

    bool is_internal(const key_type &det, const mapped_type &val) const
    {
        if constexpr (COMPACT_IP && !stores_full_hii) {
            if (ip_enabled_)
                return ip_internal_index_.find(det) != ip_internal_index_.end();
        }
        return ipentry::is_c(val);
    }

    NumericalType external_d(const mapped_type &val) const
    {
        if constexpr (COMPACT_IP) return val[0];
        else return ipentry::d(val);
    }

    template<typename H>
    void sync_ip_internal_row(const key_type &det, const mapped_type &val,
                              const H &ham, bool known_diagonal = false,
                              NumericalType diagonal = 0)
    {
        auto [it, inserted] = ip_internal_index_.emplace(det, ip_internal_rows_.size());
        if (!ip_internal_cache_enabled_) return;
        if (inserted) {
            if (!known_diagonal) {
                auto diagonal_det = det;
                if (ipentry::has_hii(val)) {
                    diagonal = ipentry::hii(val);
                } else {
                    diagonal = ham.get_diagonal(diagonal_det);
                    ++ip_diagonal_evaluations_;
                }
            }
            ip_internal_rows_.push_back({ipentry::c(val), ipentry::b(val), diagonal});
        } else {
            auto &row = ip_internal_rows_[it->second];
            row.c = ipentry::c(val);
            row.b = ipentry::b(val);
        }
    }

    void sync_ip_internal_b(const key_type &det, NumericalType b)
    {
        // The index is populated by update_x before parallel z updates start.
        // During those updates it is read-only; the main table's per-key update
        // lock also serializes writes to the corresponding dense row.
        const auto &index = ip_internal_index_;
        const auto it = index.find(det);
        if (it == index.end())
            throw std::logic_error("Internal correction row is missing");
        ip_internal_rows_[it->second].b = b;
    }

public:
    /* Constructors */
    WaveFunction(){}
    WaveFunction(size_t capacity_, NumericalType max_load_factor = 0.79)
        : Base(capacity_)
    { max_load_factor_ = max_load_factor;}

    /* Destructor */
    virtual ~WaveFunction() {}

    void clear() override
    {
        Base::clear();
        configure_ip_correction(false);
    }

    void configure_ip_correction(bool enabled)
    {
        ip_enabled_ = enabled;
        ip_internal_cache_enabled_ = false;
        ip_pt_energy_ = 0;
        ip_diagonal_evaluations_ = 0;
        ip_internal_rows_.clear();
        ip_internal_index_.clear();
    }

    const std::vector<InternalCorrectionRow> &ip_internal_rows() const
    {
        return ip_internal_rows_;
    }

    NumericalType stored_c(const key_type &det, const mapped_type &val) const
    {
        return is_internal(det, val) ? val[0] : 0.0;
    }

    data_type get_x(const key_type &det) const
    {
        if constexpr (COMPACT_IP) {
            const auto val = data_.find_val(det);
            return val ? stored_c(det, *val) : data_type{};
        }
        return Base::get_x(det);
    }

    template<typename H>
    void initialize_ip_correction(const H &ham, bool cache_internal = true)
    {
        configure_ip_correction(true);
        ip_internal_cache_enabled_ = cache_internal;
        if constexpr (NSTATES == 1 && supports_ip_correction) {
            const QUAD_PRECISION norm = this->get_xx();
            if (!(norm > 0) || !std::isfinite(static_cast<double>(norm))) {
                return;
            }
            const double energy = static_cast<double>(this->get_xz() / norm);

            // Cache diagonals and initialize the external PT2 amplitudes in one pass.
            for (auto it = this->begin(); it != this->end(); ++it) {
                auto det = it->first;
                auto &val = it->second;
                const bool owner = ipentry::is_c(val);
                const double z = ipentry::b(val);
                const double diagonal = ham.get_diagonal(det);
                ++ip_diagonal_evaluations_;
                ipentry::set_hii(val, diagonal, owner);

                if (owner) {
                    if ((COMPACT_IP && !stores_full_hii) || ip_internal_cache_enabled_)
                        sync_ip_internal_row(det, val, ham, true, diagonal);
                } else {
                    const double gap = diagonal - energy;
                    val[0] = -z / gap;
                    ip_pt_energy_ += static_cast<QUAD_PRECISION>(z) * val[0];
                }
            }

        }
    }

    bool ip_enabled() const { return ip_enabled_; }

    EnergyCorrectionResult ip_correction_result(size_t iteration) const
    {
        EnergyCorrectionResult out;
        out.iteration = iteration;
        out.diagonal_evaluations = ip_diagonal_evaluations_;
        const QUAD_PRECISION norm = this->get_xx();
        if (!(norm > 0) || !std::isfinite(static_cast<double>(norm))) {
            out.status = "invalid_norm";
            return out;
        }
        const QUAD_PRECISION energy = this->get_xz() / norm;
        out.variational_energy = static_cast<double>(energy);
        out.external_correction = static_cast<double>(ip_pt_energy_ / norm);
        out.external_valid = true;
        out.internal_correction = 0.0;
        out.internal_valid = true;
        out.valid = true;
        out.corrected_energy = out.variational_energy + out.external_correction;
        out.status = "unchecked";
        return out;
    }

    NumericalType max_load_factor() const { return max_load_factor_; }

    void reserve_capacity(size_t capacity = 16, NumericalType max_load_factor = 0.79)
    {
        data_.reserve(capacity);
        capacity_ = capacity;
        max_load_factor_ = max_load_factor;
    }

    /* Update x or z */
    template<typename H>
    void update_x(key_type &det, data_type dx, H &ham)
    {
        if constexpr (NSTATES != 1) {
            mapped_type vec_xz = {};
            mapped_type val_new = {};
            add_assign_first_half(val_new, dx);
            int new_element = 1;
            auto update_functor = [dx, &vec_xz, &new_element](mapped_type &val) {
                new_element = is_zero_on_first_half(val);
                vec_xz = val;
                add_assign_first_half(val, dx);
                return false;
            };
            upsert(data_, det, update_functor, val_new);
            data_matrix_type delta_xx = {}, delta_xz = {};
            for (size_t i = 0; i < NSTATES; i++)
                for (size_t j = 0; j < NSTATES; j++) {
                    delta_xx[IDX(i, j)] = dx[i] * vec_xz[X(j)] + vec_xz[X(i)] * dx[j] + dx[i] * dx[j];
                    delta_xz[IDX(i, j)] = dx[i] * vec_xz[Z(j)];
                }
            this->update_xx(delta_xx);
            this->update_xz(delta_xz);
            this->update_size_x(new_element);
            return;
        } else {
        mapped_type vec_xz = {};
        mapped_type val_new = {};
        mapped_type updated_owner = {};
        add_assign_first_half(val_new, dx);
        int new_element = 1;
        auto update_functor = [&, dx](mapped_type &val) {
            const bool old_c = is_internal(det, val);
            const double old_x = old_c ? val[0] : 0.0;
            const double old_z = ipentry::b(val);
            if constexpr (stores_full_hii) {
                if (ip_enabled_ && !ipentry::has_hii(val)) {
                    ipentry::set_hii(val, ham.get_diagonal(det), old_c);
                    ++ip_diagonal_evaluations_;
                }
            }
            if (ip_enabled_ && !old_c)
                ip_pt_energy_ -= static_cast<QUAD_PRECISION>(old_z) * external_d(val);
            new_element = !old_c;
            vec_xz = val;
            val[0] = old_x + dx;
            ipentry::set_owner(val, true);
            updated_owner = val;
            return false; // no deletion
        };
        // Reuse a selected row's stored diagonal when available.
        if (ip_enabled_) {
            if (!update_fn(data_, det, update_functor)) {
                if constexpr (stores_full_hii) {
                    ipentry::set_hii(val_new, ham.get_diagonal(det), true);
                    ++ip_diagonal_evaluations_;
                }
                if (upsert(data_, det, update_functor, val_new))
                    updated_owner = val_new;
            }
        } else {
            upsert(data_, det, update_functor, val_new);
        }
        if (ip_enabled_ && ((COMPACT_IP && !stores_full_hii) || ip_internal_cache_enabled_))
            sync_ip_internal_row(det, updated_owner, ham);

        // Update xx, xz and size
        // (x + dx)'*(x + dx) = x'*x + dx'*x + x'*dx +  dx'*dx
        // (x + dx)'*(z + dz) = x'*z + dx'*z + (x+dx)'*dz (the third part updated later)
        if constexpr (NSTATES == 1) {
            const auto old_x = is_internal(det, vec_xz) && !new_element
                ? ipentry::c(vec_xz) : 0.0;
            this->update_xx(2.0 * old_x * dx + dx * dx);
            this->update_xz(ipentry::b(vec_xz) * dx);
        } else {
            // update xx
            data_matrix_type delta_xx = {};
            for (size_t i = 0; i < NSTATES; i++)
                for (size_t j = 0; j < NSTATES; j++)
                    delta_xx[IDX(i, j)] = dx[i] * vec_xz[X(j)] + vec_xz[X(i)] * dx[j] + dx[i] * dx[j];
            this->update_xx(delta_xx);
            // update xz
            data_matrix_type delta_xz = {};
            for (size_t i = 0; i < NSTATES; i++)
                for (size_t j = 0; j < NSTATES; j++)
                    delta_xz[IDX(i, j)] = dx[i] * vec_xz[Z(j)];
            this->update_xz(delta_xz);
        }
        this->update_size_x(new_element);
        return;
        }
    }

    // Update z and store the updated wavefunction into sub_xz.
    // Return new_z.
    // Note:
    //   This function may be called in parallel. So dot_product_
    // and size_ are updated into sub_xz not *this, the wavefunction
    // itself, to avoid confliction. dot_product_ and size_
    // of *this need to be updated outside this function.
    //   Assume sub_xz is properly initialized.
    //   New z is inserted only if abs(z) > z_threshould.
    template<typename C, typename H>
    data_type update_z(C &column, data_type dx, H &ham, NumericalType e0,
                   wff_type &sub_xz,
                  NumericalType z_threshold = 0.0)
    {
        if constexpr (NSTATES != 1) {
            data_type new_z = {};
            data_matrix_type delta_xz = {};
            NumericalType scale = this->get_scale();
            for (auto &entry : column) {
                auto &det = entry.first;
                auto h = entry.second;
                data_type dz = multiply(dx, h);
                mapped_type vec_xz = {}, val_new = {};
                add_assign_second_half(vec_xz, dz);
                add_assign_second_half(val_new, dz);
                int new_element = 1;
                auto update_functor = [dz, &vec_xz, &new_element](mapped_type &val) {
                    new_element = is_zero_on_second_half(val);
                    add_assign_second_half(val, dz);
                    vec_xz = val;
                    return false;
                };
                auto dz_norm = max_norm(multiply(dz, scale));
                if (dz_norm > z_threshold) {
                    upsert(data_, det, update_functor, val_new);
                    sub_xz.update_size_z(new_element);
                    for (size_t i = 0; i < NSTATES; i++)
                        for (size_t j = 0; j < NSTATES; j++)
                            delta_xz[IDX(i,j)] += vec_xz[X(i)] * dz[j];
                    sub_xz.push_back(det, vec_xz);
                } else if (update_fn(data_, det, update_functor)) {
                    for (size_t i = 0; i < NSTATES; i++)
                        for (size_t j = 0; j < NSTATES; j++)
                            delta_xz[IDX(i,j)] += vec_xz[X(i)] * dz[j];
                    sub_xz.push_back(det, vec_xz);
                }
                add_assign(new_z, multiply(slice_first_half(vec_xz), h));
            }
            sub_xz.update_xz(delta_xz);
            return new_z;
        } else {
        data_type new_z = {}; // Store the recalculated z_i.
        data_matrix_type delta_xz = {};
        NumericalType scale = this->get_scale();
        double pt_delta = 0;
        size_t diagonal_evaluations = 0;
        // Loop over the column
        for (auto &entry : column)
        {
            auto &det = entry.first;
            auto h = entry.second;
            data_type dz = multiply(dx, h);

            mapped_type vec_xz = {};
            mapped_type val_new = {};
            add_assign_second_half(vec_xz, dz);
            add_assign_second_half(val_new, dz);

            int new_element = 1;
            auto update_functor = [&, dz](mapped_type &val) {
                const bool owner = is_internal(det, val);
                const double old_d = owner ? 0.0 : external_d(val);
                const double old_z = ipentry::b(val);
                if constexpr (stores_full_hii) {
                    if (ip_enabled_ && !ipentry::has_hii(val)) {
                        ipentry::set_hii(val, ham.get_diagonal(det), owner);
                        ++diagonal_evaluations;
                    }
                }
                new_element = is_zero_on_second_half(val);
                val[1] = old_z + dz;
                if (ip_internal_cache_enabled_ && owner)
                    sync_ip_internal_b(det, val[1]);
                if (ip_enabled_ && !owner) {
                    double diagonal;
                    if constexpr (stores_full_hii) {
                        diagonal = ipentry::hii(val);
                    } else {
                        diagonal = ham.get_diagonal(det);
                        ++diagonal_evaluations;
                    }
                    const double gap = diagonal - e0;
                    const double new_d = -val[1] / gap;
                    pt_delta += val[1]*new_d - old_z*old_d;
                    val[0] = new_d;
                }
                vec_xz = val;
                if constexpr (COMPACT_IP) {
                    // Candidate fragments carry c (zero for external rows),
                    // while the main map keeps c+d for incremental PT2.
                    if (!owner) vec_xz[0] = 0.0;
                }
                return false;
            };
            auto dz_norm = max_norm(multiply(dz, scale));
            // Existing external rows reuse the cached diagonal; only newly
            // stored rows evaluate one for the first time.
            const bool found = ip_enabled_ && update_fn(data_, det, update_functor);
            if (dz_norm > z_threshold)
            {
                double inserted_pt_delta = 0;
                if (!found) {
                    // The disabled path retains the original single upsert.
                    if (ip_enabled_) {
                        const double hd = ham.get_diagonal(det);
                        ++diagonal_evaluations;
                        if constexpr (stores_full_hii)
                            ipentry::set_hii(val_new, hd, false);
                        val_new[1] = dz;
                        const double gap = hd - e0;
                        val_new[0] = -dz/gap;
                        inserted_pt_delta = dz*val_new[0];
                    }
                    const bool inserted = upsert(data_, det, update_functor, val_new);
                    if (inserted) pt_delta += inserted_pt_delta;
                }

                sub_xz.update_size_z(new_element);
                if constexpr (NSTATES == 1) {
                    delta_xz += ipentry::c(vec_xz) * dz;
                } else {
                    for (size_t i = 0; i < NSTATES; i++)
                        for (size_t j = 0; j < NSTATES; j++)
                            delta_xz[IDX(i, j)] += vec_xz[X(i)] * dz[j];
                }
                sub_xz.push_back(det, vec_xz);
            }
            else
            {
                bool exist_z_flag = ip_enabled_ ? found : update_fn(data_, det, update_functor);
                if (exist_z_flag)
                {
                    if constexpr (NSTATES == 1) {
                        delta_xz += ipentry::c(vec_xz) * dz;
                    } else {
                        for (size_t i = 0; i < NSTATES; i++)
                            for (size_t j = 0; j < NSTATES; j++)
                                delta_xz[IDX(i, j)] += vec_xz[X(i)] * dz[j];
                    }

                    sub_xz.push_back(det, vec_xz);
                }
            }

            // Recall that (x + dx)'*(z + dz) = x'*z + dx'*z + (x+dx)'*dz (the third part updated here)
            // with threshold, we ensure xz = x'*z but not z = Hx.

            // Recalculate z
            add_assign(new_z, multiply(ipentry::c(vec_xz), h));
        }

        sub_xz.ip_pt_delta += pt_delta;
        sub_xz.ip_diagonal_evaluations += diagonal_evaluations;
        sub_xz.update_xz(delta_xz); // Do high-precision operations only once
        return new_z;
        }
    }

    template<typename H>
    void update_coordinate(wff_type &det_picked, H &h,
                           wff_type &sub_xz,
                           NumericalType z_threshold = 0.0)
    {
        // Initialize sub_xz
        sub_xz.clear();
        const NumericalType e0 = vecmath::index(this->get_variational_energy(), 0);

        // For (det, dx) in det_picked, update x[det] += dx
        for (const auto &keyval : det_picked)
        {
            key_type det = keyval.first;
            data_type dx = slice_first_half(keyval.second);

            update_x(det, dx, h);
        }

        // For (det, dx) in det_picked, update z[i] += h(i, det) * dx
        // and recalculate z
#ifndef CDFCI_SOLVER_SERIAL
        size_t p = omp_get_max_threads();
        size_t inner_tasks = h.nelec / 2; // so that double excitation inner task has similar workload to single excitation task

        std::shared_ptr<H> h_shared(&h, [](H*) {
            // Empty destructor to prevent std::shared_ptr from calling delete
        });

#pragma omp parallel
{
        #pragma omp for schedule(dynamic)
#endif
        for (auto &keyval : det_picked)
        {
            key_type det = keyval.first;
            data_type dx = slice_first_half(keyval.second);
            data_type new_z = {};

#ifndef CDFCI_SOLVER_SERIAL
            #pragma omp taskgroup
            {
                generate_parallel_task(det, dx, h_shared, sub_xz, new_z, z_threshold);

                for (size_t tid = 0; tid < inner_tasks; ++tid)
                    generate_parallel_task(det, dx, h_shared, sub_xz, new_z, z_threshold,
                tid, inner_tasks);
            }
#else
            auto column = h.get_column(det);
            new_z = update_z(column, dx, h, e0, sub_xz, z_threshold);
#endif
            assign_second_half(keyval.second, new_z);
        }

#ifndef CDFCI_SOLVER_SERIAL
}
#endif

        // For (det, new_z) in det_picked, reinsert new_z and update xz
        // We always ensure xz = x'*z, z could be compressed
        for (const auto &keyval : det_picked)
        {
            key_type det = keyval.first;
            data_type new_z = slice_second_half(keyval.second);

            reinsert_new_z(det, new_z);
        }

        // Reduce sub_xz to *this
        this->update_xz(sub_xz.get_xz());
        this->update_size_z(sub_xz.size_z());
        ip_pt_energy_ += sub_xz.ip_pt_delta;
        ip_diagonal_evaluations_ += sub_xz.ip_diagonal_evaluations;

        return;
    }

    template<typename H>
    void generate_parallel_task(key_type &det, data_type dx,             // read-only, defined first-private
                                std::shared_ptr<H> h_shared,                 // shared pointer
                                wff_type &sub_xz, data_type &new_z,      // shared variable that is updated in critical area
                                NumericalType z_threshold,
                                int tid = 0, int inner_tasks = 0)
    {
#ifndef CDFCI_SOLVER_SERIAL
        #pragma omp task firstprivate(det, dx) shared(sub_xz, new_z)
        {
            typename H::Column column_parallel;
            if (inner_tasks == 0)
                column_parallel = h_shared->get_column_serial_part(det);
            else
                column_parallel = h_shared->get_column_parallel_part(det, tid, inner_tasks);

            wff_type sub_xz_parallel;
            data_type new_z_parallel;
            new_z_parallel = update_z(column_parallel, dx, *h_shared,
                                      vecmath::index(this->get_variational_energy(), 0),
                                      sub_xz_parallel, z_threshold);

            // Reduce sub_xz_parallel to sub_xz
            #pragma omp critical
            {
                sub_xz.update_xz(sub_xz_parallel.get_xz());
                sub_xz.update_size_z(sub_xz_parallel.size_z());
                sub_xz.append(sub_xz_parallel);
                add_assign(new_z, new_z_parallel);
            }
        }
#endif
    }

    // Check overflow
    // The purpose is to prevent the container to expand automatically, since
    // I do not find a way to disable automatic expansion.
    void check_overflow() const {
        if (this->size_z() > max_load_factor() * this->capacity())
        {
            std::cout << this->size_z() << ' ' << max_load_factor() << ' '
                      << this->capacity() << std::endl;
            throw std::overflow_error(
                "The wave function is almost full. Please increase max_memory "
                "(max_wavefunction_size) or z_threshold.");
        }
    }

    void reinsert_new_z(key_type &det, data_type new_z)
    {
        mapped_type vec_xz = {};
        mapped_type val_new = {};
        add_assign_second_half(val_new, new_z);
        auto update_functor = [&, new_z](mapped_type &val) {
            vec_xz = val;
            assign_second_half(val, new_z);
            if (ip_internal_cache_enabled_ && is_internal(det, val))
                sync_ip_internal_b(det, ipentry::b(val));
            return false;
        };
        upsert(data_, det, update_functor, val_new);

        // Update xz
        // x*(z+dz) = xz + x*dz
        if constexpr (NSTATES == 1) {
            this->update_xz(ipentry::c(vec_xz) * (new_z - vec_xz[1]));
        } else {
            data_matrix_type delta_xz = {};
            for (size_t i = 0; i < NSTATES; i++)
                for (size_t j = 0; j < NSTATES; j++)
                    delta_xz[IDX(i, j)] = vec_xz[X(i)] * (new_z[j] - vec_xz[Z(j)]);
            this->update_xz(delta_xz);
        }
    }

    template <typename Func>
    void loop(Func&& f) {
        data_.loop(std::forward<Func>(f));
    }

    template <typename Func>
    void loop(Func&& f) const {
        data_.loop(std::forward<Func>(f));
    }
    /* I/O */
    int load_wavefunction(const std::string& path)
    {
        std::ifstream file(path, std::ios::in | std::ios::binary);
        if (!file) {
            std::cerr << "Error opening file for reading: " << path << "\n";
            return -1;
        }
        try {
            const int status = data_.load_from_file(file);
            if (status != 0) {
                return status;
            }
            configure_ip_correction(false);
        } catch (const std::exception& e) {
            std::cerr << "Exception: " << e.what() << std::endl;
            return -1; // Return error code for exception
        }
        size_t size_x, size_z;
        load(file, size_x);
        load(file, size_z);
        std::cout << "Load wavefunction with size_x = " << size_x
                    << ", size_z = " << size_z << std::endl;
        size_ = {size_x, size_z};
        data_matrix_type xx_check = {};
        data_matrix_type xz_check = {};
        load(file, xx_check);
        load(file, xz_check);
        load(file, scale_factor_);
        this->print_xx_xz(xx_check, xz_check);

        norm_square_[0] = this->to_quad_precision(xx_check);
        dot_product_ = this->to_quad_precision(xz_check);
        file.close();
        return 0; // Success
    }


    int dump_wavefunction(const std::string& path)
    {
        std::ofstream file(path, std::ios::out | std::ios::binary);
        if (!file) {
            std::cerr << "Error opening file for writing: " << path << "\n";
            return -1;
        }
        try {
            if constexpr (COMPACT_IP) {
                // The first slot holds d on external rows while IP is active.
                // A checkpoint stores only the variational c there; d is
                // reconstructed from b and H_ii when IP initializes on load.
                const uint64_t n = static_cast<uint64_t>(data_.size());
                dump(file, n);
                data_.loop([&](const auto &entry) {
                    auto value = entry.second;
                    if (!is_internal(entry.first, value)) value[0] = 0;
                    dump(file, entry.first);
                    dump(file, value);
                });
                if (!file) return -2;
            } else {
                const int status = data_.dump_to_file(file);
                if (status != 0) return status;
            }
        } catch (const std::exception& e) {
            std::cerr << "Exception: " << e.what() << std::endl;
            return -1; // Return error code for exception
        }
        size_t size_x = this->size_x();
        size_t size_z = this->size_z();
        std::cout << "Dump wavefunction with size_x = " << size_x
                    << ", size_z = " << size_z << std::endl;
        dump(file, size_x);
        dump(file, size_z);

        data_matrix_type xx = this->get_xx_double();
        data_matrix_type xz = this->get_xz_double();
        double scale = (double) scale_factor_;
        dump(file, xx);
        dump(file, xz);
        dump(file, scale);
        file.close();
        return 0; // Success
    }

};

#endif
