#include "test.h"
#include "../include/solver.h"
#include <Eigen/Dense>
#include <filesystem>

namespace {
struct DenseDiagonal {
    Eigen::MatrixXd matrix;
    double get_diagonal(int &i) const { return matrix(i, i); }
    // Deliberately no get_column / matrix-vector application API.
};

struct SparseIterate {
    static constexpr int NSTATES_val = 1;
    using key_type = int;
    using key_equal = std::equal_to<int>;
    std::vector<std::pair<int, std::array<double, 2>>> entries;
    QUAD_PRECISION norm = 0, product = 0;
    SparseIterate(const Eigen::MatrixXd &h, const Eigen::VectorXd &x) {
        Eigen::VectorXd z = h * x;
        for (int i = 0; i < x.size(); ++i) {
            if (x(i) != 0 || z(i) != 0) entries.push_back({i, {x(i), z(i)}});
            norm += static_cast<QUAD_PRECISION>(x(i)) * x(i);
            product += static_cast<QUAD_PRECISION>(x(i)) * z(i);
        }
    }
    QUAD_PRECISION get_xx() const { return norm; }
    QUAD_PRECISION get_xz() const { return product; }
    template <typename F> void loop(F &&f) const { for (const auto &entry : entries) f(entry); }
};

double direct_correction(const Eigen::MatrixXd &h, const Eigen::VectorXd &x) {
    Eigen::VectorXd v = x.normalized();
    double e = v.dot(h * v);
    Eigen::VectorXd r = e * v - h * v;
    Eigen::MatrixXd kkt = Eigen::MatrixXd::Zero(x.size() + 1, x.size() + 1);
    kkt.topLeftCorner(x.size(), x.size()).diagonal() = h.diagonal().array() - e;
    kkt.topRightCorner(x.size(), 1) = v;
    kkt.bottomLeftCorner(1, x.size()) = v.transpose();
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(x.size() + 1);
    rhs.head(x.size()) = r;
    Eigen::VectorXd t = kkt.fullPivLu().solve(rhs).head(x.size());
    CHECK(std::abs(v.dot(t)) < 1e-12);
    return -r.dot(t);
}

#ifdef CDFCI_SOLVER_SERIAL
using PendingContainer = ContainerRobinhood<int, std::array<double, 3>>;
#else
using PendingContainer = ContainerCuckoo<int, std::array<double, 3>>;
#endif
struct PendingTestWf : WaveFunction<PendingContainer, 1, true> {
    struct Diagonal {
        double get_diagonal(int &) const { return -2.0; }
    };
    void initialize(size_t rows) {
        this->reserve_capacity(65536);
        this->configure_ip_correction(true);
        this->ip_internal_cache_enabled_ = true;
        this->prepare_ip_pending_buffers();
        for (size_t i = 0; i < rows; ++i) {
            std::array<double, 3> val{{1.0, 0.0, 0.0}};
            ipentry::set_hii(val, -2.0, true);
            upsert(this->data_, static_cast<int>(i), [](auto &) { return false; }, val);
            this->sync_ip_internal_row(static_cast<int>(i), val, Diagonal{}, true, -2.0);
        }
    }
    void add_b(int key) {
        update_fn(this->data_, key, [&](auto &val) {
            val[1] += 1.0;
            // Coordinate selection must preserve a previously queued marker.
            ipentry::set_owner(val, true);
            this->queue_ip_internal_b(key, val);
            return false;
        });
    }
    const auto &unrefreshed_rows() const { return this->ip_internal_rows_; }
    auto main_value(int key) const { return this->data_.find_val(key); }
};
}

TEST_CASE("queued Olsen metadata preserves diagonals and internal ownership") {
    std::array<double, 3> val{{0.0, 7.0, 0.0}};
    ipentry::set_hii(val, -2.314159265358979, true);
    const double diagonal = ipentry::hii(val);
    CHECK(ipentry::is_c(val));
    CHECK_FALSE(ipentry::olsen_pending(val));
    ipentry::mark_olsen_pending(val);
    CHECK(ipentry::has_hii(val));
    CHECK(ipentry::is_c(val));
    CHECK(ipentry::olsen_pending(val));
    CHECK(ipentry::hii(val) == diagonal);
    CHECK(ipentry::d(val) == 0.0);
    ipentry::set_owner(val, true);
    CHECK(ipentry::olsen_pending(val));
    ipentry::clear_olsen_pending(val);
    CHECK_FALSE(ipentry::olsen_pending(val));
    CHECK(ipentry::hii(val) == diagonal);
    ipentry::set_owner(val, false);
    val[0] = 0.125;
    ipentry::mark_olsen_pending(val);
    CHECK_FALSE(ipentry::is_c(val));
    CHECK_FALSE(ipentry::olsen_pending(val));
    CHECK(ipentry::d(val) == 0.125);
    CHECK(ipentry::hii(val) == diagonal);
}

TEST_CASE("pending internal b refresh deduplicates repeated parallel updates") {
#if defined(_OPENMP) && !defined(CDFCI_SOLVER_SERIAL)
    omp_set_num_threads(4);
#endif
    constexpr size_t rows = 8192;
    PendingTestWf wf;
    wf.initialize(rows);
#if defined(_OPENMP) && !defined(CDFCI_SOLVER_SERIAL)
#pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < 8 * rows; ++i) wf.add_b(static_cast<int>(i % rows));
    CHECK(wf.pending_ip_internal_rows() == rows);
    CHECK(wf.ip_internal_b_refreshes() == 0);
    for (const auto &row : wf.unrefreshed_rows()) CHECK(row.b == 0.0);
    const auto &constant_wf = wf;
    CHECK_THROWS_AS(constant_wf.ip_internal_rows(), std::logic_error);
    const auto &updated = wf.ip_internal_rows();
    CHECK(wf.pending_ip_internal_rows() == 0);
    CHECK(wf.ip_internal_b_refreshes() == rows);
    CHECK(wf.ip_peak_pending_internal_rows() == rows);
    for (size_t i = 0; i < rows; ++i) {
        CHECK(updated[i].b == 8.0);
        const auto val = wf.main_value(static_cast<int>(i));
        REQUIRE(val.has_value());
        CHECK((*val)[0] == 1.0);
        CHECK((*val)[1] == 8.0);
        CHECK(ipentry::hii(*val) == -2.0);
        CHECK_FALSE(ipentry::olsen_pending(*val));
    }
    for (size_t i = 0; i < 1024; ++i) wf.add_b(static_cast<int>(i % 32));
    CHECK(wf.pending_ip_internal_rows() == 32);
    wf.refresh_ip_internal_rows();
    CHECK(wf.ip_internal_b_refreshes() == rows + 32);
    CHECK(wf.pending_ip_internal_rows() == 0);
    for (size_t i = 0; i < rows; ++i)
        CHECK(constant_wf.ip_internal_rows()[i].b == (i < 32 ? 40.0 : 8.0));
}

TEST_CASE("streaming correction matches the projected solve and space decomposition") {
    Eigen::Matrix4d h;
    h << -4, -.2, .1, -.15, -.2, -2, .1, .03, .1, .1, -1, -.07, -.15, .03, -.07, 1;
    Eigen::Vector4d x(1, .08, -.03, 0);
    SparseIterate wf(h, x);
    auto before = wf.entries;
    auto pe = compute_energy_correction(DenseDiagonal{h}, wf);
    REQUIRE(pe.valid);
    CHECK(pe.internal_correction + pe.external_correction == doctest::Approx(direct_correction(h, x)).epsilon(1e-11));
    auto internal = compute_energy_correction(DenseDiagonal{h}, wf, false);
    CHECK(internal.valid);
    CHECK(internal.external_correction == 0);
    CHECK(internal.internal_correction == doctest::Approx(direct_correction(h.topLeftCorner(3, 3), x.head(3))).epsilon(1e-11));
    CHECK(pe.internal_correction == internal.internal_correction);
    CHECK(pe.diagonal_evaluations == wf.entries.size());
    CHECK(wf.entries == before);
    for (double scale : {3.7, -2.3, 1e-10, 1e10}) {
        auto res = compute_energy_correction(DenseDiagonal{h}, SparseIterate(h, scale * x));
        CHECK(res.corrected_energy == doctest::Approx(pe.corrected_energy).epsilon(1e-13));
    }
    Eigen::Matrix4d shifted = h - 10 * Eigen::Matrix4d::Identity();
    auto res = compute_energy_correction(DenseDiagonal{shifted}, SparseIterate(shifted, x));
    CHECK(res.corrected_energy + 10 == doctest::Approx(pe.corrected_energy).epsilon(1e-13));
}

TEST_CASE("pivot elimination handles the single-determinant HF case") {
    Eigen::Matrix3d h;
    h << -3, -.1, .2, -.1, -1, .1, .2, .1, 2;
    Eigen::Vector3d x(1, 0, 0);
    auto pe = compute_energy_correction(DenseDiagonal{h}, SparseIterate(h, x));
    REQUIRE(pe.valid);
    CHECK(pe.internal_correction == 0);
    CHECK(pe.external_correction == doctest::Approx(-.01 / 2 - .04 / 5));
}

TEST_CASE("near diagonal energy errors have orders two three and four") {
    Eigen::Matrix3d d = Eigen::Matrix3d::Zero(), b;
    d.diagonal() << -3, -1, 1;
    b << 0, .4, -.3, .4, 0, .2, -.3, .2, 0;
    Eigen::Vector3d x(1, 0, 0);
    std::array<double, 3> previous{};
    for (double eps : {.08, .04}) {
        Eigen::Matrix3d h = d + eps * b;
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(h);
        auto pe = compute_energy_correction(DenseDiagonal{h}, SparseIterate(h, x));
        REQUIRE(pe.valid);
        Eigen::Vector3d t(0, -eps * b(1, 0) / 2, -eps * b(2, 0) / 4);
        Eigen::Vector3d corrected = x + t;
        double ev = corrected.dot(h * corrected) / corrected.squaredNorm();
        std::array<double, 3> errors{{std::abs(-3 - eig.eigenvalues()(0)),
                                     std::abs(pe.corrected_energy - eig.eigenvalues()(0)),
                                     std::abs(ev - eig.eigenvalues()(0))}};
        if (previous[0] > 0) {
            for (int i = 0; i < 3; ++i) {
                double order = std::log2(previous[i] / errors[i]);
                INFO("order " << i + 2 << " measured " << order);
                CHECK(order > i + 1.8);
                CHECK(order < i + 2.2);
            }
        }
        previous = errors;
    }
}

TEST_CASE("incremental PT and Olsen preserve the CDFCI iterate") {
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    using Det = Determinant<1>;
#ifdef CDFCI_SOLVER_SERIAL
    using Container = ContainerRobinhood<Det, std::array<double, 3>, DeterminantHash<1>, DeterminantEqual<1>>;
#else
    using Container = ContainerCuckoo<Det, std::array<double, 3>, DeterminantHashRobinhood<1>, DeterminantEqual<1>>;
#endif
    using Wf = WaveFunction<Container, 1, true>;
    using Ham = Hamiltonian<1>;
    Det::constuct_masks();
    Option hopt = {{"type", "molecule"}, {"molecule", {{"fcidump_path", "data/h2o_sto3g_psi4.FCIDUMP"}, {"verbose", 0}}}};
    auto ham = Ham::init(hopt);
    Option opts = {{"num_iterations", 100}, {"report_interval", 20}, {"z_threshold", 0.0},
                   {"stopping_dx_threshold", 0.0}, {"max_wavefunction_size", 65536}, {"verbose", 0}};
    for (int coordinates : {1, 4}) {
        opts["num_coordinates"] = coordinates;
        opts["energy_correction"] = {{"enabled", false}};
        CDFCISolver<Ham, Wf> baseline(opts);
        Wf plain;
        REQUIRE(baseline.solve(*ham, plain) == 0);

        opts["energy_correction"] = {{"enabled", true}, {"olsen_enabled", true}};
        CDFCISolver<Ham, Wf> corrected(opts);
        Wf with_correction;
        REQUIRE(corrected.solve(*ham, with_correction) == 0);
        CHECK(with_correction.ip_correction_result(100).diagonal_evaluations ==
              with_correction.size());
        size_t inspected_diagonals = 0;
        with_correction.loop([&](const auto &entry) {
            if (inspected_diagonals++ % 997 != 0) return;
            auto det = entry.first;
            CHECK(ipentry::has_hii(entry.second));
            CHECK(ipentry::hii(entry.second) ==
                  doctest::Approx(ham->get_diagonal(det)).epsilon(1e-12));
        });

        QUAD_PRECISION direct_ip_sum = 0;
        with_correction.loop([&](const auto &entry) {
            if (with_correction.stored_c(entry.first, entry.second) == 0)
                direct_ip_sum += static_cast<QUAD_PRECISION>(entry.second[0]) * entry.second[1];
        });
        CHECK(with_correction.ip_correction_result(100).external_correction ==
              doctest::Approx(static_cast<double>(direct_ip_sum / with_correction.get_xx()))
                  .epsilon(1e-11).scale(1.0));

        const auto cached_olsen = compute_cached_internal_energy_correction(with_correction);
        const auto scanned_olsen = compute_energy_correction(*ham, with_correction, false, true);
        REQUIRE(cached_olsen.valid);
        CHECK(with_correction.ip_internal_rows().size() == with_correction.size_x());
        CHECK(cached_olsen.internal_correction ==
              doctest::Approx(scanned_olsen.internal_correction).epsilon(1e-11).scale(1.0));
        CHECK(cached_olsen.internal_residual_norm ==
              doctest::Approx(scanned_olsen.internal_residual_norm).epsilon(1e-11).scale(1.0));

        const auto checkpoint = std::filesystem::temp_directory_path() /
            ("cdfci_compact_ip_" + std::to_string(coordinates) + ".dat");
        REQUIRE(with_correction.dump_wavefunction(checkpoint.string()) == 0);
        Wf reloaded;
        REQUIRE(reloaded.load_wavefunction(checkpoint.string()) == 0);
        reloaded.initialize_ip_correction(*ham, true);
        CHECK(reloaded.size_x() == with_correction.size_x());
        CHECK(reloaded.size_z() == with_correction.size_z());
        CHECK(compute_cached_internal_energy_correction(reloaded).internal_correction ==
              doctest::Approx(cached_olsen.internal_correction).epsilon(1e-11).scale(1.0));
        std::filesystem::remove(checkpoint);

        const auto a = baseline.get_result();
        const auto b = corrected.get_result();
        REQUIRE(a.energy_history.size() == b.energy_history.size());
        for (size_t i = 0; i < a.energy_history.size(); ++i)
            CHECK(b.energy_history[i] == doctest::Approx(a.energy_history[i]).epsilon(1e-13));
        CHECK(b.energy == doctest::Approx(a.energy).epsilon(1e-13));
        CHECK(with_correction.size() == plain.size());
        std::vector<std::pair<Det, std::array<double, 3>>> plain_snapshot;
        plain.loop([&](const auto &entry) {
            plain_snapshot.push_back({entry.first, entry.second});
        });
        for (const auto &entry : plain_snapshot) {
            CHECK(with_correction.get_x(entry.first) ==
                  doctest::Approx(plain.get_x(entry.first)).epsilon(1e-12).scale(1.0));
            CHECK(with_correction.get_z(entry.first) ==
                  doctest::Approx(plain.get_z(entry.first)).epsilon(1e-12).scale(1.0));
        }

        opts["energy_correction"] = {{"enabled", true}, {"olsen_enabled", false}};
        CDFCISolver<Ham, Wf> pt_only(opts);
        Wf without_olsen;
        REQUIRE(pt_only.solve(*ham, without_olsen) == 0);
        CHECK(without_olsen.ip_correction_result(100).diagonal_evaluations ==
              without_olsen.size());
        CHECK(without_olsen.ip_internal_rows().empty());
        CHECK(pt_only.get_result().energy == doctest::Approx(a.energy).epsilon(1e-13));
        QUAD_PRECISION direct_pt_only_sum = 0;
        without_olsen.loop([&](const auto &entry) {
            if (without_olsen.stored_c(entry.first, entry.second) == 0)
                direct_pt_only_sum += static_cast<QUAD_PRECISION>(entry.second[0]) * entry.second[1];
        });
        CHECK(without_olsen.ip_correction_result(100).external_correction ==
              doctest::Approx(static_cast<double>(direct_pt_only_sum / without_olsen.get_xx()))
                  .epsilon(1e-11).scale(1.0));
    }
#ifdef _OPENMP
    omp_set_num_threads(2);
    opts["num_coordinates"] = 4;
    opts["energy_correction"] = {{"enabled", true}, {"olsen_enabled", true}};
    CDFCISolver<Ham, Wf> parallel(opts);
    Wf parallel_wf;
    REQUIRE(parallel.solve(*ham, parallel_wf) == 0);
    const auto cached_parallel = compute_cached_internal_energy_correction(parallel_wf);
    const auto scanned_parallel = compute_energy_correction(*ham, parallel_wf, false, true);
    CHECK(cached_parallel.internal_correction ==
          doctest::Approx(scanned_parallel.internal_correction).epsilon(1e-11).scale(1.0));
#endif
}

TEST_CASE("H2O correction preserves the CDFCI trajectory and uses exact internal z under compression") {
#ifdef _OPENMP
    // Compare trajectories under deterministic scheduling. A separate two-
    // thread snapshot test below checks the parallel update and read-only path.
    omp_set_num_threads(1);
#endif
    using Det = Determinant<1>;
#ifdef CDFCI_SOLVER_SERIAL
    using Container = ContainerRobinhood<Det, std::array<double, 2>, DeterminantHash<1>, DeterminantEqual<1>>;
#else
    using Container = ContainerCuckoo<Det, std::array<double, 2>, DeterminantHashRobinhood<1>, DeterminantEqual<1>>;
#endif
    using Wf = WaveFunction<Container>;
    using Ham = Hamiltonian<1>;
    Det::constuct_masks();
    Option hopt = {{"type", "molecule"}, {"molecule", {{"fcidump_path", "data/h2o_sto3g_psi4.FCIDUMP"}, {"verbose", 0}}}};
    auto ham = Ham::init(hopt);
    Option opts = {{"num_iterations", 1001}, {"report_interval", 100}, {"z_threshold", 1e-3},
                   {"stopping_dx_threshold", 0}, {"max_wavefunction_size", 131072}, {"verbose", 0}};
    for (int coordinates : {1, 4}) {
        opts["num_coordinates"] = coordinates;
        opts["energy_correction"] = {{"enabled", false}};
        CDFCISolver<Ham, Wf> baseline(opts);
        Wf plain;
        REQUIRE(baseline.solve(*ham, plain) == 0);
        opts["energy_correction"] = {{"enabled", true}, {"store_history", true}};
        CDFCISolver<Ham, Wf> corrected(opts);
        Wf wf;
        REQUIRE(corrected.solve(*ham, wf) == 0);
        auto result = corrected.get_result();
        CHECK(result.energy == doctest::Approx(baseline.get_result().energy).epsilon(1e-13));
#ifdef CDFCI_SOLVER_SERIAL
        CHECK(result.energy == baseline.get_result().energy);
        CHECK(wf.size() == plain.size());
#endif
        std::vector<std::pair<Det, std::array<double, 2>>> snapshot;
        wf.loop([&](const auto &entry) { snapshot.push_back({entry.first, entry.second}); });
        // A direct snapshot evaluation must never mutate any stored entry,
        // including compressed external z in the parallel container.
        compute_energy_correction(*ham, wf);
        for (const auto &entry : snapshot) {
            CHECK(entry.second[0] == wf.get_x(entry.first));
            CHECK(entry.second[1] == wf.get_z(entry.first));
            CHECK(entry.second[0] == doctest::Approx(plain.get_x(entry.first)).epsilon(1e-11).scale(1.0));
#ifdef CDFCI_SOLVER_SERIAL
            CHECK(entry.second[1] == doctest::Approx(plain.get_z(entry.first)).epsilon(1e-11).scale(1.0));
#endif
            // Separate OpenMP runs can retain different small external z
            // contributions depending on the order of first insertions.
            if (entry.second[0] == 0) continue;
            CHECK(entry.second[1] == doctest::Approx(plain.get_z(entry.first)).epsilon(1e-11).scale(1.0));
            auto det = entry.first;
            double exact_z = 0;
            for (const auto &h : ham->get_column(det)) exact_z += h.second * wf.get_x(h.first);
            CHECK(entry.second[1] == doctest::Approx(exact_z).epsilon(1e-11).scale(1.0));
        }
        CHECK(result.energy_correction_history.size() == 11);
        CHECK(result.energy_correction_evaluations == 11);
        CHECK(result.iterations == 1001);
        CHECK(result.report_interval == 100);
        CHECK(result.energy_history.size() == 11);
        CHECK(result.hamiltonian_columns_history.size() == 11);
        CHECK(result.hamiltonian_columns > 0);
        CHECK(result.hamiltonian_columns_history.back() == result.hamiltonian_columns);
        CHECK(result.x_size_history.size() == 11);
        CHECK(result.z_size_history.size() == 11);
        CHECK(result.time_history.size() == 11);
        CHECK(result.energy_correction.iteration == 1001);
        CHECK(result.energy_correction.valid);
        CHECK(result.energy_correction.compressed_z);
        for (size_t i = 0; i + 1 < result.energy_correction_history.size(); ++i)
            CHECK(result.energy_correction_history[i].iteration == (i + 1) * 100);
        CHECK(result.energy_correction_history.front().valid);
    }
    opts["num_coordinates"] = 1;
    opts["stopping_dx_threshold"] = 100.0;
    opts["energy_correction"] = {{"enabled", true}, {"store_history", false}};
    CDFCISolver<Ham, Wf> early(opts);
    Wf early_wf;
    REQUIRE(early.solve(*ham, early_wf) == 0);
    auto stopped = early.get_result();
    CHECK(stopped.iterations == 1);
    CHECK(stopped.energy_correction.iteration == 1);
    CHECK(stopped.energy_correction_evaluations == 1);
    CHECK(stopped.energy_correction_history.empty());
#ifdef _OPENMP
    // Parallel compression and tie breaking can alter trajectories between
    // independent runs. Test the actual requirement on a single two-thread
    // trajectory: evaluation preserves all entries and internal z is exact.
    omp_set_num_threads(2);
    opts["num_coordinates"] = 4;
    opts["stopping_dx_threshold"] = 0.0;
    CDFCISolver<Ham, Wf> parallel(opts);
    Wf parallel_wf;
    REQUIRE(parallel.solve(*ham, parallel_wf) == 0);
    std::vector<std::pair<Det, std::array<double, 2>>> snapshot;
    parallel_wf.loop([&](const auto &entry) { snapshot.push_back({entry.first, entry.second}); });
    REQUIRE(compute_energy_correction(*ham, parallel_wf).valid);
    for (const auto &entry : snapshot) {
        CHECK(entry.second[0] == parallel_wf.get_x(entry.first));
        CHECK(entry.second[1] == parallel_wf.get_z(entry.first));
        if (entry.second[0] == 0) continue;
        auto det = entry.first;
        double exact_z = 0;
        for (const auto &h : ham->get_column(det)) exact_z += h.second * parallel_wf.get_x(h.first);
        CHECK(entry.second[1] == doctest::Approx(exact_z).epsilon(1e-11).scale(1.0));
    }
#endif
}
