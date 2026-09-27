#include "test.h"
#include "../include/solver.h"
#include <Eigen/Dense>

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

TEST_CASE("pivot elimination handles HF and a singular diagonal away from the largest coefficient") {
    Eigen::Matrix3d h;
    h << -3, -.1, .2, -.1, -1, .1, .2, .1, 2;
    Eigen::Vector3d x(1, 0, 0);
    auto pe = compute_energy_correction(DenseDiagonal{h}, SparseIterate(h, x));
    REQUIRE(pe.valid);
    CHECK(pe.internal_correction == 0);
    CHECK(pe.external_correction == doctest::Approx(-.01 / 2 - .04 / 5));
    h.setZero();
    h.diagonal() << 0, 2, 4;
    x << .3, std::sqrt(.91), 0;
    h(0, 1) = h(1, 0) = -2 * x(1) / (2 * x(0));
    pe = compute_energy_correction(DenseDiagonal{h}, SparseIterate(h, x));
    REQUIRE(pe.valid);
    CHECK(pe.diagonal_evaluations > 2);
    CHECK(pe.internal_correction == doctest::Approx(direct_correction(h, x)).epsilon(1e-11));
}

TEST_CASE("unsafe denominators and zero norm are reported without fabricated energies") {
    Eigen::Matrix3d h = Eigen::Matrix3d::Zero();
    Eigen::Vector3d x(1, 0, 0);
    h.diagonal() << 0, 0, 2;
    h(0, 1) = h(1, 0) = .1;
    auto pe = compute_energy_correction(DenseDiagonal{h}, SparseIterate(h, x));
    CHECK_FALSE(pe.valid);
    CHECK(pe.internal_valid);
    CHECK_FALSE(pe.external_valid);
    CHECK(pe.status == "unsafe_external_diagonal_gap");
    CHECK(std::isnan(pe.corrected_energy));
    CHECK(compute_energy_correction(DenseDiagonal{h}, SparseIterate(h, x), false).valid);
    pe = compute_energy_correction(DenseDiagonal{h}, SparseIterate(h, Eigen::Vector3d::Zero()));
    CHECK(pe.status == "invalid_norm");
    h.setZero();
    h.diagonal() << -2, -1, 3;
    x << .6, .8, 0;
    h(0, 1) = h(1, 0) = (2 * .36 + .64) / (.96);
    pe = compute_energy_correction(DenseDiagonal{h}, SparseIterate(h, x));
    CHECK_FALSE(pe.valid);
    CHECK_FALSE(pe.internal_valid);
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
        opts["energy_correction"] = {{"enabled", true}, {"interval", 70}, {"store_history", true}};
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
        CHECK(result.energy_correction_history.size() == 15);
        CHECK(result.energy_correction_evaluations == 15);
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
            CHECK(result.energy_correction_history[i].iteration == (i + 1) * 70);
        CHECK(result.energy_correction_history.front().valid);
    }
    opts["num_coordinates"] = 1;
    opts["stopping_dx_threshold"] = 100.0;
    opts["energy_correction"] = {{"enabled", true}, {"interval", 70}, {"store_history", false}};
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
