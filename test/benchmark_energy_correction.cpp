// Reproducible benchmark for the streaming energy correction. FCIDUMP loading
// is excluded from timings; solver allocation is included. Determinants up to
// two 64-bit words are supported so larger orbital spaces are not accidentally
// restricted to <= 32 spatial orbitals.
#include "../include/solver.h"
#include <fstream>
#include <sstream>

template <int N>
int run_benchmark(int argc, char **argv)
{
    using Det = Determinant<N>;
#ifdef CDFCI_SOLVER_SERIAL
    using Container = ContainerRobinhood<Det, std::array<double, 2>,
                                         DeterminantHash<N>, DeterminantEqual<N>>;
#else
    using Container = ContainerCuckoo<Det, std::array<double, 2>,
                                      DeterminantHashRobinhood<N>, DeterminantEqual<N>>;
#endif
    using Wf = WaveFunction<Container>;
    using Ham = Hamiltonian<N>;

    const size_t iterations = std::stoull(argv[2]);
    const size_t interval = std::stoull(argv[3]);
    if (iterations == 0 || interval == 0) {
        std::cerr << "iterations and interval must both be positive\n";
        return 1;
    }
    const std::string reference_argument = argv[4];
    const std::string mode = argc >= 7 ? argv[6] : "timed";
    const bool trajectory_only = mode == "trajectory-only";
    if (mode != "timed" && !trajectory_only) {
        std::cerr << "The optional mode must be timed or trajectory-only\n";
        return 1;
    }
    const double z_threshold = argc >= 8 ? std::stod(argv[7]) : 1e-3;
    const size_t max_wavefunction_size = argc >= 9 ? std::stoull(argv[8]) : 524288;
    const size_t num_coordinates = argc >= 10 ? std::stoull(argv[9]) : 1;
    const size_t timing_repeats = argc >= 11 ? std::stoull(argv[10]) : 1;
    const size_t reference_iterations = argc >= 12 ? std::stoull(argv[11]) : 5 * iterations;
    const double reference_z_threshold = argc >= 13 ? std::stod(argv[12]) : z_threshold;
    if (max_wavefunction_size == 0 || num_coordinates == 0 ||
        (!trajectory_only && timing_repeats == 0)) {
        std::cerr << "max_wavefunction_size, num_coordinates, and timing_repeats "
                     "(in timed mode) must be positive\n";
        return 1;
    }
    const bool automatic_reference = reference_argument == "auto";
    if (automatic_reference && reference_iterations <= iterations) {
        std::cerr << "automatic reference_iterations must exceed trajectory iterations\n";
        return 1;
    }

    Option hopt = {{"type", "molecule"},
                   {"molecule", {{"fcidump_path", argv[1]},
                                  {"threshold", 0.0}, {"verbose", 0}}}};
    auto ham = Ham::init(hopt);
    Det::constuct_masks();
    Option opt = {{"num_iterations", iterations}, {"report_interval", interval},
                  {"num_coordinates", num_coordinates}, {"z_threshold", z_threshold},
                  {"stopping_dx_threshold", 1e-8},
                  {"max_wavefunction_size", max_wavefunction_size}, {"verbose", 0}};

    std::ostringstream muted;
    auto original = std::cout.rdbuf(muted.rdbuf());
    auto run = [&](bool enabled, bool history, size_t steps, double run_z_threshold) {
        opt["num_iterations"] = steps;
        opt["report_interval"] = std::min(interval, steps);
        opt["z_threshold"] = run_z_threshold;
        opt["energy_correction"] = {{"enabled", enabled}, {"interval", interval},
                                    {"store_history", history}};
        CDFCISolver<Ham, Wf> solver(opt);
        Wf wf;
        const auto start = std::chrono::steady_clock::now();
        const int status = solver.solve(*ham, wf);
        const double seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        if (status != 0)
            throw std::runtime_error("CDFCI benchmark solve failed with status " +
                                     std::to_string(status));
        auto result = solver.get_result();
        Option summary = {{"seconds", seconds}, {"energy", result.energy},
                          {"stored_determinants", wf.size()},
                          {"coefficient_bytes_per_entry", sizeof(typename Wf::mapped_type)},
                          {"correction_seconds", result.energy_correction_seconds},
                          {"correction_evaluations", result.energy_correction_evaluations}};
        return std::make_pair(result, summary);
    };

    double reference = 0.0;
    Option reference_run = nullptr;
    if (automatic_reference) {
        auto computed = run(false, false, reference_iterations, reference_z_threshold);
        reference = computed.first.energy;
        reference_run = computed.second;
        reference_run["iterations"] = reference_iterations;
        reference_run["z_threshold"] = reference_z_threshold;
        if (computed.first.energy_history.size() >= 2) {
            const size_t last = computed.first.energy_history.size() - 1;
            reference_run["previous_report_energy"] = computed.first.energy_history[last - 1];
            reference_run["tail_energy_change"] =
                std::abs(computed.first.energy_history[last] -
                         computed.first.energy_history[last - 1]);
        }
    } else {
        reference = std::stod(reference_argument);
    }

    if (!trajectory_only)
        run(false, false, std::min<size_t>(iterations, 200), z_threshold); // warm caches

    // These runs provide both trajectories and the first wall-time pair. This
    // avoids two redundant production-scale solves from the old benchmark.
    auto baseline_trajectory = run(false, false, iterations, z_threshold);
    auto trajectory = run(true, true, iterations, z_threshold);
    Option runs = Option::array();
    if (!trajectory_only) {
        runs.push_back({{"baseline", baseline_trajectory.second},
                        {"enabled", trajectory.second}});
        for (size_t repeat = 1; repeat < timing_repeats; ++repeat) {
            std::pair<Result, Option> plain, pe;
            if (repeat % 2 == 0) {
                plain = run(false, false, iterations, z_threshold);
                pe = run(true, false, iterations, z_threshold);
            } else {
                pe = run(true, false, iterations, z_threshold);
                plain = run(false, false, iterations, z_threshold);
            }
            runs.push_back({{"baseline", plain.second}, {"enabled", pe.second}});
        }
    }
    std::cout.rdbuf(original);

    Option records = Option::array();
    size_t record_index = 0;
    double previous_corrected_wall_seconds = 0.0;
    double cumulative_correction_seconds = 0.0;
    for (const auto &pe : trajectory.first.energy_correction_history) {
        if (record_index >= baseline_trajectory.first.energy_history.size() ||
            record_index >= baseline_trajectory.first.time_history.size() ||
            record_index >= trajectory.first.x_size_history.size() ||
            record_index >= trajectory.first.z_size_history.size() ||
            record_index >= trajectory.first.hamiltonian_columns_history.size() ||
            record_index >= trajectory.first.time_history.size())
            throw std::runtime_error("CDFCI benchmark histories are not aligned");
        const double baseline_energy = baseline_trajectory.first.energy_history[record_index];
        const double corrected_wall_seconds = trajectory.first.time_history[record_index];
        const double interval_total_seconds =
            corrected_wall_seconds - previous_corrected_wall_seconds;
        cumulative_correction_seconds += pe.seconds;
        records.push_back({{"iteration", pe.iteration},
                           {"variational_energy", pe.variational_energy},
                           {"baseline_variational_energy", baseline_energy},
                           {"internal_correction", pe.internal_correction},
                           {"external_correction", pe.external_correction},
                           {"corrected_energy", pe.corrected_energy}, {"status", pe.status},
                           {"raw_error", std::abs(pe.variational_energy - reference)},
                           {"external_only_error",
                            std::abs(pe.variational_energy + pe.external_correction - reference)},
                           {"corrected_error", std::abs(pe.corrected_energy - reference)},
                           {"internal_residual_norm", pe.internal_residual_norm},
                           {"stored_determinants",
                            trajectory.first.x_size_history[record_index]},
                           {"stored_residual_entries",
                            trajectory.first.z_size_history[record_index]},
                           {"stored_wavefunction_entries",
                            trajectory.first.z_size_history[record_index]},
                           {"hamiltonian_columns",
                            trajectory.first.hamiltonian_columns_history[record_index]},
                           {"raw_wall_seconds",
                            baseline_trajectory.first.time_history[record_index]},
                           {"corrected_wall_seconds", corrected_wall_seconds},
                           {"correction_seconds", pe.seconds},
                           {"interval_total_seconds", interval_total_seconds},
                           {"interval_correction_overhead",
                            interval_total_seconds > 0 ? pe.seconds / interval_total_seconds : 0.0},
                           {"cumulative_correction_overhead",
                            corrected_wall_seconds > 0
                                ? cumulative_correction_seconds / corrected_wall_seconds
                                : 0.0},
                           {"diagonal_evaluations", pe.diagonal_evaluations},
                           {"seconds", pe.seconds}});
        previous_corrected_wall_seconds = corrected_wall_seconds;
        ++record_index;
    }
    Option output = {{"fcidump", argv[1]}, {"spin_orbitals", ham->norb},
                     {"electrons", ham->nelec}, {"ms2", ham->ms2},
                     {"determinant_words", N}, {"reference_energy", reference},
                     {"reference_kind", automatic_reference ? "computed" : "supplied"},
                     {"reference_run", reference_run}, {"options", opt},
                     {"baseline_trajectory_run", baseline_trajectory.second},
                     {"trajectory_run", trajectory.second}, {"timings", runs},
                     {"trajectory", records}};
    std::ofstream file(argv[5]);
    file << output.dump(2) << '\n';
    if (!file)
        return 2;
    std::cout << "Saved benchmark: " << argv[5] << '\n';
    return 0;
}

int main(int argc, char **argv)
{
    if (argc < 6 || argc > 13) {
        std::cerr
            << "Usage: benchmark_energy_correction FCIDUMP iterations interval "
               "reference_energy|auto output.json [timed|trajectory-only "
               "[z_threshold [max_wavefunction_size [num_coordinates "
               "[timing_repeats [reference_iterations [reference_z_threshold]]]]]]]]\n";
        return 1;
    }
    Option hopt = {{"type", "molecule"},
                   {"molecule", {{"fcidump_path", argv[1]}, {"verbose", 0}}}};
    const int spin_orbitals = Hamiltonian<1>::read_norb(hopt);
    const int determinant_words = 1 + (spin_orbitals - 1) / Determinant<1>::det_size;
    switch (determinant_words) {
    case 1: return run_benchmark<1>(argc, argv);
    case 2: return run_benchmark<2>(argc, argv);
    default:
        std::cerr << "Benchmark supports at most 128 spin orbitals; input has "
                  << spin_orbitals << '\n';
        return 1;
    }
}
