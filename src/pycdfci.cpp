#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "python_api.h"

namespace py = pybind11;

PYBIND11_MODULE(_cdfci, m) {
    m.doc() = "Python bindings for CDFCI";

    py::class_<EnergyCorrectionResult>(m, "EnergyCorrectionResult")
        .def_readonly("iteration", &EnergyCorrectionResult::iteration)
        .def_readonly("variational_energy", &EnergyCorrectionResult::variational_energy)
        .def_readonly("internal_correction", &EnergyCorrectionResult::internal_correction)
        .def_readonly("external_correction", &EnergyCorrectionResult::external_correction)
        .def_readonly("corrected_energy", &EnergyCorrectionResult::corrected_energy)
        .def_readonly("internal_residual_norm", &EnergyCorrectionResult::internal_residual_norm)
        .def_readonly("stored_residual_norm", &EnergyCorrectionResult::stored_residual_norm)
        .def_readonly("internal_valid", &EnergyCorrectionResult::internal_valid)
        .def_readonly("external_valid", &EnergyCorrectionResult::external_valid)
        .def_readonly("valid", &EnergyCorrectionResult::valid)
        .def_readonly("status", &EnergyCorrectionResult::status)
        .def_readonly("compressed_z", &EnergyCorrectionResult::compressed_z)
        .def_readonly("diagonal_evaluations", &EnergyCorrectionResult::diagonal_evaluations)
        .def_readonly("seconds", &EnergyCorrectionResult::seconds);

    py::class_<Result>(m, "CDFCIResult")
        .def(py::init<>())
        .def_readonly("energy", &Result::energy)
        .def_readonly("state_energies", &Result::state_energies)
        .def_readonly("iterations", &Result::iterations)
        .def_readonly("energy_correction", &Result::energy_correction)
        .def_readonly("energy_correction_history", &Result::energy_correction_history)
        .def_readonly("energy_correction_seconds", &Result::energy_correction_seconds)
        .def_readonly("energy_correction_evaluations", &Result::energy_correction_evaluations);

    py::class_<CDFCIDriverFacade>(m, "CDFCI")
        .def(py::init<const std::string &>(), py::arg("fcidump_path"))
        .def("set_num_iterations", &CDFCIDriverFacade::set_num_iterations, py::arg("n"))
        .def("set_report_interval", &CDFCIDriverFacade::set_report_interval, py::arg("n"))
        .def("set_max_memory", &CDFCIDriverFacade::set_max_memory, py::arg("gb"))
        .def("set_max_load_factor", &CDFCIDriverFacade::set_max_load_factor, py::arg("x"))
        .def("run", &CDFCIDriverFacade::run, py::call_guard<py::gil_scoped_release>()); // GIL release

        m.def("run_cdfci_json",
            &CDFCIRunnerFacade::run_cdfci_json,
            py::arg("option_json"),
            py::call_guard<py::gil_scoped_release>());

        m.def("run_xcdfci_json",
            &CDFCIRunnerFacade::run_xcdfci_json,
            py::arg("option_json"),
            py::call_guard<py::gil_scoped_release>());

        m.def("run_cdfci_rhf_integrals_json",
            &CDFCIRunnerFacade::run_cdfci_rhf_integrals_json,
            py::arg("option_json"),
            py::arg("h1e_flat"),
            py::arg("eri_flat"),
            py::arg("norb"),
            py::arg("nelec"),
            py::arg("core_energy"),
            py::arg("ms2") = 0,
            py::call_guard<py::gil_scoped_release>());

        m.def("run_optorbfci_json",
            &CDFCIRunnerFacade::run_optorbfci_json,
            py::arg("option_json"),
            py::call_guard<py::gil_scoped_release>());

        m.def("tools_frozen_orbital_json",
            &CDFCIRunnerFacade::tools_frozen_orbital_json,
            py::arg("option_json"),
            py::call_guard<py::gil_scoped_release>());

        m.def("tools_symmetry_connection_json",
            &CDFCIRunnerFacade::tools_symmetry_connection_json,
            py::arg("option_json"),
            py::call_guard<py::gil_scoped_release>());
}
