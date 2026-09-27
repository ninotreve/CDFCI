#!/usr/bin/env python3
"""Reproduce paper experiments 1--4 with deterministic dense models."""

from __future__ import annotations

import argparse
import csv
import json
import math
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write an empty data set: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def energy_postprocessor(h: np.ndarray, vector: np.ndarray) -> dict:
    """Solve the projected diagonal correction equation by a dense KKT solve."""
    v = np.asarray(vector, dtype=float)
    v = v / np.linalg.norm(v)
    hv = h @ v
    energy = float(v @ hv)
    residual = energy * v - hv
    n = len(v)
    kkt = np.zeros((n + 1, n + 1))
    kkt[:n, :n] = np.diag(np.diag(h) - energy)
    kkt[:n, n] = v
    kkt[n, :n] = v
    rhs = np.r_[residual, 0.0]
    try:
        solution = np.linalg.solve(kkt, rhs)
        status = "ok"
    except np.linalg.LinAlgError:
        solution = np.linalg.lstsq(kkt, rhs, rcond=None)[0]
        status = "least_squares"
    correction_vector = solution[:n]
    correction = -float(residual @ correction_vector)
    corrected = v + correction_vector
    corrected_rayleigh = float(corrected @ (h @ corrected) / (corrected @ corrected))
    return {
        "energy": energy,
        "correction": correction,
        "corrected_energy": energy + correction,
        "corrected_rayleigh": corrected_rayleigh,
        "residual_norm": float(np.linalg.norm(residual)),
        "orthogonality": float(abs(v @ correction_vector)),
        "status": status,
    }


def synthetic_orders(out: Path) -> dict:
    diagonal = np.diag([-3.0, -1.0, 1.0])
    coupling = np.array([[0.0, 0.4, -0.3], [0.4, 0.0, 0.2], [-0.3, 0.2, 0.0]])
    initial = np.array([1.0, 0.0, 0.0])
    eps_values = np.geomspace(0.01, 0.16, 13)
    rows = []
    for eps in eps_values:
        h = diagonal + eps * coupling
        exact = float(np.linalg.eigvalsh(h)[0])
        pe = energy_postprocessor(h, initial)
        rows.append({
            "epsilon": eps,
            "exact_energy": exact,
            "variational_energy": pe["energy"],
            "postprocessed_energy": pe["corrected_energy"],
            "corrected_vector_rayleigh": pe["corrected_rayleigh"],
            "variational_error": abs(pe["energy"] - exact),
            "postprocessed_error": abs(pe["corrected_energy"] - exact),
            "corrected_vector_error": abs(pe["corrected_rayleigh"] - exact),
        })
    fit_slice = slice(0, 8)
    slopes = {}
    for name in ("variational_error", "postprocessed_error", "corrected_vector_error"):
        slopes[name] = float(np.polyfit(
            np.log([row["epsilon"] for row in rows][fit_slice]),
            np.log([row[name] for row in rows][fit_slice]), 1)[0])
    write_csv(out / "01_synthetic_orders.csv", rows)
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    labels = {
        "variational_error": rf"$|E_0-E_*|$ (slope {slopes['variational_error']:.2f})",
        "postprocessed_error": rf"$|E_{{PE}}-E_*|$ (slope {slopes['postprocessed_error']:.2f})",
        "corrected_vector_error": rf"$|R(v+t)-E_*|$ (slope {slopes['corrected_vector_error']:.2f})",
    }
    for key, label in labels.items():
        ax.loglog(eps_values, [row[key] for row in rows], "o-", label=label)
    ax.set(xlabel=r"coupling $\epsilon$", ylabel="absolute energy error")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "01_synthetic_orders.pdf")
    plt.close(fig)
    return {"fitted_slopes": slopes, "fit_points": 8}


def fixed_weight_bits(n: int, weight: int) -> list[int]:
    return [sum(1 << i for i in occupied) for occupied in combinations(range(n), weight)]


def apply_hop(state: int, destination: int, source: int) -> tuple[int, int] | None:
    if not (state >> source) & 1 or (state >> destination) & 1:
        return None
    sign = -1 if (state & ((1 << source) - 1)).bit_count() % 2 else 1
    state ^= 1 << source
    sign *= -1 if (state & ((1 << destination) - 1)).bit_count() % 2 else 1
    state |= 1 << destination
    return state, sign


def hubbard_hamiltonian(sites: int, n_up: int, n_down: int, hopping: float,
                        interaction: float, potentials: np.ndarray) -> tuple[np.ndarray, list[int]]:
    states = [up | (down << sites)
              for up in fixed_weight_bits(sites, n_up)
              for down in fixed_weight_bits(sites, n_down)]
    index = {state: i for i, state in enumerate(states)}
    h = np.zeros((len(states), len(states)))
    for column, state in enumerate(states):
        up = state & ((1 << sites) - 1)
        down = state >> sites
        double_occupancy = (up & down).bit_count()
        occupation_energy = sum(
            potentials[site] * (((up >> site) & 1) + ((down >> site) & 1))
            for site in range(sites)
        )
        h[column, column] = interaction * double_occupancy + occupation_energy
        for spin_offset in (0, sites):
            for left in range(sites - 1):
                right = left + 1
                for destination, source in ((left, right), (right, left)):
                    hopped = apply_hop(state, destination + spin_offset, source + spin_offset)
                    if hopped is not None:
                        new_state, sign = hopped
                        h[index[new_state], column] += -hopping * sign
    if not np.allclose(h, h.T):
        raise RuntimeError("constructed Hubbard Hamiltonian is not symmetric")
    return h, states


def davidson(h: np.ndarray, initial: np.ndarray, max_iterations: int = 80,
             tolerance: float = 1e-11, max_subspace: int = 32) -> tuple[list[dict], np.ndarray]:
    diagonal = np.diag(h).copy()
    v = initial.astype(float)
    v /= np.linalg.norm(v)
    basis = [v]
    images = [h @ v]
    history = []
    rng = np.random.default_rng(90210)
    current = v
    for iteration in range(1, max_iterations + 1):
        projected = np.array([[basis[i] @ images[j] for j in range(len(basis))]
                              for i in range(len(basis))])
        values, vectors = np.linalg.eigh(projected)
        coeff = vectors[:, 0]
        current = sum(c * b for c, b in zip(coeff, basis))
        h_current = sum(c * image for c, image in zip(coeff, images))
        theta = float(values[0])
        residual = h_current - theta * current
        pe = energy_postprocessor(h, current)
        history.append({
            "iteration": iteration,
            "subspace_dimension": len(basis),
            "variational_energy": theta,
            "postprocessed_energy": pe["corrected_energy"],
            "residual_norm": float(np.linalg.norm(residual)),
            "postprocessor_status": pe["status"],
        })
        if np.linalg.norm(residual) < tolerance:
            break
        denominator = diagonal - theta
        safe = np.where(abs(denominator) < 1e-10,
                        np.where(denominator < 0, -1e-10, 1e-10), denominator)
        correction = -residual / safe
        for _ in range(2):
            for old in basis:
                correction -= old * (old @ correction)
        norm = np.linalg.norm(correction)
        if norm < 1e-12:
            correction = rng.normal(size=len(v))
            for old in basis:
                correction -= old * (old @ correction)
            norm = np.linalg.norm(correction)
        correction /= norm
        if len(basis) >= max_subspace:
            current /= np.linalg.norm(current)
            basis = [current]
            images = [h @ current]
        else:
            basis.append(correction)
            images.append(h @ correction)
    return history, current


def add_exact_diagnostics(history: list[dict], exact: float) -> None:
    for row in history:
        raw_error = abs(row["variational_energy"] - exact)
        corrected_error = abs(row["postprocessed_energy"] - exact)
        row["exact_energy"] = exact
        row["variational_error"] = raw_error
        row["postprocessed_error"] = corrected_error
        row["error_ratio"] = corrected_error / raw_error if raw_error > 0 else math.nan


def plot_trajectory(rows: list[dict], path: Path, title: str) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    iterations = [row["iteration"] for row in rows]
    ratios = [row["error_ratio"] for row in rows]
    residuals = [row["residual_norm"] for row in rows]
    axes[0].semilogy(iterations, ratios, "o-")
    axes[0].axhline(1.0, color="black", linestyle="--", linewidth=1)
    axes[0].set(xlabel="Davidson iteration", ylabel="corrected/raw error ratio")
    axes[1].loglog(residuals, ratios, "o-")
    axes[1].axhline(1.0, color="black", linestyle="--", linewidth=1)
    axes[1].set(xlabel="residual norm", ylabel="corrected/raw error ratio")
    for ax in axes:
        ax.grid(True, which="both", alpha=0.25)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def davidson_trajectories(out: Path) -> dict:
    sites = 6
    potentials = np.array([-0.23, -0.11, -0.04, 0.05, 0.13, 0.29])
    h, states = hubbard_hamiltonian(sites, 3, 3, 1.0, 4.0, potentials)
    exact = float(np.linalg.eigvalsh(h)[0])

    good = np.zeros(len(states))
    good[int(np.argmin(np.diag(h)))] = 1.0
    good_history, _ = davidson(h, good)
    add_exact_diagnostics(good_history, exact)
    write_csv(out / "02_davidson_trajectory.csv", good_history)
    plot_trajectory(good_history, out / "02_davidson_trajectory.pdf",
                    "Small Hubbard-CI Davidson trajectory")

    rng = np.random.default_rng(0)
    poor = rng.normal(size=len(states))
    poor /= np.linalg.norm(poor)
    poor_history, _ = davidson(h, poor, max_iterations=120)
    add_exact_diagnostics(poor_history, exact)
    write_csv(out / "03_poor_initial_vector.csv", poor_history)
    plot_trajectory(poor_history, out / "03_poor_initial_vector.pdf",
                    "Davidson trajectory from a deliberately poor vector")
    first_reliable = next((row["iteration"] for row in poor_history
                           if row["error_ratio"] < 1.0), None)
    return {
        "hamiltonian": {
            "model": "open six-site spinful Hubbard chain",
            "dimension": len(states),
            "n_up": 3,
            "n_down": 3,
            "t": 1.0,
            "U": 4.0,
            "site_potentials": potentials.tolist(),
        },
        "exact_energy": exact,
        "good_iterations": len(good_history),
        "poor_iterations": len(poor_history),
        "poor_first_error_ratio_below_one": first_reliable,
    }


def failure_regime(out: Path) -> dict:
    potentials = np.array([-0.23, -0.11, -0.04, 0.05, 0.13, 0.29])
    interactions = np.geomspace(0.25, 16.0, 13)
    rows = []
    for interaction in interactions:
        h, states = hubbard_hamiltonian(6, 3, 3, 1.0, float(interaction), potentials)
        exact = float(np.linalg.eigvalsh(h)[0])
        initial = np.zeros(len(states))
        initial[int(np.argmin(np.diag(h)))] = 1.0
        # The first Ritz vector deliberately probes the non-local regime where
        # the diagonal approximation need not yet be reliable.
        trajectory, vector = davidson(h, initial, max_iterations=1, tolerance=0.0)
        pe = energy_postprocessor(h, vector)
        raw_error = abs(pe["energy"] - exact)
        corrected_error = abs(pe["corrected_energy"] - exact)
        rows.append({
            "U_over_t": interaction,
            "iteration": trajectory[-1]["iteration"],
            "exact_energy": exact,
            "variational_energy": pe["energy"],
            "postprocessed_energy": pe["corrected_energy"],
            "variational_error": raw_error,
            "postprocessed_error": corrected_error,
            "error_ratio": corrected_error / raw_error,
            "overshoot": pe["corrected_energy"] < exact,
            "residual_norm": pe["residual_norm"],
            "status": pe["status"],
        })
    write_csv(out / "04_hubbard_failure_sweep.csv", rows)
    fig, ax = plt.subplots(figsize=(6.6, 4.5))
    ax.loglog(interactions, [row["error_ratio"] for row in rows], "o-")
    ax.axhline(1.0, color="black", linestyle="--", linewidth=1)
    overshoot = [row for row in rows if row["overshoot"]]
    if overshoot:
        ax.scatter([row["U_over_t"] for row in overshoot],
                   [row["error_ratio"] for row in overshoot], marker="x", s=70,
                   label="postprocessed energy below exact")
        ax.legend()
    ax.set(xlabel=r"$U/t$", ylabel="corrected/raw error ratio")
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out / "04_hubbard_failure_sweep.pdf")
    plt.close(fig)
    return {"sweep_points": len(rows), "overshoot_points": sum(row["overshoot"] for row in rows)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    summary = {
        "experiment_1": synthetic_orders(args.output),
        "experiments_2_and_3": davidson_trajectories(args.output),
        "experiment_4": failure_regime(args.output),
        "numpy_version": np.__version__,
    }
    (args.output / "small_experiments_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
