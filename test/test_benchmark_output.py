"""Integration checks for benchmark precision, stopping, and interrupted runs."""

import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time


app, fcidump = map(lambda value: str(Path(value).resolve()), sys.argv[1:])
environment = dict(os.environ, CDFCI_BENCHMARK_VERBOSE="0")


def command(output, mode="olsen-only", reference="0", steps="11", threshold=None):
    args = [app, fcidump, steps, "5", reference, str(output), mode,
            "0", "65536", "4", "1", "20", "0"]
    if threshold is not None:
        args.append(str(threshold))
    return args


def events(output):
    return [json.loads(line) for line in Path(str(output) + ".jsonl").read_text().splitlines()]


def check_errors(row, reference):
    energies = {
        "raw": row["variational_energy"],
        "external_only": row["variational_energy"] + row["external_correction"],
        "corrected": row["corrected_energy"],
    }
    for name, energy in energies.items():
        assert row[name + "_signed_error"] == energy - reference, row
        assert row[name + "_error"] == abs(row[name + "_signed_error"]), row


with tempfile.TemporaryDirectory(prefix="cdfci-benchmark-test-") as directory:
    root = Path(directory)
    for reference in ('0', '-100'):
        output = root / ('raw-only-' + reference + '.json')
        subprocess.run(command(output, 'raw-only', reference), env=environment,
                       stdout=subprocess.DEVNULL, check=True, timeout=15)
        data = json.loads(output.read_text())
        assert not data['options']['energy_correction']['enabled']
        assert data['trajectory_run']['correction_evaluations'] == 0
        assert [r['iteration'] for r in data['trajectory']] == [5, 10, 11]
        streamed = [r for r in events(output) if r['record_type'] == 'report']
        for a, b in zip(data['trajectory'], streamed):
            assert a['variational_energy'] == b['variational_energy']
            assert a['wall_seconds'] == b['wall_seconds']
            assert a['stored_determinants'] == b['stored_determinants']
            assert a['raw_signed_error'] == a['variational_energy'] - float(reference)
    for mode in ("olsen-only", "ip-only", "trajectory-only", "timed"):
        for reference in ("0", "-100"):
            output = root / (mode + reference + ".json")
            subprocess.run(command(output, mode, reference), env=environment,
                           stdout=subprocess.DEVNULL, check=True, timeout=15)
            data = json.loads(output.read_text())
            assert data["options"]["stopping_dx_threshold"] == 0
            assert [row["iteration"] for row in data["trajectory"]] == [5, 10, 11]
            for row in data["trajectory"]:
                check_errors(row, float(reference))
                for name in ("raw", "external_only", "corrected"):
                    assert (row[name + "_signed_error"] < 0) == (reference == "0")
            journal = events(output)
            run_ids = [event["run_id"] for event in journal if event["record_type"] == "run_start"]
            assert len(run_ids) == len(set(run_ids))
            reports = [event for event in journal if event["record_type"] == "report"
                       and event["phase"] == ("ip" if mode == "ip-only" else "ip+olsen")]
            assert [row["iteration"] for row in reports] == [5, 10, 11]
            for saved, streamed in zip(data["trajectory"], reports):
                check_errors(streamed, float(reference))
                for field in ("variational_energy", "external_correction", "internal_correction",
                              "corrected_energy"):
                    assert saved[field] == streamed[field], field
                time_key = "wall_seconds" if mode in ("ip-only", "olsen-only") else "corrected_wall_seconds"
                assert saved[time_key] == streamed["wall_seconds"]
            if mode == "ip-only":
                assert all(row["internal_correction"] == 0 for row in reports)

    early = root / "early.json"
    subprocess.run(command(early, threshold=100), env=environment,
                   stdout=subprocess.DEVNULL, check=True, timeout=15)
    assert [row["iteration"] for row in json.loads(early.read_text())["trajectory"]] == [1]
    assert [row["iteration"] for row in events(early) if row["record_type"] == "report"] == [1]

    fixed = root / "fixed-length.json"
    subprocess.run(command(fixed, steps="1001"), env=environment,
                   stdout=subprocess.DEVNULL, check=True, timeout=15)
    assert json.loads(fixed.read_text())["trajectory"][-1]["iteration"] == 1001

    for threshold in ("-1", "nan", "inf"):
        result = subprocess.run(command(root / "invalid.json", threshold=threshold),
                                env=environment, capture_output=True, text=True, timeout=15)
        assert result.returncode != 0
        assert "finite and nonnegative" in result.stderr

    result = subprocess.run(command(root / "missing-directory" / "output.json"),
                            env=environment, capture_output=True, text=True, timeout=15)
    assert result.returncode != 0
    assert "Cannot write benchmark journal" in result.stderr

    automatic = root / "auto.json"
    subprocess.run(command(automatic, "trajectory-only", "auto"), env=environment,
                   stdout=subprocess.DEVNULL, check=True, timeout=15)
    journal = events(automatic)
    assert journal[0]["reference_energy"] is None
    computed_reference = next(row["reference_energy"] for row in journal
                              if row["record_type"] == "reference")
    assert math.isfinite(computed_reference)
    assert all("raw_signed_error" not in row for row in journal
               if row.get("phase") == "reference" and row["record_type"] == "report")
    for row in journal:
        if row["record_type"] == "report" and row.get("phase") == "ip+olsen":
            check_errors(row, computed_reference)

    auto_early = root / "auto-early.json"
    subprocess.run(command(auto_early, "trajectory-only", "auto", threshold=100),
                   env=environment, stdout=subprocess.DEVNULL, check=True, timeout=15)
    reference_run = json.loads(auto_early.read_text())["reference_run"]
    assert reference_run["iterations"] == 1
    assert reference_run["requested_iterations"] == 20

    # Read full-precision points while the solve is still running, then kill it.
    # The final JSON must not exist, but flushed report lines remain recoverable.
    interrupted = root / "interrupted.json"
    process = subprocess.Popen(command(interrupted, steps="1000000000"), env=environment,
                               stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    try:
        deadline = time.monotonic() + 10
        while time.monotonic() < deadline:
            path = Path(str(interrupted) + ".jsonl")
            if path.exists():
                lines = path.read_text().splitlines()
                # Ignore only an in-flight last line, never an earlier record.
                complete = [json.loads(line) for line in lines[:-1]]
                if any(row["record_type"] == "report" for row in complete):
                    assert process.poll() is None
                    break
            assert process.poll() is None
            time.sleep(0.01)
        else:
            raise AssertionError("No report was persisted during the running solve")
    finally:
        process.kill()
        process.communicate(timeout=5)
    assert not interrupted.exists()
    reports = [row for row in complete if row["record_type"] == "report"]
    for row in reports:
        check_errors(row, 0)

print("Benchmark output checks passed")
