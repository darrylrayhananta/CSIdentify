"""Forensic STR profiling with regular expressions.

The program compares a crime-scene DNA sequence against a suspect database by
counting the longest consecutive run for each Short Tandem Repeat (STR) marker.
It intentionally avoids whole-sequence alignment and uses a compact STR profile
as the comparison feature.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
from pathlib import Path
from typing import Callable, Iterable


Suspect = dict[str, object]
CASE_FIELDS = ["case_id", "label", "dna_sequence", "source"]


def compact_sequence(sequence: str) -> str:
    """Uppercase a sequence and remove whitespace characters."""
    return re.sub(r"\s+", "", sequence).upper()


def validate_dna_sequence(sequence: str) -> str:
    sequence = compact_sequence(sequence)
    invalid = sorted(set(sequence) - set("ACGTN"))
    if invalid:
        raise ValueError(f"DNA sequence contains invalid characters: {', '.join(invalid)}")
    if not sequence:
        raise ValueError("DNA sequence must not be empty")
    return sequence


def read_dna_sequence(path: str | Path) -> str:
    try:
        return validate_dna_sequence(Path(path).read_text(encoding="utf-8"))
    except ValueError as exc:
        raise ValueError(str(exc).replace("DNA sequence", "DNA file")) from exc


def next_case_id(existing_rows: list[dict[str, str]]) -> str:
    highest = 0
    for row in existing_rows:
        raw_id = (row.get("case_id") or "").strip()
        match = re.fullmatch(r"CASE-(\d+)", raw_id)
        if match:
            highest = max(highest, int(match.group(1)))
    return f"CASE-{highest + 1:03d}"


def save_case_to_csv(
    path: str | Path,
    label: str,
    dna_sequence: str,
    source: str,
) -> str:
    """Persist user-provided DNA to a case CSV and return the new case id."""
    case_path = Path(path)
    case_path.parent.mkdir(parents=True, exist_ok=True)
    existing_rows: list[dict[str, str]] = []
    if case_path.exists():
        with case_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames and reader.fieldnames != CASE_FIELDS:
                raise ValueError(f"Case CSV must use columns: {', '.join(CASE_FIELDS)}")
            existing_rows = list(reader)

    case_id = next_case_id(existing_rows)
    row = {
        "case_id": case_id,
        "label": label.strip() or "DNA TKP",
        "dna_sequence": validate_dna_sequence(dna_sequence),
        "source": source.strip() or "manual input",
    }
    with case_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CASE_FIELDS)
        writer.writeheader()
        writer.writerows(existing_rows)
        writer.writerow(row)
    return case_id


def load_case_from_csv(path: str | Path, case_id: str | None = None) -> dict[str, str]:
    """Load a case by id from case CSV. If no id is supplied, use the latest row."""
    case_path = Path(path)
    with case_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != CASE_FIELDS:
            raise ValueError(f"Case CSV must use columns: {', '.join(CASE_FIELDS)}")
        rows = list(reader)

    if not rows:
        raise ValueError("Case CSV does not contain any DNA cases")
    if case_id is None:
        row = rows[-1]
    else:
        row = next((item for item in rows if item["case_id"] == case_id), None)
        if row is None:
            raise ValueError(f"Case id not found in CSV: {case_id}")

    return {
        "case_id": row["case_id"],
        "label": row["label"],
        "dna_sequence": validate_dna_sequence(row["dna_sequence"]),
        "source": row["source"],
    }


def find_repeat_runs(sequence: str, marker: str) -> list[dict[str, object]]:
    """Return every consecutive STR run with positions for an explainable trace."""
    sequence = compact_sequence(sequence)
    marker = compact_sequence(marker)
    if not marker:
        raise ValueError("STR marker must not be empty")

    pattern_text = f"(?:{re.escape(marker)})+"
    pattern = re.compile(pattern_text)
    runs: list[dict[str, object]] = []
    for match in pattern.finditer(sequence):
        matched_sequence = match.group(0)
        runs.append(
            {
                "start": match.start(),
                "end": match.end(),
                "repeat_count": len(matched_sequence) // len(marker),
                "matched_sequence": matched_sequence,
                "pattern": pattern_text,
            }
        )
    return runs


def longest_consecutive_repeats(sequence: str, marker: str) -> int:
    """Return the longest consecutive repeat count of marker inside sequence."""
    runs = find_repeat_runs(sequence, marker)
    return max((int(run["repeat_count"]) for run in runs), default=0)


def profile_sequence(sequence: str, markers: Iterable[str]) -> dict[str, int]:
    """Build an STR profile for a DNA sequence."""
    return {marker: longest_consecutive_repeats(sequence, marker) for marker in markers}


def load_suspects(path: str | Path) -> tuple[list[Suspect], list[str]]:
    """Load suspect profiles from a CSV file.

    The first column is treated as the suspect name. Remaining columns are STR
    marker names whose values must be integer repeat counts.
    """
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or len(reader.fieldnames) < 2:
            raise ValueError("CSV must contain a name column and at least one STR marker column")

        name_field = reader.fieldnames[0]
        markers = reader.fieldnames[1:]
        suspects: list[Suspect] = []
        for row_number, row in enumerate(reader, start=2):
            name = (row.get(name_field) or "").strip()
            if not name:
                raise ValueError(f"Missing suspect name at row {row_number}")

            profile: dict[str, int] = {}
            for marker in markers:
                raw_value = (row.get(marker) or "").strip()
                try:
                    profile[marker] = int(raw_value)
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid repeat count for {marker!r} at row {row_number}: {raw_value!r}"
                    ) from exc

            suspects.append({"name": name, "profile": profile})

    if not suspects:
        raise ValueError("CSV database does not contain suspect rows")
    return suspects, markers


def find_best_matches(sample_profile: dict[str, int], suspects: list[Suspect]) -> list[dict[str, object]]:
    """Rank suspects by exact STR marker matches and total profile distance."""
    markers = list(sample_profile.keys())
    ranked: list[dict[str, object]] = []

    for suspect in suspects:
        suspect_profile = suspect["profile"]
        if not isinstance(suspect_profile, dict):
            raise TypeError("suspect profile must be a dictionary")

        exact_count = sum(sample_profile[m] == int(suspect_profile[m]) for m in markers)
        distance = sum(abs(sample_profile[m] - int(suspect_profile[m])) for m in markers)
        differences = {
            marker: int(suspect_profile[marker]) - sample_profile[marker]
            for marker in markers
            if sample_profile[marker] != int(suspect_profile[marker])
        }
        ranked.append(
            {
                "name": suspect["name"],
                "profile": suspect_profile,
                "matching_markers": exact_count,
                "marker_count": len(markers),
                "score_percent": round(100 * exact_count / len(markers), 2),
                "distance": distance,
                "differences": differences,
                "is_exact_match": exact_count == len(markers),
            }
        )

    return sorted(
        ranked,
        key=lambda item: (-int(item["matching_markers"]), int(item["distance"]), str(item["name"])),
    )


def analyze_sequence(
    dna_sequence: str,
    suspects_path: str | Path,
    output_dir: str | Path,
    case_info: dict[str, str] | None = None,
) -> dict[str, object]:
    suspects, markers = load_suspects(suspects_path)
    sequence = validate_dna_sequence(dna_sequence)
    sample_profile = profile_sequence(sequence, markers)
    ranked = find_best_matches(sample_profile, suspects)
    trace = build_processing_trace(sequence, markers, sample_profile, ranked)

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    html_path = output_path / "investigation_dashboard.html"
    json_path = output_path / "analysis_result.json"

    html_path.write_text(
        build_html_dashboard(),
        encoding="utf-8",
    )
    json_path.write_text(
        json.dumps(
            {
                "case": case_info,
                "sample_profile": sample_profile,
                "ranked_suspects": ranked,
                "markers": markers,
                "processing_trace": trace,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    return {
        "sample_profile": sample_profile,
        "ranked_suspects": ranked,
        "markers": markers,
        "case": case_info,
        "processing_trace": trace,
        "html_path": html_path,
        "json_path": json_path,
    }


def analyze_case(dna_path: str | Path, suspects_path: str | Path, output_dir: str | Path) -> dict[str, object]:
    sequence = read_dna_sequence(dna_path)
    case_info = {
        "case_id": "-",
        "label": "DNA TKP",
        "dna_sequence": sequence,
        "source": f"file: {dna_path}",
    }
    return analyze_sequence(sequence, suspects_path, output_dir, case_info=case_info)


def analyze_case_csv(
    case_csv_path: str | Path,
    case_id: str | None,
    suspects_path: str | Path,
    output_dir: str | Path,
) -> dict[str, object]:
    case_info = load_case_from_csv(case_csv_path, case_id)
    return analyze_sequence(case_info["dna_sequence"], suspects_path, output_dir, case_info=case_info)


def build_processing_trace(
    sequence: str,
    markers: Iterable[str],
    sample_profile: dict[str, int],
    ranked: list[dict[str, object]],
) -> dict[str, object]:
    """Build structured narration data for interactive demos."""
    marker_traces = []
    markers = list(markers)
    for marker in markers:
        runs = find_repeat_runs(sequence, marker)
        marker_traces.append(
            {
                "marker": marker,
                "pattern": f"(?:{re.escape(marker)})+",
                "runs": runs,
                "longest_repeat": sample_profile[marker],
            }
        )

    suspect_traces = []
    for item in ranked:
        profile = item["profile"]
        if not isinstance(profile, dict):
            raise TypeError("suspect profile must be a dictionary")
        matched_markers = [marker for marker in markers if int(profile[marker]) == sample_profile[marker]]
        difference_notes = [
            f"{marker}: tersangka {int(profile[marker])} vs TKP {sample_profile[marker]}"
            for marker in markers
            if int(profile[marker]) != sample_profile[marker]
        ]
        suspect_traces.append(
            {
                "name": item["name"],
                "matched_markers": matched_markers,
                "difference_notes": difference_notes,
                "score_percent": item["score_percent"],
                "distance": item["distance"],
                "is_exact_match": item["is_exact_match"],
            }
        )

    return {"markers": marker_traces, "suspects": suspect_traces}


def make_bar(value: int, maximum: int, width: int = 18) -> str:
    if maximum <= 0:
        maximum = 1
    filled = round((value / maximum) * width)
    return "[" + "#" * filled + "-" * (width - filled) + "]"


def format_processing_trace(trace: dict[str, object]) -> list[str]:
    lines = ["", "=== Jejak Processing STR ==="]
    marker_traces = trace["markers"]
    if not isinstance(marker_traces, list):
        raise TypeError("marker traces must be a list")
    max_repeat = max(
        [1] + [int(marker_trace["longest_repeat"]) for marker_trace in marker_traces]
    )

    for marker_trace in marker_traces:
        runs = marker_trace["runs"]
        if not isinstance(runs, list):
            raise TypeError("runs must be a list")
        lines.append(f"\nScan marker {marker_trace['marker']} dengan RegEx {marker_trace['pattern']}")
        if not runs:
            lines.append("  Tidak ada run tandem yang ditemukan -> repeat = 0")
            continue
        for run in runs:
            lines.append(
                "  "
                f"posisi {run['start']}-{run['end']}: "
                f"{run['matched_sequence']} -> {run['repeat_count']} repeat"
            )
        longest = int(marker_trace["longest_repeat"])
        lines.append(f"  Longest run: {make_bar(longest, max_repeat)} {longest}")

    lines.append("\n=== Skoring Tersangka ===")
    suspect_traces = trace["suspects"]
    if not isinstance(suspect_traces, list):
        raise TypeError("suspect traces must be a list")
    for index, suspect_trace in enumerate(suspect_traces, start=1):
        status = "MATCH" if suspect_trace["is_exact_match"] else "near match"
        matched = ", ".join(suspect_trace["matched_markers"]) or "-"
        lines.append(
            f"{index}. {suspect_trace['name']} | skor {suspect_trace['score_percent']}% | "
            f"jarak {suspect_trace['distance']} | {status}"
        )
        lines.append(f"   Marker cocok: {matched}")
        differences = suspect_trace["difference_notes"]
        if differences:
            lines.append("   Selisih: " + "; ".join(differences))
    return lines


def print_lines(
    lines: Iterable[str],
    output_func: Callable[[str], None] = print,
    animate: bool = False,
    delay_seconds: float = 0.08,
) -> None:
    for line in lines:
        output_func(line)
        if animate:
            time.sleep(delay_seconds)


def print_case_summary(result: dict[str, object], output_func: Callable[[str], None] = print) -> None:
    best = result["ranked_suspects"][0]
    output_func("")
    output_func("=== Hasil Akhir Investigasi ===")
    output_func(f"Tersangka paling cocok: {best['name']}")
    output_func(f"Skor: {best['matching_markers']}/{best['marker_count']} marker ({best['score_percent']}%)")
    output_func(f"Exact match: {'YES' if best['is_exact_match'] else 'NO'}")
    output_func(f"Dashboard: {result['html_path']}")
    output_func(f"JSON: {result['json_path']}")


def run_interactive(
    args: argparse.Namespace,
    input_func: Callable[[str], str] = input,
    output_func: Callable[[str], None] = print,
) -> None:
    output_func("=== CSI STR Profiler: Mode Interaktif ===")
    output_func("1. Pakai DNA TKP demo")
    output_func("2. Ketik DNA TKP manual")
    output_func("3. Pakai file DNA sendiri")
    choice = (input_func("Pilih mode [1]: ").strip() or "1")

    dna_sequence: str | None = None
    label = "DNA TKP Demo"
    source = f"file: {args.dna}"
    if choice == "2":
        label = input_func("Label kasus [DNA TKP Manual]: ").strip() or "DNA TKP Manual"
        dna_sequence = validate_dna_sequence(input_func("Tempel sekuens DNA TKP: "))
        source = "input manual"
    elif choice == "3":
        dna_path = input_func(f"Path file DNA [{args.dna}]: ").strip() or args.dna
        dna_sequence = read_dna_sequence(dna_path)
        label = input_func("Label kasus [DNA TKP dari file]: ").strip() or "DNA TKP dari file"
        source = f"file: {dna_path}"
    else:
        dna_sequence = read_dna_sequence(args.dna)

    db_path = input_func(f"Path database tersangka [{args.db}]: ").strip() or args.db
    case_csv_path = input_func(f"File CSV kasus [{args.case_csv}]: ").strip() or args.case_csv
    output_dir = input_func(f"Folder output [{args.out}]: ").strip() or args.out
    detail_choice = (input_func("Tampilkan jejak scan RegEx? [Y/n]: ").strip().lower() or "y")
    show_detail = detail_choice != "n"

    case_id = save_case_to_csv(case_csv_path, label, dna_sequence, source)
    output_func(f"\nInput DNA tersimpan ke {case_csv_path} sebagai {case_id}.")
    output_func("Analisis membaca kembali DNA dari CSV kasus tersebut...")
    result = analyze_case_csv(case_csv_path, case_id, db_path, output_dir)
    if show_detail:
        print_lines(
            format_processing_trace(result["processing_trace"]),
            output_func=output_func,
            animate=not args.no_animation,
        )
    print_case_summary(result, output_func=output_func)


def build_html_dashboard() -> str:
    return r"""<!doctype html>
<html lang="id">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>CSI STR Profiler</title>
  <style>
    :root {
      --blue: #1F46B6;
      --magenta: #E63A6E;
      --peach: #F4C0AF;
      --cream: #F7E6D3;
      --navy: #17327D;
      --paper: rgba(255, 255, 255, 0.78);
      --paper-strong: rgba(255, 255, 255, 0.94);
      --line: rgba(23, 50, 125, 0.18);
      --text: #17327D;
      --muted: #5D6280;
      --danger: #A0183E;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      min-height: 100vh;
      font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      color: var(--text);
      background:
        linear-gradient(135deg, rgba(244, 192, 175, 0.92), rgba(247, 230, 211, 0.96) 54%, rgba(31, 70, 182, 0.16)),
        var(--cream);
      padding: 28px;
    }
    main { max-width: 1180px; margin: 0 auto; }
    header {
      display: grid;
      grid-template-columns: minmax(0, 1fr) auto;
      gap: 18px;
      align-items: center;
      margin-bottom: 16px;
    }
    h1 {
      margin: 0;
      color: var(--navy);
      font-size: clamp(32px, 5vw, 56px);
      line-height: 1;
      letter-spacing: 0;
    }
    .eyebrow {
      color: var(--magenta);
      font-size: 12px;
      font-weight: 900;
      letter-spacing: 0.1em;
      margin-bottom: 8px;
      text-transform: uppercase;
    }
    .subtitle {
      max-width: 720px;
      margin: 10px 0 0;
      color: var(--muted);
      font-size: 15px;
      line-height: 1.55;
    }
    .status-chip {
      min-width: 210px;
      padding: 12px 14px;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: var(--paper);
      text-align: right;
      backdrop-filter: blur(16px);
    }
    .status-chip span {
      display: block;
      color: var(--muted);
      font-size: 11px;
      font-weight: 800;
      letter-spacing: 0.08em;
      text-transform: uppercase;
    }
    .status-chip strong {
      display: block;
      margin-top: 4px;
      color: var(--blue);
      font-size: 16px;
    }
    .layout {
      display: grid;
      grid-template-columns: minmax(300px, 0.74fr) minmax(0, 1.26fr);
      gap: 16px;
      align-items: start;
    }
    .panel {
      border: 1px solid var(--line);
      border-radius: 8px;
      background: var(--paper);
      box-shadow: 0 18px 42px rgba(23, 50, 125, 0.13);
      backdrop-filter: blur(16px);
      padding: 18px;
    }
    .panel + .panel { margin-top: 16px; }
    .panel-title {
      margin: 0 0 14px;
      color: var(--navy);
      font-size: 18px;
      line-height: 1.2;
    }
    .upload-grid {
      display: grid;
      gap: 13px;
    }
    label {
      display: block;
      margin-bottom: 7px;
      color: var(--navy);
      font-size: 12px;
      font-weight: 900;
      letter-spacing: 0.08em;
      text-transform: uppercase;
    }
    input[type="file"] {
      width: 100%;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: rgba(255, 255, 255, 0.72);
      color: var(--text);
      padding: 11px 12px;
      font: 13px/1.45 Inter, ui-sans-serif, system-ui, sans-serif;
      outline: none;
      cursor: pointer;
    }
    input[type="file"]:focus {
      border-color: var(--blue);
      box-shadow: 0 0 0 3px rgba(31, 70, 182, 0.14);
    }
    .file-note {
      min-height: 18px;
      margin-top: 7px;
      color: var(--muted);
      font-size: 12px;
      line-height: 1.45;
      overflow-wrap: anywhere;
    }
    button {
      width: 100%;
      border: 1px solid var(--blue);
      border-radius: 8px;
      background: var(--blue);
      color: #ffffff;
      min-height: 42px;
      padding: 10px 14px;
      font-weight: 900;
      letter-spacing: 0.02em;
      cursor: pointer;
      transition: transform 160ms ease, box-shadow 160ms ease, background 160ms ease;
    }
    button:hover {
      background: var(--navy);
      box-shadow: 0 12px 24px rgba(31, 70, 182, 0.22);
      transform: translateY(-1px);
    }
    .error {
      display: none;
      margin-top: 12px;
      padding: 10px 12px;
      border: 1px solid rgba(230, 58, 110, 0.36);
      border-radius: 8px;
      background: rgba(230, 58, 110, 0.10);
      color: var(--danger);
      font-size: 13px;
      font-weight: 800;
      line-height: 1.4;
    }
    .empty-state {
      display: grid;
      min-height: 206px;
      place-items: center;
      padding: 22px;
      border: 1px dashed rgba(31, 70, 182, 0.34);
      border-radius: 8px;
      background: rgba(255, 255, 255, 0.45);
      color: var(--muted);
      text-align: center;
      line-height: 1.45;
    }
    .summary {
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 10px;
      margin-bottom: 16px;
    }
    .metric {
      min-height: 84px;
      padding: 12px;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: rgba(247, 230, 211, 0.70);
    }
    .metric span {
      display: block;
      color: var(--muted);
      font-size: 11px;
      font-weight: 900;
      letter-spacing: 0.07em;
      text-transform: uppercase;
    }
    .metric strong {
      display: block;
      margin-top: 7px;
      color: var(--navy);
      font-size: clamp(18px, 3vw, 28px);
      line-height: 1.05;
      overflow-wrap: anywhere;
    }
    .metric strong.small {
      font-size: 15px;
      line-height: 1.25;
    }
    .verdict {
      display: inline-flex;
      align-items: center;
      min-height: 28px;
      padding: 5px 10px;
      border-radius: 999px;
      background: rgba(31, 70, 182, 0.10);
      color: var(--blue);
      font-size: 12px;
      font-weight: 900;
    }
    .verdict.warn {
      background: rgba(230, 58, 110, 0.12);
      color: var(--magenta);
    }
    .chart-frame {
      overflow-x: auto;
      border-radius: 8px;
      background: var(--navy);
    }
    #profileComparisonChart {
      min-height: 360px;
    }
    .chart-empty {
      display: grid;
      min-height: 360px;
      place-items: center;
      color: rgba(247, 230, 211, 0.82);
      text-align: center;
      line-height: 1.45;
      padding: 22px;
    }
    .two-col {
      display: grid;
      grid-template-columns: minmax(0, 0.78fr) minmax(0, 1.22fr);
      gap: 16px;
      margin-top: 16px;
    }
    table {
      width: 100%;
      border-collapse: collapse;
      overflow: hidden;
      border-radius: 8px;
      background: rgba(255, 255, 255, 0.50);
    }
    th, td {
      padding: 10px 11px;
      border-bottom: 1px solid rgba(23, 50, 125, 0.10);
      text-align: left;
      vertical-align: top;
    }
    th {
      color: var(--navy);
      font-size: 11px;
      font-weight: 900;
      letter-spacing: 0.07em;
      text-transform: uppercase;
    }
    td {
      color: var(--text);
      font-size: 13px;
    }
    code {
      display: inline-block;
      max-width: 100%;
      padding: 3px 6px;
      border-radius: 6px;
      background: rgba(31, 70, 182, 0.10);
      color: var(--blue);
      overflow-wrap: anywhere;
    }
    .rank {
      display: inline-grid;
      width: 27px;
      height: 27px;
      place-items: center;
      border-radius: 999px;
      background: rgba(31, 70, 182, 0.12);
      color: var(--blue);
      font-weight: 900;
    }
    .pill {
      display: inline-flex;
      align-items: center;
      min-height: 26px;
      padding: 4px 9px;
      border-radius: 999px;
      background: rgba(31, 70, 182, 0.10);
      color: var(--blue);
      font-weight: 900;
    }
    .pill.warn {
      background: rgba(230, 58, 110, 0.10);
      color: var(--magenta);
    }
    footer {
      margin-top: 16px;
      color: var(--muted);
      font-size: 12px;
      line-height: 1.45;
    }
    .is-hidden {
      display: none !important;
    }
    @media (max-width: 900px) {
      body { padding: 18px; }
      header, .layout, .two-col {
        display: block;
      }
      .status-chip {
        margin-top: 14px;
        text-align: left;
      }
      .panel {
        margin-bottom: 16px;
      }
      .summary {
        grid-template-columns: repeat(2, minmax(0, 1fr));
      }
    }
    @media (max-width: 560px) {
      .summary {
        grid-template-columns: 1fr;
      }
      th, td {
        padding: 9px 8px;
      }
    }
  </style>
</head>
<body>
  <main>
    <header>
      <div>
        <div class="eyebrow">Forensic STR Analysis</div>
        <h1>CSI STR Profiler</h1>
        <p class="subtitle">Upload DNA TKP berbentuk file .txt dan database tersangka berbentuk file .csv. Sistem akan membaca marker STR dari header CSV, mencari repeat terpanjang dengan RegEx, lalu membandingkan profil TKP dengan semua tersangka.</p>
      </div>
      <div class="status-chip">
        <span>Status</span>
        <strong id="analysisStatus">Menunggu file</strong>
      </div>
    </header>

    <section class="layout">
      <aside>
        <div class="panel">
          <h2 class="panel-title">Input File</h2>
          <div class="upload-grid">
            <div>
              <label for="dnaFileInput">DNA TKP (.txt)</label>
              <input id="dnaFileInput" type="file" accept=".txt,text/plain">
              <div id="dnaFileStatus" class="file-note">Belum ada file DNA TKP.</div>
            </div>
            <div>
              <label for="suspectCsvFileInput">Database Tersangka (.csv)</label>
              <input id="suspectCsvFileInput" type="file" accept=".csv,text/csv">
              <div id="suspectFileStatus" class="file-note">Belum ada file CSV tersangka.</div>
            </div>
            <button id="analyzeButton" type="button">Analyze STR</button>
          </div>
          <div id="errorBox" class="error" role="alert"></div>
        </div>

        <div class="panel">
          <h2 class="panel-title">Hasil Utama</h2>
          <div id="resultEmptyState" class="empty-state">Upload file .txt DNA TKP dan file .csv tersangka terlebih dahulu.</div>
          <div id="resultSummary" class="is-hidden">
            <span id="verdictView" class="verdict">READY</span>
            <div class="summary" style="margin-top: 12px;">
              <div class="metric">
                <span>Best Match</span>
                <strong id="bestNameView" class="small">-</strong>
              </div>
              <div class="metric">
                <span>Skor</span>
                <strong id="scoreView">0%</strong>
              </div>
              <div class="metric">
                <span>Marker Cocok</span>
                <strong id="matchedView">0/0</strong>
              </div>
              <div class="metric">
                <span>Jarak</span>
                <strong id="distanceView">0</strong>
              </div>
            </div>
          </div>
        </div>
      </aside>

      <article class="panel">
        <h2 class="panel-title">Profil STR DNA TKP vs Tersangka</h2>
        <div class="chart-frame">
          <div id="profileComparisonChart">
            <div class="chart-empty">Chart akan muncul setelah kedua file diupload dan dianalisis.</div>
          </div>
        </div>
      </article>
    </section>

    <section id="detailsSection" class="two-col is-hidden">
      <article class="panel">
        <h2 class="panel-title">Profil DNA TKP</h2>
        <table>
          <thead><tr><th>STR</th><th>Repeat</th></tr></thead>
          <tbody id="profileRows"></tbody>
        </table>
      </article>

      <article class="panel">
        <h2 class="panel-title">Ranking Tersangka</h2>
        <table>
          <thead><tr><th>#</th><th>Nama</th><th>Cocok</th><th>Skor</th><th>Jarak</th><th>Status</th></tr></thead>
          <tbody id="suspectRows"></tbody>
        </table>
      </article>
    </section>

    <section id="traceSection" class="panel is-hidden">
      <h2 class="panel-title">Trace RegEx</h2>
      <table>
        <thead><tr><th>Marker</th><th>Pattern</th><th>Run Ditemukan</th><th>Repeat Terpanjang</th></tr></thead>
        <tbody id="traceRows"></tbody>
      </table>
    </section>

    <footer>Model edukatif: hasil STR ini untuk simulasi komputasi, bukan pengganti validasi laboratorium forensik.</footer>
  </main>

  <script>
    const $ = (id) => document.getElementById(id);
    const REQUIRED_FILE_MESSAGE = 'Upload file .txt DNA TKP dan file .csv tersangka terlebih dahulu.';
    const CHART_COLORS = ['#F7E6D3', '#1F46B6', '#E63A6E', '#F4C0AF', '#65D8C7', '#FFC857', '#17327D'];

    function escapeHtml(value) {
      return String(value).replace(/[&<>"']/g, (char) => ({
        '&': '&amp;',
        '<': '&lt;',
        '>': '&gt;',
        '"': '&quot;',
        "'": '&#39;',
      }[char]));
    }

    function escapeRegExp(value) {
      return String(value).replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
    }

    function compactSequence(sequence) {
      return String(sequence).replace(/\s+/g, '').toUpperCase();
    }

    function validateDna(sequence) {
      const clean = compactSequence(sequence);
      if (!clean) {
        throw new Error('File DNA TKP kosong.');
      }
      const invalid = [...new Set(clean.split('').filter((char) => !'ACGTN'.includes(char)))];
      if (invalid.length) {
        throw new Error(`DNA mengandung karakter invalid: ${invalid.join(', ')}.`);
      }
      return clean;
    }

    function parseCsvLine(line) {
      const cells = [];
      let current = '';
      let inQuotes = false;
      for (let index = 0; index < line.length; index += 1) {
        const char = line[index];
        const next = line[index + 1];
        if (char === '"' && inQuotes && next === '"') {
          current += '"';
          index += 1;
        } else if (char === '"') {
          inQuotes = !inQuotes;
        } else if (char === ',' && !inQuotes) {
          cells.push(current.trim());
          current = '';
        } else {
          current += char;
        }
      }
      cells.push(current.trim());
      return cells;
    }

    function parseSuspectCsv(csvText) {
      const lines = String(csvText).trim().split(/\r?\n/).filter((line) => line.trim());
      if (lines.length < 2) {
        throw new Error('CSV tersangka harus punya header dan minimal satu baris tersangka.');
      }
      const headers = parseCsvLine(lines[0]);
      const nameHeader = headers[0] ? headers[0].toLowerCase() : '';
      const markers = headers.slice(1).map((marker) => marker.trim()).filter(Boolean);
      if (!['nama', 'name', 'tersangka', 'suspect'].includes(nameHeader)) {
        throw new Error('Kolom pertama CSV harus berisi nama tersangka.');
      }
      if (!markers.length) {
        throw new Error('CSV harus punya minimal satu marker STR setelah kolom nama.');
      }
      const suspects = lines.slice(1).map((line, rowIndex) => {
        const rowNumber = rowIndex + 2;
        const cells = parseCsvLine(line);
        const name = cells[0] ? cells[0].trim() : '';
        if (!name) {
          throw new Error(`Nama tersangka kosong di baris ${rowNumber}.`);
        }
        const profile = {};
        markers.forEach((marker, markerIndex) => {
          const rawValue = cells[markerIndex + 1];
          const value = Number.parseInt(rawValue, 10);
          if (!Number.isInteger(value) || value < 0) {
            throw new Error(`Nilai marker ${marker} invalid di baris ${rowNumber}.`);
          }
          profile[marker] = value;
        });
        return { name, profile };
      });
      return { markers, suspects };
    }

    function readFileAsText(fileInput, label) {
      const file = fileInput.files && fileInput.files[0];
      if (!file) {
        return Promise.reject(new Error(`${label} belum dipilih.`));
      }
      return file.text();
    }

    function findRepeatRuns(sequence, marker) {
      const pattern = new RegExp(`(?:${escapeRegExp(marker)})+`, 'g');
      const runs = [];
      let match;
      while ((match = pattern.exec(sequence)) !== null) {
        runs.push({
          start: match.index,
          end: match.index + match[0].length,
          repeatCount: Math.floor(match[0].length / marker.length),
          matchedSequence: match[0],
          pattern: `(?:${marker})+`,
        });
      }
      return runs;
    }

    function profileSequence(sequence, markers) {
      const profile = {};
      const trace = markers.map((marker) => {
        const runs = findRepeatRuns(sequence, marker);
        const longestRepeat = Math.max(0, ...runs.map((run) => run.repeatCount));
        profile[marker] = longestRepeat;
        return { marker, pattern: `(?:${marker})+`, runs, longestRepeat };
      });
      return { profile, trace };
    }

    function rankSuspects(sampleProfile, suspects, markers) {
      return suspects.map((suspect) => {
        const matchedMarkers = markers.filter((marker) => suspect.profile[marker] === sampleProfile[marker]);
        const distance = markers.reduce((sum, marker) => (
          sum + Math.abs((suspect.profile[marker] || 0) - (sampleProfile[marker] || 0))
        ), 0);
        const scorePercent = Math.round((matchedMarkers.length / markers.length) * 10000) / 100;
        return {
          name: suspect.name,
          profile: suspect.profile,
          matchingMarkers: matchedMarkers.length,
          markerCount: markers.length,
          scorePercent,
          distance,
          isExactMatch: matchedMarkers.length === markers.length,
        };
      }).sort((left, right) => (
        right.matchingMarkers - left.matchingMarkers
        || left.distance - right.distance
        || left.name.localeCompare(right.name)
      ));
    }

    function showError(message) {
      $('errorBox').textContent = message;
      $('errorBox').style.display = 'block';
      $('analysisStatus').textContent = 'Input belum lengkap';
    }

    function clearError() {
      $('errorBox').textContent = '';
      $('errorBox').style.display = 'none';
    }

    function renderEmptyState() {
      $('analysisStatus').textContent = 'Menunggu file';
      $('resultEmptyState').classList.remove('is-hidden');
      $('resultSummary').classList.add('is-hidden');
      $('detailsSection').classList.add('is-hidden');
      $('traceSection').classList.add('is-hidden');
      $('profileRows').innerHTML = '';
      $('suspectRows').innerHTML = '';
      $('traceRows').innerHTML = '';
      $('profileComparisonChart').innerHTML = '<div class="chart-empty">Chart akan muncul setelah kedua file diupload dan dianalisis.</div>';
    }

    function renderGroupedProfileChart(state) {
      const markers = state.markers;
      const series = [
        { label: 'DNA TKP', profile: state.sampleProfile },
        ...state.ranked.map((suspect) => ({ label: suspect.name, profile: suspect.profile })),
      ];
      const maxValue = Math.max(1, ...markers.flatMap((marker) => series.map((item) => Number(item.profile[marker] || 0))));
      const yMax = Math.max(5, Math.ceil(maxValue / 2) * 2);
      const width = Math.max(820, markers.length * 150 + 170);
      const height = 390;
      const left = 58;
      const right = 24;
      const top = 56;
      const bottom = 82;
      const plotWidth = width - left - right;
      const plotHeight = height - top - bottom;
      const groupWidth = plotWidth / markers.length;
      const barGap = 3;
      const barWidth = Math.max(8, Math.min(20, (groupWidth - 28) / series.length - barGap));
      const tickCount = Math.min(yMax, 6);
      let svg = `<svg viewBox="0 0 ${width} ${height}" width="100%" height="${height}" role="img" aria-label="Profil STR DNA TKP vs tersangka">`;
      svg += '<rect width="100%" height="100%" rx="8" fill="#17327D"></rect>';
      svg += '<text x="18" y="30" fill="#F7E6D3" font-size="16" font-weight="800" font-family="Inter, Arial">Profil STR DNA TKP vs Tersangka</text>';
      svg += '<text x="18" y="49" fill="rgba(247,230,211,0.72)" font-size="11" font-family="Inter, Arial">Repeat count per marker</text>';

      for (let tick = 0; tick <= tickCount; tick += 1) {
        const value = Math.round((yMax / tickCount) * tick);
        const y = top + plotHeight - (value / yMax) * plotHeight;
        svg += `<line x1="${left}" y1="${y}" x2="${width - right}" y2="${y}" stroke="rgba(247,230,211,0.15)" stroke-width="1"></line>`;
        svg += `<text x="${left - 14}" y="${y + 4}" fill="#F7E6D3" font-size="11" text-anchor="end" font-family="Inter, Arial">${value}</text>`;
      }

      markers.forEach((marker, markerIndex) => {
        const groupLeft = left + markerIndex * groupWidth;
        const center = groupLeft + groupWidth / 2;
        svg += `<text x="${center}" y="${height - 44}" fill="#F7E6D3" font-size="12" font-weight="800" text-anchor="middle" font-family="Inter, Arial">${escapeHtml(marker)}</text>`;
        series.forEach((item, seriesIndex) => {
          const value = Number(item.profile[marker] || 0);
          const x = center - ((series.length * barWidth) + ((series.length - 1) * barGap)) / 2 + seriesIndex * (barWidth + barGap);
          const barHeight = (value / yMax) * plotHeight;
          const y = top + plotHeight - barHeight;
          const color = CHART_COLORS[seriesIndex % CHART_COLORS.length];
          svg += `<rect x="${x}" y="${y}" width="${barWidth}" height="${barHeight}" rx="3" fill="${color}"></rect>`;
          if (barHeight > 18) {
            svg += `<text x="${x + barWidth / 2}" y="${y + 13}" fill="#FFFFFF" font-size="9" font-weight="900" text-anchor="middle" font-family="Inter, Arial">${value}</text>`;
          } else {
            svg += `<text x="${x + barWidth / 2}" y="${Math.max(top + 11, y - 5)}" fill="#F7E6D3" font-size="9" font-weight="900" text-anchor="middle" font-family="Inter, Arial">${value}</text>`;
          }
        });
      });

      const legendY = height - 22;
      let legendX = left;
      series.forEach((item, index) => {
        const color = CHART_COLORS[index % CHART_COLORS.length];
        const label = item.label.length > 16 ? `${item.label.slice(0, 15)}...` : item.label;
        svg += `<rect x="${legendX}" y="${legendY - 10}" width="10" height="10" rx="2" fill="${color}"></rect>`;
        svg += `<text x="${legendX + 15}" y="${legendY}" fill="#F7E6D3" font-size="10" font-family="Inter, Arial">${escapeHtml(label)}</text>`;
        legendX += Math.max(92, label.length * 6 + 32);
      });

      svg += '</svg>';
      $('profileComparisonChart').innerHTML = svg;
    }

    function renderState(state) {
      const best = state.ranked[0];
      $('analysisStatus').textContent = 'Analisis selesai';
      $('resultEmptyState').classList.add('is-hidden');
      $('resultSummary').classList.remove('is-hidden');
      $('detailsSection').classList.remove('is-hidden');
      $('traceSection').classList.remove('is-hidden');

      $('verdictView').textContent = best.isExactMatch ? 'Exact match' : 'Near match';
      $('verdictView').classList.toggle('warn', !best.isExactMatch);
      $('bestNameView').textContent = best.name;
      $('scoreView').textContent = `${best.scorePercent}%`;
      $('matchedView').textContent = `${best.matchingMarkers}/${best.markerCount}`;
      $('distanceView').textContent = String(best.distance);
      $('profileRows').innerHTML = state.markers.map((marker) => (
        `<tr><td><strong>${escapeHtml(marker)}</strong></td><td>${state.sampleProfile[marker]}</td></tr>`
      )).join('');
      $('suspectRows').innerHTML = state.ranked.map((item, index) => (
        `<tr>
          <td><span class="rank">${index + 1}</span></td>
          <td><strong>${escapeHtml(item.name)}</strong></td>
          <td>${item.matchingMarkers}/${item.markerCount}</td>
          <td><span class="pill">${item.scorePercent}%</span></td>
          <td>${item.distance}</td>
          <td><span class="pill ${item.isExactMatch ? '' : 'warn'}">${item.isExactMatch ? 'Exact' : 'Selisih'}</span></td>
        </tr>`
      )).join('');
      $('traceRows').innerHTML = state.trace.map((item) => {
        const runs = item.runs.length
          ? item.runs.map((run) => `${run.repeatCount}x (${run.start}-${run.end})`).join(', ')
          : '-';
        return `<tr>
          <td><strong>${escapeHtml(item.marker)}</strong></td>
          <td><code>${escapeHtml(item.pattern)}</code></td>
          <td>${escapeHtml(runs)}</td>
          <td>${item.longestRepeat}</td>
        </tr>`;
      }).join('');
      renderGroupedProfileChart(state);
    }

    async function runInteractiveAnalysis() {
      clearError();
      const dnaFileInput = $('dnaFileInput');
      const suspectCsvFileInput = $('suspectCsvFileInput');
      if (!dnaFileInput.files.length || !suspectCsvFileInput.files.length) {
        renderEmptyState();
        showError(REQUIRED_FILE_MESSAGE);
        return;
      }

      try {
        $('analysisStatus').textContent = 'Memproses...';
        const dnaText = await readFileAsText(dnaFileInput, 'File DNA TKP');
        const suspectCsvText = await readFileAsText(suspectCsvFileInput, 'File CSV tersangka');
        const sequence = validateDna(dnaText);
        const parsed = parseSuspectCsv(suspectCsvText);
        const profiled = profileSequence(sequence, parsed.markers);
        const ranked = rankSuspects(profiled.profile, parsed.suspects, parsed.markers);
        renderState({
          markers: parsed.markers,
          sampleProfile: profiled.profile,
          trace: profiled.trace,
          ranked,
          sequence,
        });
      } catch (error) {
        renderEmptyState();
        showError(error.message);
      }
    }

    $('dnaFileInput').addEventListener('change', (event) => {
      const file = event.target.files[0];
      $('dnaFileStatus').textContent = file ? file.name : 'Belum ada file DNA TKP.';
      clearError();
    });

    $('suspectCsvFileInput').addEventListener('change', (event) => {
      const file = event.target.files[0];
      $('suspectFileStatus').textContent = file ? file.name : 'Belum ada file CSV tersangka.';
      clearError();
    });

    $('analyzeButton').addEventListener('click', runInteractiveAnalysis);
    renderEmptyState();
  </script>
</body>
</html>
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze forensic STR profiles from DNA and suspect CSV files.")
    parser.add_argument("--dna", default="data/crime_scene_dna.txt", help="Path to the crime-scene DNA sequence")
    parser.add_argument("--db", default="data/suspects.csv", help="Path to the suspect profile CSV")
    parser.add_argument("--case-csv", default="data/cases.csv", help="CSV case log for user-entered DNA samples")
    parser.add_argument("--use-case-csv", action="store_true", help="Process DNA from --case-csv instead of --dna")
    parser.add_argument("--case-id", default=None, help="Case id to process from --case-csv; defaults to latest row")
    parser.add_argument("--out", default="output", help="Directory for HTML, SVG, and JSON outputs")
    parser.add_argument("--interactive", action="store_true", help="Run an input-driven demo with step-by-step narration")
    parser.add_argument("--explain", action="store_true", help="Print the RegEx scan trace in non-interactive mode")
    parser.add_argument("--no-animation", action="store_true", help="Disable small delays between interactive narration lines")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.interactive:
        run_interactive(args)
        return

    if args.use_case_csv or args.case_id:
        result = analyze_case_csv(args.case_csv, args.case_id, args.db, args.out)
    else:
        result = analyze_case(args.dna, args.db, args.out)
    if args.explain:
        print_lines(format_processing_trace(result["processing_trace"]))
    print_case_summary(result)


if __name__ == "__main__":
    main()
