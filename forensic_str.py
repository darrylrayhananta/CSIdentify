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
DNA_SEQUENCE_COLUMNS = {"dna_sequence", "dna sequence", "sequence", "sekuens_dna", "sekuens dna"}


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


def smith_waterman_local_alignment(
    query: str,
    subject: str,
    match_score: int = 2,
    mismatch_penalty: int = -1,
    gap_penalty: int = -2,
) -> dict[str, object]:
    """Return the best local alignment between two DNA sequences."""
    query = validate_dna_sequence(query)
    subject = validate_dna_sequence(subject)
    row_count = len(query) + 1
    column_count = len(subject) + 1
    scores = [[0] * column_count for _ in range(row_count)]
    traceback = [[0] * column_count for _ in range(row_count)]
    best_score = 0
    best_position = (0, 0)

    for row in range(1, row_count):
        for column in range(1, column_count):
            diagonal = scores[row - 1][column - 1] + (
                match_score if query[row - 1] == subject[column - 1] else mismatch_penalty
            )
            up = scores[row - 1][column] + gap_penalty
            left = scores[row][column - 1] + gap_penalty
            cell_score = max(0, diagonal, up, left)
            scores[row][column] = cell_score
            if cell_score == 0:
                traceback[row][column] = 0
            elif cell_score == diagonal:
                traceback[row][column] = 1
            elif cell_score == up:
                traceback[row][column] = 2
            else:
                traceback[row][column] = 3

            if cell_score > best_score:
                best_score = cell_score
                best_position = (row, column)

    aligned_query: list[str] = []
    aligned_subject: list[str] = []
    row, column = best_position
    query_end = row
    subject_end = column
    while row > 0 and column > 0 and scores[row][column] > 0:
        direction = traceback[row][column]
        if direction == 1:
            aligned_query.append(query[row - 1])
            aligned_subject.append(subject[column - 1])
            row -= 1
            column -= 1
        elif direction == 2:
            aligned_query.append(query[row - 1])
            aligned_subject.append("-")
            row -= 1
        elif direction == 3:
            aligned_query.append("-")
            aligned_subject.append(subject[column - 1])
            column -= 1
        else:
            break

    aligned_query_text = "".join(reversed(aligned_query))
    aligned_subject_text = "".join(reversed(aligned_subject))
    match_line = "".join(
        "|" if left == right and left != "-" else " "
        for left, right in zip(aligned_query_text, aligned_subject_text)
    )
    aligned_length = len(aligned_query_text)
    matches = match_line.count("|")
    gaps = aligned_query_text.count("-") + aligned_subject_text.count("-")
    mismatches = aligned_length - matches - gaps
    query_aligned_bases = sum(1 for base in aligned_query_text if base != "-")
    subject_aligned_bases = sum(1 for base in aligned_subject_text if base != "-")

    return {
        "algorithm": "Smith-Waterman local alignment",
        "score": best_score,
        "aligned_query": aligned_query_text,
        "match_line": match_line,
        "aligned_subject": aligned_subject_text,
        "query_start": row,
        "query_end": query_end,
        "subject_start": column,
        "subject_end": subject_end,
        "aligned_length": aligned_length,
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "identity_percent": round(100 * matches / aligned_length, 2) if aligned_length else 0.0,
        "query_coverage_percent": round(100 * query_aligned_bases / len(query), 2),
        "subject_coverage_percent": round(100 * subject_aligned_bases / len(subject), 2),
        "scoring": {
            "match": match_score,
            "mismatch": mismatch_penalty,
            "gap": gap_penalty,
        },
    }


def alignment_support_percent(alignment: dict[str, object]) -> float:
    """Weight local identity by query coverage to avoid tiny high-identity matches."""
    return round(
        float(alignment["identity_percent"]) * float(alignment["query_coverage_percent"]) / 100,
        2,
    )


def load_suspects(path: str | Path) -> tuple[list[Suspect], list[str]]:
    """Load suspect profiles from a CSV file.

    The first column is treated as the suspect name. Remaining columns are STR
    marker names whose values must be integer repeat counts. A DNA_Sequence
    column is optional and is used for sequence alignment.
    """
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or len(reader.fieldnames) < 2:
            raise ValueError("CSV must contain a name column and at least one STR marker column")

        name_field = reader.fieldnames[0]
        sequence_field = next(
            (
                field
                for field in reader.fieldnames[1:]
                if field.strip().lower() in DNA_SEQUENCE_COLUMNS
            ),
            None,
        )
        markers = [field for field in reader.fieldnames[1:] if field != sequence_field]
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

            suspect: Suspect = {"name": name, "profile": profile}
            if sequence_field:
                raw_sequence = (row.get(sequence_field) or "").strip()
                if raw_sequence:
                    suspect["dna_sequence"] = validate_dna_sequence(raw_sequence)
            suspects.append(suspect)

    if not suspects:
        raise ValueError("CSV database does not contain suspect rows")
    return suspects, markers


def find_best_matches(
    sample_profile: dict[str, int],
    suspects: list[Suspect],
    sample_sequence: str | None = None,
) -> list[dict[str, object]]:
    """Rank suspects by STR profile fit, then local sequence alignment."""
    markers = list(sample_profile.keys())
    ranked: list[dict[str, object]] = []

    for suspect in suspects:
        suspect_profile = suspect["profile"]
        if not isinstance(suspect_profile, dict):
            raise TypeError("suspect profile must be a dictionary")

        exact_count = sum(sample_profile[m] == int(suspect_profile[m]) for m in markers)
        distance = sum(abs(sample_profile[m] - int(suspect_profile[m])) for m in markers)
        str_score_percent = round(100 * exact_count / len(markers), 2)
        alignment: dict[str, object] | None = None
        sequence_score_percent = 0.0
        suspect_sequence = suspect.get("dna_sequence")
        if sample_sequence and isinstance(suspect_sequence, str):
            alignment = smith_waterman_local_alignment(sample_sequence, suspect_sequence)
            sequence_score_percent = alignment_support_percent(alignment)
        combined_score_percent = round(
            (0.7 * str_score_percent) + (0.3 * sequence_score_percent)
            if alignment is not None
            else str_score_percent,
            2,
        )
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
                "score_percent": str_score_percent,
                "str_score_percent": str_score_percent,
                "distance": distance,
                "differences": differences,
                "alignment": alignment,
                "alignment_score_percent": sequence_score_percent,
                "combined_score_percent": combined_score_percent,
                "is_exact_match": exact_count == len(markers),
            }
        )

    return sorted(
        ranked,
        key=lambda item: (
            -float(item["combined_score_percent"]),
            -int(item["matching_markers"]),
            int(item["distance"]),
            -float(item["alignment_score_percent"]),
            str(item["name"]),
        ),
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
    ranked = find_best_matches(sample_profile, suspects, sequence)
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
                "combined_score_percent": item["combined_score_percent"],
                "distance": item["distance"],
                "is_exact_match": item["is_exact_match"],
                "alignment": item["alignment"],
                "alignment_score_percent": item["alignment_score_percent"],
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
            f"gabungan {suspect_trace['combined_score_percent']}% | "
            f"jarak {suspect_trace['distance']} | {status}"
        )
        lines.append(f"   Marker cocok: {matched}")
        alignment = suspect_trace["alignment"]
        if isinstance(alignment, dict):
            lines.append(
                "   Alignment: "
                f"score {alignment['score']}, identity {alignment['identity_percent']}%, "
                f"coverage TKP {alignment['query_coverage_percent']}%"
            )
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
    output_func(f"Skor gabungan STR + alignment: {best['combined_score_percent']}%")
    alignment = best.get("alignment")
    if isinstance(alignment, dict):
        output_func(
            "Alignment lokal: "
            f"score {alignment['score']}, identity {alignment['identity_percent']}%, "
            f"coverage TKP {alignment['query_coverage_percent']}%"
        )
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
  <title>CSIdentify - DNA Analysis Dashboard</title>
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700;800&family=JetBrains+Mono:wght@500;700&family=Material+Symbols+Outlined:wght,FILL@100..700,0..1&display=swap" rel="stylesheet">
  <style>
    :root {
      --background: #fff8f5;
      --surface: #fff8f5;
      --surface-lowest: #ffffff;
      --surface-low: #fff1e9;
      --surface-container: #ffeadd;
      --surface-high: #ffe3d0;
      --surface-variant: #ffdcc3;
      --primary: #002d78;
      --primary-container: #204494;
      --on-primary-container: #9db6ff;
      --secondary-container: #e0234d;
      --on-secondary: #ffffff;
      --on-background: #2a1708;
      --on-surface-variant: #444651;
      --outline: #747683;
      --outline-variant: #c4c6d3;
      --error: #ba1a1a;
      --error-container: #ffdad6;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      min-height: 100vh;
      font-family: "Inter", ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      color: var(--on-background);
      background: var(--background);
    }
    .brand-area {
      display: flex;
      justify-content: center;
      align-items: center;
      padding: 32px 16px 26px;
      background: var(--surface);
    }
    .brand-lockup {
      display: flex;
      align-items: center;
      gap: 14px;
    }
    .logo-placeholder {
      display: grid;
      place-items: center;
      width: 48px;
      height: 48px;
      border: 1px dashed var(--primary-container);
      border-radius: 4px;
      color: var(--primary-container);
      background: var(--surface-lowest);
      font: 700 10px/14px "JetBrains Mono", monospace;
      letter-spacing: 0.08em;
      text-transform: uppercase;
    }
    .logo-placeholder img {
      width: 100%;
      height: 100%;
      object-fit: cover;
      display: none;
    }
    .brand-name {
      color: var(--primary);
      font-size: 32px;
      line-height: 40px;
      font-weight: 800;
      letter-spacing: -0.02em;
    }
    main {
      width: min(1440px, calc(100% - 32px));
      margin: 0 auto;
      padding: 0 0 80px;
      display: grid;
      grid-template-columns: repeat(12, minmax(0, 1fr));
      gap: 16px;
    }
    .workspace-header,
    .action-area,
    .match-card,
    .chart-card,
    .table-card,
    .analysis-actions,
    footer {
      grid-column: 1 / -1;
    }
    .page-hidden {
      display: none !important;
    }
    .workspace-header {
      margin-bottom: 16px;
    }
    h1 {
      margin: 0;
      color: var(--primary);
      font-size: 32px;
      line-height: 40px;
      font-weight: 800;
      letter-spacing: -0.02em;
    }
    .subtitle {
      margin: 8px 0 0;
      color: var(--on-surface-variant);
      font-size: 14px;
      line-height: 20px;
    }
    .status-chip {
      display: inline-flex;
      align-items: center;
      gap: 8px;
      margin-top: 14px;
      padding: 6px 10px;
      border-radius: 2px;
      background: var(--surface-container);
      color: var(--primary-container);
      font: 500 12px/16px "JetBrains Mono", monospace;
      letter-spacing: 0.05em;
      text-transform: uppercase;
    }
    .data-tile {
      background: var(--surface-lowest);
      border: 1px solid rgba(32, 68, 148, 0.10);
      border-radius: 2px;
      box-shadow: 0 4px 18px rgba(32, 68, 148, 0.04);
    }
    .upload-zone {
      grid-column: span 6;
      min-height: 240px;
      padding: 24px;
      display: flex;
      flex-direction: column;
      align-items: center;
      justify-content: center;
      border: 2px dashed rgba(32, 68, 148, 0.22);
      cursor: pointer;
      text-align: center;
      transition: border-color 160ms ease, background 160ms ease, box-shadow 160ms ease;
    }
    .upload-zone:hover {
      border-color: var(--primary-container);
      background: rgba(32, 68, 148, 0.05);
      box-shadow: 0 4px 10px rgba(32, 68, 148, 0.07);
    }
    .upload-zone.ready {
      border-style: solid;
      background: rgba(32, 68, 148, 0.06);
    }
    .material-symbols-outlined {
      font-family: "Material Symbols Outlined";
      font-weight: normal;
      font-style: normal;
      font-size: 40px;
      line-height: 1;
      letter-spacing: normal;
      text-transform: none;
      display: inline-block;
      white-space: nowrap;
      direction: ltr;
      font-feature-settings: "liga";
      -webkit-font-feature-settings: "liga";
    }
    .upload-zone .material-symbols-outlined {
      color: var(--primary-container);
      margin-bottom: 16px;
      font-variation-settings: "FILL" 1;
    }
    .upload-title {
      margin: 0 0 8px;
      color: var(--primary);
      font-size: 20px;
      line-height: 28px;
      font-weight: 700;
    }
    .upload-copy {
      margin: 0 0 16px;
      color: var(--on-surface-variant);
      font-size: 14px;
      line-height: 20px;
    }
    .file-note {
      min-height: 18px;
      margin-top: 12px;
      color: var(--on-surface-variant);
      font: 500 12px/16px "JetBrains Mono", monospace;
      overflow-wrap: anywhere;
    }
    input[type="file"] {
      position: absolute;
      width: 1px;
      height: 1px;
      opacity: 0;
      pointer-events: none;
    }
    .btn-outline {
      display: inline-flex;
      justify-content: center;
      align-items: center;
      min-height: 42px;
      padding: 11px 22px;
      border: 1px solid var(--primary-container);
      border-radius: 2px;
      background: transparent;
      color: var(--primary-container);
      font-size: 14px;
      line-height: 20px;
      font-weight: 700;
    }
    .action-area {
      display: flex;
      justify-content: center;
      margin: 16px 0 16px;
    }
    button {
      border: 0;
      border-radius: 4px;
      background: var(--primary-container);
      color: #ffffff;
      min-height: 56px;
      padding: 14px 42px;
      display: inline-flex;
      align-items: center;
      justify-content: center;
      gap: 8px;
      font-size: 18px;
      line-height: 24px;
      font-weight: 800;
      cursor: pointer;
      box-shadow: 0 8px 18px rgba(32, 68, 148, 0.18);
      transition: background 160ms ease, transform 160ms ease;
    }
    button:hover {
      background: #001849;
      transform: translateY(-1px);
    }
    .btn-secondary {
      min-height: 44px;
      padding: 10px 18px;
      background: transparent;
      border: 1px solid var(--primary-container);
      color: var(--primary-container);
      box-shadow: none;
      font-size: 14px;
      line-height: 20px;
    }
    .btn-secondary:hover {
      background: rgba(32, 68, 148, 0.05);
      color: var(--primary-container);
    }
    .analysis-actions {
      display: flex;
      justify-content: flex-end;
      margin-bottom: -2px;
    }
    .error {
      display: none;
      grid-column: 1 / -1;
      width: 100%;
      margin: -2px 0 0;
      padding: 12px 14px;
      border-radius: 2px;
      background: var(--error-container);
      color: var(--error);
      font-size: 14px;
      line-height: 20px;
      font-weight: 700;
    }
    .match-card {
      display: none;
      padding: 24px;
      align-items: center;
      justify-content: space-between;
      gap: 24px;
      background: var(--primary-container);
      color: var(--surface-lowest);
      border: 0;
    }
    .match-card.is-visible {
      display: flex;
    }
    .match-kicker,
    .table-kicker {
      margin: 0 0 6px;
      color: var(--on-primary-container);
      font: 500 12px/16px "JetBrains Mono", monospace;
      letter-spacing: 0.05em;
      text-transform: uppercase;
    }
    .match-name-row {
      display: flex;
      flex-wrap: wrap;
      align-items: center;
      gap: 12px;
      font-size: 32px;
      line-height: 40px;
      font-weight: 800;
      letter-spacing: -0.02em;
    }
    .match-meta {
      margin: 8px 0 0;
      color: var(--on-primary-container);
      font: 500 12px/16px "JetBrains Mono", monospace;
    }
    .score-ring {
      flex: 0 0 auto;
      display: grid;
      place-items: center;
      width: 96px;
      height: 96px;
      border-radius: 999px;
      border: 4px solid var(--secondary-container);
      background: var(--surface-lowest);
      color: var(--primary-container);
      font-size: 24px;
      line-height: 32px;
      font-weight: 800;
    }
    .badge-exact,
    .badge-partial,
    .badge-low {
      display: inline-flex;
      align-items: center;
      min-height: 22px;
      padding: 3px 9px;
      border-radius: 2px;
      font: 700 10px/14px "JetBrains Mono", monospace;
      letter-spacing: 0.08em;
      text-transform: uppercase;
    }
    .badge-exact {
      background: var(--surface-lowest);
      color: var(--primary-container);
    }
    .badge-partial {
      background: var(--secondary-container);
      color: var(--on-secondary);
    }
    .badge-low {
      border: 1px solid var(--outline);
      color: var(--outline);
    }
    .chart-card {
      padding: 24px;
      margin-top: 16px;
      overflow: hidden;
    }
    .section-title {
      margin: 0 0 16px;
      padding-bottom: 10px;
      border-bottom: 1px solid var(--outline-variant);
      color: var(--primary);
      font-size: 20px;
      line-height: 28px;
      font-weight: 700;
    }
    #profileComparisonChart {
      min-height: 260px;
      overflow-x: auto;
    }
    .chart-empty {
      display: grid;
      min-height: 220px;
      place-items: center;
      border: 1px dashed var(--outline-variant);
      background: rgba(255, 220, 195, 0.28);
      color: var(--on-surface-variant);
      text-align: center;
      line-height: 20px;
      padding: 24px;
      font-size: 14px;
    }
    .details-grid {
      grid-column: 1 / -1;
      display: grid;
      grid-template-columns: repeat(12, minmax(0, 1fr));
      gap: 16px;
    }
    .profile-card {
      grid-column: span 5;
    }
    .ranking-card {
      grid-column: span 7;
    }
    .trace-card {
      grid-column: 1 / -1;
    }
    .alignment-card {
      grid-column: 1 / -1;
    }
    .table-card {
      margin-top: 16px;
      overflow: hidden;
    }
    .alignment-body {
      padding: 16px;
    }
    .alignment-metrics {
      display: grid;
      grid-template-columns: repeat(4, minmax(120px, 1fr));
      gap: 10px;
      margin-bottom: 16px;
    }
    .metric-tile {
      padding: 12px;
      border: 1px solid var(--outline-variant);
      background: var(--surface);
    }
    .metric-label {
      margin: 0 0 4px;
      color: var(--on-surface-variant);
      font: 500 11px/15px "JetBrains Mono", monospace;
      text-transform: uppercase;
    }
    .metric-value {
      margin: 0;
      color: var(--primary);
      font-size: 20px;
      line-height: 28px;
      font-weight: 800;
    }
    .alignment-block {
      margin: 0;
      padding: 14px;
      overflow-x: auto;
      border: 1px solid var(--outline-variant);
      background: #ffffff;
      color: var(--on-background);
      font: 700 13px/22px "JetBrains Mono", monospace;
      white-space: pre;
    }
    .table-header {
      padding: 16px;
      background: var(--surface-variant);
    }
    .table-header.primary {
      background: var(--primary-container);
    }
    .table-header.secondary {
      background: var(--secondary-container);
    }
    .table-header h3 {
      margin: 0;
      color: var(--on-background);
      font-size: 20px;
      line-height: 28px;
      font-weight: 700;
    }
    .table-header.primary h3,
    .table-header.secondary h3 {
      color: var(--surface-lowest);
    }
    .table-wrap {
      overflow-x: auto;
    }
    table {
      width: 100%;
      border-collapse: collapse;
      overflow: hidden;
      background: var(--surface-lowest);
    }
    tbody tr:nth-child(even) {
      background: rgba(252, 190, 167, 0.20);
    }
    th, td {
      padding: 14px 16px;
      border-bottom: 1px solid rgba(196, 198, 211, 0.65);
      text-align: left;
      vertical-align: top;
    }
    th {
      color: var(--on-surface-variant);
      background: var(--surface-lowest);
      font: 500 12px/16px "JetBrains Mono", monospace;
      letter-spacing: 0.05em;
      text-transform: uppercase;
    }
    td {
      color: var(--on-background);
      font-size: 14px;
      line-height: 20px;
    }
    code {
      display: inline-block;
      max-width: 100%;
      padding: 4px 8px;
      border: 1px solid var(--outline-variant);
      border-radius: 2px;
      background: var(--surface);
      color: var(--on-surface-variant);
      font-family: "JetBrains Mono", monospace;
      overflow-wrap: anywhere;
    }
    .rank {
      display: inline-grid;
      width: 28px;
      height: 28px;
      place-items: center;
      border-radius: 999px;
      background: #dae2ff;
      color: var(--primary-container);
      font-weight: 800;
    }
    footer {
      color: var(--on-surface-variant);
      font-size: 12px;
      line-height: 16px;
    }
    .is-hidden {
      display: none !important;
    }
    @media (max-width: 900px) {
      main {
        width: min(100% - 32px, 720px);
        grid-template-columns: repeat(4, minmax(0, 1fr));
      }
      .brand-name {
        font-size: 28px;
        line-height: 36px;
      }
      .workspace-header,
      .action-area,
      .match-card,
      .chart-card,
      .table-card,
      .analysis-actions,
      footer,
      .details-grid,
      .profile-card,
      .ranking-card,
      .alignment-card,
      .trace-card {
        grid-column: 1 / -1;
      }
      .upload-zone {
        grid-column: 1 / -1;
      }
      .match-card {
        align-items: flex-start;
        flex-direction: column;
      }
      .match-name-row {
        font-size: 26px;
        line-height: 32px;
      }
    }
    @media (max-width: 560px) {
      main {
        width: min(100% - 24px, 420px);
      }
      .brand-area {
        padding-top: 24px;
      }
      .upload-zone {
        min-height: 210px;
        padding: 20px;
      }
      button {
        width: 100%;
      }
      th, td {
        padding: 12px;
      }
      .alignment-metrics {
        grid-template-columns: 1fr;
      }
    }
  </style>
</head>
<body>
  <div class="brand-area">
    <div class="brand-lockup">
      <div class="logo-placeholder">
        <img alt="CSIdentify logo" src="">
        <span>Logo</span>
      </div>
      <span class="brand-name">CSIdentify</span>
    </div>
  </div>

  <main>
    <section id="inputHeader" class="workspace-header">
      <h1>Analysis Workspace</h1>
      <p class="subtitle">Upload forensic data to initiate matching sequences against the suspect database.</p>
      <div class="status-chip"><span class="material-symbols-outlined" style="font-size:16px;">pending</span><span id="analysisStatus">Menunggu file</span></div>
    </section>

    <label class="data-tile upload-zone input-page" id="dnaUploadZone" for="dnaFileInput">
      <span class="material-symbols-outlined">biotech</span>
      <h2 class="upload-title">Crime Scene Profile</h2>
      <p class="upload-copy">Select .txt DNA sequence file</p>
      <span class="btn-outline">Select File</span>
      <input id="dnaFileInput" type="file" accept=".txt,text/plain">
      <span id="dnaFileStatus" class="file-note">Belum ada file DNA TKP.</span>
    </label>

    <label class="data-tile upload-zone input-page" id="suspectUploadZone" for="suspectCsvFileInput">
      <span class="material-symbols-outlined">database</span>
      <h2 class="upload-title">Suspect Database</h2>
      <p class="upload-copy">Select .csv suspect database file</p>
      <span class="btn-outline">Select File</span>
      <input id="suspectCsvFileInput" type="file" accept=".csv,text/csv">
      <span id="suspectFileStatus" class="file-note">Belum ada file CSV tersangka.</span>
    </label>

    <div id="inputActions" class="action-area">
      <button id="analyzeButton" type="button">
        <span class="material-symbols-outlined">play_arrow</span>
        Execute Forensic Analysis
      </button>
    </div>
    <div id="errorBox" class="error" role="alert"></div>

    <section id="analysisHeader" class="workspace-header analysis-page page-hidden">
      <h1>Analysis Result</h1>
      <p class="subtitle">Review STR marker comparison, Smith-Waterman local alignment, suspect ranking, and RegEx trace from the uploaded files.</p>
      <div class="status-chip"><span class="material-symbols-outlined" style="font-size:16px;">verified</span><span id="analysisStatusView">Analisis selesai</span></div>
    </section>

    <div id="analysisActions" class="analysis-actions analysis-page page-hidden">
      <button id="restartButton" class="btn-secondary" type="button">
        <span class="material-symbols-outlined" style="font-size:18px;">restart_alt</span>
        New Analysis
      </button>
    </div>

    <section id="resultSummary" class="data-tile match-card analysis-page page-hidden">
      <div>
        <p class="match-kicker">Primary Match Identified</p>
        <div class="match-name-row">
          <span id="bestNameView">-</span>
          <span id="verdictView" class="badge-exact">READY</span>
        </div>
        <p class="match-meta">Marker cocok: <span id="matchedView">0/0</span> | STR: <span id="strScoreView">0%</span> | Alignment: <span id="alignmentScoreView">0%</span> | Distance: <span id="distanceView">0</span></p>
      </div>
      <div id="scoreView" class="score-ring">0%</div>
    </section>

    <section id="chartSection" class="data-tile chart-card analysis-page page-hidden">
      <h2 class="section-title">STR Marker Comparison</h2>
      <div id="profileComparisonChart">
        <div class="chart-empty">Chart akan muncul setelah kedua file diupload dan dianalisis.</div>
      </div>
    </section>

    <section id="alignmentSection" class="data-tile table-card alignment-card analysis-page page-hidden">
      <div class="table-header primary">
        <h3>Smith-Waterman Local Alignment</h3>
      </div>
      <div class="alignment-body">
        <div class="alignment-metrics">
          <div class="metric-tile">
            <p class="metric-label">Alignment score</p>
            <p id="alignmentRawScoreView" class="metric-value">0</p>
          </div>
          <div class="metric-tile">
            <p class="metric-label">Identity</p>
            <p id="alignmentIdentityView" class="metric-value">0%</p>
          </div>
          <div class="metric-tile">
            <p class="metric-label">Coverage TKP</p>
            <p id="alignmentCoverageView" class="metric-value">0%</p>
          </div>
          <div class="metric-tile">
            <p class="metric-label">Algorithm</p>
            <p id="alignmentAlgorithmView" class="metric-value">SW</p>
          </div>
        </div>
        <pre id="alignmentBlock" class="alignment-block">Upload CSV dengan kolom DNA_Sequence untuk melihat alignment.</pre>
      </div>
    </section>

    <section id="detailsSection" class="details-grid analysis-page page-hidden">
      <article class="data-tile table-card profile-card">
        <div class="table-header">
          <h3>Profil DNA TKP</h3>
        </div>
        <div class="table-wrap">
          <table>
            <thead><tr><th>STR Marker</th><th>Repeat</th></tr></thead>
            <tbody id="profileRows"></tbody>
          </table>
        </div>
      </article>

      <article class="data-tile table-card ranking-card">
        <div class="table-header primary">
          <h3>Suspect Rankings</h3>
        </div>
        <div class="table-wrap">
          <table>
            <thead><tr><th>Rank</th><th>Suspect Name</th><th>Matched</th><th>STR</th><th>Alignment</th><th>Combined</th><th>Status</th></tr></thead>
            <tbody id="suspectRows"></tbody>
          </table>
        </div>
      </article>
    </section>

    <section id="traceSection" class="data-tile table-card trace-card analysis-page page-hidden">
      <div class="table-header secondary">
        <h3>Trace RegEx Patterns</h3>
      </div>
      <div class="table-wrap">
        <table>
          <thead><tr><th>Marker</th><th>Regex Pattern</th><th>Matches Found</th><th>Longest Repeat</th></tr></thead>
          <tbody id="traceRows"></tbody>
        </table>
      </div>
    </section>

    <footer>Model edukatif: hasil STR ini untuk simulasi komputasi, bukan pengganti validasi laboratorium forensik.</footer>
  </main>

  <script>
    const $ = (id) => document.getElementById(id);
    const REQUIRED_FILE_MESSAGE = 'Upload file .txt DNA TKP dan file .csv tersangka terlebih dahulu.';
    const DNA_SEQUENCE_HEADERS = ['dna_sequence', 'dna sequence', 'sequence', 'sekuens_dna', 'sekuens dna'];
    const CHART_COLORS = ['#204494', '#e0234d', '#7c2c45', '#b3c5ff', '#ffb1c3', '#002d78', '#ffdcc3'];

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
      const sequenceColumnIndex = headers.findIndex((header, index) => (
        index > 0 && DNA_SEQUENCE_HEADERS.includes(header.trim().toLowerCase())
      ));
      const markerColumns = headers
        .map((header, index) => ({ header: header.trim(), index }))
        .filter((column) => column.index > 0 && column.index !== sequenceColumnIndex && column.header);
      const markers = markerColumns.map((column) => column.header);
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
        markerColumns.forEach((column) => {
          const rawValue = cells[column.index];
          const value = Number.parseInt(rawValue, 10);
          if (!Number.isInteger(value) || value < 0) {
            throw new Error(`Nilai marker ${column.header} invalid di baris ${rowNumber}.`);
          }
          profile[column.header] = value;
        });
        const rawSequence = sequenceColumnIndex > -1 ? cells[sequenceColumnIndex] : '';
        return {
          name,
          profile,
          dnaSequence: rawSequence ? validateDna(rawSequence) : null,
        };
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

    function smithWatermanLocalAlignment(queryInput, subjectInput, matchScore = 2, mismatchPenalty = -1, gapPenalty = -2) {
      const query = validateDna(queryInput);
      const subject = validateDna(subjectInput);
      const rowCount = query.length + 1;
      const columnCount = subject.length + 1;
      const scores = Array.from({ length: rowCount }, () => Array(columnCount).fill(0));
      const traceback = Array.from({ length: rowCount }, () => Array(columnCount).fill(0));
      let bestScore = 0;
      let bestRow = 0;
      let bestColumn = 0;

      for (let row = 1; row < rowCount; row += 1) {
        for (let column = 1; column < columnCount; column += 1) {
          const diagonal = scores[row - 1][column - 1] + (query[row - 1] === subject[column - 1] ? matchScore : mismatchPenalty);
          const up = scores[row - 1][column] + gapPenalty;
          const left = scores[row][column - 1] + gapPenalty;
          const cellScore = Math.max(0, diagonal, up, left);
          scores[row][column] = cellScore;
          if (cellScore === 0) {
            traceback[row][column] = 0;
          } else if (cellScore === diagonal) {
            traceback[row][column] = 1;
          } else if (cellScore === up) {
            traceback[row][column] = 2;
          } else {
            traceback[row][column] = 3;
          }
          if (cellScore > bestScore) {
            bestScore = cellScore;
            bestRow = row;
            bestColumn = column;
          }
        }
      }

      const alignedQuery = [];
      const alignedSubject = [];
      let row = bestRow;
      let column = bestColumn;
      const queryEnd = row;
      const subjectEnd = column;
      while (row > 0 && column > 0 && scores[row][column] > 0) {
        const direction = traceback[row][column];
        if (direction === 1) {
          alignedQuery.push(query[row - 1]);
          alignedSubject.push(subject[column - 1]);
          row -= 1;
          column -= 1;
        } else if (direction === 2) {
          alignedQuery.push(query[row - 1]);
          alignedSubject.push('-');
          row -= 1;
        } else if (direction === 3) {
          alignedQuery.push('-');
          alignedSubject.push(subject[column - 1]);
          column -= 1;
        } else {
          break;
        }
      }

      const alignedQueryText = alignedQuery.reverse().join('');
      const alignedSubjectText = alignedSubject.reverse().join('');
      const matchLine = [...alignedQueryText].map((base, index) => (
        base === alignedSubjectText[index] && base !== '-' ? '|' : ' '
      )).join('');
      const alignedLength = alignedQueryText.length;
      const matches = [...matchLine].filter((char) => char === '|').length;
      const gaps = [...alignedQueryText, ...alignedSubjectText].filter((char) => char === '-').length;
      const queryAlignedBases = [...alignedQueryText].filter((base) => base !== '-').length;
      const subjectAlignedBases = [...alignedSubjectText].filter((base) => base !== '-').length;

      return {
        algorithm: 'Smith-Waterman local alignment',
        score: bestScore,
        alignedQuery: alignedQueryText,
        matchLine,
        alignedSubject: alignedSubjectText,
        queryStart: row,
        queryEnd,
        subjectStart: column,
        subjectEnd,
        alignedLength,
        matches,
        mismatches: alignedLength - matches - gaps,
        gaps,
        identityPercent: alignedLength ? Math.round((matches / alignedLength) * 10000) / 100 : 0,
        queryCoveragePercent: Math.round((queryAlignedBases / query.length) * 10000) / 100,
        subjectCoveragePercent: Math.round((subjectAlignedBases / subject.length) * 10000) / 100,
      };
    }

    function alignmentSupportPercent(alignment) {
      return Math.round((alignment.identityPercent * alignment.queryCoveragePercent) * 100) / 10000;
    }

    function rankSuspects(sampleProfile, suspects, markers, sampleSequence) {
      return suspects.map((suspect) => {
        const matchedMarkers = markers.filter((marker) => suspect.profile[marker] === sampleProfile[marker]);
        const distance = markers.reduce((sum, marker) => (
          sum + Math.abs((suspect.profile[marker] || 0) - (sampleProfile[marker] || 0))
        ), 0);
        const scorePercent = Math.round((matchedMarkers.length / markers.length) * 10000) / 100;
        const alignment = sampleSequence && suspect.dnaSequence
          ? smithWatermanLocalAlignment(sampleSequence, suspect.dnaSequence)
          : null;
        const alignmentScorePercent = alignment ? alignmentSupportPercent(alignment) : 0;
        const combinedScorePercent = alignment
          ? Math.round(((0.7 * scorePercent) + (0.3 * alignmentScorePercent)) * 100) / 100
          : scorePercent;
        return {
          name: suspect.name,
          profile: suspect.profile,
          matchingMarkers: matchedMarkers.length,
          markerCount: markers.length,
          scorePercent,
          alignment,
          alignmentScorePercent,
          combinedScorePercent,
          distance,
          isExactMatch: matchedMarkers.length === markers.length,
        };
      }).sort((left, right) => (
        right.combinedScorePercent - left.combinedScorePercent
        || right.matchingMarkers - left.matchingMarkers
        || left.distance - right.distance
        || right.alignmentScorePercent - left.alignmentScorePercent
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

    function showInputPage() {
      document.querySelectorAll('.analysis-page').forEach((element) => element.classList.add('page-hidden'));
      document.querySelectorAll('.input-page').forEach((element) => element.classList.remove('page-hidden'));
      $('inputHeader').classList.remove('page-hidden');
      $('inputActions').classList.remove('page-hidden');
    }

    function showAnalysisPage() {
      document.querySelectorAll('.input-page').forEach((element) => element.classList.add('page-hidden'));
      document.querySelectorAll('.analysis-page').forEach((element) => element.classList.remove('page-hidden'));
      $('inputHeader').classList.add('page-hidden');
      $('inputActions').classList.add('page-hidden');
    }

    function renderEmptyState() {
      $('analysisStatus').textContent = 'Menunggu file';
      showInputPage();
      $('profileRows').innerHTML = '';
      $('suspectRows').innerHTML = '';
      $('traceRows').innerHTML = '';
      $('alignmentBlock').textContent = 'Upload CSV dengan kolom DNA_Sequence untuk melihat alignment.';
      $('profileComparisonChart').innerHTML = '<div class="chart-empty">Upload kedua file untuk melihat chart profil STR.</div>';
    }

    function renderGroupedProfileChart(state) {
      const markers = state.markers;
      const series = [
        { label: 'DNA TKP', profile: state.sampleProfile },
        ...state.ranked.map((suspect) => ({ label: suspect.name, profile: suspect.profile })),
      ];
      const maxValue = Math.max(1, ...markers.flatMap((marker) => series.map((item) => Number(item.profile[marker] || 0))));
      const yMax = Math.max(5, Math.ceil(maxValue / 2) * 2);
      const width = Math.max(820, markers.length * 142 + 160);
      const height = 300;
      const left = 58;
      const right = 24;
      const top = 34;
      const bottom = 64;
      const plotWidth = width - left - right;
      const plotHeight = height - top - bottom;
      const groupWidth = plotWidth / markers.length;
      const barGap = 3;
      const barWidth = Math.max(8, Math.min(20, (groupWidth - 28) / series.length - barGap));
      const tickCount = Math.min(yMax, 6);
      let svg = `<svg viewBox="0 0 ${width} ${height}" width="100%" height="${height}" role="img" aria-label="Profil STR DNA TKP vs tersangka">`;
      svg += '<rect width="100%" height="100%" fill="#ffffff"></rect>';

      for (let tick = 0; tick <= tickCount; tick += 1) {
        const value = Math.round((yMax / tickCount) * tick);
        const y = top + plotHeight - (value / yMax) * plotHeight;
        svg += `<line x1="${left}" y1="${y}" x2="${width - right}" y2="${y}" stroke="#c4c6d3" stroke-opacity="0.58" stroke-width="1"></line>`;
        svg += `<text x="${left - 14}" y="${y + 4}" fill="#444651" font-size="11" text-anchor="end" font-family="JetBrains Mono, monospace">${value}</text>`;
      }

      markers.forEach((marker, markerIndex) => {
        const groupLeft = left + markerIndex * groupWidth;
        const center = groupLeft + groupWidth / 2;
        svg += `<text x="${center}" y="${height - 34}" fill="#444651" font-size="11" font-weight="700" text-anchor="middle" font-family="JetBrains Mono, monospace">${escapeHtml(marker)}</text>`;
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
            svg += `<text x="${x + barWidth / 2}" y="${Math.max(top + 11, y - 5)}" fill="#2a1708" font-size="9" font-weight="900" text-anchor="middle" font-family="Inter, Arial">${value}</text>`;
          }
        });
      });

      const legendY = height - 12;
      let legendX = left;
      series.forEach((item, index) => {
        const color = CHART_COLORS[index % CHART_COLORS.length];
        const label = item.label.length > 16 ? `${item.label.slice(0, 15)}...` : item.label;
        svg += `<rect x="${legendX}" y="${legendY - 10}" width="10" height="10" rx="2" fill="${color}"></rect>`;
        svg += `<text x="${legendX + 15}" y="${legendY}" fill="#444651" font-size="10" font-family="JetBrains Mono, monospace">${escapeHtml(label)}</text>`;
        legendX += Math.max(92, label.length * 6 + 32);
      });

      svg += '</svg>';
      $('profileComparisonChart').innerHTML = svg;
    }

    function renderState(state) {
      const best = state.ranked[0];
      $('analysisStatus').textContent = 'Analisis selesai';
      $('analysisStatusView').textContent = 'Analisis selesai';
      $('resultSummary').classList.add('is-visible');
      showAnalysisPage();

      $('verdictView').textContent = best.isExactMatch ? '100% MATCH' : 'PARTIAL MATCH';
      $('verdictView').className = best.isExactMatch ? 'badge-exact' : 'badge-partial';
      $('bestNameView').textContent = best.name;
      $('scoreView').textContent = `${best.combinedScorePercent}%`;
      $('matchedView').textContent = `${best.matchingMarkers}/${best.markerCount}`;
      $('strScoreView').textContent = `${best.scorePercent}%`;
      $('alignmentScoreView').textContent = `${best.alignmentScorePercent}%`;
      $('distanceView').textContent = String(best.distance);
      if (best.alignment) {
        $('alignmentRawScoreView').textContent = String(best.alignment.score);
        $('alignmentIdentityView').textContent = `${best.alignment.identityPercent}%`;
        $('alignmentCoverageView').textContent = `${best.alignment.queryCoveragePercent}%`;
        $('alignmentAlgorithmView').textContent = 'SW local';
        $('alignmentBlock').textContent = [
          `TKP        ${best.alignment.alignedQuery}`,
          `           ${best.alignment.matchLine}`,
          `Tersangka  ${best.alignment.alignedSubject}`,
          '',
          `Region TKP: ${best.alignment.queryStart}-${best.alignment.queryEnd} | Region tersangka: ${best.alignment.subjectStart}-${best.alignment.subjectEnd}`,
        ].join('\n');
      } else {
        $('alignmentRawScoreView').textContent = '0';
        $('alignmentIdentityView').textContent = '0%';
        $('alignmentCoverageView').textContent = '0%';
        $('alignmentAlgorithmView').textContent = 'SW local';
        $('alignmentBlock').textContent = 'CSV tersangka belum memiliki kolom DNA_Sequence, sehingga ranking hanya memakai profil STR.';
      }
      $('profileRows').innerHTML = state.markers.map((marker) => (
        `<tr><td><strong>${escapeHtml(marker)}</strong></td><td>${state.sampleProfile[marker]}</td></tr>`
      )).join('');
      $('suspectRows').innerHTML = state.ranked.map((item, index) => (
        `<tr>
          <td><span class="rank">${index + 1}</span></td>
          <td><strong>${escapeHtml(item.name)}</strong></td>
          <td>${item.matchingMarkers}/${item.markerCount}</td>
          <td>${item.scorePercent}%</td>
          <td>${item.alignmentScorePercent}%</td>
          <td>${item.combinedScorePercent}%</td>
          <td><span class="${item.isExactMatch ? 'badge-exact' : 'badge-partial'}">${item.isExactMatch ? 'Exact match' : 'Partial match'}</span></td>
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

    function restartAnalysis() {
      $('dnaFileInput').value = '';
      $('suspectCsvFileInput').value = '';
      $('dnaFileStatus').textContent = 'Belum ada file DNA TKP.';
      $('suspectFileStatus').textContent = 'Belum ada file CSV tersangka.';
      $('dnaUploadZone').classList.remove('ready');
      $('suspectUploadZone').classList.remove('ready');
      clearError();
      renderEmptyState();
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
        const ranked = rankSuspects(profiled.profile, parsed.suspects, parsed.markers, sequence);
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
      $('dnaUploadZone').classList.toggle('ready', Boolean(file));
      clearError();
    });

    $('suspectCsvFileInput').addEventListener('change', (event) => {
      const file = event.target.files[0];
      $('suspectFileStatus').textContent = file ? file.name : 'Belum ada file CSV tersangka.';
      $('suspectUploadZone').classList.toggle('ready', Boolean(file));
      clearError();
    });

    $('analyzeButton').addEventListener('click', runInteractiveAnalysis);
    $('restartButton').addEventListener('click', restartAnalysis);
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
