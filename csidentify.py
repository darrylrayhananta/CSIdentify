from __future__ import annotations

import csv
import json
import re
import time
from pathlib import Path 
from typing import Callable, Iterable


Suspect = dict[str, object]
DNA_SEQUENCE_COLUMNS = {"dna_sequence", "dna sequence", "sequence", "sekuens_dna", "sekuens dna"}
DEFAULT_DNA_PATH = "data/crime_scene_dna.txt"
DEFAULT_DB_PATH = "data/suspects.csv"
DEFAULT_OUTPUT_DIR = "output"


def compact_sequence(sequence: str) -> str:
    """Ubah sekuens menjadi huruf besar dan hapus whitespace."""
    return re.sub(r"\s+", "", sequence).upper()


def validate_dna_sequence(sequence: str) -> str:
    sequence = compact_sequence(sequence)
    invalid = sorted(set(sequence) - set("ACGTN"))
    if invalid:
        raise ValueError(f"Sekuens DNA mengandung karakter tidak valid: {', '.join(invalid)}")
    if not sequence:
        raise ValueError("Sekuens DNA tidak boleh kosong")
    return sequence


def read_dna_sequence(path: str | Path) -> str:
    try:
        return validate_dna_sequence(Path(path).read_text(encoding="utf-8"))
    except ValueError as exc:
        raise ValueError(str(exc).replace("Sekuens DNA", "File DNA")) from exc


def find_repeat_runs(sequence: str, marker: str) -> list[dict[str, object]]:
    """Cari semua run STR berurutan beserta posisinya untuk jejak analisis."""
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
    """Hitung repeat berurutan terpanjang untuk satu marker."""
    runs = find_repeat_runs(sequence, marker)
    return max((int(run["repeat_count"]) for run in runs), default=0)


def profile_sequence(sequence: str, markers: Iterable[str]) -> dict[str, int]:
    """Buat profil STR dari sebuah sekuens DNA."""
    return {marker: longest_consecutive_repeats(sequence, marker) for marker in markers}


def smith_waterman_local_alignment(
    query: str,
    subject: str,
    match_score: int = 2,
    mismatch_penalty: int = -1,
    gap_penalty: int = -2,
) -> dict[str, object]:
    """Cari local alignment terbaik antara dua sekuens DNA."""
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
    """Bobot identity alignment dengan coverage agar match pendek tidak berlebihan."""
    return round(
        float(alignment["identity_percent"]) * float(alignment["query_coverage_percent"]) / 100,
        2,
    )


def load_suspects(path: str | Path) -> tuple[list[Suspect], list[str]]:
    """Baca profil tersangka dari CSV.

    Kolom pertama adalah nama tersangka. Kolom marker STR berisi jumlah repeat.
    Kolom DNA_Sequence bersifat opsional dan digunakan untuk alignment.
    """
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or len(reader.fieldnames) < 2:
            raise ValueError("CSV harus memiliki kolom nama dan minimal satu kolom marker STR")

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
                raise ValueError(f"Nama tersangka kosong pada baris {row_number}")

            profile: dict[str, int] = {}
            for marker in markers:
                raw_value = (row.get(marker) or "").strip()
                try:
                    profile[marker] = int(raw_value)
                except ValueError as exc:
                    raise ValueError(
                        f"Jumlah repeat tidak valid untuk {marker!r} pada baris {row_number}: {raw_value!r}"
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
    """Urutkan tersangka berdasarkan kecocokan profil STR dan alignment lokal."""
    markers = list(sample_profile.keys())
    ranked: list[dict[str, object]] = []

    for suspect in suspects:
        suspect_profile = suspect["profile"]
        if not isinstance(suspect_profile, dict):
            raise TypeError("Profil tersangka harus berupa dictionary")

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
    json_path = output_path / "analysis_result.json"

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
        "json_path": json_path,
    }


def build_processing_trace(
    sequence: str,
    markers: Iterable[str],
    sample_profile: dict[str, int],
    ranked: list[dict[str, object]],
) -> dict[str, object]:
    """Buat data jejak proses untuk ditampilkan di terminal."""
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
            raise TypeError("Profil tersangka harus berupa dictionary")
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
        raise TypeError("Jejak marker harus berupa list")
    max_repeat = max(
        [1] + [int(marker_trace["longest_repeat"]) for marker_trace in marker_traces]
    )

    for marker_trace in marker_traces:
        runs = marker_trace["runs"]
        if not isinstance(runs, list):
            raise TypeError("Run STR harus berupa list")
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
    output_func(f"JSON: {result['json_path']}")


def run_interactive(
    input_func: Callable[[str], str] = input,
    output_func: Callable[[str], None] = print,
    animate: bool = True,
    dna_path: str = DEFAULT_DNA_PATH,
    db_path: str = DEFAULT_DB_PATH,
    output_dir: str = DEFAULT_OUTPUT_DIR,
) -> None:
    output_func("=== CSIdentify ===")
    output_func("1. Pakai file DNA TKP")
    output_func("2. Ketik DNA TKP manual")
    output_func("3. Pakai file DNA sendiri")
    choice = (input_func("Pilih mode [1]: ").strip() or "1")

    dna_sequence: str | None = None
    label = "DNA TKP Demo"
    source = f"file: {dna_path}"
    if choice == "2":
        label = input_func("Label kasus [DNA TKP Manual]: ").strip() or "DNA TKP Manual"
        dna_sequence = validate_dna_sequence(input_func("Tempel sekuens DNA TKP: "))
        source = "input manual"
    elif choice == "3":
        selected_dna_path = input_func(f"Path file DNA [{dna_path}]: ").strip() or dna_path
        dna_sequence = read_dna_sequence(selected_dna_path)
        label = input_func("Label kasus [DNA TKP dari file]: ").strip() or "DNA TKP dari file"
        source = f"file: {selected_dna_path}"
    else:
        selected_dna_path = input_func(f"Path file DNA TKP [{dna_path}]: ").strip() or dna_path
        dna_sequence = read_dna_sequence(selected_dna_path)
        source = f"file: {selected_dna_path}"

    selected_db_path = input_func(f"Path database tersangka [{db_path}]: ").strip() or db_path

    case_info = {
        "label": label,
        "dna_sequence": dna_sequence,
        "source": source,
    }
    output_func("\nAnalisis DNA dimulai...")
    result = analyze_sequence(dna_sequence, selected_db_path, output_dir, case_info=case_info)
    print_lines(
        format_processing_trace(result["processing_trace"]),
        output_func=output_func,
        animate=animate,
    )
    print_case_summary(result, output_func=output_func)


def main(
    interactive_runner: Callable[[], None] = run_interactive,
) -> None:
    interactive_runner()


if __name__ == "__main__":
    main()
