import csv
import tempfile
import unittest
from pathlib import Path

from forensic_str import (
    analyze_case_csv,
    analyze_sequence,
    build_html_dashboard,
    build_processing_trace,
    find_best_matches,
    find_repeat_runs,
    load_case_from_csv,
    load_suspects,
    longest_consecutive_repeats,
    profile_sequence,
    save_case_to_csv,
)


class ForensicStrTests(unittest.TestCase):
    def test_longest_consecutive_repeats_uses_longest_run_not_total(self):
        sequence = "AAAGATAGATCCAGATAGATAGATAGATT"

        self.assertEqual(longest_consecutive_repeats(sequence, "AGAT"), 4)

    def test_profile_sequence_counts_each_marker(self):
        sequence = "AATGAATGAATGCCAGATAGATTTT"

        self.assertEqual(
            profile_sequence(sequence, ["AATG", "AGAT", "TATC"]),
            {"AATG": 3, "AGAT": 2, "TATC": 0},
        )

    def test_find_best_matches_identifies_exact_match_and_ranks_near_misses(self):
        sample = {"AGAT": 4, "AATG": 3, "TATC": 1}
        suspects = [
            {"name": "Alya", "profile": {"AGAT": 4, "AATG": 3, "TATC": 1}},
            {"name": "Bima", "profile": {"AGAT": 4, "AATG": 2, "TATC": 1}},
            {"name": "Citra", "profile": {"AGAT": 0, "AATG": 0, "TATC": 0}},
        ]

        ranked = find_best_matches(sample, suspects)

        self.assertEqual(ranked[0]["name"], "Alya")
        self.assertTrue(ranked[0]["is_exact_match"])
        self.assertEqual(ranked[1]["name"], "Bima")
        self.assertGreater(ranked[1]["score_percent"], ranked[2]["score_percent"])

    def test_regex_marker_is_escaped_before_matching(self):
        sequence = "A.GA.GA.GTA.G"

        self.assertEqual(longest_consecutive_repeats(sequence, "A.G"), 3)

    def test_load_suspects_reads_marker_columns_as_integers(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            csv_path = Path(tmp_dir) / "suspects.csv"
            with csv_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["Nama", "AGAT", "AATG"])
                writer.writeheader()
                writer.writerow({"Nama": "Dara", "AGAT": "5", "AATG": "2"})

            suspects, markers = load_suspects(csv_path)

        self.assertEqual(markers, ["AGAT", "AATG"])
        self.assertEqual(suspects, [{"name": "Dara", "profile": {"AGAT": 5, "AATG": 2}}])

    def test_find_repeat_runs_returns_positions_and_pattern_for_demo_trace(self):
        runs = find_repeat_runs("CCAGATAGATTTAATGAATG", "AGAT")

        self.assertEqual(
            runs,
            [
                {
                    "start": 2,
                    "end": 10,
                    "repeat_count": 2,
                    "matched_sequence": "AGATAGAT",
                    "pattern": "(?:AGAT)+",
                }
            ],
        )

    def test_analyze_sequence_accepts_user_supplied_dna_without_file(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            csv_path = Path(tmp_dir) / "suspects.csv"
            csv_path.write_text(
                "Nama,AGAT,AATG\nAlya,2,1\nBima,1,1\n",
                encoding="utf-8",
            )

            result = analyze_sequence("CCAGATAGATTAATG", csv_path, Path(tmp_dir) / "out")

        self.assertEqual(result["sample_profile"], {"AGAT": 2, "AATG": 1})
        self.assertEqual(result["ranked_suspects"][0]["name"], "Alya")

    def test_save_case_to_csv_persists_user_dna_for_later_processing(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            case_path = Path(tmp_dir) / "cases.csv"

            case_id = save_case_to_csv(case_path, "TKP Demo", "ccagatagat", "input manual")
            loaded = load_case_from_csv(case_path, case_id)

        self.assertEqual(case_id, "CASE-001")
        self.assertEqual(
            loaded,
            {
                "case_id": "CASE-001",
                "label": "TKP Demo",
                "dna_sequence": "CCAGATAGAT",
                "source": "input manual",
            },
        )

    def test_analyze_case_csv_processes_dna_from_case_csv(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            case_path = Path(tmp_dir) / "cases.csv"
            suspects_path = Path(tmp_dir) / "suspects.csv"
            suspects_path.write_text(
                "Nama,AGAT,AATG\nAlya,2,1\nBima,1,1\n",
                encoding="utf-8",
            )
            case_id = save_case_to_csv(case_path, "TKP Demo", "CCAGATAGATTAATG", "test")

            result = analyze_case_csv(case_path, case_id, suspects_path, Path(tmp_dir) / "out")

        self.assertEqual(result["case"]["case_id"], "CASE-001")
        self.assertEqual(result["sample_profile"], {"AGAT": 2, "AATG": 1})
        self.assertEqual(result["ranked_suspects"][0]["name"], "Alya")

    def test_build_processing_trace_explains_marker_scans_and_ranking(self):
        sample = {"AGAT": 2, "AATG": 1}
        suspects = [
            {"name": "Alya", "profile": {"AGAT": 2, "AATG": 1}},
            {"name": "Bima", "profile": {"AGAT": 1, "AATG": 1}},
        ]

        trace = build_processing_trace("CCAGATAGATTAATG", ["AGAT", "AATG"], sample, find_best_matches(sample, suspects))

        self.assertEqual(trace["markers"][0]["marker"], "AGAT")
        self.assertEqual(trace["markers"][0]["longest_repeat"], 2)
        self.assertEqual(trace["suspects"][0]["name"], "Alya")
        self.assertEqual(trace["suspects"][0]["matched_markers"], ["AGAT", "AATG"])
        self.assertIn("AGAT: tersangka 1 vs TKP 2", trace["suspects"][1]["difference_notes"])

    def test_html_dashboard_contains_browser_interactive_analyzer(self):
        html = build_html_dashboard()

        self.assertIn('id="dnaFileInput"', html)
        self.assertIn('accept=".txt,text/plain"', html)
        self.assertIn('id="suspectCsvFileInput"', html)
        self.assertIn('accept=".csv,text/csv"', html)
        self.assertIn('id="analyzeButton"', html)
        self.assertIn("function readFileAsText", html)
        self.assertIn("function runInteractiveAnalysis", html)
        self.assertIn("function parseSuspectCsv", html)
        self.assertIn('id="profileComparisonChart"', html)
        self.assertIn("function renderGroupedProfileChart", html)
        self.assertIn("Upload file .txt DNA TKP dan file .csv tersangka terlebih dahulu.", html)
        self.assertNotIn('id="dnaInput"', html)
        self.assertNotIn('id="suspectCsvInput"', html)
        self.assertNotIn("downloadCaseCsvButton", html)
        self.assertNotIn("resetButton", html)
        self.assertNotIn("caseLabelInput", html)
        self.assertNotIn("Alya", html)


if __name__ == "__main__":
    unittest.main()
