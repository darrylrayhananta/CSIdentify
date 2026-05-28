import csv
import tempfile
import unittest
from pathlib import Path

from csidentify import (
    analyze_sequence,
    build_processing_trace,
    find_best_matches,
    find_repeat_runs,
    load_suspects,
    longest_consecutive_repeats,
    main,
    profile_sequence,
    run_interactive,
    smith_waterman_local_alignment,
)


class CsidentifyTests(unittest.TestCase):
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

    def test_load_suspects_reads_optional_dna_sequence_column(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            csv_path = Path(tmp_dir) / "suspects.csv"
            csv_path.write_text(
                "Nama,AGAT,AATG,DNA_Sequence\nDara,2,1,ccagatagattaatg\n",
                encoding="utf-8",
            )

            suspects, markers = load_suspects(csv_path)

        self.assertEqual(markers, ["AGAT", "AATG"])
        self.assertEqual(
            suspects,
            [
                {
                    "name": "Dara",
                    "profile": {"AGAT": 2, "AATG": 1},
                    "dna_sequence": "CCAGATAGATTAATG",
                }
            ],
        )

    def test_smith_waterman_local_alignment_finds_best_matching_region(self):
        alignment = smith_waterman_local_alignment("TTAGATAGATCC", "GGAGATAGATGG")

        self.assertEqual(alignment["algorithm"], "Smith-Waterman local alignment")
        self.assertEqual(alignment["aligned_query"], "AGATAGAT")
        self.assertEqual(alignment["match_line"], "||||||||")
        self.assertEqual(alignment["aligned_subject"], "AGATAGAT")
        self.assertEqual(alignment["score"], 16)
        self.assertEqual(alignment["identity_percent"], 100.0)
        self.assertEqual(alignment["query_start"], 2)
        self.assertEqual(alignment["query_end"], 10)

    def test_analyze_sequence_uses_alignment_to_break_equal_str_profiles(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            csv_path = Path(tmp_dir) / "suspects.csv"
            csv_path.write_text(
                "Nama,AGAT,AATG,DNA_Sequence\n"
                "Bima,2,1,TTTTTTTTTTTTTTT\n"
                "Alya,2,1,CCAGATAGATTAATG\n",
                encoding="utf-8",
            )

            result = analyze_sequence("CCAGATAGATTAATG", csv_path, Path(tmp_dir) / "out")

        self.assertEqual(result["ranked_suspects"][0]["name"], "Alya")
        self.assertEqual(result["ranked_suspects"][0]["alignment"]["identity_percent"], 100.0)
        self.assertEqual(result["ranked_suspects"][0]["combined_score_percent"], 100.0)

    def test_analyze_sequence_writes_json_without_terminal_dashboard_artifact(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            csv_path = Path(tmp_dir) / "suspects.csv"
            output_dir = Path(tmp_dir) / "out"
            csv_path.write_text(
                "Nama,AGAT,AATG,DNA_Sequence\nAlya,2,1,CCAGATAGATTAATG\n",
                encoding="utf-8",
            )

            result = analyze_sequence("CCAGATAGATTAATG", csv_path, output_dir)

            self.assertTrue(result["json_path"].exists())
            self.assertFalse((output_dir / "investigation_dashboard.html").exists())
            self.assertNotIn("html_path", result)

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

    def test_interactive_mode_uses_default_output_and_always_prints_trace(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_dir = Path(tmp_dir) / "output"
            answers = iter(["1", "data/crime_scene_dna.txt", "data/suspects.csv"])
            prompts = []
            lines = []

            run_interactive(
                input_func=lambda prompt: (prompts.append(prompt), next(answers))[1],
                output_func=lines.append,
                animate=False,
                output_dir=str(output_dir),
            )

            self.assertTrue((output_dir / "analysis_result.json").exists())
            self.assertFalse(any("CSV kasus" in prompt for prompt in prompts))
            self.assertFalse(any("Folder output" in prompt for prompt in prompts))
            self.assertFalse(any("Tampilkan jejak" in prompt for prompt in prompts))
            self.assertTrue(any("Path file DNA TKP" in prompt for prompt in prompts))
            self.assertFalse(Path("data/cases.csv").exists())
            self.assertFalse(any("Mode Interaktif" in line for line in lines))
            self.assertTrue(any("Analisis DNA dimulai" in line for line in lines))
            self.assertTrue(any("=== Jejak Processing STR ===" in line for line in lines))

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

    def test_main_defaults_to_interactive_mode(self):
        calls = []

        main(interactive_runner=lambda: calls.append("interactive"))

        self.assertEqual(calls, ["interactive"])

if __name__ == "__main__":
    unittest.main()
