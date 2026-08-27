# CSIdentify

> Simulasi investigasi DNA forensik yang menggabungkan profiling Short Tandem Repeats dan Smith-Waterman local alignment.

## Demo

- **[Aplikasi web](https://csidentify.vercel.app)**
- **[Video demo](https://youtu.be/bXs0WNmRoPw)**

## Ringkasan

CSIdentify dikembangkan sebagai proyek kelompok **IF3211 — Komputasi Domain Spesifik**. Sistem menerima sekuens DNA dari tempat kejadian, menghitung repeat STR terpanjang pada setiap marker, membandingkannya dengan database tersangka, lalu menggunakan local sequence alignment sebagai bukti pendukung.

| Aspek | Detail |
|---|---|
| Domain | Bioinformatika dan digital forensics |
| Algoritma | Regex STR profiling, Smith-Waterman alignment |
| Implementasi | Python CLI dan static web dashboard |
| Quality assurance | Unit test untuk algoritma inti dan edge cases |

## Alur Analisis

```text
DNA TKP
  ├── validasi A/C/G/T/N
  ├── cari repeat STR terpanjang per marker
  └── bandingkan dengan setiap tersangka
          ├── STR similarity (70%)
          └── local alignment support (30%)
                  └── ranking kandidat
```

## Fitur

- Validasi dan normalisasi sekuens DNA.
- Jejak posisi setiap run STR yang ditemukan.
- Profil beberapa marker STR dari sampel TKP.
- Smith-Waterman dengan skor match, mismatch, dan gap yang dapat dikonfigurasi.
- Ranking tersangka berdasarkan kombinasi profil dan alignment.
- Output JSON terstruktur untuk audit atau integrasi.
- Dashboard web statis untuk demonstrasi visual.

## Menjalankan CLI

```bash
git clone https://github.com/darrylrayhananta/CSIdentify.git
cd CSIdentify
python3 csidentify.py
```

Program menawarkan tiga sumber input:

- file sampel DNA bawaan;
- sekuens yang diketik langsung; atau
- file DNA lain yang diberikan pengguna.

Hasil tersimpan di `output/analysis_result.json`.

## Menjalankan Pengujian

```bash
python3 -m unittest discover -s tests
```

Pengujian mencakup validasi DNA, STR profiling, alignment, scoring, parsing database, dan output analisis.

## Struktur

```text
.
├── csidentify.py             # Algoritma, CLI, dan ekspor hasil
├── index.html                # Dashboard web
├── data/
│   ├── crime_scene_dna.txt   # Sampel DNA TKP sintetis
│   └── suspects.csv          # Profil tersangka sintetis
├── tests/test_csidentify.py
└── docs/                     # Laporan dan presentasi
```

## Tim

| Nama | NIM |
|---|---|
| Florecita Natawirya | 18223040 |
| Darryl Rayhananta Adenan | 18223042 |
| Fhatika Adhalisman Ryanjani | 18223062 |
| Muhammad Refino Ramadhan | 18223070 |

## Batasan

Data pada repository bersifat sintetis dan sistem dibuat untuk tujuan edukasi. Hasil CSIdentify tidak boleh digunakan sebagai identifikasi forensik nyata atau dasar keputusan hukum.

