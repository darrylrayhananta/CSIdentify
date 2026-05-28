# Forensic STR Profiler

Program ini mensimulasikan investigasi DNA forensik berbasis Short Tandem Repeats (STR) dan sequence alignment. DNA dari TKP diprofilkan dengan RegEx untuk mencari jumlah pengulangan berurutan terpanjang pada setiap marker, lalu profil tersebut dibandingkan dengan database tersangka. Setelah itu, sistem melakukan Smith-Waterman local alignment antara sekuens DNA TKP dan sekuens DNA tersangka untuk memvalidasi kemiripan fragmen DNA.

## Struktur

- `csidentify.py`: modul utama dan CLI.
- `data/crime_scene_dna.txt`: sampel DNA TKP sintetis untuk demo.
- `data/suspects.csv`: database profil STR 10 tersangka beserta sekuens DNA sintetis untuk alignment.
- `tests/test_csidentify.py`: unit test untuk algoritma inti.
- `output/`: folder hasil JSON setelah program dijalankan.

## Cara Menjalankan

```bash
python3 csidentify.py
```

Secara default, program berjalan dalam mode interaktif. Program akan meminta pilihan input DNA:

- pakai DNA TKP demo,
- ketik sekuens DNA sendiri,
- atau pakai file DNA lain.

Setelah input diberikan, program langsung menganalisis DNA tersebut. Program juga akan menampilkan jejak processing: pola RegEx yang dipakai, posisi run STR yang ditemukan, repeat terpanjang per marker, skor Smith-Waterman, lalu skoring tiap tersangka.

Kalau ingin tetap non-interaktif tetapi prosesnya terlihat:

```bash
python3 csidentify.py --batch --explain --dna data/crime_scene_dna.txt --db data/suspects.csv --out output
```

Hasil program terminal:

- `output/analysis_result.json`: hasil terstruktur untuk dokumentasi/pengujian.

Untuk akses web, file `index.html` di root repo adalah dashboard statis yang siap dibuka langsung lewat Vercel. Jadi dosen bisa memilih: menjalankan analisis di terminal, atau memakai versi web di deployment.

## Cara Menjalankan Test

```bash
python3 -m unittest discover -s tests
```

## Catatan Biologi dan Komputasi

STR adalah motif DNA pendek yang berulang secara tandem. Dalam model edukatif ini, profil seseorang direpresentasikan sebagai vektor jumlah repeat maksimum untuk beberapa marker, misalnya `AGAT=7` dan `AATG=4`.

Sequence alignment ditambahkan sebagai tahap validasi karena sampel forensik sering berupa fragmen. Algoritma yang digunakan adalah Smith-Waterman local alignment dengan skor default `match=2`, `mismatch=-1`, dan `gap=-2`. Ranking akhir menggunakan kombinasi 70% skor kecocokan STR dan 30% skor dukungan alignment lokal. Skor dukungan alignment dihitung dari identity alignment yang dibobot dengan coverage sekuens TKP supaya kecocokan pendek tidak terlihat terlalu kuat.

Format CSV tersangka:

```csv
Nama,AGAT,AATG,TATC,GCTA,TCTA,DNA_Sequence
Bagas Pratama,7,4,6,3,5,GATTACAGTCCAGAT...
```
