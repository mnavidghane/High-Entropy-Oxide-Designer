# High-Entropy Oxide Single-Phase Designer

An interactive, rule-based **design and screening toolkit for single-phase
high-entropy oxides (HEOs)** across three structure families:
**Rock Salt (AO)**, **Spinel (AB₂O₄)**, and **Pyrochlore (A₂B₂O₇)**.

You pick the cation palette for each crystallographic site; the tool performs an
**exhaustive composition scan** (5 at.% grid) and keeps only the compositions that
satisfy published single-phase formation windows based on:

- **site-resolved configurational entropy** (ΔS_mix^N),
- **sublattice size-mismatch descriptors** (δr^N and δr*),
- the **cation radius ratio r_A/r_B** (pyrochlore stability).

No dependencies — pure Python 3.8+ standard library.

---

## 1. Executive Summary

| Item | Description |
|---|---|
| Language | Python 3.8+ (standard library only: `math`, `itertools`) |
| Interface | Interactive console (step-by-step prompts) |
| Structures | Rock Salt (AO), Spinel (AB₂O₄), Pyrochlore (A₂B₂O₇) |
| Scan | Exhaustive 5 at.% grid per site, strictly positive fractions, permutation-deduped |
| Criteria | ΔS_mix^N windows, δr^N / δr* limits, r_A/r_B > 1.46 (pyrochlore) |
| Database | Shannon (VI-coordinate) ionic radii embedded |
| Outputs | Per-element valid concentration ranges + a representative valid composition |

## 2. Structure families and site palettes

| Structure | Formula | A-site palette | B-site palette |
|---|---|---|---|
| Rock Salt | AO | Mg, Ni, Co, Zn, Cu, La, Ce, Pr, Nd, Sm | — |
| Spinel | AB₂O₄ | Mg, Ni, Co, Zn, Cu | Al, Cr, Fe, Ga, Ti |
| Pyrochlore | A₂B₂O₇ | Y, Gd, La, Nd, Sm, Ce, Hf | Fe, Cr, Ti, Zr, Sn, Nb |

Elements are entered per site and validated against the palette (≥ 2 elements per
occupied site, since all fractions must be strictly positive).

## 3. Single-phase formation criteria (as implemented)

| Structure | ΔS_mix^N (J·mol⁻¹·K⁻¹) | δr^N (%) | δr* (%) | r_A / r_B |
|---|---|---|---|---|
| Rock Salt | 6.5 – 8.5 | 2 – 8 | < 8 | — |
| Spinel | 4.5 – 6.75 | 0 – 12.5 | < 8 | — |
| Pyrochlore | 1.75 – 4.0 | 0 – 18 | < 5 | > 1.46 |

All three criteria (entropy, mismatch, ratio) must hold simultaneously; the
single-phase upper bound δr* < 8 % applies to every structure.

## 4. Descriptors

**Configurational entropy** — site-weighted and normalized by the atoms in the
formula unit:

```
ΔS_mix^N = −R · Σ_sites  w_site · Σ_i x_i ln x_i
w_site = C_site / (C_A + C_B + C_O)

Rock Salt:   w_A = 1/2
Spinel:      w_A = 1/7,  w_B = 2/7
Pyrochlore:  w_A = 2/11, w_B = 2/11
```

**Size mismatch:**

```
r̄_site   = Σ_i x_i r_i
δr_site  = 100 · sqrt( Σ_i x_i (1 − r_i / r̄)² )        (per sublattice)
δr*      = sqrt( δr_A² + δr_B² )
δr^N     = sqrt( w_A^Ω · δr_A² + w_B^Ω · δr_B² )

per-structure Ω-weights (n·Ω/Ω_cell):
Rock Salt:  w_A^Ω = 1.0
Spinel:     w_A^Ω = 1/24,      w_B^Ω = 1/2
Pyrochlore: w_A^Ω = 1/4,       w_B^Ω = √2/12
```

## 5. How the scan works

1. Choose the structure (1/2/3).
2. Enter the A-site (and, if applicable, B-site) elements — validated against the palette.
3. The grid generator builds all unique 5 at.% compositions per site
   (combinations + permutations, deduped, strictly positive fractions).
4. Every A×B combination is evaluated against the full criteria set (§3).
5. Valid compositions are collected and reported as:
   - **per-element concentration ranges** (min–max at.% over all valid hits),
   - **one representative valid composition** with its full descriptor printout
     (ΔS_mix^N, δr^N, δr*, r_A/r_B).

## 6. Quick sanity checks

These follow directly from the entropy normalization and confirm that the
criteria windows are self-consistent with the grid:

| Test | Expected ΔS_mix^N | In window? |
|---|---|---|
| Rock Salt, 5 equimolar A-cations | (R/2)·ln 5 ≈ **6.69** | ✅ |
| Rock Salt window ↔ equimolar count | 5–7 A-cations (n=4 → 5.76 ✗, n=8 → 8.65 ✗) | ✅ |
| Spinel, 5 + 5 equimolar | (3R/7)·ln 5 ≈ **5.73** | ✅ (4–6 per site) |
| Pyrochlore, 3 + 3 equimolar | (4R/11)·ln 3 ≈ **3.32** | ✅ (2–3 per site) |

## 7. Usage

```bash
python heo_designer.py
```

Typical session: choose structure → enter site elements (e.g. for Rock Salt:
`Mg Ni Co Zn Cu`) → scan runs → per-element ranges and one representative valid
composition are printed. For Spinel/Pyrochlore the B-site is requested after the
A-site.

## 8. Ionic radius database

Shannon (VI-coordinate) effective ionic radii are embedded in the script —
divalent 3d cations for the rock-salt/spinel A-site (Mg 0.72, Ni 0.69, Co 0.74,
Zn 0.74, Cu 0.73 Å), trivalent rare-earths (La 1.032, Ce³⁺ 1.01, Pr 0.99,
Nd 0.983, Sm 0.958, Gd 0.938 Å), and the remaining cations (Y 0.90, Ca 1.00,
Sr 1.18, Ba 1.35, Pb 1.20 Å). High-spin Fe³⁺ is used where applicable.

## 9. Limitations & scope

1. **Geometric/empirical criteria only** — no energetics (no DFT/CALPHAD);
   the formation windows are correlations from published phase-formation maps.
2. **Charge neutrality and valence are not enforced.** In particular, the
   trivalent rare-earth options on the Rock Salt A-site (La, Ce, Pr, Nd, Sm) are
   charge-imbalanced for a strict AO stoichiometry — keep to the divalent cations
   (Mg, Ni, Co, Zn, Cu) unless mixed valence is intended.
3. **Fixed stoichiometry** per structure: no oxygen non-stoichiometry, no spinel
   site inversion, no cation vacancies.
4. The 5 at.% grid can miss narrow single-phase windows; reduce `step` in
   `generate_compositions` for a finer (slower) scan.
5. Room-temperature Shannon radii; coordination/spin-state assumptions are fixed
   in the database.

## 10. Program structure

| Component | Purpose |
|---|---|
| `STRUCTURES` | registry: formula, site counts (C_A, C_B, C_O), allowed palettes, δr^N weights |
| `IONIC_RADII` | Shannon VI ionic radii database |
| `generate_compositions` | 5% grid per site, strictly positive, permutation-deduped |
| `calc_site_params` | site-average radius r̄ and site mismatch δr |
| `calc_site_entropy` | Σ x ln x per site |
| `evaluate_composition` | full descriptor set + criteria verdict (pass/fail) |
| `main` | interactive input, scan loop, ranges + representative composition output |

## 11. References

1. R. D. Shannon, *Revised effective ionic radii and systematic studies of
   interatomic distances in halides and chalcogenides*, Acta Crystallographica
   A 32, 751 (1976).
2. C. M. Rost et al., *Entropy-stabilized oxides*, Nature Communications 6, 8485 (2015).
3. Single-phase formation windows (§3): per-structure empirical criteria compiled
   from the high-entropy-oxide phase-map literature.

## License

MIT — see `LICENSE`.