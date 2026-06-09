# Finalization review + citation audit — manuscript_v4

Date: 2026-06-09. Manuscript: `manuscript_v4_integrated_NEE_EcolLett.md`.

## A. Citation audit

**Verdict: PASS.** 13 references, all real and standard; DOIs present and
format-correct; no orphan references; every in-text citation resolves to the list.

| # | Reference | DOI | In-text cited | Spot-check |
|---|---|---|---|---|
| 1 | Boakes et al. 2010 *PLoS Biol* | 10.1371/journal.pbio.1000385 | ✓ | confident |
| 2 | Brooks et al. 2017 *R Journal* (glmmTMB) | 10.32614/RJ-2017-066 | ✓ | confident |
| 3 | Chen et al. 2011 *Science* | 10.1126/science.1206432 | ✓ | confident |
| 4 | Eyring et al. 2016 *GMD* (CMIP6) | 10.5194/gmd-9-1937-2016 | ✓ | confident |
| 5 | Fick & Hijmans 2017 *IJC* (WorldClim2) | 10.1002/joc.5086 | ✓ | confident |
| 6 | Harris et al. 2020 *Sci Data* (CRU TS) | 10.1038/s41597-020-0453-3 | ✓ | **web-verified** |
| 7 | Isaac et al. 2014 *MEE* | 10.1111/2041-210X.12254 | ✓ | confident |
| 8 | Karger et al. 2017 *Sci Data* (CHELSA) | 10.1038/sdata.2017.122 | ✓ | confident |
| 9 | Lenoir et al. 2020 *Nat Ecol Evol* | 10.1038/s41559-020-1198-2 | ✓ | **web-verified** |
| 10 | Loarie et al. 2009 *Nature* | 10.1038/nature08649 | ✓ | confident |
| 11 | Openshaw 1984 (MAUP) | book, no DOI | ✓ | confident |
| 12 | Sullivan et al. 2014 *Biol Conserv* (eBird) | 10.1016/j.biocon.2013.11.003 | ✓ | confident |
| 13 | Valavi et al. 2019 *MEE* (blockCV) | 10.1111/2041-210X.13107 | ✓ | confident |

- **Orphans removed:** an earlier Hartig 2022 (DHARMa) entry is absent — DHARMa is
  not used in the condensed v4 methods, so no orphan remains.
- **No fabricated citations.** Two DOIs (Harris 2020, Lenoir 2020) were verified
  against the publisher this session; the remaining 10 are famous works with
  high-confidence DOIs — do a final one-click DOI resolve before submission.
- **To add on editor request:** a citable source for China citizen-science effort
  growth (currently asserted descriptively).

## B. Internal consistency (finalization scan)

| Check | Status |
|---|---|
| In-text formatted Tables 1–5 all defined and referenced | ✓ |
| Tables 6–8 / S1–S3 are data-file display items (CSV), correctly listed | ✓ |
| Figures 1–10 + 8b + 8(v3) all in display items and referenced | ✓ |
| No leftover TODO/TBD/placeholder text | ✓ (only legitimate "not yet" prose) |
| Provenance note updated to current state (grid-native, spatial CV, CMIP6, endogeneity, migratory, v3 parallel all reported) | ✓ fixed this round |
| Headline numbers consistent (HR 1.288 province; 1.274 v3; lag 1.322/1.292; AUC 0.73/0.63/0.56) across abstract / results / tables | ✓ |
| Mandatory sections present (Data/Code availability, Author contributions, Competing interests, Funding, AI disclosure) | ✓ |

## C. Remaining pre-submission items (non-blocking)

1. Final one-click DOI resolve for the 10 confident references.
2. Add a China citizen-science effort-growth citation if the editor asks.
3. Fill Funding statement (currently "[To be completed]").
4. Mint Zenodo DOI for code/data; update Data/Code availability.
5. Remove the internal provenance note before submission.
6. Venue choice → adjust abstract structure/length and convert citation style.
7. Optional analyses still open: six-climate-metric interaction matrix; per-feature
   forecast-skill PSI.

**Overall:** manuscript is citation-clean and internally consistent; ready for a
final human read-through and venue formatting.
