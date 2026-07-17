# Sample data

This folder holds the per-animal sample datasets used by
`../AnalysisTutorial/RunTemplate.m`, using the **post-Run2-split** data (see
`../DataStructure.md`). The sample is **swappable**: data is packaged per animal
under `SampleData/<Animal>/`, and each analysis is pointed at one by setting its
`animal_*` variable (e.g. `animal_beh`, `animal_thetaseq`, …) near the top of
`RunTemplate.m`.

All five ketamine-day-1 animals are packaged:

| Animal | pxpercm | Effective sessions |
|--------|---------|--------------------|
| URat1  | 3.82 | 3 (Run1, Run2Bi, Run3) |
| URat2  | 3.85 | 4 (Run1, Run2Bi, Run3, Run2woBi) |
| URat3  | 3.45 | 3 (Run1, Run2Bi, Run3) |
| URat4  | 3.75 | 4 (Run1, Run2Bi, Run3, Run2woBi) |
| XRat1  | 3.75 | 4 (Run1, Run2Bi, Run3, Run2woBi) |

For the 3-session animals the data still stores 4 session slots; slot 4
(Run2woBi) is empty and is skipped automatically.

## Files (per animal)

| File | Variables | Used by |
|------|-----------|---------|
| `<Animal>/Behavior.mat` | `mazel` (nsess×1 cell), `ses` (1×nsess struct), `info` (with `pxpercm`) | Analysis 1 (behavior) |
| `<Animal>/ThetaSeq.mat` | `M1`,`M2`,`PlFields`,`PlFmesh`,`ClF`, `thetapk` (nsess×1 cell of theta-peak times) | Analysis 2 (theta sequence); also loads `ses`,`mazel` from `Behavior.mat` |
| `<Animal>/WakingRep.mat` | `fratime` (nsess×1 cell of awake-rest frame times), `frainfo` (per-frame LFP features), `times` (nsess×2 cell of directional lap times) | Analysis 3 (waking maze representation); also loads `ses`,`mazel` from `Behavior.mat` and `M1`,`M2`,`PlFields`,`PlFmesh`,`ClF` from `ThetaSeq.mat` |
| (no new file) | — | Analysis 4 (spike mutual information) reuses `ThetaSeq.mat` (full `ClF`, `M2`) and `Behavior.mat` (`ses`,`mazel`) |
| (no new file) | — | Analysis 5 (spike dimensionality) reuses `ThetaSeq.mat` (`M2`,`ClF`), `Behavior.mat` (`ses`,`mazel`), and `WakingRep.mat` (`times`) |
| `<Animal>/RunLFP.mat` | `seglfp` (1×nsess cell; each a 1×nseg cell of best-tetrode LFP traces from high-velocity run segments), `fs` (sampling rate, Hz) | Analysis 6 (LFP bispectrum / theta wave asymmetry) |
| `<Animal>/LDS_REM.mat` | `behspk` (1×nSleepSession struct; `[0,s]['rem']` = per-REM-epoch `spkdata`/`shfspkdata`, time-binned smoothed spikes) | Analysis 7 (Python LDS — REM model fit) |
| `<Animal>/LDS_Run.mat` | `behspk` (nRunSession×nDir struct; `[s,d]['lap']` = per-lap `spkdata`, time-binned smoothed spikes) | Analysis 7 (Python LDS — REM→Run prediction) |

(Additional files are added as each analysis is implemented.)

## How the sample is built (swappable per animal)

Each file is extracted from the lab's preprocessed Run2-split outputs for the
chosen animal. For `Behavior.mat`:
- `mazel`, `ses` ← `KTM_Run/<Dataset>/Mazel_KetImpactRun2/Twin300Step120/{Maze.mat, ses.mat}`
- `info.pxpercm` ← the dataset entry in `ratCA_database_str.m` (URat4 = 3.75)

For `RunLFP.mat`: the best-tetrode run LFP (`Ptet`) lives at
`KTM_RunLFP/<Dataset>/BestTT_KetImpactRun2/Twin300Step120/Ptet.mat`; its time base
aligns 1:1 with `mazel`. Per session, the high-velocity run epochs (|vel| ≥ 10 cm/s,
≥ 1 s, capped at ~150 s total) are cut from the LFP and stored in `seglfp`.

To regenerate for a different animal, repackage that animal's corresponding
files under the same variable names.

To request the sample data, contact the Dragoi lab: george.dragoi@yale.edu
