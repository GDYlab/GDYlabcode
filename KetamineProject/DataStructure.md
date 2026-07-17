# Data structure — Ketamine analysis tutorial

This document describes the format of the sample dataset variables
(`SampleData/`, packaged per animal). It lets you run the tutorial on the
provided sample, or format your own data to use the same analysis functions.

Default units: cm (length), s (time), Hz (rate), cm/s (velocity).

## The Run2 split (4 effective sessions)

In each ketamine recording day the animal runs a **rest–run–rest** protocol with
three run sessions: **Run1** (pre-injection baseline), **Run2** (ketamine
injected at the start of this run), and **Run3** (recovery). Within Run2 we detect
the time window of significant behavioral impact (**Run2Bi**) and split Run2 into
the impacted window and the remainder (**Run2woBi**). The sample data is the
**post-split** version.

The number of effective sessions and their index order depend on the animal:

**4 effective sessions** — URat2, URat4, XRat1:

| index | label     | meaning |
|-------|-----------|---------|
| 1     | Run1      | pre-injection baseline run |
| 2     | Run2Bi    | Run2, ketamine behavioral-impact window |
| 3     | Run3      | recovery run |
| 4     | Run2woBi  | Run2, remainder (without behavioral impact) |

**3 effective sessions** — URat1, URat3 (no Run2woBi split element):

| index | label     | meaning |
|-------|-----------|---------|
| 1     | Run1      | pre-injection baseline run |
| 2     | Run2Bi    | Run2, ketamine behavioral-impact window |
| 3     | Run3      | recovery run |

> ⚠️ CAUTION: in the 4-session case the order is `Run1, Run2Bi, Run3, Run2woBi`
> — index 3 is **Run3**, index 4 is **Run2woBi** (Run2woBi comes AFTER Run3, not
> before it). Do not assume the two Run2 parts are adjacent. Analysis code must
> resolve sessions by their meaning, not by assuming a fixed `[1 2 3 4]` =
> `[Run1 Run2Bi Run2woBi Run3]` layout.

All per-session variables below (`ses`, `mazel`, `ClF`, place fields, …) are
indexed by the session order above. The split is parameterized by a sliding
window (`Twin`/`Step`, e.g. `Twin300Step120` = 300 s window, 120 s step).
The default sample animal URat1 has 3 effective sessions (slot 4 left empty);
the per-analysis sample animal is chosen in `RunTemplate.m`.

The animal runs a **single linear track** in every session, so `ses(i).tra_p`
is always `1` \.

---

## `ses`  — track geometry per session
`1 x 4` struct array (one element per effective session). Fields used:
- `ses(i).tra_p` — active track id(s) in session i. Always `1` here (single linear track).
- `ses(i).tralim` — linear position limits `[start end]` of the track (cm).
- `ses(i).traxlim`, `ses(i).traylim` — x / y position limits of the track.
- `ses(i).tra_plfd` — place-field sequence per running direction:
  rows `[linear distance, cell id, subfield id]`, bad/interneuron units excluded.
- `ses(i).notes` — in-file documentation of the fields above.

## `mazel` — behavioral trajectory per session
`4 x 1` cell array. Each element is an `n x 8` matrix (n = behavior samples),
columns (see `mazel_note`):

| col | content |
|-----|---------|
| 1 | x position |
| 2 | y position |
| 3 | track number |
| 4 | linear position (cm) |
| 5 | run direction (instantaneous velocity) |
| 6 | timestamp (s) |
| 7 | velocity (cm/s) |
| 8 | lap-based run direction (manually checked) |

Timestamps may repeat (track-definition overlap); call
`mazel = uniquemazeltime(mazel)` before interpolation, as in the Detour project.

## `ClF` — spikes per cell per session
`nc x 4` cell array (nc = number of cells). Each element is an `n x 5` matrix,
one row per spike, columns: `[spike time, linear position, velocity, lap-based
direction, velocity-based direction]`. Columns 2–5 come from interpolating spike
times into `mazel`. (Saved as MATLAB v7.3 / HDF5.)

## `M1`, `M2`, `PlFields`, `PlFmesh` — place fields & cell info
Same layout/conventions as the Detour project (shared codebase):
- `PlFields` — tuning curves; dims `[cell, session, direction, spatial bin]`;
  `PlFmesh` gives the spatial-bin centers for the last dim.
- `M1` — place-field parameters; dims `[cell, parameter, subfield, session]`.
- `M2` — per-cell parameters; dims `[cell, parameter, session]`; used to exclude
  interneurons via `[M1,M2,ClF,PlFields,nc,pyrid,clu2pyr,pyrintid] = CA_ExcludeInt(M1,M2,ClF,PlFields)`.

> NOTE: only the cell-type column of `M2` is used directly by the tutorial
> (`M2(:,21,1)`, 0 = interneuron); the rest of `M1`/`M2` is consumed inside
> `CA_ExcludeInt`. The full column meanings follow the Detour `DataStructure`
> reference (shared codebase).

## `ttclupyr` — tetrode / cluster / pyramidal mapping
- `ttclu` — `[tetrode, cluster]` for all sorted units.
- `ttpyr` — `[tetrode, cluster]` for pyramidal units only.

---

## Per-analysis auxiliary variables

Beyond the fundamental variables above, some analyses need extra inputs (theta
peaks, lap times, awake-rest frames, LFP segments, time-binned spikes). These are
packaged per animal in their own `.mat` files under `SampleData/<Animal>/`, and
each variable is documented in the header of the analysis function that consumes
it. The table below gives the full inventory; `SampleData/README.md` lists which
file each analysis loads.

| File | Variables | Description | Used by |
|------|-----------|-------------|---------|
| `Behavior.mat` | `mazel`, `ses`, `info` | fundamental behavior / geometry (above); `info.pxpercm` converts pixels→cm | 1–5 |
| `ThetaSeq.mat` | `M1`, `M2`, `PlFields`, `PlFmesh`, `ClF`, `thetapk` | place fields + spikes (above), plus `thetapk` = `nsess×1` cell of theta-peak times (s) per session | 2–5 |
| `WakingRep.mat` | `fratime`, `frainfo`, `times` | `fratime` = `nsess×1` cell of awake-rest frame `[start end]` times (s); `frainfo` = per-frame LFP/state features; `times` = `nsess×2` cell of directional run-lap `[start end]` times | 3, 5 |
| `RunLFP.mat` | `seglfp`, `fs` | `seglfp` = `1×nsess` cell; `seglfp{is}` is a `1×nseg` cell of best-tetrode LFP traces (row vectors), one per high-velocity run segment of session `is`. `fs` = LFP sampling rate (Hz) | 6 |
| `LDS_REM.mat` | `behspk` | `1×nSleepSession` struct; `behspk(1,s).rem` is a `1×nEpoch` struct, each REM epoch with `spkdata` (`nCells×nTimeBin`) and `shfspkdata` (time-shuffled control). Time-binned (0.2 s), Gaussian-smoothed spike matrices | 7 |
| `LDS_Run.mat` | `behspk` | `nRunSession×nDir` struct; `behspk(s,d).lap` is a `1×nLap` struct, each lap with `spkdata` (`nCells×nTimeBin`). Same binning/smoothing as `LDS_REM.mat` | 7 |

The LDS analysis (`python/`) consumes the two `LDS_*.mat` files above; their use
is also described in `python/README.md`.

To get the sample data, contact the Dragoi lab: george.dragoi@yale.edu
