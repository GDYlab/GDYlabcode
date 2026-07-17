# Ketamine project — analysis code

Code accompanying the ketamine manuscript (hippocampal CA1 spatial coding under
ketamine). The recordings use a **rest–run–rest** protocol on a single linear
track, with three run sessions — **Run1** (pre-injection baseline), **Run2**
(ketamine), **Run3** (recovery) — where Run2 is split into its behavioral-impact
window (**Run2Bi**) and remainder (**Run2woBi**). See `DataStructure.md`.

The repository has two parts.

## `FigureScripts/`
The MATLAB scripts that generate the six main figures. They read aggregated,
cross-animal source data through the lab's internal data paths and are provided
for **reproducibility / archival** — not as runnable examples. See
`FigureScripts/README.md` for the figure-to-script mapping.

## `AnalysisTutorial/`
A self-contained tutorial that runs the project's **key analyses on a
single-animal sample dataset** (post-Run2-split). Each analysis is one
documented function called from the master script, illustrating the ketamine
effect on a single example animal.

| # | Analysis | Figure | Tutorial code | Sample file(s) |
|---|----------|--------|---------------|----------------|
| 1 | Run behavior: distance/min, velocity–track angle, max velocity, end-approach/depart velocity profile | 1 | `function_Beh_RunVelocityMeasures.m` | `Behavior.mat` |
| 2 | Theta sequence by running velocity | 2 | `function_ThetaSeq_ByVelocity.m` | `Behavior.mat`, `ThetaSeq.mat` |
| 3 | Waking maze representation via population-vector cosine-similarity probability | 3 | `function_WakingMazeRep_CosineProb.m` | `Behavior.mat`, `ThetaSeq.mat`, `WakingRep.mat` |
| 4 | Pairwise mutual information of spikes during run (functional connectivity) | 4 | `function_RunSpikeMutualInfo.m` | `Behavior.mat`, `ThetaSeq.mat` |
| 5 | Intrinsic spike dimensionality during run (TwoNN) | 4 | `function_Embedding_RunSpikeDim.m` | `Behavior.mat`, `ThetaSeq.mat`, `WakingRep.mat` |
| 6 | LFP bispectrum and theta wave asymmetry | 5 | `function_LFPSpec_RunBispectrum.m` | `RunLFP.mat` |
| 7 | Linear dynamical system: REM latent dynamics and its prediction of Run | 6 | `python/LDS_FitEachREM.py`, `python/LDS_REMPredictStepErr.py` | `LDS_REM.mat`, `LDS_Run.mat` |

### Running it (MATLAB, analyses 1–6)
1. Add the `Ketamine` folder **and its subfolders** to the MATLAB path.
2. Open `AnalysisTutorial/RunTemplate.m`. Toggle individual analyses with the
   `run*` flags near the top, and pick the example animal per analysis with the
   `animal_*` variables (any of `URat1, URat2, URat3, URat4, XRat1`).
3. Run it. Approximate run times on the sample (URat1): behavior ~4 s, theta
   sequence ~3 min (Bayesian decoding), waking rep ~6 s, mutual information
   ~25 s, dimensionality ~1 s, bispectrum ~40 s.

Developed and tested on **MATLAB R2017b**. Requires the Statistics and Machine
Learning Toolbox (`ecdf`, `pdist2`, `nanmean`, …) and the Signal Processing
Toolbox (`fir1` for analysis 1, `hanning` for analysis 6). All other helper
functions are bundled in `libs/` (no further external dependencies).

### Running it (Python, analysis 7)
The LDS analysis uses the [`ssm`](https://github.com/lindermanlab/ssm) toolbox
and is run as two standalone scripts; see `AnalysisTutorial/python/README.md`
for the environment setup and run order.

## Data
- `DataStructure.md` — format of every sample variable (fundamental and
  per-analysis auxiliary).
- `SampleData/` — the per-animal sample dataset. All five animals are packaged
  (`SampleData/<Animal>/`); the sample is swappable per analysis without editing
  the analysis code. `SampleData/README.md` lists the files. Distributed on
  request (it is empty in the repository, like the Detour project).
- `libs/` — bundled MATLAB helper functions used by the tutorial.

URat1/URat3 have 3 effective run sessions; URat2/URat4/XRat1 have 4 (the
Run2woBi split element). The tutorial adapts to either automatically.

To get access to the sample data, contact the Dragoi lab: george.dragoi@yale.edu
