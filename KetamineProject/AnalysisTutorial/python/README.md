# Analysis 7 — Latent dynamics of REM sleep and its prediction of Run (Fig. 6)

This is the only analysis of the tutorial written in **Python** rather than
MATLAB, because it relies on the [`ssm`](https://github.com/lindermanlab/ssm)
linear-dynamical-system toolbox.

The analysis fits a **linear dynamical system (LDS)** to the hippocampal
population activity during each REM-sleep session, then asks how well that
REM-learned dynamics predicts the awake **Run** (maze) activity one time step
ahead. A run prediction error below a time-shuffled control indicates that REM
and Run share latent dynamical structure ("REM intrusion"); ketamine (Run2Bi)
changes this relationship.

## Scripts (run in order)

1. **`LDS_FitEachREM.py`** — fits one Gaussian LDS per REM session (latent
   dimensionality chosen by PCA at 99% variance) and saves the model. Produces
   `REM_LDS.png` (ELBO, 2-D latent flow field, eigenvalues of the dynamics
   matrix per session).
2. **`LDS_REMPredictStepErr.py`** — loads those models and computes the
   one-step-ahead RMS prediction error on every Run lap, against a shuffle.
   Produces `REM_predict_Run.png` (violin plot per session) and
   `lap_step_err.csv`.

Both default to `ANIMAL = "URat1"`; edit that variable (top of each
`__main__`) to use any of `URat1, URat2, URat3, URat4, XRat1`. Outputs are
written to `./output/<Animal>/`.

```bash
python LDS_FitEachREM.py          # fit REM models first
python LDS_REMPredictStepErr.py   # then predict Run from them
```

## Environment

The scripts need `ssm` plus `autograd`, `scikit-learn`, `scipy`, `seaborn`,
`pandas`, `matplotlib`. A minimal conda setup:

```bash
conda create -n ssm python=3.8
conda activate ssm
pip install autograd scikit-learn scipy seaborn pandas matplotlib
pip install git+https://github.com/lindermanlab/ssm.git
```

## Input data format

Shipped per animal under `../../SampleData/<Animal>/` (request from the Dragoi
lab, george.dragoi@yale.edu). Both are MATLAB `.mat` files holding a `behspk`
struct array of **time-binned, Gaussian-smoothed** spike-count matrices
(time bin 0.2 s, smoothing kernel std 2 bins).

| File | Structure |
|------|-----------|
| `LDS_REM.mat` | `behspk` is `1 × nSleepSession`. `behspk[0, s]['rem']` is `1 × nEpoch`; each REM epoch has `spkdata` (`nCells × nTimeBins`) and `shfspkdata` (time-shuffled control). |
| `LDS_Run.mat` | `behspk` is `nRunSession × nDirection`. `behspk[s, d]['lap']` is `1 × nLap`; each run lap has `spkdata` (`nCells × nTimeBins`). |

The Run sessions follow the Run2-split convention (see `../../DataStructure.md`);
`ket_run2_split_pre()` inside the prediction script maps each run session to the
pre-sleep REM session whose model should predict it (the split-off Run2woBi
reuses the Run2Bi pre-sleep model).
