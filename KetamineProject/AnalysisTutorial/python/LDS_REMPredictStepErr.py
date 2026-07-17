#!/usr/bin/env python
# coding: utf-8
"""
LDS_REMPredictStepErr.py
------------------------
Analysis 7 (manuscript Fig. 6), part 2 of 2.

Take the LDS fitted to each REM-sleep session by ``LDS_FitEachREM.py`` and test
how well its learned dynamics predict the awake *Run* (maze) population activity,
one time step ahead. For every run lap we:
  (1) invert the model's emission matrix C to read the latent state from the
      observed population vector at time t,
  (2) push it one step through the REM dynamics and emit a predicted population
      vector at t+1,
  (3) score the RMS error against the true t+1 vector, and against a
      time-shuffled control.
A prediction error below the shuffle means REM and Run share dynamical
structure ("REM intrusion"); ketamine (Run2Bi) is expected to change this.

Inputs (per animal):
    ./output/<Animal>/ses<k>_lds.pkl          (from LDS_FitEachREM.py)
    ../../SampleData/<Animal>/LDS_Run.mat      run-lap smoothed spike matrices
        behspk : nsess x ndir struct array; behspk[ses, idir]['lap'] is a
                 1 x nlap struct array, each with spkdata = ncells x ntimebins.

Output (./output/<Animal>/):
    lap_step_err.csv      per-lap RMS error (data and shuffle)
    REM_predict_Run.png   violin plot of RMS error per run session

Requires the `ssm` package and the models produced by LDS_FitEachREM.py.

Yuchen Zhou, yuchen.zhou@yale.edu, yuchenzhou93@gmail.com
"""

import os
import pickle

import autograd.numpy as np
import autograd.numpy.random as npr
npr.seed(52)

import matplotlib
matplotlib.use("Agg")              # headless; comment out to show figures
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.io import loadmat

sns.set_style("white")
sns.set_context("talk")


# ----------------------------------------------------------------------------- #
# helpers (inlined so the tutorial is self-contained)
# ----------------------------------------------------------------------------- #
def ket_run2_split_pre(runses):
    """Map a (Run2-split) run-session index to its pre-sleep REM session index
    (0-based). Run1/Run2Bi/Run3 -> their own pre-sleep; the split-off Run2woBi
    (index 3) shares Run2Bi's pre-sleep (index 1)."""
    if runses <= 2:
        return runses
    elif runses == 3:
        return 1
    return runses


def rms_step_pred_err(lds, obs):
    """One-step-ahead RMS prediction error of `obs` (T x ncells) under `lds`.
    Returns (mean RMS data, mean RMS time-shuffled control)."""
    C = lds.emissions.Cs[0]
    nobs = obs.shape[0]
    ypredict = np.copy(obs)
    for it in range(nobs - 1):
        y_t = obs[it, :]
        x_t = np.linalg.pinv(C) @ y_t                 # read latent state (no noise)
        prefix = (np.zeros(1, dtype=int), x_t[np.newaxis, :], y_t[np.newaxis, :])
        _, f_emissions = lds.sample(1, prefix=prefix)  # one step forward
        ypredict[it + 1, :] = f_emissions

    ypredict = ypredict[1:, :]
    obs = obs[1:, :]
    shfobs = obs[np.random.permutation(obs.shape[0])]
    rms = np.sqrt(np.mean((obs - ypredict) ** 2, axis=1))
    rms_shf = np.sqrt(np.mean((shfobs - ypredict) ** 2, axis=1))
    return np.mean(rms), np.mean(rms_shf)


# ----------------------------------------------------------------------------- #
def predict_animal(animal, sample_dir, out_dir):
    print(f"Processing rat {animal}")
    out = os.path.join(out_dir, animal)
    mat = loadmat(os.path.join(sample_dir, animal, "LDS_Run.mat"))
    behspk = mat["behspk"]
    nsess, ndir = behspk.shape

    rows = []
    for ses in range(nsess):
        slpses = ket_run2_split_pre(ses)
        mpath = f"{out}/ses{slpses+1}_lds.pkl"
        if not os.path.exists(mpath):
            print(f"  session {ses+1}: missing model {mpath}; run LDS_FitEachREM.py first")
            continue
        with open(mpath, "rb") as f:
            lds = pickle.load(f)
        lds.D = int(lds.D)

        nlap_tot = 0
        for idir in range(ndir):
            laps = behspk[ses, idir]["lap"]
            for il in range(laps.shape[1]):
                obs = np.transpose(laps[0, il]["spkdata"])     # (T, ncells)
                if obs.shape[0] < 3:
                    continue
                rmse, rmse_shf = rms_step_pred_err(lds, obs)
                rows.append({"ses": ses + 1, "dir": idir + 1,
                             "rmse": rmse, "rmse_shf": rmse_shf})
                nlap_tot += 1
        print(f"  run session {ses+1}/{nsess} (REM model ses{slpses+1}): {nlap_tot} laps")

    laperr = pd.DataFrame(rows)
    laperr.to_csv(f"{out}/lap_step_err.csv", index=False)

    # violin plot: RMS one-step error, data vs shuffle, per run session
    long = pd.melt(laperr, id_vars=["ses"], value_vars=["rmse", "rmse_shf"],
                   var_name="error_type", value_name="error")
    long["group"] = long["error_type"].map({"rmse": "data", "rmse_shf": "shuffle"})
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.violinplot(x="ses", y="error", hue="group", data=long, ax=ax,
                   inner="box", palette="muted")
    ax.set(title=f"{animal}: REM model -> Run one-step prediction",
           xlabel="Run session", ylabel="RMS step error")
    plt.tight_layout()
    fig.savefig(f"{out}/REM_predict_Run.png", dpi=110)
    plt.close(fig)
    print(f"  saved lap_step_err.csv + REM_predict_Run.png to {out}/")


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    SAMPLE_DIR = os.path.join(here, "..", "..", "SampleData")
    OUT_DIR = os.path.join(here, "output")
    ANIMAL = "URat1"      # must match the animal fitted by LDS_FitEachREM.py
    predict_animal(ANIMAL, SAMPLE_DIR, OUT_DIR)
