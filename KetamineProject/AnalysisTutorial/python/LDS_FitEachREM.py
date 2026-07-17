#!/usr/bin/env python
# coding: utf-8
"""
LDS_FitEachREM.py
-----------------
Analysis 7 (manuscript Fig. 6), part 1 of 2.

Fit a linear dynamical system (LDS) to the hippocampal population activity
during each REM-sleep epoch of every sleep session. The LDS learns a
low-dimensional latent dynamics matrix A whose eigenvalues describe the
rotational / expanding structure of the REM trajectory. The companion script
``LDS_REMPredictStepErr.py`` then asks how well this REM-learned model predicts
the awake Run (maze) activity.

Input (per animal, shipped in ../../SampleData/<Animal>/LDS_REM.mat):
    behspk : 1 x nsess struct array. behspk[0, ses]['rem'] is a 1 x nepoch
             struct array; each epoch has
                spkdata    : ncells x ntimebins smoothed spike matrix
                shfspkdata : a time-shuffled control of the same shape
    (time bin = 0.2 s, Gaussian-smoothed with kstd = 2; see ../../DataStructure.md)

Output (written to ./output/<Animal>/):
    ses<k>_lds.pkl       : fitted LDS model (data); consumed by part 2
    ses<k>_lds_shf.pkl   : fitted LDS model for the shuffle control
    REM_LDS.png          : ELBO, 2-D latent dynamics, eigenvalues per session

Requires the `ssm` package (https://github.com/lindermanlab/ssm) plus autograd,
scikit-learn, scipy, seaborn. See README.md for the conda environment.

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
import seaborn as sns
from scipy.io import loadmat
from sklearn.decomposition import PCA

import ssm
from ssm.plots import plot_dynamics_2d, plot_eigenvalues

sns.set_style("white")
sns.set_context("talk")
colors = sns.xkcd_palette(["windows blue", "red", "amber", "faded green"])

# ----------------------------------------------------------------------------- #
# helpers (inlined so the tutorial is self-contained)
# ----------------------------------------------------------------------------- #
def estimate_latent_dimensionality(data, var=0.99):
    """Smallest #PCs explaining `var` of the variance of the pooled trials."""
    concatenated = np.vstack(data)                     # (sum ntbins, ncells)
    pca = PCA()
    pca.fit(concatenated)
    cum = np.cumsum(pca.explained_variance_ratio_)
    return int(np.searchsorted(cum, var) + 1), cum


def reduce_dynamics_with_svd(A):
    """Best rank-2 (2x2) approximation of the dynamics matrix A via SVD, plus
    the fraction of spectral energy retained by the first two singular values."""
    U, S, Vt = np.linalg.svd(A)
    A2 = U[:2, :2] @ np.diag(S[:2]) @ Vt[:2, :2]
    explained = np.sum(S[:2] ** 2) / np.sum(S ** 2)
    return A2, explained


# ----------------------------------------------------------------------------- #
def fit_animal(animal, sample_dir, out_dir, latent_pct=0.99, num_iters=30):
    print(f"Processing rat {animal}")
    mat = loadmat(os.path.join(sample_dir, animal, "LDS_REM.mat"))
    behspk = mat["behspk"]
    nsess = behspk.shape[1]

    out = os.path.join(out_dir, animal)
    os.makedirs(out, exist_ok=True)

    fig, ax = plt.subplots(nsess, 3, figsize=(16, 5 * nsess), squeeze=False)
    for ses in range(nsess):
        epochs = behspk[0, ses]["rem"]
        nepoch = epochs.shape[1]

        ylist, yshflist = [], []
        for ie in range(nepoch):
            ylist.append(np.transpose(epochs[0, ie]["spkdata"]))        # (T, ncells)
            yshflist.append(np.transpose(epochs[0, ie]["shfspkdata"]))
        D_obs = ylist[0].shape[1]

        latent_dim, _ = estimate_latent_dimensionality(ylist, latent_pct)
        latent_dim_shf, _ = estimate_latent_dimensionality(yshflist, latent_pct)
        print(f"  session {ses+1}/{nsess}: {nepoch} REM epochs, "
              f"D_obs={D_obs}, latent_dim={latent_dim} (shuffle {latent_dim_shf})")

        # fit a Gaussian-emission LDS to the REM epochs (and to the shuffle)
        lds = ssm.LDS(D_obs, latent_dim, dynamics="gaussian", emissions="gaussian")
        elbos, q = lds.fit(ylist, num_iters=num_iters, method="laplace_em",
                           verbose=0)
        lds_shf = ssm.LDS(D_obs, latent_dim_shf, dynamics="gaussian", emissions="gaussian")
        elbos_shf, q_shf = lds_shf.fit(yshflist, num_iters=num_iters,
                                       method="laplace_em", verbose=0)

        # Save the fitted models. The variational posteriors (q, q_shf) are not
        # saved: part 2 only needs the model, and the posteriors are large.
        with open(f"{out}/ses{ses+1}_lds.pkl", "wb") as f:      pickle.dump(lds, f)
        with open(f"{out}/ses{ses+1}_lds_shf.pkl", "wb") as f:  pickle.dump(lds_shf, f)

        # (1) fitting ELBO
        ax[ses, 0].plot(elbos[1:])
        ax[ses, 0].set(xlabel="Iteration", ylabel="ELBO",
                       title=f"ses{ses+1} fitting")

        # (2) 2-D projection of the latent dynamics (flow field)
        A2, varpct = reduce_dynamics_with_svd(lds.dynamics.A)
        plot_dynamics_2d(A2, np.zeros(2), npts=10, axis=ax[ses, 1], color=colors[0])
        ax[ses, 1].set(xlabel="$x_1$", ylabel="$x_2$",
                       title=f"latent dynamics (var {varpct:.2f})")

        # (3) eigenvalues of A in the complex plane (rotation / decay structure)
        plot_eigenvalues(lds.dynamics.A, axis=ax[ses, 2])
        ax[ses, 2].set(xlabel="Real", ylabel="Imag", title="eigenvalues")

    fig.suptitle(f"{animal}: REM LDS", fontsize=14)
    plt.tight_layout()
    fig.savefig(os.path.join(out, "REM_LDS.png"), dpi=110)
    plt.close(fig)
    print(f"  saved models + REM_LDS.png to {out}/")


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    SAMPLE_DIR = os.path.join(here, "..", "..", "SampleData")
    OUT_DIR = os.path.join(here, "output")
    ANIMAL = "URat1"      # any of URat1, URat2, URat3, URat4, XRat1
    fit_animal(ANIMAL, SAMPLE_DIR, OUT_DIR)
