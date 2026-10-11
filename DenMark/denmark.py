"""
DenMark: Density-dependent Marked point process framework (Python version)
==========================================================================

User-facing functions to grid single-cell resolution spatial transcriptomics
(ST) data and fit the DenMark models with Stan through **cmdstanpy**.
This module mirrors ``DenMark-functions.R`` function by function, so results
from R and Python are directly comparable (same grid ordering, same Stan data,
same priors, same summaries).

Models (Stan files shipped alongside this module)
-------------------------------------------------
M2  ``DenMark-M2.stan``  full DenMark, two dependent latent fields::

        log lambda1(g) = log|g| + beta0 + a11 w1(g)
        log lambda2(g) = beta1 + log(Y(g) + 0.05) + a21 w1(g) + a22 w2(g)

    rho = Astar[2,1] is the density-expression correlation.
M1  ``DenMark-M1.stan``  simplified DenMark with independent fields (a21 = 0).

Quick start
-----------
>>> import denmark as dm
>>> res = dm.DenMark(df, gene="CDH1", model="M2", grid_size=10, basis=25)
>>> print(res)                     # rho, WAIC, diagnostics
>>> res.summary()                  # posterior summary of the key parameters
>>> dm.plot_fields(res)

Requirements: numpy, pandas, scipy, cmdstanpy (with CmdStan installed:
``python -m cmdstanpy.install_cmdstan``); matplotlib for the plots.

Authors: Mingchi Xu, Alexandra M. Schmidt, Qihuang Zhang
"""

from __future__ import annotations

import os
import time
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Union

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.special import logsumexp

try:
    from cmdstanpy import CmdStanModel
except ImportError:  # pragma: no cover - the module still imports for grid work
    CmdStanModel = None

__all__ = [
    "filter_genes", "buildgridpp", "grid_pts_centroid", "hsgp_indices",
    "prepare_data", "compile_model", "make_init", "DenMark", "denmark",
    "DenMarkmain", "DenMarkResult", "compare_models", "scan_genes",
    "waic", "matern32", "simulate", "save_result", "load_result", "plot_grid", "plot_fields",
    "plot_decay", "plot_scan",
]

_HERE = Path(__file__).resolve().parent
_MODEL_CACHE: Dict[str, "CmdStanModel"] = {}
STAN_DIR: Optional[Union[str, Path]] = None   # set dm.STAN_DIR to point elsewhere


# --------------------------------------------------------------------------- #
# 0. Small utilities                                                          #
# --------------------------------------------------------------------------- #

def matern32(d, sigma=1.0, ell=1.0):
    """Matern (nu = 3/2) covariance, same parameterisation as the Stan models:
    C(d) = sigma^2 (1 + sqrt(3) d / ell) exp(-sqrt(3) d / ell)."""
    u = np.sqrt(3.0) * np.asarray(d, float) / ell
    return sigma ** 2 * (1.0 + u) * np.exp(-u)


def _matern32_range_factor(threshold=0.05):
    """Distance / ell at which the Matern-3/2 correlation equals threshold
    (0.05 -> 2.7389)."""
    u = brentq(lambda u: (1 + u) * np.exp(-u) - threshold, 1e-8, 100, xtol=1e-12)
    return u / np.sqrt(3.0)


def _xyz(dataset: pd.DataFrame, gene=None, x="x", y="y") -> tuple:
    """Standardise the input to arrays x, y, gene counts."""
    df = pd.DataFrame(dataset)
    if not {x, y} <= set(df.columns):
        if df.shape[1] < 3:
            raise ValueError(f"`dataset` needs columns '{x}', '{y}' and one expression column.")
        warnings.warn(f"Columns '{x}'/'{y}' not found; using columns 1 and 2 as x and y.")
        df = df.rename(columns={df.columns[0]: x, df.columns[1]: y})
    if gene is None:
        others = [c for c in df.columns if c not in (x, y)]
        if not others:
            raise ValueError("No expression column found in `dataset`.")
        gene = others[0]
    if gene not in df.columns:
        raise ValueError(f"Gene '{gene}' is not a column of `dataset`.")
    sub = df[[x, y, gene]].astype(float)
    ok = sub.notna().all(axis=1).to_numpy()
    if not ok.all():
        warnings.warn(f"{(~ok).sum()} cell(s) with missing coordinates/expression removed.")
        sub = sub[ok]
    v = sub[gene].to_numpy()
    if (v < 0).any():
        raise ValueError("Expression values must be non-negative counts.")
    return sub[x].to_numpy(), sub[y].to_numpy(), v, gene


def _grid_dims(grid_size) -> tuple:
    gs = np.atleast_1d(np.round(np.asarray(grid_size))).astype(int)
    if gs.size == 1:
        gs = np.repeat(gs, 2)
    if gs.size != 2 or (gs < 2).any():
        raise ValueError("`grid_size` must be one integer >= 2 or (nx, ny).")
    return int(gs[0]), int(gs[1])


def _bbox(x, y, bbox=None) -> np.ndarray:
    if bbox is None:
        return np.array([x.min(), x.max(), y.min(), y.max()], float)
    bb = np.asarray(bbox, float)
    if bb.size != 4 or bb[1] <= bb[0] or bb[3] <= bb[2]:
        raise ValueError("`bbox` must be (xmin, xmax, ymin, ymax) with xmin < xmax, ymin < ymax.")
    return bb


# --------------------------------------------------------------------------- #
# 1. Data preparation                                                         #
# --------------------------------------------------------------------------- #

def filter_genes(expr: pd.DataFrame, var_quantile=0.05, min_nonzero=0.05) -> List[str]:
    """Quality-control filter: keep genes whose variance across cells is above
    the `var_quantile` quantile AND that are detected in >= `min_nonzero` of
    the cells (DenMark paper: 5% and 5%)."""
    mat = pd.DataFrame(expr).astype(float)
    var = mat.var(axis=0, ddof=1)
    keep_var = var >= np.nanquantile(var, var_quantile)
    keep_det = (mat > 0).mean(axis=0) >= min_nonzero
    return list(mat.columns[keep_var & keep_det])


def buildgridpp(grid_size, dataset: pd.DataFrame, gene=None, bbox=None) -> dict:
    """Discretise cell locations and expression on a regular grid.

    Each cell is assigned to exactly one grid cell. Vectors are ordered with x
    varying fastest (grid_index = i + j * nx, 0-based), identical to the R
    version and to ``Ymat.ravel(order="F")``.

    Returns a dict with keys
      loc_mark_mat (N x 2: point = Y, mark = M), Ymat, Mmat (nx x ny, rows = x),
      grid (DataFrame: grid_index, i, j, x, y [centroids], Y, M),
      bbox, step, grid_area, grid_size, gene.
    """
    x, y, v, gene = _xyz(dataset, gene)
    nx, ny = _grid_dims(grid_size)
    bb = _bbox(x, y, bbox)
    inside = (x >= bb[0]) & (x <= bb[1]) & (y >= bb[2]) & (y <= bb[3])
    x, y, v = x[inside], y[inside], v[inside]
    if x.size == 0:
        raise ValueError("No cells inside the bounding box.")
    dx, dy = (bb[1] - bb[0]) / nx, (bb[3] - bb[2]) / ny
    if dx <= 0 or dy <= 0:
        raise ValueError("All cells share the same x or y coordinate.")
    ix = np.clip(np.floor((x - bb[0]) / dx).astype(int), 0, nx - 1)
    iy = np.clip(np.floor((y - bb[2]) / dy).astype(int), 0, ny - 1)
    idx = ix + iy * nx
    Y = np.bincount(idx, minlength=nx * ny).astype(int)
    M = np.bincount(idx, weights=v, minlength=nx * ny)
    if np.any(np.abs(M - np.round(M)) > 1e-8):
        warnings.warn("Expression values are not integers; summed counts are rounded.")
    M = np.round(M).astype(int)
    gi = np.tile(np.arange(nx), ny)
    gj = np.repeat(np.arange(ny), nx)
    grid = pd.DataFrame({
        "grid_index": np.arange(nx * ny), "i": gi, "j": gj,
        "x": bb[0] + (gi + 0.5) * dx, "y": bb[2] + (gj + 0.5) * dy,
        "Y": Y, "M": M})
    return dict(loc_mark_mat=np.column_stack([Y, M]),
                Ymat=Y.reshape((nx, ny), order="F"),
                Mmat=M.reshape((nx, ny), order="F"),
                grid=grid, bbox=bb, step=(dx, dy), grid_area=dx * dy,
                grid_size=(nx, ny), gene=gene)


def grid_pts_centroid(grid_size, dataset: pd.DataFrame, bbox=None) -> dict:
    """Grid centroids (x varying fastest), grid area and number of grids."""
    df = pd.DataFrame(dataset)
    x = df["x"].to_numpy(float) if "x" in df else df.iloc[:, 0].to_numpy(float)
    y = df["y"].to_numpy(float) if "y" in df else df.iloc[:, 1].to_numpy(float)
    nx, ny = _grid_dims(grid_size)
    bb = _bbox(x, y, bbox)
    dx, dy = (bb[1] - bb[0]) / nx, (bb[3] - bb[2]) / ny
    gx = bb[0] + (np.arange(nx) + 0.5) * dx
    gy = bb[2] + (np.arange(ny) + 0.5) * dy
    pts = pd.DataFrame({"x": np.tile(gx, ny), "y": np.repeat(gy, nx)})
    return dict(grid_points=pts, grid_area=dx * dy, N=len(pts), step=(dx, dy), bbox=bb)


def hsgp_indices(basis) -> np.ndarray:
    """All basis-function index pairs (1..m1) x (1..m2), first index fastest."""
    b = np.atleast_1d(np.asarray(basis, int))
    if b.size == 1:
        b = np.repeat(b, 2)
    s1, s2 = np.meshgrid(np.arange(1, b[0] + 1), np.arange(1, b[1] + 1), indexing="xy")
    return np.column_stack([s1.ravel(), s2.ravel()]).astype(int)


def prepare_data(dataset: pd.DataFrame, gene=None, grid_size=10, basis=25,
                 boundfactor=1.5, bbox=None, L1=None, L2=None) -> dict:
    """Build the Stan data dict for DenMark.

    The ROI is moved to [-L1, L1] x [-L2, L2]; grid centroids are the HSGP
    coordinates and the HSGP domain is boundfactor * (L1, L2).
    Returns dict(stan_data, grid, info).
    """
    if boundfactor < 1:
        raise ValueError("`boundfactor` must be >= 1.")
    g = buildgridpp(grid_size, dataset, gene=gene, bbox=bbox)
    bb = g["bbox"]
    centre = np.array([(bb[0] + bb[1]) / 2, (bb[2] + bb[3]) / 2])
    half = np.array([(bb[1] - bb[0]) / 2, (bb[3] - bb[2]) / 2])
    for name, user, val in (("L1", L1, half[0]), ("L2", L2, half[1])):
        if user is not None and abs(user - val) > 0.01 * val:
            warnings.warn(f"{name} = {user:.4g} supplied but the data give {val:.4g}; "
                          "the data value is used.")
    grid = g["grid"].copy()
    coords = np.column_stack([grid["x"] - centre[0], grid["y"] - centre[1]])
    S = hsgp_indices(basis)
    dmax = float(np.hypot(np.ptp(coords[:, 0]), np.ptp(coords[:, 1])))
    stan_data = dict(
        n=int(coords.shape[0]),
        y1=grid["Y"].to_numpy(int),                # outcome 1: cell counts
        y2=grid["M"].to_numpy(int),                # outcome 2: summed expression
        d=2,
        mstar=int(S.shape[0]),
        coords=coords,
        log_grid_area=float(np.log(g["grid_area"])),
        mrange=dmax / (2 * 2.75),
        indices=S,
        L=(boundfactor * half).tolist(),
        is_centerted_PHI=0)
    grid["cx"], grid["cy"] = coords[:, 0], coords[:, 1]
    info = dict(gene=g["gene"], grid_size=g["grid_size"], step=g["step"],
                grid_area=g["grid_area"], bbox=bb, centre=centre,
                L1=float(half[0]), L2=float(half[1]), basis=basis,
                boundfactor=boundfactor, max_dist=dmax,
                n_cells=int(grid["Y"].sum()), Ymat=g["Ymat"], Mmat=g["Mmat"])
    return dict(stan_data=stan_data, grid=grid, info=info)


# --------------------------------------------------------------------------- #
# 2. Model compilation, initial values, fitting                               #
# --------------------------------------------------------------------------- #

def _stan_path(model: str, stan_file=None) -> Path:
    if stan_file is not None:
        p = Path(stan_file)
        if not p.exists():
            raise FileNotFoundError(f"Stan file not found: {p}")
        return p.resolve()
    fname = {"M2": "DenMark-M2.stan", "M1": "DenMark-M1.stan"}[model]
    for d in [STAN_DIR, _HERE, Path.cwd(), Path.cwd() / "DenMark"]:
        if d is not None and (Path(d) / fname).exists():
            return (Path(d) / fname).resolve()
    raise FileNotFoundError(f"Cannot find '{fname}'. Pass stan_file= or set denmark.STAN_DIR.")


def _detect_model(stan_file) -> str:
    return "M2" if "cholesky_factor_corr" in Path(stan_file).read_text() else "M1"


def compile_model(model="M2", stan_file=None, force_compile=False) -> "CmdStanModel":
    """Compile (once per session) and return a DenMark CmdStanModel."""
    if CmdStanModel is None:
        raise ImportError("cmdstanpy is required: pip install cmdstanpy")
    if model not in ("M1", "M2"):
        raise ValueError("model must be 'M2' or 'M1'")
    path = str(_stan_path(model, stan_file))
    if force_compile or path not in _MODEL_CACHE:
        _MODEL_CACHE[path] = CmdStanModel(stan_file=path,
                                          force_compile=force_compile)
    return _MODEL_CACHE[path]


def make_init(model, stan_data, rho0=0.3) -> dict:
    """Valid starting values. For M2 the Cholesky factor of the correlation
    matrix is lower triangular with unit-norm rows:
    Astar = [[1, 0], [rho0, sqrt(1 - rho0^2)]]."""
    if abs(rho0) >= 1:
        raise ValueError("`rho0` must lie in (-1, 1).")
    mstar = stan_data["mstar"]
    y1 = np.asarray(stan_data["y1"]); y2 = np.asarray(stan_data["y2"])
    beta0 = float(np.log(max(y1.mean(), 0.1)) - stan_data["log_grid_area"])
    beta1 = float(np.log(max(y2.sum(), 1) / (y1 + 0.05).sum()))
    ell = [float(stan_data["mrange"])] * 2
    base = dict(beta0=beta0, beta1=beta1, ell=ell,
                betab1=[0.0] * mstar, betab2=[0.0] * mstar)
    if model == "M2":
        base.update(Astar=[[1.0, 0.0], [rho0, float(np.sqrt(1 - rho0 ** 2))]],
                    sigma=[1.0, 0.5])
    else:
        base.update(a11=1.0, a22=0.5)
    return base


def waic(log_lik: np.ndarray) -> dict:
    """WAIC from a draws x observations log-likelihood matrix (same
    definitions as loo::waic)."""
    ll = np.asarray(log_lik, float)
    S = ll.shape[0]
    lppd_i = logsumexp(ll, axis=0) - np.log(S)
    p_i = ll.var(axis=0, ddof=1)
    elpd_i = lppd_i - p_i
    n = ll.shape[1]
    se = lambda v: float(np.sqrt(n * np.var(v, ddof=1)))
    return dict(elpd_waic=float(elpd_i.sum()), elpd_waic_se=se(elpd_i),
                p_waic=float(p_i.sum()), p_waic_se=se(p_i),
                waic=float(-2 * elpd_i.sum()), waic_se=2 * se(elpd_i),
                pointwise=elpd_i, n_high_p=int((p_i > 0.4).sum()))


_KEY_VARS = {
    "M2": ["beta0", "beta1", "ell[1]", "ell[2]", "Astar[2,1]", "sigma[1]", "sigma[2]"],
    "M1": ["beta0", "beta1", "ell[1]", "ell[2]", "a11", "a22"],
}
_LABELS = {
    "beta0": "intercept, cell density", "beta1": "intercept, expression",
    "ell[1]": "length-scale, density field", "ell[2]": "length-scale, expression field",
    "Astar[2,1]": "rho: density-expression correlation",
    "sigma[1]": "SD, density field", "sigma[2]": "SD, expression field",
    "a11": "SD, density field", "a22": "SD, expression field",
}


def _rhat_ess(x: np.ndarray) -> tuple:
    """Rank-normalised split R-hat and bulk/tail ESS (Vehtari et al. 2021)
    for one parameter; x has shape (draws, chains)."""
    from scipy.stats import norm, rankdata

    def split(a):
        h = a.shape[0] // 2
        return np.concatenate([a[:h], a[-h:]], axis=1) if a.shape[0] >= 4 else a

    def rhat_basic(a):
        n = a.shape[0]
        cm, cv = a.mean(0), a.var(0, ddof=1)
        W, B = cv.mean(), n * cm.var(ddof=1)
        return np.sqrt(((n - 1) / n * W + B / n) / W) if W > 0 else np.nan

    def ess_basic(a):
        n, m = a.shape
        if np.var(a) == 0:
            return np.nan
        acov = []
        for c in range(m):
            z = a[:, c] - a[:, c].mean()
            f = np.fft.rfft(z, 2 * n)
            ac = np.fft.irfft(f * np.conj(f))[:n] / n
            acov.append(ac)
        acov = np.array(acov)
        cv, cm = acov[:, 0] * n / (n - 1), a.mean(0)
        W = cv.mean()
        var_plus = W * (n - 1) / n + (cm.var(ddof=1) if m > 1 else 0)
        rho = 1 - (W - acov.mean(0)) / var_plus
        rho[0] = 1
        t, s = 0.0, 1
        prev = rho[0] + rho[1]
        tau = -1 + 2 * prev  # Geyer initial positive + monotone sequence
        s = 2
        while s + 1 < n:
            cur = rho[s] + rho[s + 1]
            if cur < 0:
                break
            cur = min(cur, prev)
            tau += 2 * cur
            prev = cur
            s += 2
        return n * m / max(tau, 1 / np.log10(n * m))

    def z(a):
        r = rankdata(a, method="average").reshape(a.shape)
        return norm.ppf((r - 3 / 8) / (a.size - 2 * 3 / 8 + 1))

    xs = split(x)
    rb = rhat_basic(z(xs))
    med = np.median(xs)
    rf = rhat_basic(z(np.abs(xs - med)))
    rhat = np.nanmax([rb, rf])
    ess_bulk = ess_basic(z(xs))
    q05, q95 = np.quantile(xs, [0.05, 0.95])
    ess_tail = np.nanmin([ess_basic((xs <= q05).astype(float)),
                          ess_basic((xs <= q95).astype(float))])
    return float(rhat), float(ess_bulk), float(ess_tail)


@dataclass
class DenMarkResult:
    """Result of a DenMark fit (Python counterpart of the R 'denmark' object)."""
    fit: object
    model: str
    gene: str
    stan_data: dict
    grid: pd.DataFrame
    info: dict
    settings: dict
    time: float = np.nan                       # minutes
    waic: Optional[dict] = None
    diagnostics: dict = field(default_factory=dict)
    hsgp: Optional[pd.DataFrame] = None
    _summary: Optional[pd.DataFrame] = None
    _key_draws: Optional[np.ndarray] = None    # (draws, chains, vars) when fit dropped

    # --- draws ------------------------------------------------------------- #
    def _key_array(self) -> np.ndarray:
        """(draws, chains, n_key_vars) array of the key scalar parameters."""
        if self._key_draws is not None:
            return self._key_draws
        arr = self.fit.draws(concat_chains=False)       # (draws, chains, columns)
        cols = list(self.fit.column_names)
        idx = [cols.index(v) for v in _KEY_VARS[self.model]]
        return arr[:, :, idx]

    def summary(self, probs=(0.025, 0.25, 0.5, 0.75, 0.975)) -> pd.DataFrame:
        """Posterior summary of the key parameters (mean, sd, quantiles,
        R-hat, bulk/tail ESS)."""
        if self._summary is not None and tuple(probs) == (0.025, 0.25, 0.5, 0.75, 0.975):
            return self._summary
        arr = self._key_array()
        rows = []
        for k, v in enumerate(_KEY_VARS[self.model]):
            x = arr[:, :, k]
            r, eb, et = _rhat_ess(x)
            row = dict(variable=v, mean=x.mean(), sd=x.std(ddof=1))
            for p in probs:
                row[f"q{p * 100:g}"] = np.quantile(x, p)
            row.update(rhat=r, ess_bulk=eb, ess_tail=et, label=_LABELS[v])
            rows.append(row)
        return pd.DataFrame(rows)

    def draws(self) -> pd.DataFrame:
        """Harmonised draws: beta0, beta1, ell1, ell2, rho, sigma1, sigma2 and
        a11 = sigma1, a21 = sigma2*rho, a22 = sigma2*sqrt(1 - rho^2).
        For M1: rho = a21 = 0, sigma1 = a11, sigma2 = a22."""
        arr = self._key_array()
        flat = arr.transpose(1, 0, 2).reshape(-1, arr.shape[2])   # chain-major, like R
        d = dict(zip(_KEY_VARS[self.model], flat.T))
        if self.model == "M2":
            rho, s1, s2 = d["Astar[2,1]"], d["sigma[1]"], d["sigma[2]"]
            return pd.DataFrame(dict(beta0=d["beta0"], beta1=d["beta1"],
                                     ell1=d["ell[1]"], ell2=d["ell[2]"], rho=rho,
                                     sigma1=s1, sigma2=s2, a11=s1, a21=s2 * rho,
                                     a22=s2 * np.sqrt(np.clip(1 - rho ** 2, 0, None))))
        a11, a22 = d["a11"], d["a22"]
        z = np.zeros_like(a11)
        return pd.DataFrame(dict(beta0=d["beta0"], beta1=d["beta1"],
                                 ell1=d["ell[1]"], ell2=d["ell[2]"], rho=z,
                                 sigma1=a11, sigma2=a22, a11=a11, a21=z, a22=a22))

    @property
    def rho(self) -> np.ndarray:
        """Posterior draws of rho (zeros for M1)."""
        return self.draws()["rho"].to_numpy()

    # --- derived quantities ------------------------------------------------ #
    def fields(self, probs: Optional[Sequence[float]] = None) -> pd.DataFrame:
        """Latent-field decomposition per grid cell (posterior means):
        eta1 = beta0 + a11 w1, a11w1, eta2 = beta1 + a21 w1 + a22 w2, a21w1,
        a22w2; plus log_lambda1 = log|g| + eta1, log_lambda2 = log(Y+0.05) + eta2."""
        if self.fit is None:
            raise RuntimeError("fields() needs the fit (use keep_fit=True).")
        W1 = self.fit.stan_variable("w1")          # (draws, n), chain-major
        W2 = self.fit.stan_variable("w2")
        d = self.draws()
        a11, a21, a22 = (d[c].to_numpy()[:, None] for c in ("a11", "a21", "a22"))
        comp = dict(a11w1=a11 * W1, a21w1=a21 * W1, a22w2=a22 * W2)
        comp["eta1"] = d["beta0"].to_numpy()[:, None] + comp["a11w1"]
        comp["eta2"] = d["beta1"].to_numpy()[:, None] + comp["a21w1"] + comp["a22w2"]
        out = self.grid[["grid_index", "i", "j", "x", "y", "Y", "M"]].copy()
        for nm in ("eta1", "a11w1", "eta2", "a21w1", "a22w2"):
            out[nm] = comp[nm].mean(axis=0)
            for p in (probs or []):
                out[f"{nm}_q{p * 100:g}"] = np.quantile(comp[nm], p, axis=0)
        out["log_lambda1"] = self.stan_data["log_grid_area"] + out["eta1"]
        out["log_lambda2"] = np.log(out["Y"] + 0.05) + out["eta2"]
        return out

    def decay(self, distances=None, probs=(0.025, 0.975), threshold=0.05) -> dict:
        """Spatial correlation decay curves and practical ranges.
        r1(d) = M(d; ell1);  r2(d) = [a21^2 M(d; ell1) + a22^2 M(d; ell2)]/(a21^2 + a22^2)."""
        d = self.draws()
        if distances is None:
            distances = np.linspace(0, self.info["max_dist"], 200)
        distances = np.asarray(distances, float)
        e1, e2 = d["ell1"].to_numpy(), d["ell2"].to_numpy()
        w21, w22 = d["a21"].to_numpy() ** 2, d["a22"].to_numpy() ** 2
        m1 = matern32(distances[None, :], 1, e1[:, None])
        m2 = matern32(distances[None, :], 1, e2[:, None])
        r2 = (w21[:, None] * m1 + w22[:, None] * m2) / (w21 + w22)[:, None]

        def summ(mat, name):
            return pd.DataFrame(dict(field=name, distance=distances, mean=mat.mean(0),
                                     lower=np.quantile(mat, probs[0], axis=0),
                                     upper=np.quantile(mat, probs[1], axis=0)))
        k = _matern32_range_factor(threshold)
        range1 = k * e1
        range2 = np.empty_like(e1)
        for s in range(e1.size):
            if w21[s] == 0:
                range2[s] = k * e2[s]
                continue
            f = lambda r: (w21[s] * matern32(r, 1, e1[s]) + w22[s] * matern32(r, 1, e2[s])) \
                / (w21[s] + w22[s]) - threshold
            range2[s] = brentq(f, 0, k * max(e1[s], e2[s]) * 1.01, xtol=1e-10)

        def rsum(v, name):
            return dict(field=name, mean=v.mean(), median=np.median(v),
                        lower=np.quantile(v, probs[0]), upper=np.quantile(v, probs[1]))
        return dict(curves=pd.concat([summ(m1, "density"), summ(r2, "expression")],
                                     ignore_index=True),
                    ranges=pd.DataFrame([rsum(range1, "density"), rsum(range2, "expression")]),
                    range_draws=pd.DataFrame(dict(density=range1, expression=range2)),
                    threshold=threshold)

    def check_hsgp(self, warn=True) -> pd.DataFrame:
        """Rule-of-thumb HSGP check (Riutort-Mayol et al. 2023, Matern-3/2):
        ell_min = 3.42 c S / m and c >= max(1.2, 4.5 ell / S), per axis.
        Conservative for estimating rho; treat as guidance."""
        d = self.draws()
        m = np.atleast_1d(self.info["basis"]).astype(float)
        m = np.repeat(m, 2) if m.size == 1 else m
        S = np.array([self.info["L1"], self.info["L2"]])
        c = self.info["boundfactor"]
        ell = {"density": np.median(d["ell1"]), "expression": np.median(d["ell2"])}
        rows = []
        for a, ax in enumerate(("x", "y")):
            for f, e in ell.items():
                emin = 3.42 * c * S[a] / m[a]
                cneed = max(1.2, 4.5 * e / S[a])
                rows.append(dict(field=f, axis=ax, ell_median=e, ell_min=emin,
                                 c_needed=cneed, basis_ok=e >= emin, boundary_ok=c >= cneed))
        tab = pd.DataFrame(rows)
        if warn and not (tab.basis_ok.all() and tab.boundary_ok.all()):
            msg = []
            if not tab.basis_ok.all():
                msg.append(f"posterior length-scale below {tab.ell_min[~tab.basis_ok].max():.3g}: "
                           "consider increasing `basis`")
            if not tab.boundary_ok.all():
                msg.append(f"consider boundfactor >= {tab.c_needed[~tab.boundary_ok].max():.2f}")
            print(f"DenMark HSGP check ({self.gene}, rule of thumb): " + "; ".join(msg))
        return tab

    # --- display ----------------------------------------------------------- #
    def __repr__(self) -> str:
        nx, ny = self.info["grid_size"]
        b = self.info["basis"]
        lines = [f"DenMark fit  |  model {self.model}  |  gene {self.gene}",
                 f"  grid {nx} x {ny} ({self.stan_data['n']} grids, {self.info['n_cells']} cells)"
                 f"  |  basis {b}  |  boundfactor {self.info['boundfactor']:.2f}"]
        s = self.summary()
        if self.model == "M2":
            r = s.loc[s.variable == "Astar[2,1]"].iloc[0]
            lines.append(f"  rho (density-expression correlation): {r['mean']:.3f}  "
                         f"[95% CrI {r['q2.5']:.3f}, {r['q97.5']:.3f}]")
        if self.waic is not None:
            lines.append(f"  WAIC: {self.waic['waic']:.2f} (SE {self.waic['waic_se']:.2f})")
        dg = self.diagnostics
        lines.append(f"  max R-hat {dg.get('max_rhat', np.nan):.3f} | min ESS bulk "
                     f"{dg.get('min_ess_bulk', np.nan):.0f} | divergences "
                     f"{dg.get('num_divergent', 'NA')} | time {self.time:.2f} min\n")
        cols = ["variable", "mean", "sd", "q2.5", "q50", "q97.5", "rhat", "ess_bulk"]
        lines.append(s[cols].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        return "\n".join(lines)

    # backward-compatible accessors (match the R list elements)
    @property
    def Ymat(self):
        return self.info["Ymat"]

    @property
    def Mmat(self):
        return self.info["Mmat"]


def save_result(res: DenMarkResult, path) -> str:
    """Pickle a DenMarkResult after loading all draws into memory, so it can
    be reloaded with load_result() after the CmdStan CSV files are gone."""
    import pickle
    if res.fit is not None:
        res.fit.draws()                       # forces the CSVs into memory
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as fh:
        pickle.dump(res, fh)
    return str(path)


def load_result(path) -> DenMarkResult:
    """Load a DenMarkResult written by save_result()."""
    import pickle
    with open(path, "rb") as fh:
        return pickle.load(fh)


def _diagnostics(fit, summ: pd.DataFrame, max_treedepth: int) -> dict:
    ndiv = ntd = np.nan
    try:
        ndiv = int(np.sum(fit.divergences))
        ntd = int(np.sum(fit.max_treedepths))
    except Exception:
        try:
            mv = fit.method_variables()
            ndiv = int(np.sum(mv["divergent__"]))
            ntd = int(np.sum(mv["treedepth__"] >= max_treedepth))
        except Exception:
            pass
    return dict(num_divergent=ndiv, num_max_treedepth=ntd,
                max_rhat=float(summ["rhat"].max()),
                min_ess_bulk=float(summ["ess_bulk"].min()),
                min_ess_tail=float(summ["ess_tail"].min()))


def DenMark(dataset: pd.DataFrame, gene=None, model="M2", grid_size=10, basis=25,
            boundfactor=1.5, bbox=None, chains=4, parallel_chains=None,
            iter_warmup=1000, iter_sampling=1000, adapt_delta=0.95,
            max_treedepth=10, seed=2024, refresh=400, inits=None, rho0=0.3,
            stan_file=None, output_dir=None, compute_waic=True, keep_fit=True,
            quiet=False, L1=None, L2=None, **sample_kwargs) -> DenMarkResult:
    """Fit DenMark to one gene.

    Parameters mirror the R function DenMark(): `dataset` has columns x, y and
    expression counts (one row per cell); `model` is "M2" (full, estimates
    rho) or "M1" (simplified); grid_size, basis may be scalars or (nx, ny);
    remaining arguments go to CmdStanModel.sample().
    """
    if stan_file is not None:
        model = _detect_model(stan_file)
    if model not in ("M1", "M2"):
        raise ValueError("model must be 'M2' or 'M1'")
    prep = prepare_data(dataset, gene=gene, grid_size=grid_size, basis=basis,
                        boundfactor=boundfactor, bbox=bbox, L1=L1, L2=L2)
    sd = prep["stan_data"]
    if sd["y2"].sum() == 0:
        raise ValueError(f"Gene '{prep['info']['gene']}' has zero counts in the ROI.")
    mod = compile_model(model, stan_file)
    if inits is None:
        inits = make_init(model, sd, rho0=rho0)
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)

    import logging
    log = logging.getLogger("cmdstanpy")
    old_level = log.level
    if quiet:
        log.setLevel(logging.ERROR)
    t0 = time.perf_counter()
    try:
        fit = mod.sample(data=sd, chains=chains,
                         parallel_chains=parallel_chains or chains,
                         iter_warmup=iter_warmup, iter_sampling=iter_sampling,
                         adapt_delta=adapt_delta, max_treedepth=max_treedepth,
                         seed=seed, inits=inits, refresh=refresh,
                         show_progress=not quiet, show_console=False,
                         output_dir=output_dir, **sample_kwargs)
    finally:
        log.setLevel(old_level)
    elapsed = (time.perf_counter() - t0) / 60

    res = DenMarkResult(fit=fit, model=model, gene=prep["info"]["gene"],
                        stan_data=sd, grid=prep["grid"], info=prep["info"],
                        settings=dict(chains=chains, iter_warmup=iter_warmup,
                                      iter_sampling=iter_sampling, adapt_delta=adapt_delta,
                                      max_treedepth=max_treedepth, seed=seed),
                        time=elapsed)
    res._summary = res.summary()
    res.diagnostics = _diagnostics(fit, res._summary, max_treedepth)
    if compute_waic:
        res.waic = waic(fit.stan_variable("log_lik"))
    res.hsgp = res.check_hsgp(warn=not quiet)
    if not keep_fit:
        res._key_draws = res._key_array()
        res.fit = None
    return res


denmark = DenMark   # snake_case alias


def DenMarkmain(grid_size, dataset, basis, boundfactor, L1=None, L2=None,
                denmark_version="M2", **kwargs) -> DenMarkResult:
    """Backward-compatible wrapper of the original interface; `denmark_version`
    is "M2"/"M1" or a path to a DenMark Stan file."""
    if denmark_version in ("M1", "M2"):
        return DenMark(dataset, model=denmark_version, grid_size=grid_size, basis=basis,
                       boundfactor=boundfactor, L1=L1, L2=L2, **kwargs)
    return DenMark(dataset, stan_file=denmark_version, grid_size=grid_size, basis=basis,
                   boundfactor=boundfactor, L1=L1, L2=L2, **kwargs)


# --------------------------------------------------------------------------- #
# 3. Model comparison and DCG discovery                                       #
# --------------------------------------------------------------------------- #

def compare_models(**fits: DenMarkResult) -> pd.DataFrame:
    """WAIC comparison, e.g. compare_models(M2=f2, M1=f1). Lower WAIC is better;
    elpd_diff / se_diff are relative to the best model (as loo::loo_compare)."""
    rows = []
    for name, f in fits.items():
        if f.waic is None:
            raise ValueError(f"{name}: fit without WAIC")
        w = f.waic
        rows.append(dict(model=name, waic=w["waic"], waic_se=w["waic_se"],
                         elpd_waic=w["elpd_waic"], p_waic=w["p_waic"]))
    tab = pd.DataFrame(rows).sort_values("waic").reset_index(drop=True)
    if len(tab) > 1:
        best = fits[tab.model[0]].waic["pointwise"]
        tab["elpd_diff"] = [tab.elpd_waic[k] - tab.elpd_waic[0] for k in range(len(tab))]
        tab["se_diff"] = [float(np.sqrt(len(best) * np.var(fits[m].waic["pointwise"] - best,
                                                          ddof=1))) for m in tab.model]
    return tab


def scan_genes(data: pd.DataFrame, genes: Optional[Iterable[str]] = None,
               fit_simple=False, save_dir=None, cred_level=0.95, quiet=True,
               verbose=True, **kwargs) -> pd.DataFrame:
    """Fit DenMark gene by gene and call density-correlated genes:
    PCG if the credible interval of rho is above 0, NCG if below 0, else
    "not DCG". With save_dir, one CSV row per gene is cached (resumable)."""
    data = pd.DataFrame(data)
    if genes is None:
        genes = [c for c in data.columns if c not in ("x", "y")]
    genes = [g for g in genes if g in data.columns]
    if not genes:
        raise ValueError("None of the requested genes are columns of `data`.")
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    a = (1 - cred_level) / 2
    rows = []
    for g in genes:
        fcache = Path(save_dir) / f"denmark_{g}.csv" if save_dir else None
        if fcache is not None and fcache.exists():
            rows.append(pd.read_csv(fcache).iloc[0].to_dict())
            continue
        try:
            f2 = DenMark(data[["x", "y", g]], gene=g, model="M2", keep_fit=False,
                         quiet=quiet, **kwargs)
            rho = f2.rho
            row = dict(gene=g, rho_mean=rho.mean(), rho_median=np.median(rho),
                       rho_lower=np.quantile(rho, a), rho_upper=np.quantile(rho, 1 - a),
                       prob_positive=(rho > 0).mean(), total_counts=int(f2.stan_data["y2"].sum()),
                       max_rhat=f2.diagnostics["max_rhat"],
                       min_ess_bulk=f2.diagnostics["min_ess_bulk"],
                       divergences=f2.diagnostics["num_divergent"],
                       waic_M2=f2.waic["waic"] if f2.waic else np.nan,
                       time_min=f2.time, error=None)
            if fit_simple:
                f1 = DenMark(data[["x", "y", g]], gene=g, model="M1", keep_fit=False,
                             quiet=quiet, **kwargs)
                row["waic_M1"] = f1.waic["waic"]
                row["delta_waic"] = row["waic_M1"] - row["waic_M2"]   # > 0 favours M2
                row["time_min"] += f1.time
            if fcache is not None:
                pd.DataFrame([row]).to_csv(fcache, index=False)
            if verbose:
                print(f"[{g}] rho = {row['rho_mean']:+.3f} [{row['rho_lower']:+.3f}, "
                      f"{row['rho_upper']:+.3f}]  R-hat {row['max_rhat']:.3f}  "
                      f"({row['time_min']:.1f} min)")
        except Exception as ex:   # keep scanning
            row = dict(gene=g, error=str(ex))
            if verbose:
                print(f"[{g}] failed: {ex}")
        rows.append(row)
    res = pd.DataFrame(rows)
    res["class"] = np.where(res["rho_lower"] > 0, "PCG",
                            np.where(res["rho_upper"] < 0, "NCG", "not DCG"))
    res.loc[res["rho_lower"].isna(), "class"] = None
    res["converged"] = (res["max_rhat"] < 1.01) & (res["divergences"] == 0)
    return res.sort_values("rho_mean", ascending=False).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# 4. Visualisation (matplotlib)                                               #
# --------------------------------------------------------------------------- #

def _tile(ax, grid, values, title, vmin=None, vmax=None, cmap="magma"):
    import matplotlib.pyplot as plt
    nx, ny = int(grid["i"].max()) + 1, int(grid["j"].max()) + 1
    Z = np.asarray(values, float).reshape((ny, nx))         # rows = y, x fastest
    dx = grid["x"].iloc[1] - grid["x"].iloc[0] if nx > 1 else 1.0
    dy = grid["y"].iloc[nx] - grid["y"].iloc[0] if ny > 1 else 1.0
    ext = [grid["x"].min() - dx / 2, grid["x"].max() + dx / 2,
           grid["y"].min() - dy / 2, grid["y"].max() + dy / 2]
    im = ax.imshow(Z, origin="lower", extent=ext, cmap=cmap, vmin=vmin, vmax=vmax,
                   interpolation="nearest")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("x"); ax.set_ylabel("y")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    return im


def plot_grid(res, figsize=(9, 3.8)):
    """Gridded cell counts (Y) and summed expression (M)."""
    import matplotlib.pyplot as plt
    grid = res.grid if isinstance(res, DenMarkResult) else res["grid"]
    gene = res.gene if isinstance(res, DenMarkResult) else res["gene"]
    fig, axes = plt.subplots(1, 2, figsize=figsize, constrained_layout=True)
    _tile(axes[0], grid, grid["Y"], "Gridded cell counts")
    _tile(axes[1], grid, grid["M"], f"Gridded expression: {gene}")
    return fig


def plot_fields(res: DenMarkResult, fields: Optional[pd.DataFrame] = None, figsize=(12, 7)):
    """Latent-field decomposition maps; each row shares one colour scale."""
    import matplotlib.pyplot as plt
    f = res.fields() if fields is None else fields
    lim1 = (min(f.eta1.min(), f.a11w1.min()), max(f.eta1.max(), f.a11w1.max()))
    lim2 = (min(f.eta2.min(), f.a21w1.min(), f.a22w2.min()),
            max(f.eta2.max(), f.a21w1.max(), f.a22w2.max()))
    fig, axes = plt.subplots(2, 3, figsize=figsize, constrained_layout=True)
    _tile(axes[0, 0], f, f.eta1, r"cell density: $\beta_0 + a_{11}w_1$", *lim1)
    _tile(axes[0, 1], f, f.a11w1, r"$a_{11}w_1$", *lim1)
    axes[0, 2].axis("off")
    _tile(axes[1, 0], f, f.eta2, r"expression: $\beta_1 + a_{21}w_1 + a_{22}w_2$", *lim2)
    _tile(axes[1, 1], f, f.a21w1, r"$a_{21}w_1$", *lim2)
    _tile(axes[1, 2], f, f.a22w2, r"$a_{22}w_2$", *lim2)
    return fig


def plot_decay(res: DenMarkResult, decay: Optional[dict] = None, xlabel="Distance",
               figsize=(9, 3.5)):
    """Spatial correlation decay curves with posterior bands."""
    import matplotlib.pyplot as plt
    dc = res.decay() if decay is None else decay
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True, constrained_layout=True)
    for ax, (name, title) in zip(axes, [("density", r"$r_1$: cell density field"),
                                        ("expression", r"$r_2$: expression field")]):
        c = dc["curves"][dc["curves"].field == name]
        ax.fill_between(c.distance, c.lower, c.upper, color="0.8")
        ax.plot(c.distance, c["mean"], color="k", lw=1.5)
        ax.axhline(dc["threshold"], ls="--", color="0.4", lw=1)
        ax.set_title(title, fontsize=9); ax.set_xlabel(xlabel)
    axes[0].set_ylabel("Posterior spatial correlation")
    return fig


def plot_scan(scan: pd.DataFrame, top: Optional[int] = None, figsize=None):
    """Forest plot of rho across genes from scan_genes()."""
    import matplotlib.pyplot as plt
    s = scan.dropna(subset=["rho_mean"]).copy()
    if top is not None and len(s) > top:
        s = s.reindex(s.rho_mean.abs().sort_values(ascending=False).index[:top])
    s = s.sort_values("rho_mean")
    col = s["class"].map({"PCG": "#B2182B", "NCG": "#2166AC"}).fillna("0.5")
    fig, ax = plt.subplots(figsize=figsize or (6, 0.35 * len(s) + 1.2), constrained_layout=True)
    y = np.arange(len(s))
    ax.axvline(0, ls="--", color="0.5", lw=1)
    ax.hlines(y, s.rho_lower, s.rho_upper, colors=col, lw=1.5)
    ax.scatter(s.rho_mean, y, c=col, s=25, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(s.gene)
    ax.set_xlim(-1, 1)
    ax.set_xlabel(r"Density-expression correlation $\rho$ (posterior mean, credible interval)")
    return fig


# --------------------------------------------------------------------------- #
# 5. Simulation                                                               #
# --------------------------------------------------------------------------- #

def _rff_field(ell, rng, n_features=500) -> Callable:
    """Zero-mean unit-variance field with ~Matern-3/2 covariance via random
    Fourier features (2-D Matern-3/2 spectral density = bivariate t, 3 df)."""
    z = rng.standard_normal((n_features, 2))
    u = rng.chisquare(3, n_features)
    omega = (np.sqrt(3) / ell) * z / np.sqrt(u / 3)[:, None]
    phase = rng.uniform(0, 2 * np.pi, n_features)
    amp = rng.standard_normal(n_features)

    def f(x, y):
        arg = np.outer(x, omega[:, 0]) + np.outer(y, omega[:, 1]) + phase
        return np.cos(arg) @ amp * np.sqrt(2 / n_features)
    return f


def simulate(rho: Optional[Dict[str, float]] = None, width=6.0, height=4.0,
             n_cells=5000, ell1=1.0, ell2=1.0, a11=1.0, sigma2=0.8, beta1=0.5,
             seed=1) -> pd.DataFrame:
    """Simulate a marked point pattern from the DenMark data-generating
    process with known rho per gene (see denmark_simulate() in R). The true
    values are in ``df.attrs["truth"]``."""
    rho = rho or {"geneA": 0.6, "geneB": 0.0, "geneC": -0.6}
    rng = np.random.default_rng(seed)
    w1 = _rff_field(ell1, rng)
    xs = rng.uniform(0, width, 60000); ys = rng.uniform(0, height, 60000)
    f1 = w1(xs, ys)
    lam = np.exp(a11 * f1)
    keep = rng.uniform(size=lam.size) < lam / lam.max()
    xs, ys, f1 = xs[keep], ys[keep], f1[keep]
    if xs.size > n_cells:
        s = rng.choice(xs.size, n_cells, replace=False)
        xs, ys, f1 = xs[s], ys[s], f1[s]
    out = pd.DataFrame({"x": xs, "y": ys})
    z1 = (f1 - f1.mean()) / f1.std(ddof=1)
    for g, r in rho.items():
        w2 = _rff_field(ell2, rng)(xs, ys)
        # orthogonalise w2 against w1 over the observed cells so that the
        # realised field correlation equals rho
        w2 = w2 - w2.mean() - (np.dot(w2 - w2.mean(), z1) / np.dot(z1, z1)) * z1
        w2 = w2 / w2.std(ddof=1)
        eta = beta1 + sigma2 * (r * z1 + np.sqrt(1 - r ** 2) * w2)
        out[g] = rng.poisson(np.exp(eta))
    out.attrs["truth"] = dict(rho=rho, ell1=ell1, ell2=ell2, a11=a11, sigma2=sigma2,
                              beta1=beta1, window=(width, height))
    return out
