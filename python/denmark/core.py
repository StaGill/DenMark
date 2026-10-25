"""
DenMark: Density-dependent Marked point process framework (Python version)
==========================================================================

User-facing functions to grid single-cell resolution spatial transcriptomics
(ST) data and fit the DenMark model with Stan through **cmdstanpy**.
This module mirrors the R package DenMark function by function, so results
from R and Python are directly comparable (same grid ordering, same Stan data,
same priors, same summaries).

Model (Stan file shipped inside the package)
--------------------------------------------
``DenMark.stan`` treats the cell count Y(g) and the summed expression M(g) of
grid cell g as Poisson counts with means lambda1(g) and lambda2(g), linked by
two dependent latent fields::

        log lambda1(g) = log|g| + beta0 + a11 w1(g)
        log lambda2(g) = beta1 + log(Y(g) + 0.05) + a21 w1(g) + a22 w2(g)

rho = Astar[2,1] is the density-expression correlation. With a tissue
boundary (alpha shape), ``DenMark-vararea.stan`` is used: the same model with
a per-grid area offset log|g cap tissue|.

Quick start
-----------
>>> import denmark as dm
>>> res = dm.DenMark(df, gene="CDH1", grid_size=10, basis=25)
>>> print(res)                     # rho and sampler diagnostics
>>> res.summary()                  # posterior summary of the key parameters
>>> dm.plot_fields(res)

Requirements: numpy, pandas, scipy, cmdstanpy (with CmdStan installed:
``python -m cmdstanpy.install_cmdstan``); matplotlib for the plots; shapely
for the tissue boundary.

Authors: Mingchi Xu, Alexandra M. Schmidt, Qihuang Zhang
"""

from __future__ import annotations

import hashlib
import os
import re
import time
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Union

import numpy as np
import pandas as pd
from scipy.optimize import brentq

try:
    from cmdstanpy import CmdStanModel
except ImportError:  # pragma: no cover - the module still imports for grid work
    CmdStanModel = None

__all__ = [
    "filter_genes", "alpha_shape", "TissueBoundary", "buildgridpp", "grid_pts_centroid",
    "hsgp_indices", "prepare_data", "compile_model", "make_init", "DenMark", "denmark",
    "DenMarkmain", "DenMarkResult", "scan_genes",
    "matern32", "simulate", "save_result", "load_result", "plot_boundary", "plot_grid",
    "plot_fields", "plot_decay", "plot_scan",
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
    if not pd.api.types.is_numeric_dtype(df[gene]):
        raise ValueError(f"Column '{gene}' is not numeric: pass `gene=`, and drop non-numeric "
                         "columns such as a cell-type label.")
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


def _shapely(what="A tissue boundary"):
    """shapely, which only the tissue boundary needs (R: package sf)."""
    try:
        import shapely
    except ImportError:
        raise ImportError(f"{what} needs the 'shapely' package: pip install shapely") from None
    if int(shapely.__version__.split(".")[0]) < 2:
        raise ImportError(f"{what} needs shapely 2.0 or later: pip install --upgrade shapely")
    return shapely


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


@dataclass(frozen=True)
class TissueBoundary:
    """Tissue region from alpha_shape(): the shapely ``polygon`` and the
    circumradius threshold ``alpha`` that produced it (the R package returns
    an sf polygon with the attribute "alpha").

    :param polygon: shapely Polygon or MultiPolygon.
    :param alpha: the circumradius threshold.
    """
    polygon: object
    alpha: float

    @property
    def area(self) -> float:
        """Area of the tissue region."""
        return float(self.polygon.area)


def _tissue_polygon(boundary):
    """The shapely geometry of a tissue boundary, given as a TissueBoundary or
    as a shapely Polygon or MultiPolygon."""
    poly = boundary.polygon if isinstance(boundary, TissueBoundary) else boundary
    if getattr(poly, "geom_type", None) not in ("Polygon", "MultiPolygon"):
        raise TypeError("`boundary` must be a TissueBoundary from alpha_shape() or a "
                        "shapely Polygon or MultiPolygon.")
    return poly


def alpha_shape(x, y, alpha=None, alpha_mult=3, keep="largest",
                fill_holes=True) -> TissueBoundary:
    """Tissue boundary as an alpha shape of the cell locations.

    The Delaunay triangulation of the cells is computed and every triangle
    whose circumradius is smaller than `alpha` is kept; the union of the kept
    triangles is the tissue region. This follows the boundary step of the
    real-data analyses: compute it from all cells of the section (not only the
    cell type analyzed), so that the density of one cell type is modeled over
    the whole tissue rather than over its own convex hull. The triangulation
    and the unions are computed by GEOS, as in the R package (sf), so both
    give the same polygon.

    :param x: x coordinates of the cells (all cells of the section).
    :param y: y coordinates of the cells.
    :param alpha: circumradius threshold; the default is `alpha_mult` times the
        median Delaunay edge length.
    :param alpha_mult: multiplier for the default `alpha` (3, as in the
        analyses).
    :param keep: "largest" keeps the largest connected piece of tissue, "all"
        keeps every piece.
    :param fill_holes: fill holes inside the tissue.
    :return: a TissueBoundary with the shapely polygon (``.polygon``), the
        threshold (``.alpha``) and the tissue area (``.area``). Needs shapely.
    """
    shapely = _shapely("alpha_shape()")
    if keep not in ("largest", "all"):
        raise ValueError("`keep` must be 'largest' or 'all'.")
    xy = np.column_stack([np.asarray(x, float), np.asarray(y, float)])
    first = np.unique(xy, axis=0, return_index=True)[1]
    xy = xy[np.sort(first)]                            # distinct locations, input order
    if xy.shape[0] < 3:
        raise ValueError("At least 3 distinct cell locations are needed.")
    tri = shapely.get_parts(shapely.delaunay_triangles(shapely.multipoints(xy)))
    cc = shapely.get_coordinates(tri).reshape(-1, 4, 2)     # 4 rows per triangle
    A, B, C = cc[:, 0], cc[:, 1], cc[:, 2]
    la = np.sqrt(((A - B) ** 2).sum(axis=1)); lb = np.sqrt(((B - C) ** 2).sum(axis=1))
    lc = np.sqrt(((C - A) ** 2).sum(axis=1))
    area2 = np.abs((B[:, 0] - A[:, 0]) * (C[:, 1] - A[:, 1])
                   - (C[:, 0] - A[:, 0]) * (B[:, 1] - A[:, 1]))
    with np.errstate(divide="ignore", invalid="ignore"):
        R = np.where(area2 > 0, la * lb * lc / (2 * area2), np.inf)   # circumradius
    if alpha is None:
        alpha = alpha_mult * np.median(np.concatenate([la, lb, lc]))
    kept = tri[R < alpha]
    if kept.size == 0:
        raise ValueError("No triangle passes the alpha threshold; increase `alpha`.")
    parts = shapely.get_parts(shapely.union_all(kept))
    if keep == "largest" and parts.size > 1:
        parts = parts[[np.argmax(shapely.area(parts))]]
    if fill_holes:
        parts = shapely.polygons(shapely.get_exterior_ring(parts))
    return TissueBoundary(polygon=shapely.union_all(parts), alpha=float(alpha))


def _grid_area(grid: pd.DataFrame, step, boundary) -> pd.DataFrame:
    """Area of every grid rectangle that lies inside the tissue polygon, and
    the grid type: "interior" (fully inside), "boundary" (partly) or "exterior"."""
    shapely = _shapely()
    dx, dy = step
    full = dx * dy
    x0 = grid["x"].to_numpy(float) - dx / 2
    y0 = grid["y"].to_numpy(float) - dy / 2
    corners = [(x0, y0), (x0 + dx, y0), (x0 + dx, y0 + dy), (x0, y0 + dy), (x0, y0)]
    rects = shapely.polygons(np.stack([np.column_stack(c) for c in corners], axis=1))
    tissue = _tissue_polygon(boundary)
    hit = shapely.intersects(rects, tissue)
    area = np.zeros(len(grid))
    if hit.any():
        area[hit] = shapely.area(shapely.intersection(rects[hit], tissue))
    tol = 1e-6 * full
    gtype = np.where(area <= tol, "exterior", np.where(area >= full - tol, "interior", "boundary"))
    area[gtype == "interior"] = full
    return pd.DataFrame({"area": area, "type": gtype})


def buildgridpp(grid_size, dataset: pd.DataFrame, gene=None, bbox=None,
                boundary=None) -> dict:
    """Discretise cell locations and expression on a regular grid.

    Each cell is assigned to exactly one grid cell. Vectors are ordered with x
    varying fastest (grid_index = i + j * nx, 0-based), identical to the R
    version and to ``Ymat.ravel(order="F")``. The bounding box `bbox`
    defaults to the data range, or to the extent of `boundary` when one is
    given; cells outside it are dropped. With a tissue `boundary` (from
    alpha_shape(), or a shapely Polygon or MultiPolygon), every grid cell gets
    its area inside the tissue (`area`) and a `type`: interior, boundary or
    exterior.

    Returns a dict with keys
      loc_mark_mat (N x 2: point = Y, mark = M), Ymat, Mmat (nx x ny, rows = x),
      grid (DataFrame: grid_index, i, j, x, y [centroids], Y, M, area [area
      inside the tissue; the full grid cell area without a boundary], type),
      bbox, step, grid_area, grid_size, gene, boundary (None without one).
    """
    x, y, v, gene = _xyz(dataset, gene)
    nx, ny = _grid_dims(grid_size)
    bbox_from_boundary = bbox is None and boundary is not None
    if bbox_from_boundary:
        b = _tissue_polygon(boundary).bounds          # (xmin, ymin, xmax, ymax)
        bbox = (b[0], b[2], b[1], b[3])
    bb = _bbox(x, y, bbox)
    inside = (x >= bb[0]) & (x <= bb[1]) & (y >= bb[2]) & (y <= bb[3])
    if bbox_from_boundary and not inside.all():
        print(f"{int((~inside).sum())} cell(s) lie outside the extent of the tissue "
              "boundary and are dropped.")
    x, y, v = x[inside], y[inside], v[inside]
    if x.size == 0:
        raise ValueError("No cells inside the extent of the tissue boundary; are the cells "
                         "and the boundary in the same units?" if bbox_from_boundary
                         else "No cells inside the bounding box.")
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
        "Y": Y, "M": M, "area": dx * dy, "type": "interior"})
    if boundary is not None:
        ga = _grid_area(grid, (dx, dy), boundary)
        grid["area"] = ga["area"]
        grid["type"] = ga["type"]
    return dict(loc_mark_mat=np.column_stack([Y, M]),
                Ymat=Y.reshape((nx, ny), order="F"),
                Mmat=M.reshape((nx, ny), order="F"),
                grid=grid, bbox=bb, step=(dx, dy), grid_area=dx * dy,
                grid_size=(nx, ny), gene=gene, boundary=boundary)


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
                 boundfactor=1.5, bbox=None, L1=None, L2=None, boundary=None,
                 min_area_frac=0) -> dict:
    """Build the Stan data dict for DenMark.

    The grid is centered at the origin, so the ROI becomes [-L1, L1] x
    [-L2, L2]; grid centroids are the HSGP coordinates and the HSGP domain is
    boundfactor * (L1, L2). With a tissue `boundary` (from alpha_shape()),
    the grid cells outside the tissue are removed first and the ROI is that of
    the retained grid; each remaining grid cell enters with its own area
    inside the tissue, so `log_grid_area` becomes a vector and DenMark() uses
    the variable-area model DenMark-vararea.stan.

    :param min_area_frac: with a boundary, drop the grid cells whose area
        inside the tissue is below this fraction of a full grid cell (default
        0: keep all).
    :return: dict(stan_data, grid, info). `grid` lacks the grid cells outside
        the tissue; `info` also holds vararea, boundary, n_grid_total,
        n_interior and n_boundary.
    """
    if boundfactor < 1:
        raise ValueError("`boundfactor` must be >= 1.")
    g = buildgridpp(grid_size, dataset, gene=gene, bbox=bbox, boundary=boundary)
    bb = g["bbox"]
    vararea = boundary is not None
    grid = g["grid"].copy()
    if vararea:
        keep = (grid["type"] != "exterior") & (grid["area"] >= min_area_frac * g["grid_area"])
        lost = int(grid["Y"][~keep].sum())
        if lost > 0:
            print(f"{lost} cell(s) fall in grid cells outside the tissue and are not modeled.")
        grid = grid[keep].reset_index(drop=True)
    # center the retained grid at the origin: ROI = [-L1, L1] x [-L2, L2]
    centre = np.array([(grid["x"].min() + grid["x"].max()) / 2,
                       (grid["y"].min() + grid["y"].max()) / 2])
    half = np.array([(grid["x"] - centre[0]).abs().max() + g["step"][0] / 2,
                     (grid["y"] - centre[1]).abs().max() + g["step"][1] / 2])
    for name, user, val in (("L1", L1, half[0]), ("L2", L2, half[1])):
        if user is not None and abs(user - val) > 0.01 * val:
            warnings.warn(f"{name} = {user:.4g} supplied but the data give {val:.4g}; "
                          "the data value is used.")
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
        log_grid_area=(np.log(grid["area"].to_numpy(float)) if vararea
                       else float(np.log(g["grid_area"]))),     # area offset
        mrange=dmax / (2 * 2.75),
        indices=S,
        L=(boundfactor * half).tolist(),
        is_centerted_PHI=0)
    grid["cx"], grid["cy"] = coords[:, 0], coords[:, 1]
    info = dict(gene=g["gene"], grid_size=g["grid_size"], step=g["step"],
                grid_area=g["grid_area"], bbox=bb, centre=centre,
                L1=float(half[0]), L2=float(half[1]), basis=basis,
                boundfactor=boundfactor, max_dist=dmax,
                n_cells=int(grid["Y"].sum()), Ymat=g["Ymat"], Mmat=g["Mmat"],
                vararea=vararea, boundary=boundary, n_grid_total=len(g["grid"]),
                n_interior=int((grid["type"] == "interior").sum()),
                n_boundary=int((grid["type"] == "boundary").sum()))
    return dict(stan_data=stan_data, grid=grid, info=info)


# --------------------------------------------------------------------------- #
# 2. Model compilation, initial values, fitting                               #
# --------------------------------------------------------------------------- #

def _stan_path(stan_file=None, vararea=False) -> Path:
    if stan_file is not None:
        p = Path(stan_file)
        if not p.exists():
            raise FileNotFoundError(f"Stan file not found: {p}")
        return p.resolve()
    fname = "DenMark" + ("-vararea" if vararea else "") + ".stan"
    for d in [STAN_DIR, _HERE / "stan"]:
        if d is not None and (Path(d) / fname).exists():
            return (Path(d) / fname).resolve()
    raise FileNotFoundError(f"Cannot find '{fname}'. Pass stan_file= or set denmark.STAN_DIR.")


def _cached_stan_file(path: Path) -> Path:
    """Copy of a Stan file in the writable cache dir, where CmdStan then
    writes the executable (the installed package may be read-only). Each
    source folder gets its own subfolder, so two files with the same name
    never share an executable; the copy is rewritten only when its contents
    differ, so an unchanged model is not recompiled."""
    cache = Path(os.environ.get("DENMARK_CACHE_DIR", Path.home() / ".cache" / "denmark"))
    dest = cache / hashlib.sha256(os.fsencode(path.parent)).hexdigest()[:12] / path.name
    dest.parent.mkdir(parents=True, exist_ok=True)
    code = path.read_bytes()
    if not dest.exists() or dest.read_bytes() != code:
        dest.write_bytes(code)
    return dest


def compile_model(stan_file=None, vararea=False, force_compile=False) -> "CmdStanModel":
    """Compile (once per session) and return a DenMark CmdStanModel. The Stan
    file is copied to DENMARK_CACHE_DIR (default ~/.cache/denmark) and
    compiled there.

    :param stan_file: optional path to a Stan file, compiled instead of the
        packaged model.
    :param vararea: use the variable-area model DenMark-vararea.stan, the same
        model with one area offset `log_grid_area` per grid cell (used with a
        tissue boundary).
    :param force_compile: recompile even if a compiled model is available.
    """
    if CmdStanModel is None:
        raise ImportError("cmdstanpy is required: pip install cmdstanpy")
    path = str(_stan_path(stan_file, vararea))
    if force_compile or path not in _MODEL_CACHE:
        _MODEL_CACHE[path] = CmdStanModel(stan_file=str(_cached_stan_file(Path(path))),
                                          force_compile=force_compile)
    return _MODEL_CACHE[path]


def make_init(stan_data, rho0=0.3) -> dict:
    """Valid starting values. The Cholesky factor of the correlation matrix
    is lower triangular with unit-norm rows:
    Astar = [[1, 0], [rho0, sqrt(1 - rho0^2)]]. The length-scales start at
    the prior mean (`mrange`) and the intercepts at their empirical values
    (for beta0, the number of cells per unit of tissue area).

    :param stan_data: the Stan data dict, as in ``prepare_data()["stan_data"]``.
    :param rho0: starting value of rho, in (-1, 1).
    :return: a dict of initial values for CmdStanModel.sample().
    """
    if abs(rho0) >= 1:
        raise ValueError("`rho0` must lie in (-1, 1).")
    mstar = stan_data["mstar"]
    y1 = np.asarray(stan_data["y1"]); y2 = np.asarray(stan_data["y2"])
    area = np.exp(np.resize(np.asarray(stan_data["log_grid_area"], float), y1.size))
    beta0 = float(np.log(max(y1.sum(), 1) / area.sum()))
    beta1 = float(np.log(max(y2.sum(), 1) / (y1 + 0.05).sum()))
    ell = [float(stan_data["mrange"])] * 2
    base = dict(beta0=beta0, beta1=beta1, ell=ell,
                betab1=[0.0] * mstar, betab2=[0.0] * mstar)
    base.update(Astar=[[1.0, 0.0], [rho0, float(np.sqrt(1 - rho0 ** 2))]],
                sigma=[1.0, 0.5])
    return base


_KEY_VARS = ["beta0", "beta1", "ell[1]", "ell[2]", "Astar[2,1]", "sigma[1]", "sigma[2]"]
_LABELS = {
    "beta0": "intercept, cell density", "beta1": "intercept, expression",
    "ell[1]": "length-scale, density field", "ell[2]": "length-scale, expression field",
    "Astar[2,1]": "rho: density-expression correlation",
    "sigma[1]": "SD, density field", "sigma[2]": "SD, expression field",
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
    gene: str
    stan_data: dict
    grid: pd.DataFrame
    info: dict
    settings: dict
    time: float = np.nan                       # minutes
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
        idx = [cols.index(v) for v in _KEY_VARS]
        return arr[:, :, idx]

    def summary(self, probs=(0.025, 0.25, 0.5, 0.75, 0.975)) -> pd.DataFrame:
        """Posterior summary of the key parameters (mean, sd, quantiles,
        R-hat, bulk/tail ESS)."""
        if self._summary is not None and tuple(probs) == (0.025, 0.25, 0.5, 0.75, 0.975):
            return self._summary
        arr = self._key_array()
        rows = []
        for k, v in enumerate(_KEY_VARS):
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
        a11 = sigma1, a21 = sigma2*rho, a22 = sigma2*sqrt(1 - rho^2)."""
        arr = self._key_array()
        flat = arr.transpose(1, 0, 2).reshape(-1, arr.shape[2])   # chain-major, like R
        d = dict(zip(_KEY_VARS, flat.T))
        rho, s1, s2 = d["Astar[2,1]"], d["sigma[1]"], d["sigma[2]"]
        return pd.DataFrame(dict(beta0=d["beta0"], beta1=d["beta1"],
                                 ell1=d["ell[1]"], ell2=d["ell[2]"], rho=rho,
                                 sigma1=s1, sigma2=s2, a11=s1, a21=s2 * rho,
                                 a22=s2 * np.sqrt(np.clip(1 - rho ** 2, 0, None))))

    @property
    def rho(self) -> np.ndarray:
        """Posterior draws of rho, the density-expression correlation."""
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
        lines = [f"DenMark fit  |  gene {self.gene}",
                 f"  grid {nx} x {ny} ({self.stan_data['n']} grids, {self.info['n_cells']} cells)"
                 f"  |  basis {b}  |  boundfactor {self.info['boundfactor']:.2f}"]
        s = self.summary()
        r = s.loc[s.variable == "Astar[2,1]"].iloc[0]
        lines.append(f"  rho (density-expression correlation): {r['mean']:.3f}  "
                     f"[95% CrI {r['q2.5']:.3f}, {r['q97.5']:.3f}]")
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


def _removed_args(kwargs) -> None:
    """Fail before any gridding or compilation when a call passes an argument
    removed in 0.9.0, when the package kept only the DenMark model."""
    bad = [k for k in ("model", "compute_waic", "fit_simple") if k in kwargs]
    if bad:
        raise TypeError(f"`{bad[0]}` was removed in denmark 0.9.0: the package fits "
                        "only the DenMark model and no longer computes WAIC.")


def DenMark(dataset: pd.DataFrame, gene=None, grid_size=10, basis=25,
            boundfactor=1.5, bbox=None, chains=4, parallel_chains=None,
            iter_warmup=1000, iter_sampling=1000, adapt_delta=0.95,
            max_treedepth=10, seed=2024, refresh=400, inits=None, rho0=0.3,
            stan_file=None, output_dir=None, keep_fit=True,
            quiet=False, L1=None, L2=None, boundary=None, min_area_frac=0,
            **sample_kwargs) -> DenMarkResult:
    """Fit DenMark to one gene.

    Parameters mirror the R function DenMark(): `dataset` has columns x, y and
    expression counts (one row per cell); grid_size, basis may be scalars or
    (nx, ny); `stan_file` replaces DenMark.stan and must declare the same
    parameters; remaining arguments go to CmdStanModel.sample().

    A tissue section is rarely a rectangle. With `boundary` (an alpha shape of
    all cells of the section, from alpha_shape()), the grid cells outside the
    tissue are dropped and every remaining grid cell enters with its area
    inside the tissue, through the variable-area model DenMark-vararea.stan
    (the same model with one `log_grid_area` per grid cell); a custom
    `stan_file` must then declare ``vector[n] log_grid_area``.
    `min_area_frac` is passed to prepare_data(). Without `boundary`, the
    whole bounding box is gridded with a common grid area.
    """
    _removed_args(sample_kwargs)
    prep = prepare_data(dataset, gene=gene, grid_size=grid_size, basis=basis,
                        boundfactor=boundfactor, bbox=bbox, L1=L1, L2=L2,
                        boundary=boundary, min_area_frac=min_area_frac)
    sd = prep["stan_data"]
    if sd["y2"].sum() == 0:
        raise ValueError(f"Gene '{prep['info']['gene']}' has zero counts in the ROI.")
    if stan_file is not None and prep["info"]["vararea"] and not any(
            re.search(r"vector\[n\]\s*log_grid_area", line)
            for line in Path(stan_file).read_text(errors="replace").splitlines()):
        raise ValueError("With a tissue boundary the Stan model needs `vector[n] log_grid_area` "
                         "(use the '-vararea' models).")
    mod = compile_model(stan_file, vararea=prep["info"]["vararea"])
    if inits is None:
        inits = make_init(sd, rho0=rho0)
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

    res = DenMarkResult(fit=fit, gene=prep["info"]["gene"],
                        stan_data=sd, grid=prep["grid"], info=prep["info"],
                        settings=dict(chains=chains, iter_warmup=iter_warmup,
                                      iter_sampling=iter_sampling, adapt_delta=adapt_delta,
                                      max_treedepth=max_treedepth, seed=seed),
                        time=elapsed)
    res._summary = res.summary()
    res.diagnostics = _diagnostics(fit, res._summary, max_treedepth)
    res.hsgp = res.check_hsgp(warn=not quiet)
    if not keep_fit:
        res._key_draws = res._key_array()
        res.fit = None
    return res


denmark = DenMark   # snake_case alias


def DenMarkmain(grid_size, dataset, basis, boundfactor, L1=None, L2=None,
                denmark_version="M2", **kwargs) -> DenMarkResult:
    """Backward-compatible wrapper of the original interface; the other
    arguments are passed to DenMark().

    :param denmark_version: "M2", the label of the DenMark model before
        0.9.0, or a path to a DenMark Stan file; "M1" raises ValueError.
    :return: a DenMarkResult.
    """
    if denmark_version == "M1":
        raise ValueError("The simplified model (M1) was removed in denmark 0.9.0; "
                         "DenMark() fits the main model.")
    if denmark_version == "M2":
        return DenMark(dataset, grid_size=grid_size, basis=basis,
                       boundfactor=boundfactor, L1=L1, L2=L2, **kwargs)
    return DenMark(dataset, stan_file=denmark_version, grid_size=grid_size, basis=basis,
                   boundfactor=boundfactor, L1=L1, L2=L2, **kwargs)


# --------------------------------------------------------------------------- #
# 3. DCG discovery                                                            #
# --------------------------------------------------------------------------- #

def scan_genes(data: pd.DataFrame, genes: Optional[Iterable[str]] = None,
               save_dir=None, cred_level=0.95, quiet=True,
               verbose=True, **kwargs) -> pd.DataFrame:
    """Fit DenMark gene by gene and call density-correlated genes:
    PCG if the credible interval of rho is above 0, NCG if below 0, else
    "not DCG". With save_dir, one CSV row per gene is cached (resumable).
    Further arguments (grid_size, basis, boundary, ...) go to DenMark()."""
    _removed_args(kwargs)
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
            f2 = DenMark(data[["x", "y", g]], gene=g, keep_fit=False,
                         quiet=quiet, **kwargs)
            rho = f2.rho
            row = dict(gene=g, rho_mean=rho.mean(), rho_median=np.median(rho),
                       rho_lower=np.quantile(rho, a), rho_upper=np.quantile(rho, 1 - a),
                       prob_positive=(rho > 0).mean(), total_counts=int(f2.stan_data["y2"].sum()),
                       max_rhat=f2.diagnostics["max_rhat"],
                       min_ess_bulk=f2.diagnostics["min_ess_bulk"],
                       divergences=f2.diagnostics["num_divergent"],
                       time_min=f2.time, error=None)
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

def _add_boundary(ax, boundary, color="black"):
    """Draw every ring of the tissue boundary, exterior and holes."""
    if boundary is None:
        return
    poly = _tissue_polygon(boundary)
    for part in getattr(poly, "geoms", [poly]):
        for ring in (part.exterior, *part.interiors):
            xy = np.asarray(ring.coords)
            ax.plot(xy[:, 0], xy[:, 1], color=color, lw=1)


def _tile(ax, grid, values, title, vmin=None, vmax=None, cmap="magma", boundary=None):
    import matplotlib.pyplot as plt
    i, j = grid["i"].to_numpy(int), grid["j"].to_numpy(int)
    x, y = grid["x"].to_numpy(float), grid["y"].to_numpy(float)
    # grid cells may be missing (outside the tissue); they stay empty
    Z = np.full((np.ptp(j) + 1, np.ptp(i) + 1), np.nan)    # rows = y, x fastest
    Z[j - j.min(), i - i.min()] = np.asarray(values, float)
    dx = np.ptp(x) / np.ptp(i) if np.ptp(i) > 0 else 1.0
    dy = np.ptp(y) / np.ptp(j) if np.ptp(j) > 0 else 1.0
    ext = [x.min() - dx / 2, x.max() + dx / 2, y.min() - dy / 2, y.max() + dy / 2]
    im = ax.imshow(Z, origin="lower", extent=ext, cmap=cmap, vmin=vmin, vmax=vmax,
                   interpolation="nearest")
    _add_boundary(ax, boundary, color="0.6")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("x"); ax.set_ylabel("y")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    return im


def plot_grid(res, figsize=(9, 3.8)):
    """Gridded cell counts (Y) and summed expression (M). With a tissue
    boundary, only the grid cells inside the tissue are drawn, together with
    the tissue outline."""
    import matplotlib.pyplot as plt
    grid = res.grid if isinstance(res, DenMarkResult) else res["grid"]
    gene = res.gene if isinstance(res, DenMarkResult) else res["gene"]
    boundary = res.info.get("boundary") if isinstance(res, DenMarkResult) else res.get("boundary")
    if "type" in grid:
        grid = grid[grid["type"] != "exterior"]
    fig, axes = plt.subplots(1, 2, figsize=figsize, constrained_layout=True)
    _tile(axes[0], grid, grid["Y"], "Gridded cell counts", boundary=boundary)
    _tile(axes[1], grid, grid["M"], f"Gridded expression: {gene}", boundary=boundary)
    return fig


def plot_boundary(cells: pd.DataFrame, boundary, grid: Optional[dict] = None,
                  highlight=None, xlabel="x", ylabel="y", point_size=1.0, figsize=None):
    """Tissue boundary and grid classification.

    Left: all cells (gray), the cells of the analyzed type (colored) and the
    alpha-shape tissue boundary. Right (when `grid` is given): the grid cells
    colored by type (interior, boundary, exterior) with the boundary outline.

    :param cells: DataFrame with the columns x and y (all cells).
    :param boundary: the tissue boundary from alpha_shape(), or a shapely
        Polygon or MultiPolygon.
    :param grid: optional output of ``buildgridpp(..., boundary=boundary)``.
    :param highlight: optional boolean array, one element per row of `cells`:
        True for the cells of the analyzed type.
    :param xlabel: label of the x axis.
    :param ylabel: label of the y axis.
    :param point_size: marker area of the cells, in points squared.
    :param figsize: figure size; the default depends on the number of panels.
    :return: a matplotlib Figure with one panel, or two with `grid`.
    """
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch, Rectangle
    from matplotlib.transforms import ScaledTranslation
    cells = pd.DataFrame(cells)
    group = (np.full(len(cells), "cells") if highlight is None else
             np.where(np.asarray(highlight, bool), "analyzed cell type", "other cells"))
    colors = {"cells": "0.4", "other cells": "0.75", "analyzed cell type": "#B2182B"}
    npanel = 1 if grid is None else 2
    fig, axes = plt.subplots(1, npanel, figsize=figsize or (5 * npanel, 4.2),
                             constrained_layout=True, squeeze=False)

    def legend_below(ax, handles):     # below the x axis label, like R's legend.position = "bottom"
        ax.legend(handles=handles, loc="upper center", ncol=len(handles), frameon=False,
                  fontsize=8, bbox_to_anchor=(0.5, 0), bbox_transform=ax.transAxes
                  + ScaledTranslation(0, -0.45, fig.dpi_scale_trans))

    ax = axes[0, 0]
    shown = [k for k in colors if (group == k).any()]    # analyzed cells drawn last, on top
    for k in shown:
        sel = group == k
        ax.scatter(cells["x"].to_numpy()[sel], cells["y"].to_numpy()[sel], s=point_size,
                   color=colors[k], linewidths=0)
    _add_boundary(ax, boundary)
    a = getattr(boundary, "alpha", None)
    ax.set_title(f"Alpha-shape tissue boundary (alpha = {'NA' if a is None else f'{a:.3g}'})",
                 fontsize=9)
    legend_below(ax, [Line2D([], [], ls="", marker="o", ms=5, color=colors[k], label=k)
                      for k in shown])
    ax.set_aspect("equal"); ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    if grid is None:
        return fig
    g, (dx, dy) = grid["grid"], grid["step"]
    fills = {"interior": "#9ECAE1", "boundary": "#FDAE6B", "exterior": "white"}
    ax = axes[0, 1]
    rects = [Rectangle((xc - dx / 2, yc - dy / 2), dx, dy) for xc, yc in zip(g["x"], g["y"])]
    ax.add_collection(PatchCollection(rects, facecolors=g["type"].map(fills).tolist(),
                                      edgecolors="0.4", linewidths=0.3))
    _add_boundary(ax, boundary)
    n = g["type"].value_counts()
    ax.set_title(f"Grid cells: {n.get('interior', 0)} interior, {n.get('boundary', 0)} "
                 f"boundary, {n.get('exterior', 0)} exterior", fontsize=9)
    legend_below(ax, [Patch(facecolor=c, edgecolor="0.4", label=t) for t, c in fills.items()])
    ax.autoscale_view()
    ax.set_aspect("equal"); ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    return fig


def plot_fields(res: DenMarkResult, fields: Optional[pd.DataFrame] = None, figsize=(12, 7)):
    """Latent-field decomposition maps; each row shares one color scale.
    With a tissue boundary, the tissue outline is drawn on every map."""
    import matplotlib.pyplot as plt
    f = res.fields() if fields is None else fields
    bd = res.info.get("boundary")
    lim1 = (min(f.eta1.min(), f.a11w1.min()), max(f.eta1.max(), f.a11w1.max()))
    lim2 = (min(f.eta2.min(), f.a21w1.min(), f.a22w2.min()),
            max(f.eta2.max(), f.a21w1.max(), f.a22w2.max()))
    fig, axes = plt.subplots(2, 3, figsize=figsize, constrained_layout=True)
    _tile(axes[0, 0], f, f.eta1, r"cell density: $\beta_0 + a_{11}w_1$", *lim1, boundary=bd)
    _tile(axes[0, 1], f, f.a11w1, r"$a_{11}w_1$", *lim1, boundary=bd)
    axes[0, 2].axis("off")
    _tile(axes[1, 0], f, f.eta2, r"expression: $\beta_1 + a_{21}w_1 + a_{22}w_2$", *lim2,
          boundary=bd)
    _tile(axes[1, 1], f, f.a21w1, r"$a_{21}w_1$", *lim2, boundary=bd)
    _tile(axes[1, 2], f, f.a22w2, r"$a_{22}w_2$", *lim2, boundary=bd)
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


def _tissue_radius(theta):
    """Irregular, lobed tissue shape inside [0, width] x [0, height]."""
    return 0.95 * (1 + 0.22 * np.sin(3 * theta + 0.5) + 0.12 * np.cos(2 * theta)
                   + 0.08 * np.sin(5 * theta + 1)) / 1.42


def _in_tissue(x, y, width, height):
    u = (x - width / 2) / (width / 2); v = (y - height / 2) / (height / 2)
    return np.sqrt(u ** 2 + v ** 2) < _tissue_radius(np.arctan2(v, u))


def simulate(rho: Optional[Dict[str, float]] = None, width=6.0, height=4.0,
             n_cells=5000, ell1=1.0, ell2=1.0, a11=1.0, sigma2=0.8, beta1=0.5,
             seed=1, tissue="rectangle", n_other=0, celltype="Target") -> pd.DataFrame:
    """Simulate a marked point pattern from the DenMark data-generating
    process with known rho per gene (see denmark_simulate() in R). The true
    values are in ``df.attrs["truth"]``.

    The fields are standardized, and each expression field is orthogonalized
    against the density field, on a regular lattice over the tissue, so that
    their spatial correlation equals rho. With ``tissue="irregular"``, the
    cells lie in a lobed tissue section inside the window, so the analysis
    needs a tissue boundary (see alpha_shape()). The cells of other types are
    drawn from the candidate locations that the thinning rejected, so they are
    densest where the analyzed cells are sparse; their gene counts are Poisson
    with mean 0.3. The draws use numpy's generator, so a seed does not give
    the cells that R's denmark_simulate() gives.

    :param tissue: "rectangle" (cells fill the window) or "irregular" (cells
        lie in a lobed tissue section).
    :param n_other: number of cells of other types; together with the analyzed
        cells they define the tissue boundary.
    :param celltype: label of the analyzed cell type (used when n_other > 0;
        the other cells are labeled "Other").
    :return: a DataFrame with the columns x, y, celltype (when n_other > 0)
        and one count column per gene. The truth includes `tissue` and, for
        ``tissue="irregular"``, the true tissue `outline` as a dict of lists x
        and y (``pd.DataFrame(outline)`` makes it a table).
    """
    if tissue not in ("rectangle", "irregular"):
        raise ValueError("`tissue` must be 'rectangle' or 'irregular'.")
    rho = rho or {"geneA": 0.6, "geneB": 0.0, "geneC": -0.6}
    rng = np.random.default_rng(seed)
    w1 = _rff_field(ell1, rng)
    # candidate locations: uniform over the window, restricted to the tissue
    xs = rng.uniform(0, width, 80000); ys = rng.uniform(0, height, 80000)
    if tissue == "irregular":
        inside = _in_tissue(xs, ys, width, height)
        xs, ys = xs[inside], ys[inside]
    # thinning: keep a candidate with probability lambda / max(lambda)
    f1 = w1(xs, ys)
    lam = np.exp(a11 * f1)
    keep = rng.uniform(size=lam.size) < lam / lam.max()
    xs_o, ys_o = xs[~keep], ys[~keep]                  # pool for other cells
    xs, ys, f1 = xs[keep], ys[keep], f1[keep]
    if xs.size > n_cells:
        s = rng.choice(xs.size, n_cells, replace=False)
        xs, ys, f1 = xs[s], ys[s], f1[s]
    out = pd.DataFrame({"x": xs, "y": ys})
    # regular lattice over the tissue: the fields are standardized and made
    # orthogonal over SPACE (area-weighted), which is what rho measures in the
    # model; two independent smooth fields can be correlated by chance on a
    # finite window, so w2 is orthogonalized against w1 on this lattice
    lx = np.tile(np.linspace(0, width, 150), 100)
    ly = np.repeat(np.linspace(0, height, 100), 150)
    if tissue == "irregular":
        inside = _in_tissue(lx, ly, width, height)
        lx, ly = lx[inside], ly[inside]
    l1 = w1(lx, ly)
    m1, s1 = l1.mean(), l1.std(ddof=1)
    z1 = (f1 - m1) / s1; zl1 = (l1 - m1) / s1
    for g, r in rho.items():
        w2f = _rff_field(ell2, rng)
        l2 = w2f(lx, ly)
        b = np.cov(l2, zl1)[0, 1] / zl1.var(ddof=1)
        r2 = l2 - b * zl1
        w2 = (w2f(xs, ys) - b * z1 - r2.mean()) / r2.std(ddof=1)
        eta = beta1 + sigma2 * (r * z1 + np.sqrt(1 - r ** 2) * w2)
        out[g] = rng.poisson(np.exp(eta))
    if n_other > 0:
        s = rng.choice(xs_o.size, min(n_other, xs_o.size), replace=False)
        other = pd.DataFrame({"x": xs_o[s], "y": ys_o[s]})
        for g in rho:
            other[g] = rng.poisson(0.3, len(other))
        out["celltype"] = celltype; other["celltype"] = "Other"
        out = pd.concat([out, other], ignore_index=True)[["x", "y", "celltype", *rho]]
    outline = None
    if tissue == "irregular":
        th = np.linspace(-np.pi, np.pi, 361)
        r = _tissue_radius(th)
        outline = dict(x=(width / 2 * (1 + r * np.cos(th))).tolist(),
                       y=(height / 2 * (1 + r * np.sin(th))).tolist())
    out.attrs["truth"] = dict(rho=rho, ell1=ell1, ell2=ell2, a11=a11, sigma2=sigma2,
                              beta1=beta1, window=(width, height), tissue=tissue,
                              outline=outline)
    return out
