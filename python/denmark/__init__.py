"""DenMark for Python: density-dependent marked point process models for
single-cell resolution spatial transcriptomics.

Implementation accompanying Xu, Schmidt and Zhang (under review at Cell
Reports Methods). It builds the same Stan data as the DenMark R package and
fits the same Stan model through cmdstanpy. The workflow maps to:

  tissue boundary from all cells         alpha_shape()
  grid the cells, build the Stan data    buildgridpp(), prepare_data()
  fit one gene                           DenMark()
  scan genes and call DCGs               scan_genes()
  summaries, latent fields, decay        DenMarkResult.summary(), .fields(), .decay()
  figures                                plot_boundary(), plot_grid(), plot_fields(),
                                         plot_decay(), plot_scan()

denmark.datasets.load_toy() returns simulated toy data with known answers.
"""

import sys
import types

__version__ = "0.9.0"

from . import core, datasets
from .core import *  # noqa: F401,F403

__all__ = ["__version__", *core.__all__, "datasets"]


class _Package(types.ModuleType):
    """Forwards ``denmark.STAN_DIR`` to ``denmark.core``, where the Stan files
    are looked up."""

    @property
    def STAN_DIR(self):
        return core.STAN_DIR

    @STAN_DIR.setter
    def STAN_DIR(self, value):
        core.STAN_DIR = value


sys.modules[__name__].__class__ = _Package
