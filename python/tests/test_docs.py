"""Run the Python examples of the README files and of the Get started in Python page as
written, on the bundled toy data."""

import re
import sys

import pytest

from .conftest import PKG_ROOT, REPO_ROOT, cmdstan_available

PAGES = {"README.md": REPO_ROOT / "README.md", "python/README.md": PKG_ROOT / "README.md",
         "Get started in Python": REPO_ROOT / "vignettes" / "articles" / "python.Rmd"}
# Examples that read the reader's own files cannot run here.
NEEDS_FILES = ('read_csv("',)
# Examples that fit a model need CmdStan, and the examples after them use the fits.
NEEDS_CMDSTAN = ("DenMark(", "denmark(", "DenMarkmain(", "scan_genes(", "compile_model(")


def _blocks(path):
    """The ```python blocks of a page that can run here, in order."""
    text = path.read_text(encoding="utf-8")
    blocks = [b for b in re.findall(r"```python[ \t]*\n(.*?)```", text, flags=re.S)
              if not any(s in b for s in NEEDS_FILES)]
    if not cmdstan_available():
        for k, block in enumerate(blocks):
            if any(s in block for s in NEEDS_CMDSTAN):
                return blocks[:k]
    return blocks


@pytest.mark.usefixtures("toy")   # skips until the toy data file is in place
@pytest.mark.parametrize("page", [p for p, path in PAGES.items() if path.exists()])
def test_examples_run_as_written(page, tmp_path, monkeypatch):
    blocks = _blocks(PAGES[page])
    if not blocks:
        pytest.skip(f"{page} has no Python example that can run here")
    if any(re.search(r"\bplot_\w+\(|matplotlib|plt\.", b) for b in blocks):
        pytest.importorskip("matplotlib").use("Agg")
    if any(re.search(r"\balpha_shape\(|shapely", b) for b in blocks):
        pytest.importorskip("shapely")
    monkeypatch.chdir(tmp_path)   # files the examples write land here
    namespace = {}
    for block in blocks:
        exec(compile(block, page, "exec"), namespace)
    if "matplotlib.pyplot" in sys.modules:
        sys.modules["matplotlib.pyplot"].close("all")
