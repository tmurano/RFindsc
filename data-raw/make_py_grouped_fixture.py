"""Fixture for tests/testthat/test-baseline-groups.R: per-group baseline scores from the Python package rfindsc.

The R test scores the same counts with score_cells(baseline_groups = ) and requires agreement with these Python
scores to 1e-8 (and identical signatures). Input: the simulated dataset shared with rfindsc-py
(rfindsc-py/tests/data/toy_*), 30 cells of each of its 3 groups; panels are written into the fixture so that the test
does not depend on the bundled references. Not run at install time. Run from the RFindsc folder:

    python data-raw/make_py_grouped_fixture.py [path/to/rfindsc-py]
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import scipy.io as sio
import scipy.sparse as sp

PY_REPO = Path(sys.argv[1] if len(sys.argv) > 1 else "../rfindsc-py").resolve()
sys.path.insert(0, str(PY_REPO / "src"))
import rfindsc  # noqa: E402

SRC = PY_REPO / "tests" / "data"
OUT = Path("tests/testthat/fixtures/py_grouped")
OUT.mkdir(parents=True, exist_ok=True)

with gzip.open(SRC / "toy_counts.mtx.gz", "rb") as fh:
    X = sp.csc_matrix(sio.mmread(fh), dtype=np.float64)
genes = (SRC / "toy_genes.txt").read_text(encoding="utf-8").splitlines()
cells = pd.read_csv(SRC / "toy_cells.tsv", sep="\t")

# 30 cells per group, evenly spaced (group A: half with a DAM programme; B: half with a glucocorticoid programme)
pick = np.concatenate([np.flatnonzero(cells["group"].to_numpy() == g)[np.linspace(0, n - 1, 30).round().astype(int)]
                       for g, n in cells["group"].value_counts(sort=False).items()])
pick.sort()
X, cells = X[:, pick], cells.iloc[pick].reset_index(drop=True)   # every gene is kept (zero rows too)

hb, dd = rfindsc.default_biosets(), rfindsc.load_dedonno_biosets()
panels = {"DAM": hb["DAM"], "Glucocorticoid": dd["Glucocorticoid"], "G2M": hb["G2M"]}
# the panels are written in full so that R scores exactly the same thing: both lists, `eligible` (pair universe, 0.3.5)
# and whether the panel is ranked (G2M is a gene group, used whole)
rows = [(nm, side, g) for nm, p in panels.items() for side in ("up", "down") for g in p[side]["gene"]]
rows += [(nm, "eligible", g) for nm, p in panels.items() for g in (p.get("eligible") or [])]
pd.DataFrame(rows, columns=["panel", "side", "gene"]).to_csv(OUT / "panels.tsv", sep="\t", index=False)
pd.DataFrame({"panel": list(panels), "ranked": [bool(p.get("ranked", True)) for p in panels.values()]}).to_csv(
    OUT / "panels_ranked.tsv", sep="\t", index=False)

cfg = rfindsc.rfindsc_config(n_cell=100, n_cores=1)
res = rfindsc.score_cells(X, panels, cfg, return_signatures=True, genes=genes, cells=cells["cell"],
                          baseline_groups=cells["group"])

with gzip.open(OUT / "counts.mtx.gz", "wb") as fh:
    sio.mmwrite(fh, sp.coo_matrix(X.astype(np.int64)), field="integer")
(OUT / "genes.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")
cells.to_csv(OUT / "cells.tsv", sep="\t", index=False)
S = res["scores"]
with open(OUT / "scores.tsv", "w", encoding="utf-8") as fh:
    fh.write("\t".join(["cell"] + list(S.columns)) + "\n")
    for cell, row in zip(S.index, S.to_numpy()):
        fh.write("\t".join([str(cell)] + [f"{v:.17g}" for v in row]) + "\n")
(OUT / "universe.txt").write_text("\n".join(res["universe"]) + "\n", encoding="utf-8")
with gzip.open(OUT / "signatures.tsv.gz", "wt", encoding="utf-8") as fh:   # 1-based, as in R
    fh.write("cell\tside\tidx\n")
    for side in ("up_idx", "down_idx"):
        for cell, v in zip(S.index, res["signatures"][side]):
            fh.write(f"{cell}\t{side}\t{','.join(str(i + 1) for i in v)}\n")
(OUT / "README.txt").write_text(
    f"Written by data-raw/make_py_grouped_fixture.py with rfindsc {rfindsc.__version__} (Python), "
    f"numpy {np.__version__}, scipy {scipy.__version__}, pandas {pd.__version__}.\n"
    "rfindsc_config(n_cell=100, n_cores=1); baseline_groups = cells.tsv$group.\n", encoding="utf-8")
print(f"{X.shape[0]} genes x {X.shape[1]} cells, {len(res['universe'])} universe genes, panels {list(S.columns)}")
