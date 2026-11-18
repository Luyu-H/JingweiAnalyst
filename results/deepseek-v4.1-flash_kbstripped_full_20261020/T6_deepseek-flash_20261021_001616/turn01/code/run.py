
"""Re-run this turn's computation exactly as the agent executed it.

    cd /home/lorraine/projects/2025-jingwei_QA_GMD && python <bundle>/code/run.py

The agent runs generated code with every registered tool already in scope and
the selected dataset ids bound to `selected_datasets`, so this script rebuilds
that namespace instead of importing the code as a module. Figures are written to
figures/rerun/ so they can be compared against the originals in figures/.
"""
import sys
from pathlib import Path

BUNDLE = Path(__file__).resolve().parent.parent
REPO = Path(r"/home/lorraine/projects/2025-jingwei_QA_GMD")
sys.path.insert(0, str(REPO))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

# The same rcParams the agent's sandbox sets: different font metrics shift the
# figures by a pixel or two and an exact reproduction then looks like a mismatch.
plt.rcParams["font.sans-serif"] = ["Noto Sans CJK SC"]
plt.rcParams["axes.unicode_minus"] = False

from src.tools.registry import TOOLS, auto_import_tools
auto_import_tools()

namespace = {name: meta["func"] for name, meta in TOOLS.items()}
namespace.update({"np": np, "pd": pd, "xr": xr, "plt": plt})
namespace["selected_datasets"] = ['woa23', 'cmip6_historical']


with xr.set_options(keep_attrs=True):
    exec(compile((BUNDLE / "code" / "compute.py").read_text(), "compute.py", "exec"), namespace)
    visualize = BUNDLE / "code" / "visualize.py"
    if visualize.exists():
        exec(compile(visualize.read_text(), "visualize.py", "exec"), namespace)

out = BUNDLE / "figures" / "rerun"
out.mkdir(parents=True, exist_ok=True)
for index, number in enumerate(plt.get_fignums(), 1):
    plt.figure(number).savefig(out / ("figure_%d.png" % index), dpi=300, bbox_inches="tight")
print("re-ran %s; %d figure(s) written to %s" % (BUNDLE.name, len(plt.get_fignums()), out))
