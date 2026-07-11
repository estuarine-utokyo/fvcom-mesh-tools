import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from fvcom_mesh_tools.io import read_fort14

m = read_fort14('outputs/pipeline_v6s/tokyo_bay_v6s_raw.14')
fig, ax = plt.subplots(figsize=(10, 12))
ax.triplot(m.nodes[:, 0], m.nodes[:, 1], m.elements, lw=0.25,
           color="steelblue")
ax.set_aspect(1 / np.cos(np.deg2rad(35.4)))
ax.set_title(f"v6s RAW NP={len(m.nodes):,}")
fig.savefig('outputs/figures/v6s_raw_full.png', dpi=170,
            bbox_inches="tight")
print("saved", flush=True)
