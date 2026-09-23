# Does the base depth file meet the r-factor its name claims?
#
# The finish step reported 374 frozen-pair edges over r = 0.2 on a patch that
# may not touch any of them.  A frozen pair is an edge between two RETAINED
# nodes, which is a base edge unless the patch invented it, so either the
# base does not meet 0.2 or the patch created 374 new edges between retained
# nodes.  This says which.
import sys
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.io.fvcom_native import read_fvcom_case

G = Path(sys.argv[1])
for dep in sys.argv[2:]:
    m = read_fvcom_case(G, Path(dep), None)
    e = np.unique(np.sort(np.vstack([m.elements[:, [0, 1]], m.elements[:, [1, 2]],
                                     m.elements[:, [2, 0]]]), axis=1), axis=0)
    h = m.depths
    r = np.abs(h[e[:, 0]] - h[e[:, 1]]) / (h[e[:, 0]] + h[e[:, 1]])
    print(f"{Path(dep).name}: depth {h.min():.3f}..{h.max():.3f} m, "
          f"{len(e):,} edges, max r {r.max():.4f}, "
          f"over 0.2: {int((r > 0.2).sum()):,} ({100 * (r > 0.2).mean():.2f} %), "
          f"p99 r {np.percentile(r, 99):.4f}")
