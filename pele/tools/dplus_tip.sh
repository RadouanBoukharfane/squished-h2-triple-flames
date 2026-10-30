#!/bin/bash
# Delta+ at the tip height (active_control.height of the reactive case) for finished cold runs COLDF-<case>,
# from the fit Delta+^2 = a + b Y+ of delta_profile.py. Usage: tools/dplus_tip.sh COLDF-A [COLDF-B ...]
# environment providing fextract (AMReX Tools/Plotfile) must be loaded
cd "${RUNS:?set RUNS to the directory holding the run directories}"
for cf in "$@"; do
  c=${cf#COLDF-}
  fit=$(python "$(dirname "$0")"/delta_profile.py $cf 2>&1 | grep "fit:")
  a=$(echo "$fit" | awk '{print $4}'); b=$(echo "$fit" | awk '{print $6}')
  h=$(grep -oP '^active_control.height = \K\S+' $c/input.case)
  w=$(grep -oP '^geometry.prob_hi\s*=\s*\K\S+' $c/input.case)
  python3 -c "import math;dL=317e-6*$w/0.075;Y=$h/dL;print('%-22s a=%-7s b=%-6s Y+=%6.1f Delta+=%.2f'%('$c','$a','$b',Y,math.sqrt($a+$b*Y)))"
done
