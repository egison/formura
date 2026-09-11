#!/bin/sh
# Grid coordinates read inside step must equal to_pos_* for every cell,
# with and without temporal blocking, and after the shifting frame has
# wrapped around the periodic domain.
#
# Before the fix, the kernel saw a coordinate larger than the true one by
# sleeve * sub-step under temporal blocking, and the raw index was never
# reduced modulo the grid size, so cells saw y + L once n->offset_* wrapped.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"
compiler=${FORMURA:-$(cabal list-bin exe:formura)}
work=$(mktemp -d "${TMPDIR:-/tmp}/formura-coordinate-probe.XXXXXX")
trap 'rm -rf "$work"' EXIT HUP INT TERM
cat > "$work/model.fmr" <<'FMR'
dimension :: 3
axes :: x,y,z
double :: dt = 0.1 * dx
begin function (u,cx,cy,cz) = init()
  double [] :: u, cx, cy, cz
  u[i,j,k] = 0
  cx[i,j,k] = 0
  cy[i,j,k] = 0
  cz[i,j,k] = 0
end function
begin function (u',cx',cy',cz') = step(u,cx,cy,cz)
  double [] :: u', cx', cy', cz'
  u'[i,j,k] = u[i,j,k] + dt * (u[i-1,j,k] + u[i+1,j,k] + u[i,j-1,k] + u[i,j+1,k] + u[i,j,k-1] + u[i,j,k+1] + (-6) * u[i,j,k]) / dx**2
  cx'[i,j,k] = dx * i
  cy'[i,j,k] = dy * j
  cz'[i,j,k] = dz * k
end function
FMR
cat > "$work/check.c" <<'C'
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "model.h"
/* After every Formura_Forward, the coordinates the kernel used at its last
   sub-step must equal to_pos_* of the cell that now holds the result. */
int main(int argc, char **argv) {
  Formura_Navi n;
  Formura_Init(&argc, &argv, &n);
  int forwards = atoi(argv[1]);
  long bad = 0;
  for (int f = 0; f < forwards; f++) {
    Formura_Forward(&n);
    for (int i = n.lower_x; i < n.upper_x; i++)
      for (int j = n.lower_y; j < n.upper_y; j++)
        for (int k = n.lower_z; k < n.upper_z; k++) {
          double used[3] = {formura_data.cx[i][j][k], formura_data.cy[i][j][k], formura_data.cz[i][j][k]};
          double truth[3] = {to_pos_x(i, n), to_pos_y(j, n), to_pos_z(k, n)};
          for (int a = 0; a < 3; a++)
            if (fabs(used[a] - truth[a]) > 1e-12) {
              if (bad < 5)
                fprintf(stderr, "step=%d cell=(%d,%d,%d) axis=%d kernel=%.17g to_pos=%.17g\n",
                        n.time_step, i, j, k, a, used[a], truth[a]);
              bad++;
            }
        }
  }
  Formura_Finalize();
  if (bad) { fprintf(stderr, "%ld mismatching coordinates\n", bad); return 1; }
  return 0;
}
C
for interval in 0 1 2 3 4; do
  cat > "$work/model.yaml" <<'YAML'
length_per_node: [1.0, 1.0, 1.0]
grid_per_node: [16, 16, 16]
YAML
  if [ "$interval" != 0 ]; then
    block=$((16 + 2 * interval))
    cat >> "$work/model.yaml" <<YAML
grid_per_block: [$block, $block, $block]
temporal_blocking_interval: $interval
YAML
  fi
  (cd "$work" && "$compiler" model.fmr)
  ${CC:-cc} -std=c11 -O1 -fsanitize=undefined -I"$work" "$work/check.c" "$work/model.c" -lm -o "$work/check"
  # 40 steps: the frame offset wraps around the 16-cell domain more than twice.
  "$work/check" $((40 / (interval > 0 ? interval : 1)))
done
printf 'Formura grid coordinates inside step with and without temporal blocking: ok\n'
