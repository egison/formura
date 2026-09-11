#!/bin/sh
# Grid coordinates read inside step must equal to_pos_* for every cell at
# every sub-step, with and without temporal blocking, read directly or
# through a shifted intermediate array, and after the shifting frame has
# wrapped around the periodic domain.
#
# Before the fix, the kernel saw a coordinate off by a multiple of the
# sleeve under temporal blocking (up to 2*sleeve*(interval-1) cells at the
# first sub-step), and the raw index was never reduced modulo the grid size,
# so cells saw y + L once n->offset_* wrapped.
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
begin function (u,cx,cy,cz,sx,sy,sz) = init()
  double [] :: u, cx, cy, cz, sx, sy, sz
  u[i,j,k] = 0
  cx[i,j,k] = 0
  cy[i,j,k] = 0
  cz[i,j,k] = 0
  sx[i,j,k] = 0
  sy[i,j,k] = 0
  sz[i,j,k] = 0
end function
begin function (u',cx',cy',cz',sx',sy',sz') = step(u,cx,cy,cz,sx,sy,sz)
  double [] :: u', cx', cy', cz', sx', sy', sz'
  u'[i,j,k] = u[i,j,k] + dt * (u[i-1,j,k] + u[i+1,j,k] + u[i,j-1,k] + u[i,j+1,k] + u[i,j,k-1] + u[i,j,k+1] + (-6) * u[i,j,k]) / dx**2
  # the coordinate of the updated cell, accumulated over every step
  cx'[i,j,k] = cx[i,j,k] + dx * i
  cy'[i,j,k] = cy[i,j,k] + dy * j
  cz'[i,j,k] = cz[i,j,k] + dz * k
  # the same coordinate read through a shifted intermediate array
  mx[i,j,k] = dx * i + 0 * u[i,j,k]
  my[i,j,k] = dy * j + 0 * u[i,j,k]
  mz[i,j,k] = dz * k + 0 * u[i,j,k]
  sx'[i,j,k] = sx[i,j,k] + mx[i+1,j,k] - dx
  sy'[i,j,k] = sy[i,j,k] + my[i,j-1,k] + dy
  sz'[i,j,k] = sz[i,j,k] + (mz[i,j,k+1] + mz[i,j,k-1]) / 2
end function
FMR
cat > "$work/check.c" <<'C'
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "model.h"
/* A cell keeps its physical position while its index moves with the frame,
   so after every Formura_Forward the accumulated coordinate of a cell must
   equal time_step times to_pos_* of the cell; this checks every sub-step of
   a temporally blocked step, not only the last one.  The neighbor of a cell
   across the periodic boundary sits at position +-L, which the shifted
   reads of sx..sz see; the expected sums account for it. */
int main(int argc, char **argv) {
  Formura_Navi n;
  Formura_Init(&argc, &argv, &n);
  int forwards = atoi(argv[1]);
  long bad = 0;
  double length[3] = {n.length_x, n.length_y, n.length_z};
  double d[3] = {n.space_interval_x, n.space_interval_y, n.space_interval_z};
  for (int f = 0; f < forwards; f++) {
    Formura_Forward(&n);
    for (int i = n.lower_x; i < n.upper_x; i++)
      for (int j = n.lower_y; j < n.upper_y; j++)
        for (int k = n.lower_z; k < n.upper_z; k++) {
          double truth[3] = {to_pos_x(i, n), to_pos_y(j, n), to_pos_z(k, n)};
          double direct[3] = {formura_data.cx[i][j][k], formura_data.cy[i][j][k], formura_data.cz[i][j][k]};
          double shifted[3] = {formura_data.sx[i][j][k], formura_data.sy[i][j][k], formura_data.sz[i][j][k]};
          for (int a = 0; a < 3; a++) {
            double expect = n.time_step * truth[a];
            /* the shifted reads see the wrapped neighbor at the domain edge */
            double edge = 0;
            if (a == 0 && truth[0] + d[0] >= length[0]) edge = -length[0];
            if (a == 1 && truth[1] - d[1] < 0) edge = length[1];
            if (a == 2 && (truth[2] + d[2] >= length[2])) edge = -length[2] / 2;
            if (a == 2 && (truth[2] - d[2] < 0)) edge = length[2] / 2;
            double expect_shifted = expect + n.time_step * edge;
            if (fabs(direct[a] - expect) > 1e-9 || fabs(shifted[a] - expect_shifted) > 1e-9) {
              if (bad < 5)
                fprintf(stderr, "step=%d cell=(%d,%d,%d) axis=%d direct=%.17g shifted=%.17g expected=%.17g/%.17g\n",
                        n.time_step, i, j, k, a, direct[a], shifted[a], expect, expect_shifted);
              bad++;
            }
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
