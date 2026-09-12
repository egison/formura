#!/bin/sh
# Temporal blocking on domains with walls must reproduce the non-blocked
# program cell for cell, bit for bit, and the grid coordinates read inside
# step must stay exact.
#
# The model mixes stencil radii (a pass-through variable next to radius-1
# updates), reads diagonal neighbors (so corner ghosts matter), and adds a
# coordinate-dependent source.  It runs with mirror, fixed, and periodic
# axes in two arrangements, without blocking and with several blocking
# intervals and block sizes (including blocks of exactly 2*sleeve*interval,
# where the halo occupies a whole block).  Every state variable of every
# cell is compared with the non-blocked run after every Formura_Forward;
# periodic axes drift under blocking, so cells are matched through
# n.offset_*.  Compiled with FP contraction off so that identical
# arithmetic yields identical bits.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"
compiler=${FORMURA:-$(cabal list-bin exe:formura)}
work=$(mktemp -d "${TMPDIR:-/tmp}/formura-tb-boundary.XXXXXX")
trap 'rm -rf "$work"' EXIT HUP INT TERM
cat > "$work/model.fmr" <<'FMR'
dimension :: 3
axes :: x,y,z
double :: dt = 0.05
double :: kappa = 0.3
extern function :: sin
begin function (u,v,w,cx,cy,cz) = init()
  double [] :: u, v, w, cx, cy, cz
  u[i,j,k] = sin(0.7*i + 0.3*j) + 0.1*k
  v[i,j,k] = 0.01*(i + 2*j + 3*k)
  w[i,j,k] = 1.5 + 0.02*i - 0.01*k
  cx[i,j,k] = 0
  cy[i,j,k] = 0
  cz[i,j,k] = 0
end function
begin function (u',v',w',cx',cy',cz') = step(u,v,w,cx,cy,cz)
  double [] :: u', v', w', cx', cy', cz'
  lap[i,j,k] = u[i-1,j,k] + u[i+1,j,k] + u[i,j-1,k] + u[i,j+1,k] + u[i,j,k-1] + u[i,j,k+1] + (-6)*u[i,j,k]
  u'[i,j,k] = u[i,j,k] + dt*kappa*lap[i,j,k] + dt*0.01*v[i,j,k] + dt*0.001*(dx*i + dy*j)
  v'[i,j,k] = 0.999*v[i,j,k] + 0.001*(u[i+1,j+1,k] + u[i-1,j-1,k+1] + w[i,j,k-1])
  w'[i,j,k] = w[i,j,k]
  cx'[i,j,k] = cx[i,j,k] + dx*i
  cy'[i,j,k] = cy[i,j,k] + dy*j
  cz'[i,j,k] = cz[i,j,k] + dz*k
end function
FMR
cat > "$work/dump.c" <<'C'
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "model.h"
/* Print every state variable of every cell in physical cell order after
   each of the requested forwards, and check the accumulated coordinates
   against time_step * to_pos_* for every cell. */
int main(int argc, char **argv) {
  Formura_Navi n;
  Formura_Init(&argc, &argv, &n);
  int forwards = atoi(argv[1]);
  long bad = 0;
  for (int f = 0; f < forwards; f++) {
    Formura_Forward(&n);
    printf("time_step %d\n", n.time_step);
    /* physical cell order: slot = cell - offset on a drifting axis */
    for (int ci = 0; ci < n.total_grid_x; ci++)
      for (int cj = 0; cj < n.total_grid_y; cj++)
        for (int ck = 0; ck < n.total_grid_z; ck++) {
          int i = (ci - n.offset_x + n.total_grid_x) % n.total_grid_x;
          int j = (cj - n.offset_y + n.total_grid_y) % n.total_grid_y;
          int k = (ck - n.offset_z + n.total_grid_z) % n.total_grid_z;
          printf("%d %d %d %a %a %a\n", ci, cj, ck,
                 formura_data.u[i][j][k], formura_data.v[i][j][k], formura_data.w[i][j][k]);
          double c[3] = {formura_data.cx[i][j][k], formura_data.cy[i][j][k], formura_data.cz[i][j][k]};
          double p[3] = {to_pos_x(i, n), to_pos_y(j, n), to_pos_z(k, n)};
          for (int a = 0; a < 3; a++)
            if (fabs(c[a] - n.time_step * p[a]) > 1e-9) {
              if (bad < 5)
                fprintf(stderr, "step=%d cell=(%d,%d,%d) axis=%d coordinate=%.17g expected=%.17g\n",
                        n.time_step, ci, cj, ck, a, c[a], n.time_step * p[a]);
              bad++;
            }
        }
  }
  Formura_Finalize();
  if (bad) { fprintf(stderr, "%ld mismatching coordinates\n", bad); return 1; }
  return 0;
}
C
# the same program with a radius-2 term, so that the sleeve is 2
sed 's/^  w.\[i,j,k\] = w\[i,j,k\]$/  w'"'"'[i,j,k] = w[i,j,k] + 0.001*(u[i-2,j,k] + u[i+2,j,k] + u[i,j,k+2] + u[i,j-2,k])/' \
  "$work/model.fmr" > "$work/model2.fmr"
grep -q 'u\[i+2,j,k\]' "$work/model2.fmr"
run_case() {
  # $1: boundary line, $2: label, $3: interval (0 = no blocking), $4: block sizes, $5: model file
  boundary=$1; label=$2; interval=$3; blocks=$4; model=$5
  dir="$work/$label"
  mkdir -p "$dir"
  cp "$work/$model" "$dir/model.fmr"
  cp "$work/dump.c" "$dir/"
  {
    echo "length_per_node: [1.2, 1.6, 2.0]"
    echo "grid_per_node: [12, 16, 20]"
    echo "boundary: $boundary"
    if [ "$interval" != 0 ]; then
      echo "grid_per_block: $blocks"
      echo "temporal_blocking_interval: $interval"
    fi
  } > "$dir/model.yaml"
  (cd "$dir" && "$compiler" model.fmr > generate.log)
  ${CC:-cc} -std=c11 -O1 -ffp-contract=off -I"$dir" "$dir/dump.c" "$dir/model.c" -lm -o "$dir/dump"
}
total=24
cells=$((12 * 16 * 20))
check_pair() {
  # $1: tag of the plain run, $2: label of the blocked run, $3: interval
  # keep only the records the blocked run also produced (its steps are
  # multiples of the interval), then require bit equality
  awk -v nt="$3" '/^time_step/ { keep = ($2 % nt == 0); if (keep) print; next } { if (keep) print }' \
    "$work/$1.txt" > "$work/$2-expected.txt"
  expected_lines=$(( (total / $3) * (cells + 1) ))
  [ "$(wc -l < "$work/$2-expected.txt")" -eq "$expected_lines" ]
  [ "$(wc -l < "$work/$2.txt")" -eq "$expected_lines" ]
  if ! cmp -s "$work/$2-expected.txt" "$work/$2.txt"; then
    echo "temporal blocking differs from the plain run: $2" >&2
    diff "$work/$2-expected.txt" "$work/$2.txt" | head -20 >&2
    exit 1
  fi
}
for boundary in "[mirror, periodic, fixed 0.5]" "[fixed 0.0, mirror, mirror]"; do
  case "$boundary" in
    "[mirror, periodic, fixed 0.5]") tag=mpf ;;
    *) tag=fmm ;;
  esac
  run_case "$boundary" "$tag-plain" 0 "" model.fmr
  "$work/$tag-plain/dump" $total > "$work/$tag-plain.txt"
  # interval:blocks; the floor is grid + 2*interval per axis (sleeve 1)
  for spec in "1:[7, 9, 11]" "2:[8, 5, 6]" "2:[4, 4, 4]" "3:[9, 11, 13]" "4:[10, 12, 14]" "4:[20, 8, 28]"; do
    interval=${spec%%:*}
    blocks=${spec#*:}
    label="$tag-tb$interval-$(echo "$blocks" | tr -d '[], ')"
    run_case "$boundary" "$label" "$interval" "$blocks" model.fmr
    "$work/$label/dump" $((total / interval)) > "$work/$label.txt"
    check_pair "$tag-plain" "$label" "$interval"
  done
  # sleeve 2: the floor is grid + 4*interval per axis
  run_case "$boundary" "$tag-s2-plain" 0 "" model2.fmr
  "$work/$tag-s2-plain/dump" $total > "$work/$tag-s2-plain.txt"
  for spec in "1:[8, 10, 12]" "1:[4, 4, 4]" "2:[10, 12, 14]" "3:[12, 14, 16]"; do
    interval=${spec%%:*}
    blocks=${spec#*:}
    label="$tag-s2-tb$interval-$(echo "$blocks" | tr -d '[], ')"
    run_case "$boundary" "$label" "$interval" "$blocks" model2.fmr
    "$work/$label/dump" $((total / interval)) > "$work/$label.txt"
    check_pair "$tag-s2-plain" "$label" "$interval"
  done
done
printf 'Formura temporal blocking with walls: ok\n'
