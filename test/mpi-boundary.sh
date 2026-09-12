#!/bin/sh
# A decomposed run (mpi_shape with more than one rank) on a domain with
# walls must reproduce the single-rank program cell for cell, bit for bit,
# with and without temporal blocking.  Every rank dumps its cells with
# global indices (n.offset_* is the rank's origin on a walled axis and
# includes the drift on a periodic one); the gathered records are sorted
# and compared with the single-rank plain run.  Skipped when mpicc or
# mpirun is unavailable.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"
if ! command -v mpicc >/dev/null 2>&1 || ! command -v mpirun >/dev/null 2>&1; then
  printf 'Formura decomposed runs with walls: skipped (no mpicc/mpirun)\n'
  exit 0
fi
compiler=${FORMURA:-$(cabal list-bin exe:formura)}
mpirun_args=${MPIRUN_ARGS:---bind-to none}
work=$(mktemp -d "${TMPDIR:-/tmp}/formura-mpi-boundary.XXXXXX")
trap 'rm -rf "$work"' EXIT HUP INT TERM
# the model of tb-boundary.sh: mixed radii, diagonal reads, coordinates
sed -n '/^cat > "$work\/model.fmr" <<.FMR.$/,/^FMR$/p' "$root/test/tb-boundary.sh" | sed '1d;$d' > "$work/model.fmr"
grep -q '^dimension :: 3' "$work/model.fmr"
cat > "$work/dump.c" <<'C'
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "model.h"
/* One record per cell and forward: global cell, every state variable, and
   the coordinate check against time_step * to_pos_*. */
int main(int argc, char **argv) {
  Formura_Navi n;
  Formura_Init(&argc, &argv, &n);
  int forwards = atoi(argv[1]);
  /* one file per rank: the standard output of the ranks would interleave */
  char path[512];
  snprintf(path, sizeof path, "%s.%d", argv[2], n.my_rank);
  FILE *out = fopen(path, "w");
  if (!out) { perror(path); return 2; }
  long bad = 0;
  for (int f = 0; f < forwards; f++) {
    Formura_Forward(&n);
    for (int i = n.lower_x; i < n.upper_x; i++)
      for (int j = n.lower_y; j < n.upper_y; j++)
        for (int k = n.lower_z; k < n.upper_z; k++) {
          int ci = (i + n.offset_x) % n.total_grid_x;
          int cj = (j + n.offset_y) % n.total_grid_y;
          int ck = (k + n.offset_z) % n.total_grid_z;
          fprintf(out, "%d %d %d %d %a %a %a\n", n.time_step, ci, cj, ck,
                  formura_data.u[i][j][k], formura_data.v[i][j][k], formura_data.w[i][j][k]);
          double c[3] = {formura_data.cx[i][j][k], formura_data.cy[i][j][k], formura_data.cz[i][j][k]};
          double p[3] = {to_pos_x(i, n), to_pos_y(j, n), to_pos_z(k, n)};
          for (int a = 0; a < 3; a++)
            if (fabs(c[a] - n.time_step * p[a]) > 1e-9) {
              if (bad < 5)
                fprintf(stderr, "rank=%d step=%d cell=(%d,%d,%d) axis=%d coordinate=%.17g expected=%.17g\n",
                        n.my_rank, n.time_step, ci, cj, ck, a, c[a], n.time_step * p[a]);
              bad++;
            }
        }
  }
  if (fclose(out)) { perror(path); return 2; }
  Formura_Finalize();
  if (bad) { fprintf(stderr, "%ld mismatching coordinates\n", bad); return 1; }
  return 0;
}
C
# grid 12 x 16 x 20 in total; per-rank grid and length follow the shape
build() {
  # $1 label, $2 boundary, $3 mpi shape (px py pz), $4 interval, $5 block sizes
  label=$1; boundary=$2; shape=$3; interval=$4; blocks=$5
  set -- $shape; px=$1; py=$2; pz=$3
  dir="$work/$label"; mkdir -p "$dir"
  cp "$work/model.fmr" "$work/dump.c" "$dir/"
  {
    printf 'length_per_node: [%s, %s, %s]\n' "$(python3 -c "print(1.2/$px)")" "$(python3 -c "print(1.6/$py)")" "$(python3 -c "print(2.0/$pz)")"
    printf 'grid_per_node: [%d, %d, %d]\n' $((12 / px)) $((16 / py)) $((20 / pz))
    printf 'boundary: %s\n' "$boundary"
    printf 'mpi_shape: [%d, %d, %d]\n' "$px" "$py" "$pz"
    if [ "$interval" != 0 ]; then
      printf 'grid_per_block: %s\n' "$blocks"
      printf 'temporal_blocking_interval: %s\n' "$interval"
    fi
  } > "$dir/model.yaml"
  (cd "$dir" && "$compiler" model.fmr > generate.log)
  mpicc -std=c11 -O1 -ffp-contract=off -I"$dir" "$dir/dump.c" "$dir/model.c" -lm -o "$dir/dump"
}
run() {
  # $1 label, $2 ranks, $3 forwards
  mpirun $mpirun_args -n "$2" "$work/$1/dump" "$3" "$work/$1/cells"
  cat "$work/$1"/cells.* | sort -n -k1,1 -k2,2 -k3,3 -k4,4 > "$work/$1.txt"
}
total=24
cells=$((12 * 16 * 20))
for boundary in "[mirror, periodic, fixed 0.5]" "[fixed 0.0, mirror, mirror]"; do
  case "$boundary" in
    "[mirror, periodic, fixed 0.5]") tag=mpf ;;
    *) tag=fmm ;;
  esac
  build "$tag-ref" "$boundary" "1 1 1" 0 ""
  run "$tag-ref" 1 $total
  [ "$(wc -l < "$work/$tag-ref.txt")" -eq $((total * cells)) ]
  # label:shape:interval:blocks (blocks divide the per-rank floor grid + 2*interval)
  for spec in "nb-x:2 1 1:0:" "nb-yz:1 2 2:0:" "nb-xyz:2 2 2:0:" \
              "tb2-x:2 1 1:2:[5, 10, 12]" "tb2-y:1 2 1:2:[8, 6, 12]" "tb2-z:1 1 2:2:[8, 5, 7]" \
              "tb2-xyz:2 2 2:2:[5, 6, 7]" "tb4-xz:2 1 2:4:[14, 12, 9]" "tb3-xy:2 2 1:3:[12, 14, 26]"; do
    label="$tag-${spec%%:*}"; rest=${spec#*:}
    shape=${rest%%:*}; rest=${rest#*:}
    interval=${rest%%:*}; blocks=${rest#*:}
    set -- $shape; ranks=$(( $1 * $2 * $3 ))
    build "$label" "$boundary" "$shape" "$interval" "$blocks"
    forwards=$total
    [ "$interval" != 0 ] && forwards=$((total / interval))
    run "$label" "$ranks" "$forwards"
    if [ "$interval" = 0 ]; then
      cp "$work/$tag-ref.txt" "$work/$label-expected.txt"
    else
      awk -v nt="$interval" '$1 % nt == 0' "$work/$tag-ref.txt" > "$work/$label-expected.txt"
    fi
    [ "$(wc -l < "$work/$label.txt")" -eq "$(wc -l < "$work/$label-expected.txt")" ]
    if ! cmp -s "$work/$label-expected.txt" "$work/$label.txt"; then
      echo "decomposed run differs from the single-rank run: $label" >&2
      diff "$work/$label-expected.txt" "$work/$label.txt" | head -20 >&2
      exit 1
    fi
  done
done
printf 'Formura decomposed runs with walls: ok\n'
