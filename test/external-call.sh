#!/bin/sh
# Multi-argument scalar external functions in ordinary Formura C output.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"
compiler=${FORMURA:-$(cabal list-bin exe:formura)}
work=$(mktemp -d "${TMPDIR:-/tmp}/formura-external-call.XXXXXX")
trap 'rm -rf "$work"' EXIT HUP INT TERM
cat > "$work/model.fmr" <<'FMR'
dimension :: 1
axes :: x
extern function :: atan2
extern function :: pow
begin function q = init()
  double [] :: q
  q[i] = atan2(1,1) + pow(2,3)
end function
begin function q' = step(q)
  double [] :: q'
  q'[i] = (q[i-1] + q[i+1])/2 + atan2(1,1) + pow(2,3)
end function
FMR
cat > "$work/check.c" <<'C'
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "model.h"
int main(int argc, char **argv) {
  Formura_Navi n;
  Formura_Init(&argc, &argv, &n);
  const double increment = atan2(1.0, 1.0) + pow(2.0, 3.0);
  for (int pass = 0; pass < 4; ++pass) {
    for (int i = n.lower_x; i < n.upper_x; ++i) {
      if (fabs(formura_data.q[i] - (n.time_step + 1)*increment) >= 1e-12)
        fprintf(stderr, "step=%d i=%d actual=%.17g expected=%.17g\n", n.time_step, i, formura_data.q[i], (n.time_step + 1)*increment);
      assert(fabs(formura_data.q[i] - (n.time_step + 1)*increment) < 1e-12);
    }
    Formura_Forward(&n);
  }
  Formura_Finalize();
  return 0;
}
C
for blocking in no yes; do
  cat > "$work/model.yaml" <<'YAML'
grid_per_node: [32]
length_per_node: [1]
YAML
  if [ "$blocking" = yes ]; then
    cat >> "$work/model.yaml" <<'YAML'
grid_per_block: [12]
temporal_blocking_interval: 2
YAML
  fi
  (cd "$work" && "$compiler" model.fmr)
  ${CC:-cc} -std=c11 -O1 -fsanitize=undefined -I"$work" "$work/check.c" "$work/model.c" -lm -o "$work/check"
  "$work/check"
done
printf 'Formura scalar external calls with and without temporal blocking: ok\n'
