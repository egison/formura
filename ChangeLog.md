# Changelog for formura
## unreleased

- Temporal blocking on domains with walls: `grid_per_block` and
  `temporal_blocking_interval` are accepted together with `mirror` and
  `fixed` boundaries.  A walled axis stays anchored: its interior is
  placed at `sleeve*interval` in the floor with that many cells of halo on
  either side, every block overwrites the ghost cells that fall into its
  buffer with the boundary values before each sub-step kernel, and the
  interior is back at the start of the floor when it is copied back.
  Periodic axes keep the shifting frame.
  `test/tb-boundary.sh` checks bit equality with the non-blocked program
  for mixed boundaries, several intervals and block sizes, and sleeves 1
  and 2, including the grid coordinates read inside `step`.
- Decomposed runs (`mpi_shape` with more than one rank) on domains with
  walls, with and without temporal blocking.  A walled axis exchanges
  halos on both sides with its neighbor ranks; the rank across a wall of
  the domain is `MPI_PROC_NULL`, and the ranks at the walls fill the ghost
  cells from the boundary condition (`pos_<axis>` in `Formura_Navi` gives
  the rank position).  `test/mpi-boundary.sh` checks bit equality with the
  single-rank program for decompositions of one to three axes, with and
  without blocking; it needs `mpicc` and `mpirun`.
- The blocked step with walls overlaps communication with computation as
  the all-periodic step does: the blocks whose reads stay above the low
  halos (every axis in its upper range) run after only the slabs from
  beside and above have been placed, while the slabs from below are still
  in flight; the remaining blocks wait for them.
- Reject temporal blocking when a periodic axis is shorter than the
  one-sided halo `2*sleeve*temporal_blocking_interval` of a blocked step:
  the halo is copied from the neighbor's grid, and a shorter axis silently
  produced wrong values before.

## version 2.3.2

- Add install.sh

## version 2.3.1

- Improve compilation performance

## version 2.3

- Run also without MPI
- Implement LoadIndex
- Fix `to_pos` functions
- `Formura_Init` function initializes MPI and global data
- Add `Formura_Finalize` function
- Add `Formura_Custom_Init` function
- Add `space_interval_x` field to `Formura_Navi` struct
- Add `total_grid_x` field to `Formura_Navi` struct
- Support `first_step`
- Support `filter`

## version 2.2

- Support OpenMP

## version 2.1

- Add no blocking mode
- Add the name of the global data structure

## version 2.0

- Change the temporal blocking form
- Change the config format in yaml
- Fix bug on MPI
- Fix bug on temporal blocking
