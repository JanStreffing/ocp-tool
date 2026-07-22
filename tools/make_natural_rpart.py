#!/usr/bin/env python3
"""Emit a synthetic 'natural order' FESOM rpart.out (identity permutation).

The oasis_reorder_tool maps one FESOM partition's PE-contiguous node ordering
to another's, using dist_<NPES>/rpart.out at each end. ocp-tool's rmp files are
in *natural* global order (feom built straight from nod2d.out), which is exactly
the single-PE identity partition: npes=1, PE0 owns all nodes in order 1..nod2D.

Feeding this file as DIST_OLD (with a real dist_<N> as DIST_NEW) makes the
reorder tool sort rmp addresses from natural order into that partition's order.

Format (see reorder_oasis.F90 read_rpart_mapping):
  line 1:      npes                     -> 1
  line 2:      ncount(1:npes)           -> nod2D
  next nod2D:  mapping(n), one per line -> n   (identity; one value PER LINE,
               because the reader does `read(u,*) mapping(n)` per record)

Usage:
    python make_natural_rpart.py <mesh_dir> [out_dir]
      mesh_dir : dir containing nod2d.out (to read nod2D)
      out_dir  : where to write rpart.out (default: <mesh_dir>/dist_natural)
"""
import sys
from pathlib import Path


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    mesh_dir = Path(sys.argv[1])
    out_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else mesh_dir / "dist_natural"

    nod2d_file = mesh_dir / "nod2d.out"
    with open(nod2d_file) as f:
        nod2d = int(f.readline().split()[0])
    print(f"nod2D = {nod2d} (from {nod2d_file})")

    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "rpart.out"

    # Buffered write of 1..nod2D, one per line, in chunks to keep memory flat.
    chunk = 1_000_000
    with open(out_file, "w") as f:
        f.write(f"{1}\n{nod2d}\n")
        n = 1
        while n <= nod2d:
            hi = min(n + chunk - 1, nod2d)
            f.write("\n".join(map(str, range(n, hi + 1))))
            f.write("\n")
            n = hi + 1
    print(f"Wrote identity rpart.out ({nod2d} nodes) -> {out_file}")


if __name__ == "__main__":
    main()
