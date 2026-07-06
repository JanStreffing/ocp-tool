"""Reusable scipy-style interpolation for a fixed source -> target grid.

``scipy.interpolate.griddata`` rebuilds, on *every* call, both the Delaunay
triangulation of the SOURCE points (``method='linear'``) and the point-location
of the TARGET points within it. When many value arrays are interpolated between
the SAME source and target grids -- e.g. the ~60 vertical levels of a 3D CO2
field, or a series of 2D fields on one grid -- *both* are redundant, and the
point-location dominates (each ``LinearNDInterpolator.__call__`` re-runs
``find_simplex`` over all target points).

``ReusableGridInterp`` builds the triangulation once AND the target-point
barycentric weights once, so each subsequent field is just a gather + weighted
sum (an einsum) -- fast at every resolution, and with no process/cluster
parallelism. The linear weights are the same barycentric weights scipy's
``LinearNDInterpolator`` uses (``tri.transform``); results match ``griddata`` to
interpolation/round-off precision (and bit-for-bit after GRIB packing).
"""
import numpy as np
from scipy.spatial import Delaunay, cKDTree


class ReusableGridInterp:
    """Delaunay + KD-tree over a fixed source-point set, plus a cached
    barycentric plan for a fixed target-point set, reused across value arrays."""

    def __init__(self, source_points):
        self.source_points = np.asarray(source_points, dtype=float)
        self._tri = None
        self._tree = None
        self._plan_key = None      # id() of the target array the plan was built for
        self._plan = None          # (vertices, weights, outside, outside_nearest_idx)

    def _triangulation(self):
        if self._tri is None:
            self._tri = Delaunay(self.source_points)
        return self._tri

    def _kdtree(self):
        if self._tree is None:
            self._tree = cKDTree(self.source_points)
        return self._tree

    def _linear_plan(self, target_points):
        """Locate the target points in the source triangulation ONCE and cache
        the per-target simplex vertices + barycentric weights. Also caches, for
        target points outside the source convex hull, their nearest source index
        (that set is value-independent, so it too is computed once)."""
        tp = np.asarray(target_points, dtype=float)
        if self._plan_key != id(tp):
            tri = self._triangulation()
            ndim = tp.shape[1]
            simplex = tri.find_simplex(tp)
            outside = simplex < 0
            # Barycentric coords via the affine transform scipy stores per simplex
            # (same construction LinearNDInterpolator uses internally).
            T = tri.transform[simplex]
            bary = np.einsum('nij,nj->ni', T[:, :ndim, :], tp - T[:, ndim, :])
            weights = np.concatenate([bary, (1.0 - bary.sum(axis=1))[:, None]], axis=1)
            vertices = tri.simplices[simplex]        # (n_target, ndim+1)
            near = self._kdtree().query(tp[outside])[1] if outside.any() else None
            self._plan = (vertices, weights, outside, near)
            self._plan_key = id(tp)
        return self._plan

    def linear(self, values, target_points, fill_value=np.nan):
        """Linear interpolation; points outside the source hull get ``fill_value``
        (== griddata's ``fill_value``, default NaN)."""
        vertices, weights, outside, _ = self._linear_plan(target_points)
        values = np.asarray(values)
        out = np.einsum('nj,nj->n', values[vertices], weights)
        out[outside] = fill_value
        return out

    def nearest(self, values, target_points):
        """== griddata(source_points, values, target_points, method='nearest')."""
        _, idx = self._kdtree().query(np.asarray(target_points, dtype=float))
        return np.asarray(values)[idx]

    def linear_with_nearest_fill(self, values, target_points):
        """Linear, with the out-of-hull points filled by nearest neighbour --
        the common ``griddata(linear)`` + ``griddata(nearest, nan_mask)`` idiom,
        with both the triangulation and the point-location computed once."""
        vertices, weights, outside, near = self._linear_plan(target_points)
        values = np.asarray(values)
        out = np.einsum('nj,nj->n', values[vertices], weights)
        if near is not None:
            out[outside] = values[near]
        return out


def _same_points(a, b):
    """Cheap guard: are two source-point arrays identical (so a cached plan/
    triangulation stays valid)? O(N) but tiny next to an interpolation."""
    a = np.asarray(a); b = np.asarray(b)
    return a.shape == b.shape and np.array_equal(a, b)
