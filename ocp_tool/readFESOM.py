import numpy as np
import os
import warnings
import sys
import time
import math
np.set_printoptions(threshold=sys.maxsize)

def read_fesom_grid(griddir, rot=False, rot_invert=False, rot_abg=None, threeD=True, remove_empty_lev=False, read_boundary=True,
                    reorder_ccw=True, maxmaxneigh=12, findneighbours_maxiter=10, repeatlastpoint=True, onlybaryc=False,
                    omitcoastnds=False, calcpolyareas=True, Rearth=6371000, basicreadonly=False, fesom2=True, verbose=True):


    def convert_elements(arr):
        new_arr = []
        for elem in arr:
            if math.isnan(elem):
                new_arr.append(-1)
            else:
                new_arr.append(int(elem))
        return new_arr


    def barycenter(lon=None, lat=None, x=None, y=None, z=None, weights=None, rm_na=True):
        rad = np.pi / 180

        if lon is None:
            if x is None or y is None or z is None:
                raise ValueError("Please provide lon, lat, and optionally weights.")

            if rm_na:
                if weights is None:
                    inds = ~np.isnan(x) & ~np.isnan(y) & ~np.isnan(z)
                    x, y, z = np.array(x)[inds], np.array(y)[inds], np.array(z)[inds]
                else:
                    inds = ~np.isnan(x) & ~np.isnan(y) & ~np.isnan(z) & ~np.isnan(weights)
                    x, y, z, weights = np.array(x)[inds], np.array(y)[inds], np.array(z)[inds], np.array(weights)[inds]

            lon = np.arctan2(y, x) * 180 / np.pi
            lat = np.arcsin(z) * 180 / np.pi

        else:
            if lat is None:
                raise ValueError("Please provide lon and lat.")

            if rm_na:
                if weights is None:
                    inds = ~np.isnan(lon) & ~np.isnan(lat)
                    lon, lat = np.array(lon)[inds], np.array(lat)[inds]
                else:
                    inds = ~np.isnan(lon) & ~np.isnan(lat) & ~np.isnan(weights)
                    lon, lat, weights = np.array(lon)[inds], np.array(lat)[inds], np.array(weights)[inds]

            lon = np.array(lon) * rad
            lat = np.array(lat) * rad
            x = np.cos(lat) * np.cos(lon)
            y = np.cos(lat) * np.sin(lon)
            z = np.sin(lat)

        if weights is None:
            x_mean = np.mean(x)
            y_mean = np.mean(y)
            z_mean = np.mean(z)
        else:
            weights = weights / np.sum(weights)
            x_mean = np.sum(x * weights)
            y_mean = np.sum(y * weights)
            z_mean = np.sum(z * weights)

        dist = np.sqrt(x_mean ** 2 + y_mean ** 2 + z_mean ** 2)
        x_mean = x_mean / dist
        y_mean = y_mean / dist
        z_mean = z_mean / dist
        lon_result = np.arctan2(y_mean, x_mean) * 180 / np.pi
        lat_result = np.arcsin(z_mean) * 180 / np.pi

        return lon_result, lat_result



    def checkposition(a, b, c, ccw_defined_as=1):
        if len(a) != 2 or len(b) != 2 or len(c) != 2:
            raise ValueError("a, b, and c must be lon-lat vectors of length 2!")

        if abs(ccw_defined_as) != 1:
            raise ValueError("ccw_defined_as must be one of -1 and 1!")

        def lonlat2xyz(lon_lat):
            lon, lat = np.deg2rad(lon_lat[0]), np.deg2rad(lon_lat[1])
            x = np.cos(lat) * np.cos(lon)
            y = np.cos(lat) * np.sin(lon)
            z = np.sin(lat)
            return np.array([x, y, z])

        def crossvec(u, v):
            x = u[1] * v[2] - u[2] * v[1]
            y = u[2] * v[0] - u[0] * v[2]
            z = u[0] * v[1] - u[1] * v[0]
            return np.array([x, y, z])

        a_xyz = lonlat2xyz(a)
        b_xyz = lonlat2xyz(b)
        c_xyz = lonlat2xyz(c)

        alpha = b_xyz - a_xyz
        beta = c_xyz - a_xyz

        alpha_cross_beta = crossvec(alpha, beta)
        alpha_cross_beta_dot_a = np.dot(alpha_cross_beta, a_xyz)

        return int(np.sign(alpha_cross_beta_dot_a) * ccw_defined_as)

    def rotate(lon, lat, abg, invert=False):
        abg_rad = np.deg2rad(abg)
        rot_matrix = np.array([[np.cos(abg_rad[0]) * np.cos(abg_rad[1]), -np.sin(abg_rad[0]), np.cos(abg_rad[0]) * np.sin(abg_rad[1])],
                               [np.sin(abg_rad[0]) * np.cos(abg_rad[1]), np.cos(abg_rad[0]), np.sin(abg_rad[0]) * np.sin(abg_rad[1])],
                               [-np.sin(abg_rad[1]), 0, np.cos(abg_rad[1])]])
        if invert:
            rot_matrix = np.linalg.inv(rot_matrix)
        xyz = np.column_stack((np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(lon)),
                              np.cos(np.deg2rad(lat)) * np.sin(np.deg2rad(lon)),
                              np.sin(np.deg2rad(lat))))
        xyz_rotated = np.dot(xyz, rot_matrix.T)
        lon_rotated = np.arctan2(xyz_rotated[:, 1], xyz_rotated[:, 0]) * 180 / np.pi
        lat_rotated = np.arcsin(xyz_rotated[:, 2]) * 180 / np.pi
        return lon_rotated, lat_rotated, xyz_rotated[:, 0], xyz_rotated[:, 1], xyz_rotated[:, 2]

    def fill_equidist(lon, lat, nfill):
        lon_endpoints = np.array([lon[-1], lon[0]])
        lat_endpoints = np.array([lat[-1], lat[0]])
        lon_fill = np.interp(np.linspace(0, 1, nfill + 2), [0, 1], lon_endpoints)
        lat_fill = np.interp(np.linspace(0, 1, nfill + 2), [0, 1], lat_endpoints)
        return lon_fill[1:-1], lat_fill[1:-1]


    def triag_area(lon, lat):
        if len(lon) != 3 or len(lat) != 3:
            raise ValueError("lon and lat must be vectors of length 3!")

        ang = np.zeros(3)

        for i in range(3):
            alpha = lon[i] + 90
            beta = 90 - lat[i]
            gamma = 0
            rot_res = rotate(lon, lat, alpha, beta, gamma, return_xyz=True)
            x = rot_res['x']
            y = rot_res['y']
            a = np.array([x[i], y[i]])
            b = np.array([x[(i + 1) % 3], y[(i + 1) % 3]])
            c = np.array([x[(i + 2) % 3], y[(i + 2) % 3]])
            ab = b - a
            ac = c - a

            with np.errstate(invalid='ignore'):
                ang[i] = np.arccos((np.dot(ab / np.linalg.norm(ab), ac / np.linalg.norm(ac))))

        if np.sum(np.isnan(ang)) > 0:
            return 0
        else:
            return np.sum(ang) - np.pi


    def rotate(lon, lat, alpha, beta, gamma, return_xyz=False, invert=False):
        rad = np.pi / 180
        a = alpha * rad
        b = beta * rad
        c = gamma * rad
        lon = lon * rad
        lat = lat * rad

        rotmat = np.zeros((3, 3))
        rotmat[0, 0] = np.cos(c) * np.cos(a) - np.sin(c) * np.cos(b) * np.sin(a)
        rotmat[0, 1] = np.cos(c) * np.sin(a) + np.sin(c) * np.cos(b) * np.cos(a)
        rotmat[0, 2] = np.sin(c) * np.sin(b)
        rotmat[1, 0] = -np.sin(c) * np.cos(a) - np.cos(c) * np.cos(b) * np.sin(a)
        rotmat[1, 1] = -np.sin(c) * np.sin(a) + np.cos(c) * np.cos(b) * np.cos(a)
        rotmat[1, 2] = np.cos(c) * np.sin(b)
        rotmat[2, 0] = np.sin(b) * np.sin(a)
        rotmat[2, 1] = -np.sin(b) * np.cos(a)
        rotmat[2, 2] = np.cos(b)

        if invert:
            rotmat = np.linalg.inv(rotmat)

        x = np.cos(lat) * np.cos(lon)
        y = np.cos(lat) * np.sin(lon)
        z = np.sin(lat)

        x_rot = rotmat[0, 0] * x + rotmat[0, 1] * y + rotmat[0, 2] * z
        y_rot = rotmat[1, 0] * x + rotmat[1, 1] * y + rotmat[1, 2] * z
        z_rot = rotmat[2, 0] * x + rotmat[2, 1] * y + rotmat[2, 2] * z
        z_rot[z_rot > 1] = 1
        z_rot[z_rot < -1] = -1

        lon_rot = np.arctan2(y_rot, x_rot) / rad
        np.seterr(invalid='ignore')  # Suppress warnings for invalid values
        lat_rot = np.arcsin(z_rot) / rad
        np.seterr(invalid='warn')  # Restore warnings for invalid values

        if return_xyz:
            return {'lon': lon_rot, 'lat': lat_rot, 'x': x_rot, 'y': y_rot, 'z': z_rot}
        else:
            return {'lon': lon_rot, 'lat': lat_rot}



    def read_aux3d_out(file_path):
        with open(file_path, 'r') as file:
            lines = file.readlines()
        Nlev = int(lines[0])
        depth_bounds = np.array([float(line.strip()) for line in lines[1:Nlev+2]]) * -1
        depth = (depth_bounds[:-1] + depth_bounds[1:]) / 2
        return Nlev, depth, depth_bounds

    def read_nlvls_out(file_path):
        with open(file_path, 'r') as file:
            lines = file.readlines()
        depth_lev = np.array([int(line.strip()) for line in lines]) - 1
        return depth_lev

    def read_nod3d_out(file_path):
        with open(file_path, 'r') as file:
            lines = file.readlines()
        depth = np.array([float(line.strip()) for line in lines[4::5]]) * -1
        return depth

    def read_elem2d_out(file_path):
        with open(file_path, 'r') as file:
            lines = file.readlines()
        Ne = int(lines[0])
        elem = np.array([list(map(int, line.strip().split())) for line in lines[1:]])
        return Ne, elem

    def read_nod2d_out(file_path):
        with open(file_path, 'r') as file:
            lines = file.readlines()
        N = int(lines[0])
        data = np.array([line.strip().split() for line in lines[1:]], dtype=float)
        lon_orig, lat_orig, coast = data[:, 1], data[:, 2], data[:, 3]
        return N, lon_orig, lat_orig, coast


    def find_neighbors(elem, maxmaxneigh=12, reverse=True, verbose=False, max_iter=10):
        if np.any(np.isnan(elem)):
            raise ValueError("'elem' must not contain NaNs.")

        N = np.max(elem)
        Ne = elem.shape[0]
        neighmat = np.full((N, maxmaxneigh), np.nan)
        Nneigh = np.zeros(N, dtype=int)
        barmat = np.full((N, maxmaxneigh), np.nan)
        iscomplete = np.full(N, False)
        iekdone = np.full((Ne, 3), False)
        niter = 0
        completed = True

        while np.sum(iekdone) < Ne * 3:
            niter += 1
            if niter > max_iter:
                warnings.warn("Some elements could not be arranged in order due to multi-domain nodes! Returned neighbourhood information is incomplete.")
                completed = False
                break

            if verbose:
                print(f"Starting iteration {niter}...")
            
            for ie in range(Ne):
                if np.all(iekdone[ie, :]):
                    continue
            
                for k in range(3):
                    if iekdone[ie, k]:
                        continue

                    i = elem[ie, k]-1

                    if iscomplete[i]:
                        raise ValueError("Ups! Trying to add neighbors to a node labeled complete!")

                    neigh1 = elem[ie, (k + 1) % 3]
                    neigh2 = elem[ie, (k + 2) % 3]

                    if np.isnan(neighmat[i, 0]):
                        barmat[i, 0] = ie
                        neighmat[i, 0] = neigh1
                        neighmat[i, 1] = neigh2
                        Nneigh[i] = 2
                        iekdone[ie, k] = True
                    else:
                        found1 = np.any(neighmat[i, :Nneigh[i]] == neigh1)
                        found2 = np.any(neighmat[i, :Nneigh[i]] == neigh2)

                        if found1 and found2:
                            if verbose:
                                print("Found both, node complete.")
                            barmat[i, Nneigh[i]-1] = ie
                            iscomplete[i] = True
                            iekdone[ie, k] = True
                        else:
                            if Nneigh[i] == maxmaxneigh:
                                raise ValueError("Ups! maxmaxneigh is insufficient!")

                            if found1:
                                if verbose:
                                    print("Found 1.")
                                neighmat[i, Nneigh[i]] = neigh2
                                barmat[i, Nneigh[i]-1] = ie
                                Nneigh[i] += 1
                                iekdone[ie, k] = True

                            elif found2:
                                if verbose:
                                    print("Found 2.")
                                neighmat[i, 1:Nneigh[i] + 1] = neighmat[i, :Nneigh[i]]
                                neighmat[i, 0] = neigh1
                                barmat[i, 1:Nneigh[i] + 1] = barmat[i, :Nneigh[i]]
                                barmat[i, 0] = ie
                                Nneigh[i] += 1
                                iekdone[ie, k] = True
                            else:
                                if verbose:
                                    print("Found none, retry element in next iteration.")

        maxneigh = max(Nneigh)
        neighmat = neighmat[:, :maxneigh]
        barmat = barmat[:, :maxneigh]

        if reverse:
            print("Reversing order of neighbors")
            for i in range(N):
                if Nneigh[i] > 1:
                    neighmat[i, :Nneigh[i]] = neighmat[i, Nneigh[i] - 1::-1]
                    barmat[i, :Nneigh[i] - 1] = barmat[i, Nneigh[i] - 2::-1]

        # Calculate the average number of neighbors
        avg_num_neighbors = np.mean(Nneigh)

        return neighmat, barmat, iscomplete, Nneigh, completed, avg_num_neighbors

    ##########################################
    # END of function defs, START of program #
    ##########################################


    fun_call = None
    if basicreadonly:
        print("Reading only basic grid data without further computation of neighborhood etc.")
        if rot or threeD or reorder_ccw:
            print("Reading would be even faster with rot, reorder_ccw, and threeD all set to False")

    if not os.path.exists(os.path.join(griddir, "nod2d.out")) or not os.path.exists(os.path.join(griddir, "elem2d.out")):
        raise FileNotFoundError(f"Files nod2d.out and/or elem2d.out not found in {griddir}.")

    if threeD:
        if not os.path.exists(os.path.join(griddir, "aux3d.out")):
            raise FileNotFoundError("3D information (file aux3d.out) missing. To read 2D only, rerun with threeD=False")
        if not fesom2:
            if not os.path.exists(os.path.join(griddir, "nod3d.out")):
                raise FileNotFoundError("3D information (file nod3d.out) missing. To read 2D only, rerun with threeD=False. To read FESOM2 grid, set fesom2=True.")
        else:
            use_nlvls_out = False
            if not os.path.exists(os.path.join(griddir, "fesom.mesh.diag.nc")):
                if not os.path.exists(os.path.join(griddir, "nlvls.out")):
                    raise FileNotFoundError("3D information (file nlvls.out and fesom.mesh.diag.nc) missing. To read 2D only, rerun with threeD=False. To read FESOM1 grid, set fesom2=False.")
                use_nlvls_out = True
                warnings.warn("File fesom.mesh.diag.nc missing, nlvls.out is used to derive depths only for nodes, not for elements. To read 2D only, rerun with threeD=False. To read FESOM1 grid, set fesom2=False.")
            else:
                try:
                    import netCDF4
                except ImportError:
                    raise ImportError("Package 'netCDF4' is required to read FESOM2 grid data.")
    if verbose:
        print("reading node (grid point) coordinates and coast information ...")
        start_time = time.time()
    N, lon_orig, lat_orig, coast = read_nod2d_out(os.path.join(griddir, "nod2d.out"))
    if verbose:
        print(f"... done. grid contains {N} nodes of which {np.sum(coast)} are coastal (according to info in nod2d.out).")
        end_time = time.time()
        print(f"... execution Time:", round(end_time - start_time, 2), "seconds")

    if rot:
        if verbose:
            print("rotating grid ...")
            start_time = time.time()
        lon, lat, x, y, z = rotate(lon_orig, lat_orig, rot_abg, invert=rot_invert)
        if verbose:
            print("... done.")
            end_time = time.time()
            print(f"... execution Time:", round(end_time - start_time, 2), "seconds")

    else:
        lon, lat, x, y, z = lon_orig, lat_orig, np.cos(np.deg2rad(lat_orig)) * np.cos(np.deg2rad(lon_orig)), np.cos(np.deg2rad(lat_orig)) * np.sin(np.deg2rad(lon_orig)), np.sin(np.deg2rad(lat_orig))
    lon[lon > 180] -= 360
    lon[lon <= -180] += 360

    if verbose:
        print("reading neighbourhood (triangular elements) information ...")
        start_time = time.time()
    Ne, elem = read_elem2d_out(os.path.join(griddir, "elem2d.out"))
    if verbose:
        print(f"... done. grid contains {Ne} triangular elements.")
        end_time = time.time()
        print(f"... execution Time:", round(end_time - start_time, 2), "seconds")

    # Reorder clockwise triangular elements counterclockwise if specified
    if reorder_ccw:
        if verbose:
            start_time = time.time()
            print("reordering clockwise triangular elements counterclockwise ...")
        ord_c = 0
        for ie in range(Ne):
            a = np.array([lon_orig[elem[ie, 0] - 1], lat_orig[elem[ie, 0] - 1]])
            b = np.array([lon_orig[elem[ie, 1] - 1], lat_orig[elem[ie, 1] - 1]])
            c = np.array([lon_orig[elem[ie, 2] - 1], lat_orig[elem[ie, 2] - 1]])
            if checkposition(a, b, c) == -1:
                elem[ie, :] = elem[ie, ::-1]
                ord_c += 1
        if verbose:
            print(f"... done. {ord_c} of {Ne} elements reordered.")
            end_time = time.time()
            print(f"... execution Time:", round(end_time - start_time, 2), "seconds")

    N3D = None
    Nlev = None
    depth = None
    depth_bounds = None
    depth_lev = None
    elemdepth_lev = None
    boundary = None
    if threeD:
        if verbose:
            start_time = time.time()
            print("reading 3D information ...")
        Nlev, depth, depth_bounds = read_aux3d_out(os.path.join(griddir, "aux3d.out"))
        if fesom2:
            Nlev -= 1
            if use_nlvls_out:
                depth_lev = read_nlvls_out(os.path.join(griddir, "nlvls.out"))
            else:
                mesh_diag_fl = netCDF4.Dataset(os.path.join(griddir, "fesom.mesh.diag.nc"))
                depth_lev = mesh_diag_fl.variables["nlevels_nod2D"][:] - 1
                elemdepth_lev = mesh_diag_fl.variables["nlevels"][:] - 1
                mesh_diag_fl.close()
            if remove_empty_lev and np.max(depth_lev) < Nlev:
                if verbose:
                    print(f"removing {Nlev - np.max(depth_lev)} empty levels from data")
                Nlev = np.max(depth_lev)
                depth_bounds = depth_bounds[:Nlev + 1]
                depth = depth[:Nlev]
            N3D = np.sum(depth_lev)
        else:
            aux3d_mat = np.genfromtxt(os.path.join(griddir, "aux3d.out"), skip_header=1, dtype=int, missing_values="-999", usemask=True)
            depth_lev = np.repeat(Nlev, N)
            for lev in range(Nlev, 0, -1):
                isna_lev = aux3d_mat[:, lev - 1].mask
                if remove_empty_lev and np.sum(isna_lev) == N:
                    if verbose:
                        print(f"removing empty level {lev} from data")
                    Nlev -= 1
                    aux3d_mat = aux3d_mat[:, :Nlev]
                depth_lev[isna_lev] -= 1
            N3D = np.sum(~aux3d_mat.mask)

            if verbose:
                print("retrieving depth information from nod3d.out ...")
            depth = read_nod3d_out(os.path.join(griddir, "nod3d.out"))
            if len(depth) != Nlev:
                raise ValueError("data in aux3d.out is inconsistent with the number of depth levels; consider trying with 'remove_empty_lev=True'")
            depth_bounds = np.concatenate(([depth[0]], (depth[1:] + depth[:-1]) / 2, [depth[-1]]))
            if read_boundary:
                if verbose:
                    print("retrieving 'coast/bottom' information from nod3d.out ...")
                boundary = read_nod3d_out(os.path.join(griddir, "nod3d.out"))[4::5]
        if verbose:
            end_time = time.time()
            print(f"... done. Grid over all levels contains {N3D} elements.")
            print(f"... execution Time:", round(end_time - start_time, 2), "seconds")
    if basicreadonly:
        return {
            'N': N, 'Nlev': Nlev, 'N3D': N3D, 'lon': lon, 'lat': lat, 'elem': elem, 'coast': coast,
            'neighnodes': None, 'neighelems': None, 'stamppoly_lon': None, 'stamppoly_lat': None,
            'baryc_lon': None, 'baryc_lat': None, 'cellareas': None, 'elemareas': None,
            'depth': depth, 'depth_lev': depth_lev, 'boundary': boundary
        }


    if verbose:
        start_time = time.time()
        print("searching all neighbors of each node based on the triangular elements ...")
    neighnodes, neighelems, internal_nodes, Nneighs, elems_completed, all_elements_arranged = find_neighbors(elem, maxmaxneigh, reverse=False, max_iter=findneighbours_maxiter)

    if np.any(coast == internal_nodes):
        warnings.warn("coast information from nod2d.out seems to be corrupt, using diagnosed coast flag instead.")
        coast = ~internal_nodes
    badnodes = None if all_elements_arranged else np.arange(1, N + 1)[~elems_completed]
    if verbose:
        print(f"... done. number of neighbors ranges from {np.min(Nneighs)} to {np.max(Nneighs)} nodes and is {np.mean(Nneighs):.4f} on average.")
        end_time = time.time()
        print(f"... execution Time:", round(end_time - start_time, 2), "seconds")
    if badnodes is not None:
        warnings.warn(f"if 'findneighbours_maxiter' was not set too low, the grid contains {len(badnodes)} 'bad nodes'. consider increasing 'findneighbours_maxiter'. if the problem remains, the grid indeed contains bad nodes that should not exist in the first place. for such nodes only one part of the corresponding ocean patches will be returned by this function (which introduces a slight grid inaccuracy).")

    if verbose:
        start_time = time.time()
        print("determining which elements include coastal nodes ...")
    elemcoast = np.array([np.sum(coast[elem[ie] - 1]) > 1 for ie in range(Ne)])
    Nelemcoast = np.sum(elemcoast)
    if verbose:
        print(f"... done. grid features {Nelemcoast} elements that contain coastal nodes.")
        end_time = time.time()
        print(f"... execution Time:", round(end_time - start_time, 2), "seconds")

    if verbose:
        start_time = time.time()
        print("computing barycenters (centroids) for all triangular elements ...")
    baryc_lon = np.zeros(Ne)
    baryc_lat = np.zeros(Ne)
    for ie in range(Ne):
        elem_ie = elem[ie, :] - 1  # Adjust the indices
        lon_ie, lat_ie = barycenter(lon[elem_ie], lat[elem_ie], z[elem_ie])
        baryc_lon[ie] = lon_ie
        baryc_lat[ie] = lat_ie
    baryc_lon[baryc_lon > 180] -= 360
    baryc_lon[baryc_lon <= -180] += 360
    if verbose:
        print("... done.")
        end_time = time.time()
        print(f"... execution Time:", round(end_time - start_time, 2), "seconds")


    if verbose:
        start_time = time.time()
        print("generate 'stamp polygons' around each node ...")
    maxneighs = neighnodes.shape[1]
    maxNstamp = 2 * maxneighs
    stampmat_lon = np.full((N, maxNstamp), np.nan)
    stampmat_lat = np.full((N, maxNstamp), np.nan)
    Nstamp = np.full(N, np.nan)
    for i in range(N):
        Nstamp_i = 0
        for j in range(maxneighs):
            nn = neighnodes[i, j]
            for bla in range(18):
                value = stampmat_lat[0, bla]
                if np.isnan(value):
                    value = "nan"
            if np.isnan(nn):
                break
            if not onlybaryc or (coast[i] and (j == 0 or j == Nneighs[i] - 1)):
                # compute median of central node and neighbor node
                nn_index = int(nn) - 1  # Subtract 1 to correct the index
                lon_ij, lat_ij = barycenter([lon[i], lon[nn_index]], [lat[i], lat[nn_index]], [z[i], z[nn_index]])
                stampmat_lon[i, Nstamp_i] = lon_ij
                stampmat_lat[i, Nstamp_i] = lat_ij
                Nstamp_i += 1
            ne = neighelems[i, j]
            if np.isnan(ne):
                break
            stampmat_lon[i, Nstamp_i] = baryc_lon[int(ne)]
            stampmat_lat[i, Nstamp_i] = baryc_lat[int(ne)]
            Nstamp_i += 1
        if coast[i] and not omitcoastnds:
            stampmat_lon[i, Nstamp_i] = lon[i]
            stampmat_lat[i, Nstamp_i] = lat[i]
            Nstamp_i += 1
        Nstamp[i] = Nstamp_i


    if maxNstamp > int(np.max(Nstamp)):
        maxNstamp = int(np.max(Nstamp))
        stampmat_lon = stampmat_lon[:, :maxNstamp]
        stampmat_lat = stampmat_lat[:, :maxNstamp]

    for i in range(N):
        nstamp_i_int = int(Nstamp[i])
        if nstamp_i_int < maxNstamp:
            if repeatlastpoint:
                stampmat_lon[i, nstamp_i_int:maxNstamp] = stampmat_lon[i, nstamp_i_int - 1]
                stampmat_lat[i, nstamp_i_int:maxNstamp] = stampmat_lat[i, nstamp_i_int - 1]
            else:
                lon_endpoints = np.array([stampmat_lon[i, nstamp_i_int - 1], stampmat_lon[i, 0]])
                lat_endpoints = np.array([stampmat_lat[i, nstamp_i_int - 1], stampmat_lat[i, 0]])
                nfill = maxNstamp - nstamp_i_int
                lon_fill, lat_fill = fill_equidist(lon_endpoints, lat_endpoints, nfill)
                stampmat_lon[i, nstamp_i_int:maxNstamp] = lon_fill
                stampmat_lat[i, nstamp_i_int:maxNstamp] = lat_fill

    stampmat_lon[stampmat_lon > 180] -= 360
    stampmat_lon[stampmat_lon <= -180] += 360
    if verbose:
        print(f"... done. number of 'stamp polygon' vertices per node ranges from {int(np.min(Nstamp))} (before padding) to {maxNstamp} and is {np.mean(Nstamp):.4f} on average (before padding).")
        end_time = time.time()
        print(f"... execution Time:", round(end_time - start_time, 2), "seconds")

    cellareas, elemareas = None, None
    if calcpolyareas:
        if verbose:
            start_time = time.time()
            print("computing element and 'stamp polygon' areas ...")
        elemareas = np.zeros(Ne)
        cellareas = np.zeros(N)
        for ie in range(Ne):
            elemareas[ie] = triag_area(lon[elem[ie, :]-1], lat[elem[ie, :]-1])
        elemareas *= Rearth ** 2
        neighelems_masked = np.ma.masked_invalid(neighelems)
        for i in range(N):
            cellareas[i] = np.sum(elemareas[neighelems_masked[i, :]].astype(int), where=neighelems_masked)
        cellareas /= 3
        if verbose:
            print("... done.")
            end_time = time.time()
            print(f"... execution Time:", round(end_time - start_time, 2), "seconds")
    
    return {
        'N': N, 'Nelem': Ne, 'Nlev': Nlev, 'N3D': N3D, 'lon': lon, 'lat': lat, 'elem': elem, 'elemcoast': elemcoast, 'coast': coast,
        'neighnodes': neighnodes, 'neighelems': neighelems, 'Nneighs': Nneighs, 'Nstamp': Nstamp, 'stampmat_lon': stampmat_lon, 'stampmat_lat': stampmat_lat,
        'baryc_lon': baryc_lon, 'baryc_lat': baryc_lat, 'cellareas': cellareas, 'elemareas': elemareas,
        'depth': depth, 'depth_bounds': depth_bounds, 'depth_lev': depth_lev, 'boundary': boundary, 'elemdepth_lev': elemdepth_lev
    }

out = read_fesom_grid(griddir='/work/ab0246/a270092/input/fesom2/core2/', basicreadonly=False)
print(out)
