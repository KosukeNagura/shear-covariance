"""Spherical, signal-only shear covariance. NumPy vectorized pair-pair blocks.

All components use the west/north basis, g=(healpy Q, -healpy U), matching
existing TreeCorr simulations. Spin-2 rotations follow the handwritten notes.
No extra global factor 2 is applied to unordered i<j pairs.
"""
import hashlib
import os
from pathlib import Path
import numpy as np

T = np.array([[1., 1.], [1., -1.]])


def exclusive_save(path, **data):
    """Create a final result, never replace one."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as f:
        np.savez(f, **data)


def checkpoint_save(path, **data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    with tmp.open('wb') as f:
        np.savez(f, **data)
    os.replace(tmp, path)


def file_hash(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def local_bases(vectors):
    v = np.asarray(vectors, float)
    longitude = np.arctan2(v[:, 1], v[:, 0])
    east = np.column_stack((-np.sin(longitude), np.cos(longitude), np.zeros(len(v))))
    north = np.cross(v, east)
    return -east, north


def rotations(angle):
    c, s = np.cos(2 * angle), np.sin(2 * angle)
    return np.stack((-c, -s, s, -c), axis=-1).reshape(angle.shape + (2, 2))


def connection_rotations(vectors, u, v):
    """Broadcast endpoint indices; return rotations into their great-circle basis."""
    east, north = local_bases(vectors)
    pu, pv = vectors[u], vectors[v]
    phi_u = np.arctan2(np.sum(pv * north[u], axis=-1), np.sum(pv * east[u], axis=-1))
    phi_v = np.arctan2(np.sum(pu * north[v], axis=-1), np.sum(pu * east[v], axis=-1))
    return rotations(phi_u), rotations(phi_v)


def xi_from_cl(theta, cl_ee, cl_bb=None):
    """Exact spherical d^ell_(2,+/-2) sums via stable Jacobi recurrence.

    d22 = cos(theta/2)^4 P_(ell-2)^(0,4)(cos theta)
    d2-2 = sin(theta/2)^4 P_(ell-2)^(4,0)(cos theta).
    C_ell is the actual spin-2 E/B spectrum passed to healpy (no conversion).
    """
    theta = np.asarray(theta, float)
    ee = np.asarray(cl_ee, float)
    bb = np.zeros_like(ee) if cl_bb is None else np.asarray(cl_bb, float)
    if ee.ndim != 1 or bb.ndim != 1 or ee.shape != bb.shape or len(ee) < 3:
        raise ValueError('E/B spectra must have identical length >= 3')
    if np.any(~np.isfinite(theta)) or np.any((theta < 0) | (theta > np.pi)):
        raise ValueError("theta must be finite and in [0, pi]")
    ee, bb = ee.copy(), bb.copy()
    ee[:2] = 0.; bb[:2] = 0.  # array index is ell; never shift the spectrum
    x = np.cos(theta)
    result = []
    for alpha, beta, spectrum, factor in [
        (0., 4., ee + bb, np.cos(theta / 2)**4),
        (4., 0., ee - bb, np.sin(theta / 2)**4),
    ]:
        older = np.ones_like(x)
        total = 5 * spectrum[2] * older
        if len(ee) > 3:
            previous = .5 * ((alpha - beta) + (alpha + beta + 2) * x)
            total += 7 * spectrum[3] * previous
            for n in range(2, len(ee) - 2):
                z = 2*n + alpha + beta
                current = ((z - 1) * (z * (z - 2) * x + alpha**2 - beta**2) * previous
                           - 2 * (n + alpha - 1) * (n + beta - 1) * z * older
                           ) / (2 * n * (n + alpha + beta) * (z - 2))
                total += (2*(n+2) + 1) * spectrum[n+2] * current
                older, previous = previous, current
        result.append(factor * total / (4 * np.pi))
    ell = np.arange(2, len(ee))
    zero_plus = np.sum((2*ell + 1)*(ee[2:] + bb[2:]))/(4*np.pi)
    result[0] = np.where(theta == 0, zero_plus, result[0])
    result[1] = np.where(theta == 0, 0., result[1])
    return tuple(result)


def signal_covariance(vectors, theta_grid, xip, xim, xi_zero, chunk=128):
    """S[u,v,component_u,component_v] in fixed local west/north bases."""
    n = len(vectors)
    out = np.empty((n, n, 2, 2), dtype=float)
    for start in range(0, n, chunk):
        u = np.arange(start, min(n, start + chunk))[:, None]
        v = np.arange(n)[None, :]
        dot = np.sum(vectors[u] * vectors[v], axis=-1)
        theta = np.arccos(np.clip(dot, -1, 1))
        if np.any((dot < -1 + 1e-12) & (u != v)):
            raise ValueError('Antipodal pixels need a separately specified transport convention')
        ru, rv = connection_rotations(vectors, u, v)
        pp = np.interp(theta, theta_grid, xip)
        mm = np.interp(theta, theta_grid, xim)
        diag = np.stack(((pp + mm)/2, (pp - mm)/2), axis=-1)
        block = np.einsum('...ai,...a,...aj->...ij', ru, diag, rv, optimize=True)
        same = np.broadcast_to(u == v, theta.shape)
        # At zero separation the direction is undefined; use the isotropic limit.
        block[same] = np.eye(2) * xi_zero / 2
        out[start:start + len(u)] = block
    return out


def pair_geometry(pairs, vectors, weights):
    pairs = np.asarray(pairs, dtype=np.int64)
    i, j = pairs.T
    ri, rj = connection_rotations(vectors, i, j)
    return pairs, weights[i] * weights[j], ri, rj


def rotated_cross(signal, u, v, ru, rv):
    return np.einsum('...ac,...cd,...bd->...ab', ru, signal[u, v], rv, optimize=True)


def kernel_block(left, right, signal):
    """Return (Nleft,Nright,2,2), last axes are X,Y = +,-."""
    pairs, _, ri, rj = left
    other, _, rk, rl = right
    i, j = pairs[:, 0, None], pairs[:, 1, None]
    k, l = other[None, :, 0], other[None, :, 1]
    hik = rotated_cross(signal, i, k, ri[:, None], rk[None])
    hjl = rotated_cross(signal, j, l, rj[:, None], rl[None])
    w = hik * hjl
    del hik, hjl
    hil = rotated_cross(signal, i, l, ri[:, None], rl[None])
    hjk = rotated_cross(signal, j, k, rj[:, None], rk[None])
    w += hil * hjk
    return np.einsum('xa,...ab,yb->...xy', T, w, T, optimize=True)


def geometry_slice(g, sl):
    return tuple(x[sl] for x in g)


def raw_pair_sum(left, right, signal, right_chunk=2000, left_chunk=0):
    """Same pair-by-pair matrix strategy as convergence; chunk RHS by default.

    left_chunk=0 keeps all LHS pairs (legacy strategy, can require much RAM).
    Set a positive left_chunk to tile both axes without changing the estimator.
    """
    n, m = len(left[0]), len(right[0])
    result = np.zeros((2, 2))
    if not n or not m:
        return result
    if right_chunk < 1 or left_chunk < 0:
        raise ValueError('Invalid chunk size')
    for a in range(0, n, left_chunk or n):
        ga = geometry_slice(left, slice(a, a + (left_chunk or n)))
        for b in range(0, m, right_chunk):
            gb = geometry_slice(right, slice(b, b + right_chunk))
            kernel = kernel_block(ga, gb, signal)
            result += np.einsum('i,j,ijxy->xy', ga[1], gb[1], kernel, optimize=True)
    return result


def identical_sum(geometry, signal, chunk=100000):
    """Compute q_r^2 K_rr for each sampled edge, never r-by-r matrices."""
    total = np.zeros((2, 2))
    for start in range(0, len(geometry[0]), chunk):
        pairs, q, ri, rj = geometry_slice(geometry, slice(start, start + chunk))
        i, j = pairs.T
        hii = rotated_cross(signal, i, i, ri, ri)
        hjj = rotated_cross(signal, j, j, rj, rj)
        hij = rotated_cross(signal, i, j, ri, rj)
        hji = rotated_cross(signal, j, i, rj, ri)
        w = hii*hjj + hij*hji
        total += np.einsum('n,xa,nab,yb->xy', q*q, T, w, T, optimize=True)
    return total


def measure_grid(g, geometries, denominators):
    result = np.empty((2, len(geometries)))
    for b, (pairs, q, ri, rj) in enumerate(geometries):
        u = np.einsum('nab,nb->na', ri, g[pairs[:, 0]])
        v = np.einsum('nab,nb->na', rj, g[pairs[:, 1]])
        sums = np.sum(q[:, None] * u*v, axis=0) / denominators[b]
        result[:, b] = T @ sums
    return result.reshape(-1)  # all + bins, then all - bins


def bernoulli_pairs(pairs, probability, rng, chunk=1000000):
    if probability == 1:
        return np.asarray(pairs)
    selected = []
    for start in range(0, len(pairs), chunk):
        piece = pairs[start:start + chunk]
        selected.append(piece[rng.random(len(piece)) < probability])
    return np.concatenate(selected) if selected else np.empty((0, 2), np.int32)


def assemble_blocks(blocks):
    # blocks[a,b,X,Y] -> joint[X*n + a, Y*n + b], X,Y = 0 (+), 1 (-)
    n = blocks.shape[0]
    return blocks.transpose(2, 0, 3, 1).reshape(2*n, 2*n)
