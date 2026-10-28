import os
import sys
import numpy as np
import healpy as hp

from config_custom_10bins import (
    NSIDE, WORK_ROOT, THETA_EDGES_RAD, NBINS,
    pair_meta_path, pair_chunk_dir, ensure_dirs
)

def angular_distance(vec1, vec2):
    c = np.dot(vec1, vec2)
    c = np.clip(c, -1.0, 1.0)
    return np.arccos(c)

def find_bin_index(angle, theta_bin_edges):
    b = np.searchsorted(theta_bin_edges, angle, side="right") - 1
    if 0 <= b < len(theta_bin_edges) - 1:
        return b
    return -1

def main():
    if len(sys.argv) != 4:
        raise ValueError("Usage: python pair_chunk_worker.py <chunk_id> <i_start> <i_end>")

    chunk_id = int(sys.argv[1])
    i_start = int(sys.argv[2])
    i_end = int(sys.argv[3])

    ensure_dirs()
    os.makedirs(pair_chunk_dir(), exist_ok=True)

    with np.load(pair_meta_path()) as d:
        pix_idx = d["pix_idx"]
        weights = d["weights"]
        theta_vecs = d["theta_vecs"]
        pix_to_local = d["pix_to_local"]
        nside_loaded = int(d["nside"])

    assert nside_loaded == NSIDE

    pair_lists = [[] for _ in range(NBINS)]
    n_active = len(pix_idx)
    i_end = min(i_end, n_active)
    theta_max = THETA_EDGES_RAD[-1]

    print(f"chunk_id={chunk_id}, i_start={i_start}, i_end={i_end}, n_active={n_active}", flush=True)

    for i in range(i_start, i_end):
        if (i - i_start) % 500 == 0:
            print(f"[chunk {chunk_id}] i={i}/{i_end}", flush=True)

        vec_i = theta_vecs[i]
        neighbor_pix = hp.query_disc(NSIDE, vec_i, theta_max, inclusive=False)

        for pix_j in neighbor_pix:
            j = pix_to_local[pix_j]
            if j <= i:
                continue

            ang = angular_distance(theta_vecs[i], theta_vecs[j])
            b = find_bin_index(ang, THETA_EDGES_RAD)
            if b >= 0:
                pair_lists[b].append((i, j))

    for b in range(NBINS):
        arr = np.asarray(pair_lists[b], dtype=np.int32)
        outpath = os.path.join(pair_chunk_dir(), f"pairs_bin_{b:02d}_chunk_{chunk_id:04d}.npy")
        np.save(outpath, arr)
        print(f"saved {outpath}, shape={arr.shape}", flush=True)

    print(f"chunk {chunk_id} finished.", flush=True)

if __name__ == "__main__":
    main()