import os
import glob
import numpy as np

from config_custom_10bins import (
    NBINS, pair_meta_path, pair_chunk_dir, pair_merged_dir,
    pair_counts_path, ensure_dirs
)

def main():
    ensure_dirs()
    os.makedirs(pair_merged_dir(), exist_ok=True)

    pair_counts = np.zeros(NBINS, dtype=np.int64)

    for b in range(NBINS):
        files = sorted(glob.glob(os.path.join(pair_chunk_dir(), f"pairs_bin_{b:02d}_chunk_*.npy")))
        arrays = []

        for f in files:
            arr = np.load(f)
            if arr.size > 0:
                arrays.append(arr)

        if len(arrays) == 0:
            merged = np.empty((0, 2), dtype=np.int32)
        else:
            merged = np.vstack(arrays).astype(np.int32, copy=False)

        outpath = os.path.join(pair_merged_dir(), f"pairs_bin_{b:02d}.npy")
        np.save(outpath, merged)
        pair_counts[b] = len(merged)

        print(f"bin {b}: merged shape = {merged.shape}, saved to {outpath}", flush=True)

    np.save(pair_counts_path(), pair_counts)
    print(f"saved pair_counts to {pair_counts_path()}", flush=True)
    print("merge finished.", flush=True)

if __name__ == "__main__":
    main()