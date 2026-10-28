"""Portable path/settings adapter for the established XMM Nside=512 scripts.

The pair worker and merger are copied from the existing XMM/Nside512 workflow.
Set the three environment variables below before running the metadata stage.
"""
import os
from pathlib import Path
import numpy as np

IZ = 4
NSIDE = 512
NSIMS = 1000
AREA = "XMM"
HERE = Path(__file__).resolve().parent
BASE_ROOT = str(HERE)
WORK_ROOT = os.environ.get("SHEAR_PAIR_WORK_ROOT", str(HERE))
CATALOG_PATH = os.environ.get("SHEAR_PAIR_CATALOG")
THETA_NPZ_PATH = os.environ.get("SHEAR_PAIR_THETA_NPZ")
if not CATALOG_PATH or not THETA_NPZ_PATH:
    raise RuntimeError("Set SHEAR_PAIR_CATALOG and SHEAR_PAIR_THETA_NPZ")

with np.load(THETA_NPZ_PATH) as _data:
    theta = _data["theta"]
    theta_edges = _data["theta_edges"]
THETA_EDGES_DEG = theta_edges
THETA_EDGES_RAD = np.deg2rad(THETA_EDGES_DEG)
THETA_DEG = np.sqrt(THETA_EDGES_DEG[:-1] * THETA_EDGES_DEG[1:])
THETA_RAD = np.deg2rad(THETA_DEG)
NBINS = len(THETA_DEG)


def ensure_dirs():
    for relative in (
        "log", "jobs_pair_chunks", "jobs_cov_blocks", "outputs/precompute",
        "outputs/precompute/pair_chunks", "outputs/precompute/pairs_merged",
        "outputs/cov_blocks",
    ):
        (Path(WORK_ROOT) / relative).mkdir(parents=True, exist_ok=True)


def pair_meta_path():
    return str(Path(WORK_ROOT) / "outputs/precompute/pair_meta.npz")


def pair_counts_path():
    return str(Path(WORK_ROOT) / "outputs/precompute/pair_counts.npy")


def pair_chunk_dir():
    return str(Path(WORK_ROOT) / "outputs/precompute/pair_chunks")


def pair_merged_dir():
    return str(Path(WORK_ROOT) / "outputs/precompute/pairs_merged")


def cov_block_dir():
    return str(Path(WORK_ROOT) / "outputs/cov_blocks")


def final_cov_path():
    return str(Path(WORK_ROOT) / f"outputs/covariance_XMM_Nside{NSIDE}_bin{IZ}_10bins.npz")
