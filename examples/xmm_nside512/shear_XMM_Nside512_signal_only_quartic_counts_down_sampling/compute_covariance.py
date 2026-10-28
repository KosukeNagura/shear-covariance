#!/usr/bin/env python3
"""Thinned entry point; share the validated full covariance implementation."""
import sys
from pathlib import Path
import runpy
common=Path(__file__).resolve().parents[1]/'shear_XMM_Nside512_signal_only_quartic_counts'
sys.path.insert(0,str(common))
if '--mode' not in sys.argv: sys.argv += ['--mode','thinned']
if '--mode' in sys.argv and sys.argv[sys.argv.index('--mode')+1]!='thinned':
    raise SystemExit('This entry point is for thinned mode only.')
runpy.run_path(str(common/'compute_covariance.py'),run_name='__main__')
