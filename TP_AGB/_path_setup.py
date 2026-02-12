"""
Auto-generated path setup for shared python_helpers.
Import this at the top of any script that needs shared helpers:

    import _path_setup  # noqa: F401

This adds the repo root and python_helpers/ to sys.path.
"""
import os
import sys

_this_dir = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.normpath(os.path.join(_this_dir, os.pardir, ""))
_helpers = os.path.join(_repo_root, "python_helpers")

for _p in [_repo_root, _helpers]:
    if _p not in sys.path:
        sys.path.insert(0, _p)
