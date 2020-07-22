#!/usr/bin/env python

import numpy as np


def check_symmetric(a, tol=1e-8):
    """Checks if a matrix is symmetric."""
    return np.allclose(a, a.T, atol=tol)
