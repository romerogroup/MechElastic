#!/usr/bin/env python

import numpy as np


def positive_evals(cnew, verbose=True):
    """This method checks the postivity of the eigenvalues
    of a matrix."""
    to_print = ""
    to_print += ("\nEigen Values of the matrix CNEW:\n")
    evals = np.linalg.eigvals(cnew)
    to_print += ("%s" % list(np.around(np.array(evals), 3)))+"\n"
    check = 0
    for i in range(len(evals)):
        if evals[i] > 0.0:
            pass
        else:
            check = 1
    if check == 1:
        to_print += (
            "ATTENTION: One or more eigen values are negative indicating elastic instability.\n"
        )
    if check == 0:
        to_print += ("All eigen values are positive indicating elastic stability.\n")
    if verbose:
        print(to_print)
    return not(bool(check))
