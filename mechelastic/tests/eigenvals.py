#!/usr/bin/env python

from numpy import linalg as LA


def positive_evals(cnew):
    """This method checks the postivity of the eigenvalues
    of a matrix."""

    print("\nEigen Values of the matrix CNEW:")
    evals = LA.eigvals(cnew)
    print(evals)
    check = 0
    for i in range(len(evals)):
        if evals[i] > 0.0:
            pass
        else:
            check = 1
    if check == 1:
        print(
            "ATTENTION: One or more eigen values are negative indicating elastic instability."
        )
    if check == 0:
        print("All eigen values are positive indicating elastic stability.")

    return not(bool(check))
