#!/usr/bin/env python

import numpy as np
from numpy import linalg as LA


def stability_test(matrix, crystal_type, verbose=True):
    """This methods tests for the stability of a structure, then returns a boolean if ther structure is stable.

        References
        [1] Necessary and Sufficient Elastic Stability Conditions in Various Crystal Systems. Félix Mouhat and François-Xavier Coudert. Phys. Rev. B (2014)
        [2] Crystal Structures and Elastic Properties of Superhard IrN2 and IrN3 from First Principles. Zhi-jian Wu et al. Phys. Rev. B (2007)
    """

    c = np.copy(matrix)
    stable = True
    to_print = ""
    if crystal_type == "cubic":
        to_print += ("Cubic crystal system \n")
        to_print += (
            "Born stability criteria for the stability of cubic system are: Ref.[1]  \n"
        )
        to_print += ("(i) C11 - C12 > 0;    (ii) C11 + 2C12 > 0;   (iii) C44 > 0 \n ")

        # check (i)   keep in mind list starts with 0, so c11 is stored as c00
        if c[0][0] - c[0][1] > 0.0:
            to_print += ("Condition (i) satisfied.")
            condition1 = True
        else:
            to_print += ("Condition (i) NOT satisfied.")
            condition1 = False

        if c[0][0] + 2 * c[0][1] > 0.0:
            to_print += ("Condition (ii) satisfied.")
            condition2 = True
        else:
            to_print += ("Condition (ii) NOT satisfied.")
            condition2 = False

        if c[3][3] > 0.0:
            to_print += ("Condition (iii) satisfied.")
            condition3 = True
        else:
            to_print += ("Condition (iii) NOT satisfied.")
            condition3 = False

        conditions = (condition1, condition2, condition3)

        for condition in conditions:
            if condition == False:
                stable = False

        return stable

    if crystal_type == "hexagonal":
        to_print += ("Hexagonal crystal system \n")
        to_print += (
            "Born stability criteria for the stability of hexagonal system are: Ref.[1]  \n"
        )
        to_print += (
            "(i) C11 - C12 > 0;    (ii) 2*C13^2 < C33(C11 + C12);   (iii) C44 > 0 \n "
        )

        # check (i)   keep in mind list starts with 0, so c11 is stored as c00
        if c[0][0] - c[0][1] > 0.0:
            to_print += ("Condition (i) satisfied.")
            condition1 = True
        else:
            to_print += ("Condition (i) NOT satisfied.")
            condition1 = False

        if 2 * (c[0][2] * c[0][2]) < c[2][2] * (c[0][0] + c[0][1]):
            to_print += ("Condition (ii) satified.")
            condition2 = True
        else:
            to_print += ("Condition (ii) NOT satisfied.")
            condition2 = False

        if c[3][3] > 0.0:
            to_print += ("Condition (iii) satisfied.")
            condition3 = True
        else:
            to_print += ("Condition (iii) NOT satisfied.")
            condition3 = False
        conditions = (condition1, condition2, condition3)

        for condition in conditions:
            if condition == False:
                stable = False

    if crystal_type == "tetragonal":
        to_print += ("Tetragonal crystal system \n")
        to_print += (
            "Born stability criteria for the stability of Tetragonal system are: Ref.[1]  \n"
        )
        to_print += (
            "(i) C11 - C12 > 0;    (ii) 2*C13^2 < C33(C11 + C12);   (iii) C44 > 0; \n(iv) C66 > 0;    (v) 2C16^2 < C66*(C11-C12) \n "
        )

        # check (i)   keep in mind list starts with 0, so c11 is stored as c00
        if c[0][0] - c[0][1] > 0.0:
            to_print += ("Condition (i) is satisfied.")
            condition1 = True
        else:
            to_print += ("Condition (i) is NOT satisfied.")
            condition1 = False

        if 2 * (c[0][2] * c[0][2]) < c[2][2] * (c[0][0] + c[0][1]):
            to_print += ("Condition (ii) is satisfied.")
            condition2 = True
        else:
            to_print += ("Condition (ii) is NOT satisfied.")
            condition2 = False

        if c[3][3] > 0.0:
            to_print += ("Condition (iii) is satisfied.")
            condition3 = True
        else:
            to_print += ("Condition (iii) is NOT satisfied.")
            condition3 = False

        if c[5][5] > 0.0:
            to_print += ("Condition (iv) is satisfied.")
            condition4 = True
        else:
            to_print += ("Condition (iv) is NOT satisfied.")
            condition4 = False

        if 2 * c[0][5] * c[0][5] < c[5][5] * (c[0][0] - c[0][1]):
            to_print += ("Condition (v) is satisfied.")
            condition5 = True
        else:
            to_print += ("Condition (v) is NOT satisfied.")
            condition5 = False

        conditions = (condition1, condition2,
                      condition3, condition4, condition5)

        for condition in conditions:
            if condition == False:
                stable = False

    if crystal_type == "rhombohedral-1":
        to_print += ("Rhombohedral (class-1): point group: 3m, -3m and 32 \n")
        to_print += ("Born stability criteria for this class are: Ref.[1]  \n")
        to_print += (
            "(i) C11 - C12 > 0;    (ii) C13^2 < (1/2)*C33(C11 + C12);   (iii) C14^2 < (1/2)*C44*(C11-C12) = C44*C66; \n(iv)  C44 > 0; \n "
        )

        if c[0][0] - c[0][1] > 0.0:
            to_print += ("Condition (i) is satisfied.")
            condition1 = True
        else:
            to_print += ("Condition (i) is NOT satisfied.")
            condition1 = False

        if (c[0][2] * c[0][2]) < (0.5) * c[2][2] * (c[0][0] + c[0][1]):
            to_print += ("Condition (ii) is satisfied.")
            condition2 = True
        else:
            to_print += ("Condition (ii) is NOT satisfied.")
            condition2 = False

        if c[0][3] * c[0][3] < 0.5 * c[3][3] * (c[0][0] - c[0][1]):
            to_print += ("Condition (iii) is satisfied.")
            condition3 = True
        else:
            to_print += ("Condition (iii) is NOT satisfied.")
            condition3 = False

        if c[3][3] > 0.0:
            to_print += ("Condition (iv) is satisfied.")
            condition4 = True
        else:
            to_print += ("Condition (iv) is NOT satisfied.")
            condition4 = False

        conditions = (condition1, condition2, condition3, condition4)

        for condition in conditions:
            if condition == False:
                stable = False

    if crystal_type == "rhombohedral-2":
        to_print += ("Rhombohedral (class-2): i.e structures with point group: 3, -3 \n")
        to_print += (
            "Born stability criteria for the stability of Rhombohedral-1 class system are: Ref.[1]  \n"
        )
        to_print += (
            "(i) C11 - C12 > 0;    (ii) C13^2 < (1/2)*C33(C11 + C12);   (iii) C14^2 + C15^2 < (1/2)*C44*(C11-C12) = C44*C66; \n(iv)  C44 > 0;  Note: C15 is added.. \n "
        )

        if c[0][0] - c[0][1] > 0.0:
            to_print += ("Condition (i) is satisfied.")
            condition1 = True
        else:
            to_print += ("Condition (i) is NOT satisfied.")
            condition1 = False

        if c[0][2] * c[0][2] < (0.5) * c[2][2] * (c[0][0] + c[0][1]):
            to_print += ("Condition (ii) is satisfied.")
            condition2 = True
        else:
            to_print += ("Condition (ii) is NOT satisfied.")
            condition2 = False

        if c[0][3] * c[0][3] + c[0][4] * c[0][4] < 0.5 * c[3][3] * (c[0][0] - c[0][1]):
            to_print += ("Condition (iii) is satified.")
            condition3 = True
        else:
            to_print += ("Condition (iii) is NOT satisfied.")
            condition3 = False

        if c[3][3] > 0.0:
            to_print += ("Condition (iv) is satisfied.")
            condition4 = True
        else:
            to_print += ("Condition (iv) is NOT satisfied.")
            condition4 = False

        conditions = (condition1, condition2, condition3, condition4)

        for condition in conditions:
            if condition == False:
                stable = False

    if crystal_type == "orthorhombic":
        to_print += ("Orthorhombic crystal system.... \n")
        to_print += (
            "Born stability criteria for the stability of Orthorhombic systems are: Ref.[1]  \n"
        )
        to_print += (
            "(i) C11 > 0;   (ii) C11*C22 > C12^2;   (iii) C11*C22*C33 + 2C12*C13*C23 - C11*C23^2 - C22*C13^2 - C33*C12^2 > 0; \n(iv)  C44 > 0;   (v)  C55 > 0 ;   (vi)  C66 > 0 \n"
        )

        # check (i)   keep in mind list starts with 0, so c11 is stored as c00
        if c[0][0] > 0.0:
            to_print += ("Condition (i) is satisfied.")
            condition1 = True
        else:
            to_print += ("Condition (i) is NOT satisfied.")
            condition1 = False

        if c[0][0] * c[1][1] > c[0][1] * (c[0][1]):
            to_print += ("Condition (ii) is satisfied.")
            condition2 = True
        else:
            to_print += ("Condition (ii) is NOT satisfied.")
            condition2 = False

        if (
            c[0][0] * c[1][1] * c[2][2]
            + 2 * c[0][1] * c[0][2] * c[1][2]
            - c[0][0] * c[1][2] * c[1][2]
            - c[1][1] * c[0][2] * c[0][2]
            - c[2][2] * c[0][1] * c[0][1]
            > 0
        ):
            to_print += ("Condition (iii) is satisfied.")
            condition3 = True
        else:
            to_print += ("Condition (iii) is NOT satisfied.")
            condition3 = False

        if c[3][3] > 0.0:
            to_print += ("Condition (iv) is satisfied.")
            condition4 = True
        else:
            to_print += ("Condition (iv) is NOT satisfied.")
            condition4 = False

        if c[4][4] > 0.0:
            to_print += ("Condition (v) is satisfied.")
            condition5 = True
        else:
            to_print += ("Condition (v) is NOT satisfied.")
            condition5 = False

        if c[5][5] > 0.0:
            to_print += ("Condition (vi) is satisfied.")
            condition6 = True
        else:
            to_print += ("Condition (vi) is NOT satisfied.")
            condition6 = False

        conditions = (
            condition1,
            condition2,
            condition3,
            condition4,
            condition5,
            condition6,
        )

        for condition in conditions:
            if condition == False:
                stable = False

    if crystal_type == "monoclinic":
        to_print += ("Monoclinic crystal system.... \n")
        to_print += (
            "Born stability criteria for the stability of monoclinic systems are: Ref.[1,2]  \n"
        )
        to_print += (
            "(i) C11 > 0;  (ii)  C22 > 0; (iii)  C33 > 0; \n(iv)  C44 > 0;   (v)  C55 > 0 ;   (vi)  C66 > 0  "
        )
        to_print += (
            " (vii) [C11 + C22 + C33 + 2*(C12 + C13 + C23)] > 0;    (viii)  C33*C55 - C35^2 > 0; \n(ix)  C44*C66 - C46^2 > 0;   (x) C22 + C33 - 2*C23  > 0 "
        )
        to_print += (
            " (xi) C22*(C33*C55 - C35^2) + 2*C23*C25*C35 - (C23^2)*C55 - (C25^2)*C33   > 0  "
        )
        to_print += (
            " (xii)  2*[C15*C25*(C33*C12 - C13*C23) + C15*C35*(C22*C13 - C12*C23) + C25*C35*(C11*C23 - C12*C13)] - [C15*C15*(C22*C33 - C23^2) + C25*C25*(C11*C33 - C13^2) + C35*C35*(C11*C22 - C12^2)] + C55*g > 0  "
        )
        to_print += (
            "          where, g = [C11*C22*C33 - C11*C23*C23 - C22*C13*C13 - C33*C12*C12 + 2*C12*C13*C23 ] "
        )

        if c[0][0] > 0.0:
            to_print += ("Condition (i) is satisfied.")
            condition1 = True
        else:
            to_print += ("Condition (i) is NOT satisfied.")
            condition1 = False

        if c[1][1] > 0.0:
            to_print += ("Condition (ii) is satisfied.")
            condition2 = True
        else:
            to_print += ("Condition (ii) is NOT satisfied.")
            condition2 = False

        if c[2][2] > 0.0:
            to_print += ("Condition (iii) is satisfied.")
            condition3 = True
        else:
            to_print += ("Condition (iii) is NOT satisfied.")
            condition3 = False

        if c[3][3] > 0.0:
            to_print += ("Condition (iv) is satisfied.")
            condition4 = True
        else:
            to_print += ("Condition (iv) is NOT satisfied.")
            condition4 = False

        if c[4][4] > 0.0:
            to_print += ("Condition (v) is satisfied.")
            condition5 = True
        else:
            to_print += ("Condition (v) is NOT satisfied.")
            condition5 = False

        if c[5][5] > 0.0:
            to_print += ("Condition (vi) is satisfied.")
            condition6 = True
        else:
            to_print += ("Condition (vi) is NOT satisfied.")
            condition6 = False

        # for i in range(0, 6):
        #    if(c[i][i]  > 0.0):
        #        to_print += ("Condition (%2d) is satified." % (i+1))
        #    else:
        #        to_print += ("Condition (%2d) is NOT satified." % (i+1))

        if c[0][0] + c[1][1] + c[2][2] + 2 * (c[0][1] + c[0][2] + c[1][2]) > 0:
            to_print += ("Condition (vii) is satisfied.")
            condition7 = True
        else:
            to_print += ("Condition (vii) is NOT satisfied.")
            condition7 = False

        if c[2][2] * c[4][4] - c[2][4] * c[2][4] > 0:
            to_print += ("Condition (viii) is satisfied.")
            condition8 = True
        else:
            to_print += ("Condition (viii) is NOT satisfied.")
            condition8 = False

        if c[3][3] * c[5][5] - c[3][5] * c[3][5] > 0.0:
            to_print += ("Condition (ix) is satisfied.")
            condition9 = True
        else:
            to_print += ("Condition (ix) is NOT satisfied.")
            condition9 = False

        if c[1][1] + c[2][2] - 2 * c[1][2] > 0.0:
            to_print += ("Condition (x) is satisfied.")
            condition10 = True
        else:
            to_print += ("Condition (x) is NOT satisfied.")
            condition10 = False

        if (
            c[1][1] * (c[2][2] * c[4][4] - c[2][4] * c[2][4])
            + 2 * c[1][2] * c[1][4] * c[2][4]
            - c[1][2] * c[1][2] * c[4][4]
            - c[1][4] * c[1][4] * c[2][2]
            > 0.0
        ):
            to_print += ("Condition (xi) is satisfied.")
            condition11 = True
        else:
            to_print += ("Condition (xi) is NOT satisfied.")
            condition11 = False

        g = (
            (c[0][0] * c[1][1] * c[2][2])
            - (c[0][0] * c[1][2] * c[1][2])
            - (c[1][1] * c[0][2] * c[0][2])
            - (c[2][2] * c[0][1] * c[0][1])
            + 2.0 * (c[0][1] * c[0][2] * c[1][2])
        )

        h1 = 2 * (
            c[0][4] * c[1][4] * (c[2][2] * c[0][1] - c[0][2] * c[1][2])
            + c[0][4] * c[2][4] * (c[1][1] * c[0][2] - c[0][1] * c[1][2])
            + c[1][4] * c[2][4] * (c[0][0] * c[1][2] - c[0][1] * c[0][2])
        )

        h2 = (
            c[0][4] * c[0][4] * (c[1][1] * c[2][2] - c[1][2] * c[1][2])
            + c[1][4] * c[1][4] * (c[0][0] * c[2][2] - c[0][2] * c[0][2])
            + c[2][4] * c[2][4] * (c[0][0] * c[1][1] - c[0][1] * c[0][1])
        )

        x = h1 - h2 + c[4][4] * g

        #  to_print +=  'x = %10.4f' % x
        #  to_print +=  'g = %10.4f' % g

        if x > 0.0:
            to_print += ("Condition (xii) is satisfied.")
            condition12 = True
        else:
            to_print += ("Condition (xii) is NOT satisfied.")
            condition12 = False
        conditions = (
            condition1,
            condition2,
            condition3,
            condition4,
            condition5,
            condition6,
            condition7,
            condition8,
            condition9,
            condition10,
            condition11,
            condition12,
        )

        for condition in conditions:
            if condition == False:
                stable = False

    if crystal_type == "triclinic":
        to_print += ("Triclinic crystal system.... \n")
        to_print += (
            "Triclinic systems only criteria for stability are positive eigen values of the elastic tensor  \n"
        )
        evals = list(LA.eigvals(c))
        evals_print = list(np.around(np.array(evals), 3))
        to_print += ("%s" % evals_print)
        check = 0
        for i in range(len(evals)):
            if evals[i] > 0.0:
                stable = True
                pass
            else:
                check = 1
                stable = False

    if verbose:
        print(to_print)
    return stable


def stability_test_2d(matrix, lattice_type, verbose):
    c = np.zeros(shape=(3, 3))

    c[0, 0] = matrix[0, 0]
    c[0, 1] = matrix[0, 1]
    c[0, 2] = matrix[0, 5]

    c[1, 0] = matrix[1, 0]
    c[1, 1] = matrix[1, 1]
    c[1, 2] = matrix[1, 5]

    c[2, 0] = matrix[5, 0]
    c[2, 1] = matrix[5, 1]
    c[2, 2] = matrix[5, 5]

    stable = True
    to_print = ""
    if lattice_type == "hexagonal":
        to_print += ("Hexagonal lattice \n")
        to_print += ("Stability criteria for the stability of hexagonal system are: \n")
        to_print += ("(i) C11 + C12 > 0;    (ii) C11 - C12 > 0;  \n ")

        if c[0][0] + c[0][1] > 0.0:
            to_print += ("Condition (i) satisfied.")
            condition1 = True
        else:
            to_print += ("Condition (i) NOT satisfied.")
            condition1 = False

        if c[0][0] - c[0][1] > 0.0:
            to_print += ("Condition (ii) satisfied.")
            condition2 = True
        else:
            to_print += ("Condition (ii) NOT satisfied.")
            condition2 = False

        conditions = (condition1, condition2)

        for condition in conditions:
            if condition == False:
                stable = False

    elif lattice_type == "square":
        to_print += ("Square lattice \n")
        to_print += ("Stability criteria for the stability of square system are: \n")
        to_print += ("(i) C11 + C12 > 0;    (ii) C11 - C12 > 0;  (iii) C33 > 0     \n ")

        if c[0][0] + c[0][1] > 0.0:
            to_print += ("Condition (i) satisfied.")
            condition1 = True
        else:
            to_print += ("Condition (i) NOT satisfied.")
            condition1 = False

        if c[0][0] - c[0][1] > 0.0:
            to_print += ("Condition (ii) satisfied.")
            condition2 = True
        else:
            to_print += ("Condition (ii) NOT satisfied.")
            condition2 = False

        if c[2][2] > 0.0:
            to_print += ("Condition (iii) satisfied.")
            condition3 = True
        else:
            to_print += ("Condition (iii) NOT satisfied.")
            condition3 = False

        conditions = (condition1, condition2, condition3)

        for condition in conditions:
            if condition == False:
                stable = False

    elif lattice_type == "rectangular" or lattice_type == "rectangular-center":
        to_print += ("Rectangular lattice or Rectangular Center lattice \n")
        to_print += (
            "Stability criteria for the stability of rectangular lattice or rectangular Center lattic are: \n"
        )
        to_print += (
            "(i) 1/2*(C11 + C22 + (4(C12)**2 - (C11 -C22))**0.5)> 0;    (ii) 1/2*(C11 + C22 - (4(C12)**2 - (C11 -C22))**0.5)> 0 ;  (iii) C33 > 0     \n "
        )

        if (
            0.5
            * (
                c[0][0]
                + c[1][1]
                + (4 * (c[0][1]) ** 2 + (c[0][0] - c[1][1]) ** 2) ** 0.5
            )
            > 0.0
        ):
            to_print += ("Condition (i) satisfied.")
            condition1 = True
        else:
            to_print += ("Condition (i) NOT satisfied.")
            condition1 = False

        if (
            0.5
            * (
                c[0][0]
                + c[1][1]
                - (4 * (c[0][1]) ** 2 + (c[0][0] - c[1][1]) ** 2) ** 0.5
            )
            > 0.0
        ):
            to_print += ("Condition (ii) satisfied.")
            condition2 = True
        else:
            to_print += ("Condition (ii) NOT satisfied.")
            condition2 = False

        if c[2][5] > 0.0:
            to_print += ("Condition (iii) satisfied.")
            condition3 = True
        else:
            to_print += ("Condition (iii) NOT satisfied.")
            condition3 = False

        conditions = (condition1, condition2, condition3)

        for condition in conditions:
            if condition == False:
                stable = False

    elif lattice_type == "oblique":
        to_print += ("Oblique lattice \n")
        to_print += ("Stability criteria for the stability of oblique lattice are: \n")
        to_print += ("(i) C11 > 0 ;    (ii) C11 *C22 > (C12)*2 ;  (iii) det(Cij) > 0     \n ")

        if c[0][0] > 0.0:
            to_print += ("Condition (i) satisfied.")
            condition1 = True
        else:
            to_print += ("Condition (i) NOT satisfied.")
            condition1 = False

        if c[0][0] * c[1][1] > (c[0][1]) ** 2:
            to_print += ("Condition (ii) satisfied.")
            condition2 = True
        else:
            to_print += ("Condition (ii) NOT satisfied.")
            condition2 = False

        if np.linalg.det(c) > 0.0:
            to_print += ("Condition (iii) satisfied.")
            condition3 = True
        else:
            to_print += ("Condition (iii) NOT satisfied.")
            condition3 = False

        conditions = (condition1, condition2, condition3)

        for condition in conditions:
            if condition == False:
                stable = False
    else:
        raise NameError("lattice type %s not in recognized" % lattice_type)

    if verbose:
        print(to_print)
    return stable
