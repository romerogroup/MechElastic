#!/usr/bin/env python

import numpy as np


def stability_test(matrix, crystal_type):

    """This methods tests for the stability of a structure."""

    c = np.copy(matrix)

    if crystal_type == "cubic":
        print("Cubic crystal system \n")
        print(
            "Born stability criteria for the stability of cubic system are: [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n"
        )
        print("(i) C11 - C12 > 0;    (ii) C11 + 2C12 > 0;   (iii) C44 > 0 \n ")

        # check (i)   keep in mind list starts with 0, so c11 is stored as c00
        if c[0][0] - c[0][1] > 0.0:
            print("Condition (i) satified.")
        else:
            print("Condition (i) NOT satisfied.")

        if c[0][0] + 2 * c[0][1] > 0.0:
            print("Condition (ii) satified.")
        else:
            print("Condition (ii) NOT satisfied.")

        if c[3][3] > 0.0:
            print("Condition (iii) satified.")
        else:
            print("Condition (iii) NOT satisfied.")

    if crystal_type == "hexagonal":
        print("Hexagonal crystal system \n")
        print(
            "Born stability criteria for the stability of hexagonal system are: [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n"
        )
        print(
            "(i) C11 - C12 > 0;    (ii) 2*C13^2 < C33(C11 + C12);   (iii) C44 > 0 \n "
        )

        # check (i)   keep in mind list starts with 0, so c11 is stored as c00
        if c[0][0] - c[0][1] > 0.0:
            print("Condition (i) satified.")
        else:
            print("Condition (i) NOT satisfied.")

        if 2 * (c[0][2] * c[0][2]) < c[2][2] * (c[0][0] + c[0][1]):
            print("Condition (ii) satified.")
        else:
            print("Condition (ii) NOT satisfied.")

        if c[3][3] > 0.0:
            print("Condition (iii) satified.")
        else:
            print("Condition (iii) NOT satisfied.")

    if crystal_type == "tetragonal":
        print("Tetragonal crystal system \n")
        print(
            "Born stability criteria for the stability of Tetragonal system are: [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n"
        )
        print(
            "(i) C11 - C12 > 0;    (ii) 2*C13^2 < C33(C11 + C12);   (iii) C44 > 0;   (iv) C66 > 0;    (v) 2C16^2 < C66*(C11-C12) \n "
        )

        # check (i)   keep in mind list starts with 0, so c11 is stored as c00
        if c[0][0] - c[0][1] > 0.0:
            print("Condition (i) is satified.")
        else:
            print("Condition (i) is NOT satisfied.")

        if 2 * (c[0][2] * c[0][2]) < c[2][2] * (c[0][0] + c[0][1]):
            print("Condition (ii) is satified.")
        else:
            print("Condition (ii) is NOT satisfied.")

        if c[3][3] > 0.0:
            print("Condition (iii) is satified.")
        else:
            print("Condition (iii) is NOT satisfied.")

        if c[5][5] > 0.0:
            print("Condition (iv) is satified.")
        else:
            print("Condition (iv) is NOT satisfied.")

        if 2 * c[0][5] * c[0][5] < c[5][5] * (c[0][0] - c[0][1]):
            print("Condition (v) is satified.")
        else:
            print("Condition (v) is NOT satisfied.")

    if crystal_type == "rhombohedral-1":
        print(
            "Rhombohedral (class-1): i.e. structures with point group: 3m, -3m and 32 \n"
        )
        print(
            "Born stability criteria for the stability of Rhombohedral-1 class system are: [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n"
        )
        print(
            "(i) C11 - C12 > 0;    (ii) C13^2 < (1/2)*C33(C11 + C12);   (iii) C14^2 < (1/2)*C44*(C11-C12) = C44*C66;   (iv)  C44 > 0; \n "
        )

        if c[0][0] - c[0][1] > 0.0:
            print("Condition (i) is satified.")
        else:
            print("Condition (i) is NOT satisfied.")

        if (c[0][2] * c[0][2]) < (0.5) * c[2][2] * (c[0][0] + c[0][1]):
            print("Condition (ii) is satified.")
        else:
            print("Condition (ii) is NOT satisfied.")

        if c[0][3] * c[0][3] < 0.5 * c[3][3] * (c[0][0] - c[0][1]):
            print("Condition (iii) is satified.")
        else:
            print("Condition (iii) is NOT satisfied.")

        if c[3][3] > 0.0:
            print("Condition (iv) is satified.")
        else:
            print("Condition (iv) is NOT satisfied.")

    if crystal_type == "rhombohedral-2":
        print("Rhombohedral (class-2): i.e structures with point group: 3, -3 \n")
        print(
            "Born stability criteria for the stability of Rhombohedral-1 class system are: [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n"
        )
        print(
            "(i) C11 - C12 > 0;    (ii) C13^2 < (1/2)*C33(C11 + C12);   (iii) C14^2 + C15^2 < (1/2)*C44*(C11-C12) = C44*C66;   (iv)  C44 > 0;  Note: C15 is added.. \n "
        )

        if c[0][0] - c[0][1] > 0.0:
            print("Condition (i) is satified.")
        else:
            print("Condition (i) is NOT satisfied.")

        if c[0][2] * c[0][2] < (0.5) * c[2][2] * (c[0][0] + c[0][1]):
            print("Condition (ii) is satified.")
        else:
            print("Condition (ii) is NOT satisfied.")

        if c[0][3] * c[0][3] + c[0][4] * c[0][4] < 0.5 * c[3][3] * (c[0][0] - c[0][1]):
            print("Condition (iii) is satified.")
        else:
            print("Condition (iii) is NOT satisfied.")

        if c[3][3] > 0.0:
            print("Condition (iv) is satified.")
        else:
            print("Condition (iv) is NOT satisfied.")

    if crystal_type == "orthorhombic":
        print("Orthorhombic crystal system.... \n")
        print(
            "Born stability criteria for the stability of Orthorhombic systems are: [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n"
        )
        print(
            "(i) C11 > 0;   (ii) C11*C22 > C12^2;   (iii) C11*C22*C33 + 2C12*C13*C23 - C11*C23^2 - C22*C13^2 - C33*C12^2 > 0;   (iv)  C44 > 0;   (v)  C55 > 0 ;   (vi)  C66 > 0 \n"
        )

        # check (i)   keep in mind list starts with 0, so c11 is stored as c00
        if c[0][0] > 0.0:
            print("Condition (i) is satified.")
        else:
            print("Condition (i) is NOT satisfied.")

        if c[0][0] * c[1][1] > c[0][1] * (c[0][1]):
            print("Condition (ii) is satified.")
        else:
            print("Condition (ii) is NOT satisfied.")

        if (
            c[0][0] * c[1][1] * c[2][2]
            + 2 * c[0][1] * c[0][2] * c[1][2]
            - c[0][0] * c[1][2] * c[1][2]
            - c[1][1] * c[0][2] * c[0][2]
            - c[2][2] * c[0][1] * c[0][1]
            > 0
        ):
            print("Condition (iii) is satified.")
        else:
            print("Condition (iii) is NOT satisfied.")

        if c[3][3] > 0.0:
            print("Condition (iv) is satified.")
        else:
            print("Condition (iv) is NOT satisfied.")

        if c[4][4] > 0.0:
            print("Condition (iv) is satified.")
        else:
            print("Condition (iv) is NOT satisfied.")

        if c[5][5] > 0.0:
            print("Condition (iv) is satified.")
        else:
            print("Condition (iv) is NOT satisfied.")

    if crystal_type == "monoclinic":
        print("Monoclinic crystal system.... \n")
        print(
            "Born stability criteria for the stability of monoclinic systems are: [Ref- Mouhat and Coudert, PRB 90, 224104 (2014), and Wu et al. PRB 76, 054115 (2007)]  \n"
        )
        print(
            "(i) C11 > 0;  (ii)  C22 > 0; (iii)  C33 > 0; (iv)  C44 > 0;   (v)  C55 > 0 ;   (vi)  C66 > 0  "
        )
        print(
            " (vii) [C11 + C22 + C33 + 2*(C12 + C13 + C23)] > 0;    (viii)  C33*C55 - C35^2 > 0;   (ix)  C44*C66 - C46^2 > 0;   (x) C22 + C33 - 2*C23  > 0 "
        )
        print(
            " (xi) C22*(C33*C55 - C35^2) + 2*C23*C25*C35 - (C23^2)*C55 - (C25^2)*C33   > 0  "
        )
        print(
            " (xii)  2*[C15*C25*(C33*C12 - C13*C23) + C15*C35*(C22*C13 - C12*C23) + C25*C35*(C11*C23 - C12*C13)] - [C15*C15*(C22*C33 - C23^2) + C25*C25*(C11*C33 - C13^2) + C35*C35*(C11*C22 - C12^2)] + C55*g > 0  "
        )
        print(
            "          where, g = [C11*C22*C33 - C11*C23*C23 - C22*C13*C13 - C33*C12*C12 + 2*C12*C13*C23 ] "
        )

        if c[0][0] > 0.0:
            print("Condition (i) is satified.")
        else:
            print("Condition (i) is NOT satified.")

        if c[1][1] > 0.0:
            print("Condition (ii) is satified.")
        else:
            print("Condition (ii) is NOT satified.")

        if c[2][2] > 0.0:
            print("Condition (iii) is satified.")
        else:
            print("Condition (iii) is NOT satified.")

        if c[3][3] > 0.0:
            print("Condition (iv) is satified.")
        else:
            print("Condition (iv) is NOT satified.")

        if c[4][4] > 0.0:
            print("Condition (v) is satified.")
        else:
            print("Condition (v) is NOT satified.")

        if c[5][5] > 0.0:
            print("Condition (vi) is satified.")
        else:
            print("Condition (vi) is NOT satified.")

        # for i in range(0, 6):
        #    if(c[i][i]  > 0.0):
        #        print("Condition (%2d) is satified." % (i+1))
        #    else:
        #        print("Condition (%2d) is NOT satified." % (i+1))

        if c[0][0] + c[1][1] + c[2][2] + 2 * (c[0][1] + c[0][2] + c[1][2]) > 0:
            print("Condition (vii) is satified.")
        else:
            print("Condition (vii) is NOT satisfied.")

        if c[2][2] * c[4][4] - c[2][4] * c[2][4] > 0:
            print("Condition (viii) is satified.")
        else:
            print("Condition (viii) is NOT satisfied.")

        if c[3][3] * c[5][5] - c[3][5] * c[3][5] > 0.0:
            print("Condition (ix) is satified.")
        else:
            print("Condition (ix) is NOT satisfied.")

        if c[1][1] + c[2][2] - 2 * c[1][2] > 0.0:
            print("Condition (x) is satified.")
        else:
            print("Condition (x) is NOT satisfied.")

        if (
            c[1][1] * (c[2][2] * c[4][4] - c[2][4] * c[2][4])
            + 2 * c[1][2] * c[1][4] * c[2][4]
            - c[1][2] * c[1][2] * c[4][4]
            - c[1][4] * c[1][4] * c[2][2]
            > 0.0
        ):
            print("Condition (xi) is satified.")
        else:
            print("Condition (xi) is NOT satisfied.")

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

        #  print 'x = %10.4f' % x
        #  print 'g = %10.4f' % g

        if x > 0.0:
            print("Condition (xii) is satified.")
        else:
            print("Condition (xii) is NOT satisfied.")

    return
