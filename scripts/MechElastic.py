#!/usr/bin/env python

"""This script reads elastic constants from OUTCAR, calculates elastic moduli, and performs mechanical stability test:
   Authors: Sobhit Singh (1,2) and Aldo Romero (2)
   Email: smsingh@mix.wvu.edu, alromero@mail.wvu.edu
   (1) Rutgers University, Piscataway, NJ, USA
   (2) West Virginia University, Morgantown, WV, USA
Version: 12.2019   #Dec. 2019

Please cite the below paper if you use this code for your research:
     Sobhit Singh, Irais Valencia-Jaime, Olivia Pavlic, and Aldo H. Romero; Phys. Rev. B 97, 054108 (2018).


To RUN this code for a 3D system assuming you know the crystal type before hand:
 >> python MechElastic.py -i OUTCAR-file -c hexagonal --dim 3D

Note: The crystal type is needed only to perform the elastic stabilty test.
      This test is not yet implemented for 2D systems, so you can ignore the crystal type for 2D systems.

If the crystal symmetry is not provided by user, then the code will rely on spglib to find it
 >> python MechElastic.py -i OUTCAR-file --dim 3D

To RUN this code for a 2D monolayer system
 >> python MechElastic.py -i OUTCAR-file --dim 2D
      (please pay attention to the warning for 2D systems)

OUTCAR file, if present in the current directory, is read as default unless a filename is specified by the user.


Disclaimer: Please check the authenticity of your results before publishing.
            AUTHORS of this script do not guarantee the quality and/or accuracy
            of results generated using this script.

"""

import argparse
import sys
import mechelastic

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i", "--input", type=str, help="input the OUTCAR", default="OUTCAR"
)
parser.add_argument(
    "-c",
    "--crystal",
    type=str,
    default=None,
    help="Provide the crystal type. Otherwise it would be determined from OUTCAR",
)
parser.add_argument(
    "-d",
    "--dim",
    type=str,
    help="Enter the dimension, 2D or 3D: For example: '-d 2D' or '--dim 2D'. Default is '3D' ",
    default="3D",
)
parser.add_argument(
    "-co", "--code", type=str, help="DFT code", default="vasp", choices=["vasp"]
)
args = parser.parse_args()

print("-----------------------------")
print("List of arguments entered:")
print("Input file name:", args.input)
print("DFT code:", args.code)
print("Crystal type: ", args.crystal)
print("Dimensions:", args.dim)
print("-----------------------------")


# calculate elastic properties
def main():
    mechelastic.calculate_elastic(
        code=args.code, dim=args.dim, infile=args.input, crystal=args.crystal
    )


if __name__ == "__main__":
    main()
