This example shows the user how to use mechelastic with qe_thermo_pw code and implement it in python.

Input units must be in atomic units, mechelastic converts everything to Angstrom at the end

Do not use the quantum espresso inputs A , B , C, COSAB, COSBC, or COSAC, use [celldm(1)-celldm(2)] as defined in the input descriptions of quantum espresso.

Description of the Input files:
outfile = 'si.elastic.out'    : This is the outputput you feed into thermo_pw. 
infile  = 'si.elastic.in'     : This is the input you feed into thermo_pw. 

Make sure you run the python script in the same directory as these files