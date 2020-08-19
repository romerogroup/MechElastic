This example shows the user how to use mechelastic and implement it in python.

Input units must be in atomic units, mechelastic converts everything to Angstrom at the end

Do not use the quantum espresso inputs A , B , C, COSAB, COSBC, or COSAC, use [celldm(1)-celldm(2)] as defined in the input descriptions of quantum espresso.

Description of the Input files:
outfile = 'ElaStic_2nd.out' or 'ElaStic_3rd.out'  : This is a generated output of ElaStic for either 2nd order or 3rd order
infile  = 'scf.in'     : This is the starting relaxed structure input you feed into Elastic. The name doesn't have to be 'scf.in'

Make sure you run the python script in the same directory as these files