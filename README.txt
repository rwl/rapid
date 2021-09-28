Example for sequential simulation using Fine integration
Changed files are

1) pdefaults.py
added option --nsave, to define how many increments between tstart and tend to save in the output, default=10
- added option --fault to initiate faults, e.g. --fault fault.json

2) fault.json
Reading the fault.json file creates dictionary that can be used in the program. You can add more levels and data as needed.
The fault data is stored md.fault
See access to data in para_real_proto.py and pdefault.py

3) Fault.py
Creates object for fault and assigns specific functions from ffault.py to the object based on the fault type
Once created, the object can be later invoked with .run and will be executed based on the logic of run function on Fault.py

4) ffault.py
Fault functions that operate on main program data, e.g. Y bus matrix
Similar purpose as fn.py
The functions get assigned to Fault objects so that when Fault.run is invoked, the proper function is executed.

5) PowerModel.py
Read the original data and solve the power flow problem to calculate intial values for all other variables. 
Save the Ybus matrix and pre-factorized LU object to be used in "fn.py". 

Example
mpiexec -n 19 python para_real.py 0.0 10. --nCoarse 25 --nFine 250 --tol 0.1 --tolcheck maxabs --debug 1 --fault fault.json -o ofile.csv

Saves 5 increments of solution for times
0.0000e+00
4.0000e-02
8.0000e-02
1.2000e-01
1.6000e-01
2.0000e-01

ofile.csv has the same format as before.
In this example ofile.txt has 6 lines, where first entry in each line is time stamp, and the rest of the line is solution vector
