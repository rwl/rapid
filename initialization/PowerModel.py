import pandapower as pp
import pandas as pd
import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix, csc_matrix, lil_matrix, identity
from scipy.sparse.linalg import splu
from scipy.io import loadmat
from initialization.initial import initialization
from initialization.linearization import linearization_all
from initialization.bustype import bustype
from initialization.makeYbus import makeYbus
from initialization.makeSbus import makeSbus
from initialization.newtonpf import newtonpf
from initialization.pfsoln import pfsoln
import pickle
import csv
########################################################################################################################
import opendssdirect as dss
import sys
import os
import json
from timeit import default_timer as timer

########################################################################################################################
solveYbus = None
Lower = None
Upper = None

def updateYbus1(temp):
    global Ybus0
    global Ybusbf
    Ybus0 = Ybusbf + temp

def updateYbus2():
    global Ybus0
    global Ybusbf
    Ybus0 = Ybusbf

def solve_Ybus(stage):
    global solveYbus
    global Ybus0
    global Lower
    global Upper

    solveYbus = None
    solveYbus = splu(Ybus0.tocsc())  # factorize the Ybus matrix
    #print("solve stage", stage)

def info_Ybus(stage):
    global solveYbus
    global Ybus0

def get_Ybus():
    global solveYbus
    return solveYbus

# read the original data
#ReadFile = 'C:/.../v6_sequential_dummy_complete/Data_EI.mat'
temp = loadmat('./data/Data_NE.mat', struct_as_record=True)
##
ExcData = temp['ExcData']      # dynamic data
GenData = temp['GenData']
satData = temp['satdata']
TurbData = temp['TurbData']
mpc = temp['mpc']              # network data

baseMVA = 1.0*mpc['baseMVA'][0,0][0,0]     # 1.0 is added to make this number float64
SYSMVA = baseMVA
bus = csc_matrix(mpc['bus'][0,0])          # bus data of network data
gen = csc_matrix(mpc['gen'][0,0])          # gen data of network data
branch = csc_matrix(mpc['branch'][0,0])    # branch data of network data
gencost = csc_matrix(mpc['gencost'][0,0])  # gencost data of network data (may not need for the dynamic simulation)

##
nb = len(bus.toarray())
ng = len(gen.toarray())

## get bus index lists of each type of bus
ref, pv, pq = bustype(bus, gen)
ref_pv = list(set(ref) | set(pv))
ref_pv = np.sort(ref_pv)

## construct the Ybus matrix
Ybus, Yf, Yt = makeYbus(baseMVA, bus, branch)

## construct the complex bus power (Sbus) injecton vector
Sbus = makeSbus(baseMVA, bus, gen)

## initialize voltages for the power flow problem
V0 = bus[:, 7].toarray()*np.exp(1j*bus[:, 8].toarray()*(np.pi/180))
vcb = np.ones(len(V0))
on = sp.find(gen[:, 7] > 0)[0]
gbus = gen[on, 0].toarray().astype('int64')-1

row = gen[:, 0].toarray().reshape(-1)-1
col = np.arange(ng)
data = (gen[:, 7] > 0).toarray().reshape(-1)
Cg = csc_matrix((data, (row, col)), shape=(nb, ng))
bus_gen_status = Cg*np.ones(ng)
vcb[pq] = 0
k = sp.find(vcb[gbus] > 0 )[0]
V0[gbus[k].reshape(-1)] = gen[on[k], 5].toarray() / abs(V0[gbus[k].reshape(-1)])*V0[gbus[k].reshape(-1)]
V0 = V0.reshape(-1)

# (Co-simulation option started)
########################################################################################################################
Cosimulation = '0'

if Cosimulation == '1':
    idxDistBus = [3]
else:
    idxDistBus = []

First_Iteration = True
I_Dist_Dynamic = np.zeros((len(idxDistBus), 1), dtype=np.complex)  # Distribution dynamic current initialization
num_distribution = len(idxDistBus)
Loads_P_Dist = np.zeros(num_distribution)  # initialize OpenDSS circuits' Loads_P_Dist (real power)
Loads_Q_Dist = np.zeros(num_distribution)  # initialize OpenDSS circuits' Loads_Q_Dist (reactive power)
owd = os.getcwd()
Voltages_from_transmission = np.ones(len(idxDistBus))  # Transmission voltage initialization

if len(idxDistBus) > 0:

    start = timer()  # for execution time calculation
    Dist_fileName = './ckt24/master_ckt24.dss'  # Distribution System File directory
    I_Dist_Static = np.zeros((len(idxDistBus), 1), dtype=np.complex)  # Distribution steady state current initialization

    OpenDSSfileName = Dist_fileName
    num_distribution = len(idxDistBus)         # number of distribution system solved in parallel
    Loads_P_Dist = np.zeros(num_distribution)  # initialize OpenDSS circuits' Loads_P_Dist (real power)
    Loads_Q_Dist = np.zeros(num_distribution)  # initialize OpenDSS circuits' Loads_Q_Dist (reactive power)
    S = np.zeros(num_distribution)             # initialize OpenDSS circuits' S (apparent) power
    S0 = np.zeros(num_distribution)            # initialize OpenDSS circuits' S0 (apparent) power
    deltaS = np.zeros(num_distribution)        # initialize OpenDSS circuits' deltaS (difference in apparent power)
    DistConv = np.zeros(num_distribution)
    tol_d = 1e-8  # set tolerance for convergence if convergence checked outside of OpenDSS

    V_Mags_trans = np.ones(num_distribution)  # initialize Voltages of transmission system

    dss.run_command('Compile (' + OpenDSSfileName + ')')
    dss.run_command('VSource.source.pu=' + str(V_Mags_trans))  # set the slack voltage from transmission bus voltage
    dss.run_command('set mode = snap')  # Set mode at the start of the iteration
    dss.run_command('Solve')

    Load_Dist = dss.Circuit.TotalPower()    # Total Power of the system in KW and KVar (-ve sign for absorbed power)
    Loads_P_Dist[0] = -Load_Dist[0] / 1000  # converting to MW and +ve
    Loads_Q_Dist[0] = -Load_Dist[1] / 1000  # converting to MVar and +Ve
    print('iterations : ', dss.Solution.Iterations())

    #
    if dss.Solution.Converged() == True:
        DistConv[0] = 1
    else:
        DistConv[0] = 0

    Pload_from_Tran = np.zeros(len(idxDistBus))  # initialize transmission system real load at the buses where distribution system are connected
    QLoad_from_Tran = np.zeros(len(idxDistBus))  # initialize transmission system reactive load at the buses where distribution system are connected

    for i in range(1, len(idxDistBus) + 1):
        Pload_from_Tran[i-1] = bus[idxDistBus[i-1],2] - Loads_P_Dist[i-1]  # transmission system real load at the buses where distribution system are connected
        QLoad_from_Tran[i-1] = bus[idxDistBus[i-1],3] - Loads_Q_Dist[i-1]  # transmission system reactive load at the buses where distribution system are connected

    DistConv = np.zeros(len(idxDistBus))
    TransConv = 0
    iteration = np.zeros(len(idxDistBus))

    #### start iterative transmission - distribution system steady state power flow solution
    itera = 0
    while sum(DistConv) < len(idxDistBus) or TransConv < 1:

        itera = itera + 1
        # Transmission system power flow solution
        for i in range(1, len(idxDistBus) + 1):
            bus[idxDistBus[i-1],2] = Pload_from_Tran[i-1] + Loads_P_Dist[i-1]  # add distibution system load to transmission system
            bus[idxDistBus[i-1],3] = QLoad_from_Tran[i-1] + Loads_Q_Dist[i-1]  # add distibution system load to transmission system

        V, success, iterations = newtonpf(baseMVA, Ybus, Sbus, V0, ref, pv, pq, 1)
        V0 = V

        if success:
            TransConv = 1
        else:
            TransConv = 0

        # Distribution system power flow solution
        Voltages_from_transmission = np.array(V[idxDistBus])
        V_Mags_trans = np.abs(Voltages_from_transmission)
        if itera == 1:
            dss.run_command('set Maxiterations = 2')

        dss.run_command('VSource.source.pu=' + str(V_Mags_trans))  # set the slack voltage from transmission bus voltage
        dss.run_command('Solve')
        Load_Dist = dss.Circuit.TotalPower()    # Total Power of the sytem in KW and KVar (-ve sign for absorbed power)
        Loads_P_Dist[0] = -Load_Dist[0]/1000  # converting to MW and +ve
        Loads_Q_Dist[0] = -Load_Dist[1]/1000  # converting to MVar and +Ve
        #I_Dist_Static[0] = -(Loads_P_Dist[0] - 1j*Loads_Q_Dist[0]) / (SYSMVA*V_Mags_trans[0])  # complex number? S=VI*
        I_Dist_Static[0] = -(Loads_P_Dist[0] - 1j*Loads_Q_Dist[0]) / (SYSMVA*np.conjugate(Voltages_from_transmission[0]))  # complex number? S=VI*
        S[0] = np.abs(Loads_P_Dist[0] + 1j*Loads_Q_Dist[0])
        deltaS[0] = S[0] - S0[0]
        S0[0] = S[0]

        if dss.Solution.Converged() == True:
            DistConv[0] = 1
        else:
            DistConv[0] = 0

    # Distribution system power flow end
    os.chdir(owd)
    end = timer()
    # print('time: ', (end - start))
    #### end iterative transmission- distribution system power flow solution
    print('Co-simulation Converged')
    print('S Steady State:', Loads_P_Dist + 1j*Loads_Q_Dist)
    I_Dist_Dynamic = I_Dist_Static[0]
    print(I_Dist_Static)

# if there is no distribution system only transmission system power flow solution
########################################################################################################################
else:
    V, success, iterations = newtonpf(baseMVA, Ybus, Sbus, V0, ref, pv, pq, 20)

Sbase =  baseMVA # 100 #
## organize the solution
pg_sol, qg_sol = pfsoln(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V, ref, pv, pq)
v_sol = V

#########################################################################################################################
# function for Dynamic simulation of Distribution System with the voltage defined by Transmission System
def DistLoad_Dynamic(idxDistBus, Voltages_from_transmission, Sbase, Step):

    V_trans = Voltages_from_transmission[idxDistBus]
    V_Mags_trans = np.abs(V_trans)

    if V_Mags_trans < 0.001: # To avoid  problem in case of Fault (V is near 0)
        V_trans.real = 0.001*np.sign(V_trans.real)  # Byungkwon: set the real value of the complex voltage to 0.001 to avoid the numerical instability
        V_trans.imag = 0.001*np.sign(V_trans.imag)  # Byungkwon: set the imaginary value of the complex voltage to 0.001 to avoid the numerical instability

    dss.run_command('VSource.source.pu=' + str(V_Mags_trans))  # set the slack voltage from transmission bus voltage
    dss.run_command('set mode = dynamic')
    dss.run_command('set control = time')
    dss.run_command('Set stepsize =' + str(Step))
    dss.run_command('Set Number = 1')
    dss.run_command('Solve')
    Load_Dist_Dynamic = dss.Circuit.TotalPower()  # Total Power of the sytem in KW and KVar (-ve sign for absorbed power)
    Loads_P_Dist[0] = -Load_Dist_Dynamic[0] / 1000  # converting to MW and +ve
    Loads_Q_Dist[0] = -Load_Dist_Dynamic[1] / 1000  # converting to MVar and +Ve
    I_Dist_Dynamic = -(Loads_P_Dist[0] - 1j*Loads_Q_Dist[0])/(SYSMVA*np.conjugate(V_trans))  # Byungkwon: fixed the issue so that it now takes the complex voltage not just voltage magnitude to calculate the current injection

    return I_Dist_Dynamic

########################################################################################################################
Ybus = lil_matrix(Ybus)

## Reorganize bus information
if len(idxDistBus) > 0:
    for i in range(1, len(idxDistBus) + 1):
        bus[idxDistBus[i-1], 2] = Pload_from_Tran[i-1] # bus[idxDistBus[i - 1] - 1, 2] - Loads_P_Dist[i - 1] /100 # add distibution system load to transmission system
        bus[idxDistBus[i-1], 3] = QLoad_from_Tran[i-1] # bus[idxDistBus[i - 1] - 1, 3] - Loads_Q_Dist[i - 1] /100

########################################################################################################################
## Define indices for variables
# Number of Generator States For IEEE Model 2.2 = 9, IEEE Model 1.1 = 5 etc.
nGenST  = 9
nExcST = 4   # Number of Excitation System States
nGovST = 1   # Number of Governor States
nTurbST = 1  # Number of Turbine States
nLoadST = 2  # Number of Load States
nBusIn = 2   # Number of Inputs

nogen = len(GenData)
nobus = len(v_sol)

LoadBuses = sp.find(bus.tocsc()[:, 2] > 0)[0]
noload = len(LoadBuses)

nSTGEP = (nGenST + nExcST + nTurbST + nGovST)  # Total States for each Generator
nTotDE = nSTGEP*nogen + nLoadST*noload    # Total No of States in the System
nTotIn = nBusIn*nogen   # Total No of Inputs
fB  =  GenData[:,16]    # BaseFrequency

# indices for each device
idxGov = np.arange(nogen*nGovST, dtype=np.int64)   #### Governor states indices
idxTurb = np.arange(idxGov[-1]+1, nogen*nTurbST + idxGov[-1]+1, 1, dtype=np.int64)  #### Turbine States
idxExc = np.arange(idxTurb[-1]+1, nogen*nExcST + idxTurb[-1]+1, 1, dtype=np.int64)  #### Excitation System States
idxGen = np.arange(idxExc[-1]+1, nogen*nGenST + idxExc[-1]+1, 1, dtype=np.int64)    #### Generator states indices
idxLoad = np.arange(idxGen[-1]+1, noload*2 + idxGen[-1]+1, 1, dtype=np.int64)       #### Load dummy states indices
idxV = np.arange(idxLoad[-1]+1, nobus + idxLoad[-1]+1, 1, dtype=np.int64)       #### terminal voltage
idxI = np.arange(idxV[-1]+1, nobus + idxV[-1]+1, 1, dtype=np.int64)       #### Currents
idxTe = np.arange(idxI[-1]+1, nogen + idxI[-1]+1, 1, dtype=np.int64)      #### Torque
idxs = np.concatenate((idxGov, idxTurb, idxExc, idxGen, idxLoad, idxV, idxI, idxTe))

# indices for each variable
idxEfd = (idxTurb[-1] + 1) + np.arange(nogen, dtype=np.int64)
idxV2 = (idxEfd[-1] + 1) + np.arange(nogen, dtype=np.int64)
idxV1 = (idxV2[-1] + 1) + np.arange(nogen, dtype=np.int64)
idxVR = (idxV1[-1] + 1) + np.arange(nogen, dtype=np.int64)
idxdelta = (idxVR[-1] + 1) + np.arange(nogen, dtype=np.int64)
idxsm = (idxdelta[-1] + 1) + np.arange(nogen, dtype=np.int64)
idxsif = (idxsm[-1] + 1) + np.arange(nogen, dtype=np.int64)
idxsih = (idxsif[-1] + 1) + np.arange(nogen, dtype=np.int64)
idxsig = (idxsih[-1] + 1) + np.arange(nogen, dtype=np.int64)
idxsik = (idxsig[-1] + 1) + np.arange(nogen, dtype=np.int64)
idxEdc = (idxsik[-1] + 1) + np.arange(nogen, dtype=np.int64)
idxXadspp = (idxEdc[-1] + 1) + np.arange(nogen, dtype=np.int64)
idxXaqspp = (idxXadspp[-1] + 1) + np.arange(nogen, dtype=np.int64)
idxILr = (idxXaqspp[-1] + 1) + np.arange(noload, dtype=np.int64)
idxILi = (idxILr[-1] + 1) + np.arange(noload, dtype=np.int64)

# indices for the network
idxGenBus = (GenData[:, 0] - 1).astype('int64')
idxLoadBus = LoadBuses
idxGenBus_new = np.unique(idxGenBus)

########################################################################################################################
## Obtain initial values for all other variables using the power flow solution
print("initializing during PowerModel ")
X0, Vbus0, Ibus0, Ybus0, GenData, Pc, Vref = initialization(GenData, ExcData, satData, TurbData, bus, v_sol, pg_sol, qg_sol, Ybus=Ybus)
Ibus = Ibus0
if len(idxDistBus) > 0: 
    for i in range(1, len(idxDistBus) + 1): 
        Ibus[idxDistBus[i-1]] = I_Dist_Dynamic[i-1] 

Ybusbf = Ybus0.copy()     # copy the modified Ybus0 as before fault (Ybusbf), the lil_matrix

solve_Ybus("init")
info_Ybus("init")

## Linearization for the reduced model
d_ref = np.array([1]).reshape(-1)[0]-1
A, B, states, inputs, d_delta = linearization_all(GenData, ExcData, satData, TurbData, bus, X0, Vbus0, Ibus0, Pc, idxV, d_ref)

########################################################################################################################
# Parameters
Vref = (Vref.real).reshape(-1,1)
Pc = (Pc.real).reshape(-1,1)
GenStatus = np.ones([nogen, 1], dtype=np.int64)
GENMVA = GenData[:, 17].reshape(-1,1)
Xfl = GenData[:,26].reshape(-1,1)
Xhl = GenData[:,27].reshape(-1,1)
Xkl = GenData[:,28].reshape(-1,1)
Xgl = GenData[:,29].reshape(-1,1)
Xl = GenData[:,14].reshape(-1,1)
Ra = GenData[:,13].reshape(-1,1)
Asd = GenData[:,18].reshape(-1,1)
Bsd = GenData[:,19].reshape(-1,1)
siTd = GenData[:,20].reshape(-1,1)
Asq = GenData[:,21].reshape(-1,1)
Bsq = GenData[:,22].reshape(-1,1)
siTq = GenData[:,23].reshape(-1,1)
Xadu = GenData[:,24].reshape(-1,1)
Xaqu = GenData[:,25].reshape(-1,1)
Rf = GenData[:,30].reshape(-1,1)
Rh = GenData[:,31].reshape(-1,1)
Rk = GenData[:,32].reshape(-1,1)
Rg = GenData[:,33].reshape(-1,1)
H = GenData[:,11].reshape(-1,1)
D = GenData[:,12].reshape(-1,1)
Tc = GenData[:,15].reshape(-1,1)
wB = 2*np.pi*GenData[:,16].reshape(-1,1)
Txd = np.array([0.01]).reshape(-1,1)
Txq = np.array([0.01]).reshape(-1,1)

RD = TurbData[:,1].reshape(-1,1)
TSV = TurbData[:,2].reshape(-1,1)
Psvmax = TurbData[:,3].reshape(-1,1)
Psvmin = TurbData[:,4].reshape(-1,1)
TCH = TurbData[:,0].reshape(-1,1)

KA = ExcData[:,0].reshape(-1,1)
TA = ExcData[:,1].reshape(-1,1)
KE = ExcData[:,2].reshape(-1,1)
TE = ExcData[:,3].reshape(-1,1)
KF = ExcData[:,4].reshape(-1,1)
TF = ExcData[:,5].reshape(-1,1)
AE = ExcData[:,6].reshape(-1,1)
BE = ExcData[:,7].reshape(-1,1)
VRmax = ExcData[:,8].reshape(-1,1)
VRmin = ExcData[:,9].reshape(-1,1)
TR = ExcData[:,10].reshape(-1,1)

VL0 = np.abs(v_sol[LoadBuses]).reshape(-1,1)
PL0 = (bus.tocsr()[LoadBuses, 2]/baseMVA).toarray()
QL0 = (bus.tocsr()[LoadBuses, 3]/baseMVA).toarray()
a1 = np.array([0]).reshape(-1,1)
a2 = np.array([0]).reshape(-1,1)
a3 = np.array([1]).reshape(-1,1)
b1 = np.array([0]).reshape(-1,1)
b2 = np.array([0]).reshape(-1,1)
b3 = np.array([1]).reshape(-1,1)
TLr = np.array([0.1]).reshape(-1,1)
TLi = np.array([0.1]).reshape(-1,1)

ui = inputs.reshape(-1,1)
xi = states.reshape(-1,1)
d_delta = d_delta
d_thres = np.array([0.005]).reshape(-1)[0]
d_ref = d_ref
A = A
B = B

##
Xdpp = GenData[:,3].reshape(-1,1)

##

print("Read PowerModel.py")

