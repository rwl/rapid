# Fault function module
# functions to be assigned to faults
# operate on main program data, e.g. Y bus matrix
# similar to fn.py

import numpy as np
from scipy.sparse import csr_matrix, csc_matrix, lil_matrix, identity
import scipy.sparse as sp
from initialization.PowerModel import Ybus0, Ybusbf, branch, nobus, idxGenBus, Ra, Xdpp, SYSMVA, GENMVA, GenStatus, solve_Ybus, info_Ybus, solveYbus, updateYbus1, updateYbus2

################### bus fault functions
# execute during the fault
# probably a dummy function
def bus(t,parm,linetrip):
    # just for illustration
    result = parm[0]
    return result

# execute at the start of the fault
def bus_start(t,parm,linetrip):
    global solveYbus
    # This is where one would modify Y bus matrix that is included from Data_init
    # no need to return Ybus by the function because it is modified in place by the command below

    result = 0
    # Ybus during fault
    Fbus = parm[0]  #### Bus having this fault
    Ybus0[Fbus-1, Fbus-1] = Ybus0[Fbus-1, Fbus-1] + 1/(0.000001) # Modify the Ybus

    solve_Ybus("start")
    print("startfault",solveYbus)
    # solveYbus = splu(Ybus0.tocsc())

    return result

# execute at the end of the fault
def bus_end(t,parm,linetrip):
    global solveYbus
    # This is where you would modify Y bus matrix to what you want it to be after the fault
    # no need to return Ybus by the function because it is modified in place by the command below

    result = 0
    # Ybus after fault
    if linetrip == 1:
        Fbus = parm[0]   #### Fault bus having this fault
        FrBus = Fbus     #### From Bus Number; not valid for this fault so set to default = 1
        tmp = sp.find(branch[:,0] == Fbus)[0]
        if not tmp.size == 0:
            fault_br = tmp[0]
            ToBus = branch[fault_br, 1]  #### to bus number of that faulted branch
        else:
            tmp = sp.find(branch[:, 1] == Fbus)[0] #### find the line indices that have Fbus as the to bus of that line.
            fault_br = tmp[0]
            ToBus = branch[fault_br, 0]  #### from bus number of that faulted branch

        tag = 0
        k = 1
        if not fault_br.size == 0:
            k = fault_br
            Yf = 1/(branch[k,2] + 1j*branch[k,3])  ### Self admittance of the line switched
            Ysf = 1j*branch[k, 4]/2                ### The shunt reactance of the line switched
            tag = 1

        while tag != 1:
            if ([sp.find(branch[k,0] == FrBus)[0] and sp.find(branch[k,1] == ToBus)[0]]) or ([sp.find(branch[k,0] == ToBus)[0] and sp.find(branch[k,1] == FrBus)[0]]):
                Yf = 1/(branch[k, 2] + 1j*branch[k, 3])  ### Self admittance of the line switched
                Ysf = 1j*branch[k, 4] / 2                ### The shunt reactance of the line switched
                tag = 1
            k = k + 1

        #### Remove the self and shunt admittances from the Ybus
        #### diagonal elements corresponding to From bus and To bus diagonal
        #### Make the off - diagonal elements zero.
        Ybus0[FrBus-1, FrBus-1] = Ybusbf[FrBus-1, FrBus-1] - Yf - Ysf
        Ybus0[ToBus-1, ToBus-1] = Ybusbf[ToBus-1, ToBus-1] - Yf - Ysf
        Ybus0[FrBus-1, ToBus-1] = 0
        Ybus0[ToBus-1, FrBus-1] = 0

    elif linetrip == 0:
        Fbus = parm[0]  #### Fault bus having this fault
        Ybus0[Fbus-1, Fbus-1] = Ybusbf[Fbus-1, Fbus-1]

    solve_Ybus("end")
    print("endfault", solveYbus)

    return result

################### line fault functions
# execute during the fault
# probably a dummy function
def line(t,parm,linetrip):
    # just for illustration
    result = parm[0]
    return result

# execute at the start of the fault
def line_start(t,parm,linetrip):

    fault_dist = 0.000001
    fault_br = parm[0]   #### faulted line number
    FrBus = branch[fault_br-1, 0].astype('int64')    #### From Bus Number of this line;
    ToBus = branch[fault_br-1, 1].astype('int64')    #### To Bus Number of this line;

    # Ybus during fault
    U = lil_matrix((nobus, 2), dtype=np.complex64)
    Uk = lil_matrix((nobus-1, 1), dtype=np.complex64)
    kk = fault_br
    frBusOut = FrBus
    toBusOut = ToBus

    Yij = 1/(branch[kk-1, 2] + 1j*branch[kk-1, 3])
    Ysh = 1j*branch[kk-1, 4] / 2

    d = fault_dist/100

    D1 = sp.hstack([-((1 - 1/d)*Yij + (1 - d)*Ysh), Yij])
    D2 = sp.hstack([Yij, -((1 - 1/(1 - d))*Yij + d*Ysh)])
    D = sp.vstack([D1,D2])

    U[frBusOut-1, 0] = 1
    U[toBusOut-1, 1] = 1

    temp = U.tocsc()*D*U.tocsc().transpose()
    #Ybus0 = (Ybusbf + temp)                     # qs9 : this is not correct
    Ybus0[:,:] = (Ybusbf + temp).tolil()[:,:]   # to change the values of Ybus0 permanently. qs9 : this is correct
    #updateYbus1(temp)

    Uk[frBusOut-1, 0] = -1/d*Yij
    Uk[toBusOut-1, 0] = -1/(1 - d)*Yij
    Vk = Uk.tocsc().getH().conjugate()
    ak = (1/d * 1/(1 - d)*Yij + Ysh)

    Ybus0[0:nobus-1:1, nobus-1] = Uk
    Ybus0[nobus-1, 0:nobus-1:1] = Vk
    Ybus0[nobus-1, nobus-1] = ak
    Ybus0[nobus-1, nobus-1] = Ybus0[nobus-1, nobus-1] + 1 / (0.000001)

    solve_Ybus("start")
    print("startfault", solveYbus)

    result = 0
    return result

# execute at the end of the fault
def line_end(t,parm,linetrip):

    #Ybus0 = Ybusbf              # qs9 : this is not correct
    Ybus0[:, :] = Ybusbf[:, :]  # qs9 : this is correct
    #updateYbus2()

    fault_br = parm[0]  #### faulted line number
    #### Ybus After Fault
    if linetrip == 1:   # If Line is tripped
        FrBus = branch[fault_br-1, 0].astype('int64')  #### From Bus Number of this line;
        ToBus = branch[fault_br-1, 1].astype('int64')  #### To Bus Number of this line;
        From = FrBus
        To = ToBus

        k = fault_br-1
        Yf = 1/(branch[k, 2] + 1j*branch[k, 3])  ### Self admittance of the line switched
        Ysf = 1j*branch[k, 4] / 2  ### The shunt reactance of the line switched

        #### Remove the self and shunt admittances from the Ybus
        #### diagonal elements corresponding to From bus and Tobus diagonal
        #### Make the off-diagonal elements zero.
        Ybus0[FrBus-1, FrBus-1] = Ybusbf[FrBus-1, FrBus-1] - Yf - Ysf
        Ybus0[ToBus-1, ToBus-1] = Ybusbf[ToBus-1, ToBus-1] - Yf - Ysf
        Ybus0[FrBus-1, ToBus-1] = 0
        Ybus0[ToBus-1, FrBus-1] = 0

    elif linetrip == 0:
        #Ybus0 = Ybusbf              # qs9 : this is not correct
        Ybus0[:, :] = Ybusbf[:, :]  # qs9 : this is correct
        #updateYbus2()

    solve_Ybus("start")
    print("endfault",solveYbus)

    result = 0
    return result

################### generator fault functions
# execute during the fault
# probably a dummy function
def generator(t,parm,linetrip):
    # just for illustration
    result = parm[0]
    return result

# execute at the start of the fault
def generator_start(t,parm,linetrip):
    ### Generator Admittance is removed at the Tripped Generator
    ### Bus, correspondong to that generator, current will be made zero
    ### while solving the algebraic equations.
    GenNoOut = parm[0]-1  #### Gen number of that faulted gen
    BusGenNoOut = idxGenBus[GenNoOut]  ## disconnected generator bus number
    Ybus0[BusGenNoOut, BusGenNoOut] = Ybusbf[BusGenNoOut, BusGenNoOut] - 1/(1./((Ra[GenNoOut] + 1j*Xdpp[GenNoOut])*SYSMVA/GENMVA[GenNoOut]))
    GenStatus[GenNoOut] = 0  ### the status of the disconnected generator becomes 0

    result = 0
    return result

# execute at the end of the fault
def generator_end(t,parm,linetrip):
    ### Do not need to do anything for the generator_end
    result = parm[0]
    return result

# dictionary of fault functions

# functions that execute during the fault
def fault_dict(var):
    func_dict = {'bus':bus,'line':line,'generator':generator}
    return func_dict.get(var)

# functions that execute at the start of the fault
def fault_start_dict(var):
    func_dict_start = {'bus':bus_start,'line':line_start,'generator':generator_start}
    return func_dict_start.get(var)

# functions that execute at the end of the fault
def fault_end_dict(var):
    func_dict_end = {'bus':bus_end,'line':line_end,'generator':generator_end}
    return func_dict_end.get(var)   # returns the corresponding operator (string) which is defined above
