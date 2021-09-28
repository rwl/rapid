import numpy as np
import scipy.sparse as sp

def Initial_Gen(GenData, satData, v_sol, bus, pg_sol, qg_sol, Ybus=[]):
    SYSMVA = 100
    GENMVA = GenData[:,17]
    GenBuses = (GenData[:,0]-1).astype('int64')
    nogen = len(GENMVA)
    nobus = len(v_sol)

    Xd = GenData[:, 1]
    Xdp = GenData[:, 2]
    Xdpp = GenData[:, 3]
    Td0p = GenData[:, 4]
    Td0pp = GenData[:, 5]
    Xq = GenData[:, 6]
    Xqp = GenData[:, 7]
    Xqpp = GenData[:, 8]
    Tq0p = GenData[:, 9]
    Tq0pp = GenData[:, 10]
    H = GenData[:, 11]
    D = GenData[:, 12]
    Ra = GenData[:, 13]
    Xl = GenData[:, 14]
    Tc = GenData[:, 15]
    wB = GenData[:, 16]*2*np.pi

    ### Generator currents and Voltages
    Vg0bar = v_sol[GenBuses]
    Ig0bar = ((pg_sol - 1j*qg_sol)/np.conjugate(Vg0bar))*(SYSMVA/GENMVA)

    # --------------UnsaturatedParameters Calculations - ----------------
    Xadu = Xd - Xl
    Xaqu = Xq - Xl

    Xfl = Xadu*(Xdp - Xl)/(Xd - Xdp)
    Xhl = 1.0*(((Xdpp - Xl)*Xadu*Xfl)/((Xadu*Xfl - (Xfl + Xadu)*(Xdpp - Xl)) + 1e-15))
    RF = 1.0*((Xadu + Xfl)/(Td0p*wB))
    RH = (Xhl + (Xadu*Xfl)/(Xadu + Xfl))/(Td0pp*wB)

    Xgl = 1.0*(Xaqu*(Xqp - Xl)/((Xq - Xqp) + 1e-15))
    Xkl = 1.0*(((Xqpp - Xl)*Xaqu*Xgl)/((Xaqu*Xgl - (Xgl + Xaqu)*(Xqpp - Xl)) + 1e-15))
    RG = (Xaqu + Xgl)/(Tq0p*wB)
    RK = (Xkl + (Xaqu*Xgl)/(Xaqu + Xgl))/(Tq0pp*wB)

    # -----------------  Saturation Calculations - -----------------
    Ksd = np.ones((nogen,), dtype=np.float64)
    Ksq = np.ones((nogen,), dtype=np.float64)

    Ea = Vg0bar + 1j*Xl*Ig0bar
    siaT0 = np.abs(Ea)

    sidT1 = satData[:, 1]
    siaT1 = satData[:, 2]
    siaTu1 = satData[:, 3]
    siaT2 = satData[:, 4]
    siaTu2 = satData[:, 5]

    siId1 = siaTu1 - siaT1
    siId2 = siaTu2 - siaT2
    Bsd = (np.log(siId1) - np.log(siId2))/(siaT1 - siaT2)
    Asd = np.exp(np.log(siId2) - Bsd*(siaT2 - sidT1))

    Asd_sat = np.zeros((nobus,), dtype=np.float64)
    Asd_sat[(satData[:,0]-1).astype('int64')] = Asd    # why i need int here but not above
    Asd = Asd_sat[GenBuses]

    Bsd_sat = np.zeros((nobus,), dtype=np.float64)
    Bsd_sat[(satData[:,0]-1).astype('int64')] = Bsd
    Bsd = Bsd_sat[GenBuses]

    siqT1 = satData[:, 6]
    siaT1 = satData[:, 7]
    siaTu1 = satData[:, 8]
    siaT2 = satData[:, 9]
    siaTu2 = satData[:, 10]

    siIq1 = siaTu1 - siaT1
    siIq2 = siaTu2 - siaT2
    Bsq = (np.log(siIq1) - np.log(siIq2))/(siaT1 - siaT2)
    Asq = np.exp(np.log(siIq2) - Bsq*(siaT2 - siqT1))

    Asq_sat = np.zeros((nobus,), dtype=np.float64)
    Asq_sat[(satData[:,0]-1).astype('int64')] = Asq
    Asq = Asq_sat[GenBuses]

    Bsq_sat = np.zeros((nobus,), dtype=np.float64)
    Bsq_sat[(satData[:,0]-1).astype('int64')] = Bsq
    Bsq = Bsq_sat[GenBuses]

    sidT1_sat = np.zeros((nobus,), dtype=np.float64)
    sidT1_sat[(satData[:,0]-1).astype('int64')] = sidT1
    sidT1 = sidT1_sat[GenBuses]

    siqT1_sat = np.zeros((nobus,), dtype=np.float64)
    siqT1_sat[(satData[:,0]-1).astype('int64')] = siqT1
    siqT1 = siqT1_sat[GenBuses]

    # ----------------Including Saturation Effect - ----------------
    f_sat = np.ones((nobus,), dtype=np.float64)
    f_nosat = GenBuses  # Enter the Gen  no. in [] for no saturation
    f_sat[f_nosat] = 0
    f_sat = f_sat[GenBuses]

    # ----------------Option to make ksq = ksd for salient pole machines----------------
    f_sil = np.ones((nobus,), dtype=np.float64)
    f_qed = []   # Enter the Gen no. in [] for ksq=ksd.
    f_sil[f_qed] = 0
    f_sil = f_sil[GenBuses] > 0   # f_sil = f_sil(GenData(:, 1))
    f_sil = f_sil*1

    Asd = Asd*f_sat
    Asq = Asq*f_sat
    Asq = Asq*f_sil + Asd*(1-f_sil)
    Bsq = Bsq*f_sil + Bsd*(1-f_sil)
    siqT1 = siqT1*f_sil + sidT1*(1-f_sil)

    Asd1 = Asd*((siaT0 > sidT1)*1)
    Asq1 = Asq*((siaT0 > siqT1)*1)

    siId = Asd1*np.exp(Bsd*(siaT0 - sidT1))
    siIq = Asq1*np.exp(Bsq*(siaT0 - siqT1))

    Ksd = siaT0/(siaT0 + siId)
    Ksq = siaT0/(siaT0 + siIq)

    Xad = Ksd*Xadu
    Xaq = Ksq*Xaqu

    Xd = Xad + Xl
    Xq = Xaq + Xl

    Xadpp = 1/((1/Xad) + (1/Xfl) + (1/Xhl))
    Xaqpp = 1/((1/Xaq) + (1/Xgl) + (1/Xkl))

    Xdpp = Xadpp + Xl
    Xqpp = Xaqpp + Xl
    # -------------End of Saturation Calculations - ---------------

    Eqbar = (Vg0bar + 1j*Xq*Ig0bar)
    Eq = np.abs(Eqbar)
    delta0 = np.angle(Eqbar)


    VQg0 = Vg0bar.real
    VDg0 = Vg0bar.imag

    iDg0 = Ig0bar.imag
    iQg0 = Ig0bar.real

    Ig0bardq = Ig0bar*(np.cos(delta0) - 1j*np.sin(delta0))
    iq0 = Ig0bardq.real
    id0 = Ig0bardq.imag

    Vg0bardq = Vg0bar*(np.cos(delta0) - 1j*np.sin(delta0))
    vqg0 = Vg0bardq.real
    vdg0 = Vg0bardq.imag

    # -----------------d and q axes States - -------------------
    sid0 = vqg0
    siq0 = -vdg0

    ifd0 = (sid0 - Xd*id0)/Xad
    Efd0 = Xadu*ifd0
    IFD0 = Xadu*ifd0

    siad0 = sid0 - Xl*id0
    siaq0 = siq0 - Xl*iq0

    sif0 = siad0 + (Xfl/Xadu)*Efd0
    sih0 = siad0

    sig0 = siaq0
    sik0 = siaq0

    # ----------------------Torque - ---------------------------
    Eqdd0 = Xadpp*((sif0/Xfl) + (sih0/Xhl))
    Eddd0 = -Xaqpp*((sig0/Xgl) + (sik0/Xkl))
    Tm0 = Eqdd0*iq0 + Eddd0*id0 + id0*iq0*(Xadpp - Xaqpp)

    # -------------------Dummy coil - --------------------------
    Edummydd0 = -(Xqpp - Xdpp)*iq0

    # generator impedance embedded in the Ybus matrix
    for ii in range(0,nogen):
        Ybus[GenBuses[ii], GenBuses[ii]] = Ybus[GenBuses[ii], GenBuses[ii]] + 1.0/((Ra[ii] + 1j*Xdpp[ii])*SYSMVA/GENMVA[ii])

    # Load Modelled as Impedance and added to the YBUS Matrix
    for ii in range(0,nobus):
        Ybus[ii, ii] = Ybus[ii, ii] + (bus[ii,2]/SYSMVA - 1j*bus[ii,3]/SYSMVA) / (np.abs(v_sol[ii]))**2

    GenPrev = GenData
    GenData = np.zeros((nogen, 34), dtype=np.float64)
    GenData[:, 0:18] = GenPrev     # supplement Gendata
    GenData[:, 18] = Asd
    GenData[:, 19] = Bsd
    GenData[:, 20] = sidT1
    GenData[:, 21] = Asq
    GenData[:, 22] = Bsq
    GenData[:, 23] = siqT1
    GenData[:, 24] = Xadu
    GenData[:, 25] = Xaqu
    GenData[:, 26] = Xfl
    GenData[:, 27] = Xhl
    GenData[:, 28] = Xkl
    GenData[:, 29] = Xgl
    GenData[:, 30] = RF
    GenData[:, 31] = RH
    GenData[:, 32] = RK
    GenData[:, 33] = RG

    Smo = np.zeros((nogen,1), dtype=np.float64)
    Xgen0 = np.concatenate((delta0.reshape(-1,1), Smo, sif0.reshape(-1,1), sih0.reshape(-1,1), sig0.reshape(-1,1), sik0.reshape(-1,1), Edummydd0.reshape(-1,1), \
                            Xadpp.reshape(-1,1), Xaqpp.reshape(-1,1), Efd0.reshape(-1,1), Tm0.reshape(-1,1), iq0.reshape(-1,1), id0.reshape(-1,1), \
                            vqg0.reshape(-1,1), vdg0.reshape(-1,1), Ig0bar.reshape(-1,1)))

    return Xgen0, Ybus, GenData

def Initial_Exc(ExcData, v_gen, Efd0):

    KA = ExcData[:, 0]
    TA = ExcData[:, 1]
    KE = ExcData[:, 2]
    TE = ExcData[:, 3]
    KF = ExcData[:, 4]
    TF = ExcData[:, 5]
    AE = ExcData[:, 6]
    BE = ExcData[:, 7]
    VRmax = ExcData[:, 8]
    VRmin = ExcData[:, 9]

    VR1 = (KE + AE*np.exp(BE*Efd0))*Efd0
    V2 = KF/TF*Efd0
    V1 = np.abs(v_gen)
    VF = KF/TF*Efd0 - V2
    Vrefo = VR1/KA + V1
    VR = KA*(Vrefo - V1 - VF)

    XExc0 = np.concatenate((Efd0, V2, V1, VR, Vrefo)).reshape(-1,1)

    chkmax = np.where(VR > VRmax)[0]
    if not chkmax.size == 0:
        VR[chkmax] = VRmax[chkmax]

    chkmin = np.where(VR < VRmax)[0]
    if not chkmin.size == 0:
        VR[chkmin] = VRmin[chkmin]

    return XExc0

def Initial_Gov(Pc):
    Psv = Pc
    XGov0 = Psv.reshape(-1,1)

    return XGov0

def Initial_Turb(XGov0):
    Tm0 = XGov0
    XTurb0 = Tm0.reshape(-1,1)

    return XTurb0

def Initial_Load(noload):
    XLoad0 = np.zeros((2*noload,1), dtype=np.float64)

    return XLoad0

def Initial_Bus(v_sol, Ybus):
    Ibus = (Ybus.dot(v_sol)).reshape(-1,1)

    return Ibus

def initialization(GenData, ExcData, satData, TurbData, bus, v_sol, pg_sol, qg_sol, Ybus=[]):
    GenBuses = (GenData[:, 0] - 1).astype('int64')
    LoadBuses = sp.find(bus.tocsc()[:,2] > 0)[0]
    nGenST = 9
    nExcST = 4
    nogen = len(GenBuses)
    noload = len(LoadBuses)

    ##
    Xgen0, Ybus, GenData = Initial_Gen(GenData, satData, v_sol, bus, pg_sol, qg_sol, Ybus)
    Pc = Xgen0[(nGenST+1)*nogen:(nGenST+2)*nogen:1, 0]
    Te0 = Pc  # initial torque

    XExc0 = Initial_Exc(ExcData, v_sol[GenBuses], Xgen0[nGenST*nogen:nogen*(nGenST+1):1, 0])
    Vref = XExc0[(nogen*nExcST):nogen*(nExcST+1):1, 0]

    XGov0 = Initial_Gov(Pc)
    XTurb0 = Initial_Turb(XGov0)
    XLoad0 = Initial_Load(noload)
    Ibus = Initial_Bus(v_sol, Ybus)

    ##
    X0 = np.concatenate((XGov0, XTurb0, XExc0[0:nogen*nExcST:1,0].reshape(-1,1), Xgen0[0:nogen*nGenST,0].reshape(-1,1), XLoad0)).real
    Vbus0 = v_sol.reshape(-1,1)
    Ibus0 = Ibus
    Ybus0 = Ybus

    return X0, Vbus0, Ibus0, Ybus0, GenData, Pc, Vref


