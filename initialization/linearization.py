import numpy as np
from scipy import sparse
import scipy.sparse as sp
from scipy.sparse import csr_matrix, csc_matrix, lil_matrix, identity

def linearization_all(GenData, ExcData, satData, TurbData, bus, X0, Vbus0, Ibus0, Pc, idxV, d_ref):

    TSV = TurbData[:,2]
    RD = TurbData[:,1]
    TCH = TurbData[:,0]
    AE = ExcData[:,6]
    BE = ExcData[:,7]
    KE = ExcData[:,2]
    TE = ExcData[:,3]
    KF = ExcData[:,4]
    TF = ExcData[:,5]
    KA = ExcData[:,0]
    TA = ExcData[:,1]
    TR = ExcData[:,10]
    Ra = GenData[:,13]
    Xadu = GenData[:,24]
    Xfl = GenData[:,26]
    Xhl = GenData[:,27]
    Xkl = GenData[:,28]
    Xgl = GenData[:,29]
    Xl = GenData[:,14]
    Rf = GenData[:,30]
    Rh = GenData[:,31]
    Rk = GenData[:,32]
    Rg = GenData[:,33]
    wB = 2*np.pi*GenData[:,16]
    H = GenData[:,11]
    D = GenData[:,12]
    Tc = GenData[:,15]
    Txd = 0.01
    Txq = 0.01
    TLr = 0.1
    TLi = 0.1

    nogen = len(GenData)
    LoadBuses = sp.find(bus.tocsc()[:, 2] > 0)[0]
    noload = len(LoadBuses)

    Psv = X0[0:nogen:1, 0]
    Tm = X0[1*nogen:2*nogen, 0]
    Efd = X0[2*nogen:3*nogen, 0]
    V2 = X0[3*nogen:4*nogen, 0]
    V1 = X0[4*nogen:5*nogen, 0]
    VR = X0[5*nogen:6*nogen, 0]
    delta = X0[6*nogen:7*nogen, 0]
    sm = X0[7*nogen:8*nogen, 0]
    sif = X0[8*nogen:9*nogen, 0]
    sih = X0[9*nogen:10*nogen, 0]
    sig = X0[10*nogen:11*nogen, 0]
    sik = X0[11*nogen:12*nogen, 0]
    Edc = X0[12*nogen:13*nogen, 0]
    Xadspp = X0[13*nogen:14*nogen, 0]
    Xaqspp = X0[14*nogen:15*nogen, 0]
    XLoad = X0[15*nogen:18*noload, 0]
    states = np.concatenate((Psv, Tm, Efd, V2, V1, VR, delta, sm, sif, sih, sig, sik, Edc, Xadspp, Xaqspp, XLoad))

    idxGenBus = (GenData[:, 0] - 1).astype('int64')
    Vbus = Vbus0
    Vgen = Vbus[idxGenBus, 0]
    V = np.abs(Vgen)
    theta = np.angle(Vgen)

    inputs = np.concatenate((theta, V))

    A11 = sp.dia_matrix((-1.0/TSV, 0), shape=(len(TSV), len(TSV)), dtype=np.float64).tocsc()
    A12 = csc_matrix((nogen,nogen), dtype=np.float64)
    A13 = csc_matrix((nogen,nogen), dtype=np.float64)
    A14 = csc_matrix((nogen,nogen), dtype=np.float64)
    A15 = csc_matrix((nogen,nogen), dtype=np.float64)
    A16 = csc_matrix((nogen,nogen), dtype=np.float64)
    A17 = csc_matrix((nogen,nogen), dtype=np.float64)
    A18 = sp.dia_matrix((-1.0/(TSV*RD), 0), shape=(len(TSV), len(TSV)), dtype=np.float64).tocsc()
    A19 = csc_matrix((nogen,nogen), dtype=np.float64)
    A1A = csc_matrix((nogen,nogen), dtype=np.float64)
    A1B = csc_matrix((nogen,nogen), dtype=np.float64)
    A1C = csc_matrix((nogen,nogen), dtype=np.float64)
    A1D = csc_matrix((nogen,nogen), dtype=np.float64)
    A1E = csc_matrix((nogen,nogen), dtype=np.float64)
    A1F = csc_matrix((nogen,nogen), dtype=np.float64)

    A21 = sp.dia_matrix((1.0/TCH, 0), shape=(len(TCH), len(TCH)), dtype=np.float64).tocsc()
    A22 = sp.dia_matrix((-1.0/TCH, 0), shape=(len(TCH), len(TCH)), dtype=np.float64).tocsc()
    A23 = csc_matrix((nogen,nogen), dtype=np.float64)
    A24 = csc_matrix((nogen,nogen), dtype=np.float64)
    A25 = csc_matrix((nogen,nogen), dtype=np.float64)
    A26 = csc_matrix((nogen,nogen), dtype=np.float64)
    A27 = csc_matrix((nogen,nogen), dtype=np.float64)
    A28 = csc_matrix((nogen,nogen), dtype=np.float64)
    A29 = csc_matrix((nogen,nogen), dtype=np.float64)
    A2A = csc_matrix((nogen,nogen), dtype=np.float64)
    A2B = csc_matrix((nogen,nogen), dtype=np.float64)
    A2C = csc_matrix((nogen,nogen), dtype=np.float64)
    A2D = csc_matrix((nogen,nogen), dtype=np.float64)
    A2E = csc_matrix((nogen,nogen), dtype=np.float64)
    A2F = csc_matrix((nogen,nogen), dtype=np.float64)

    A31 = csc_matrix((nogen,nogen), dtype=np.float64)
    A32 = csc_matrix((nogen,nogen), dtype=np.float64)
    A33 = sp.dia_matrix( ((-AE*BE*np.exp(BE*Efd)*Efd-KE-AE*np.exp(BE*Efd))/TE, 0), shape=(len(AE), len(AE)), dtype=np.float64).tocsc()
    A34 = csc_matrix((nogen,nogen), dtype=np.float64)
    A35 = csc_matrix((nogen,nogen), dtype=np.float64)
    A36 = sp.dia_matrix((1.0/TE, 0), shape=(len(TE), len(TE)), dtype=np.float64).tocsc()
    A37 = csc_matrix((nogen,nogen), dtype=np.float64)
    A38 = csc_matrix((nogen,nogen), dtype=np.float64)
    A39 = csc_matrix((nogen,nogen), dtype=np.float64)
    A3A = csc_matrix((nogen,nogen), dtype=np.float64)
    A3B = csc_matrix((nogen,nogen), dtype=np.float64)
    A3C = csc_matrix((nogen,nogen), dtype=np.float64)
    A3D = csc_matrix((nogen,nogen), dtype=np.float64)
    A3E = csc_matrix((nogen,nogen), dtype=np.float64)
    A3F = csc_matrix((nogen,nogen), dtype=np.float64)

    A41 = csc_matrix((nogen,nogen), dtype=np.float64)
    A42 = csc_matrix((nogen,nogen), dtype=np.float64)
    A43 = sp.dia_matrix((KF/TF**2, 0), shape=(len(KF),len(KF)), dtype=np.float64).tocsc()
    A44 = sp.dia_matrix((-1.0/TF, 0), shape=(len(TF),len(TF)), dtype=np.float64).tocsc()
    A45 = csc_matrix((nogen,nogen), dtype=np.float64)
    A46 = csc_matrix((nogen,nogen), dtype=np.float64)
    A47 = csc_matrix((nogen,nogen), dtype=np.float64)
    A48 = csc_matrix((nogen,nogen), dtype=np.float64)
    A49 = csc_matrix((nogen,nogen), dtype=np.float64)
    A4A = csc_matrix((nogen,nogen), dtype=np.float64)
    A4B = csc_matrix((nogen,nogen), dtype=np.float64)
    A4C = csc_matrix((nogen,nogen), dtype=np.float64)
    A4D = csc_matrix((nogen,nogen), dtype=np.float64)
    A4E = csc_matrix((nogen,nogen), dtype=np.float64)
    A4F = csc_matrix((nogen,nogen), dtype=np.float64)

    A51 = csc_matrix((nogen,nogen), dtype=np.float64)
    A52 = csc_matrix((nogen,nogen), dtype=np.float64)
    A53 = csc_matrix((nogen,nogen), dtype=np.float64)
    A54 = csc_matrix((nogen,nogen), dtype=np.float64)
    A55 = sp.dia_matrix((-1.0/TR, 0), shape=(len(TR),len(TR)), dtype=np.float64).tocsc()
    A56 = csc_matrix((nogen,nogen), dtype=np.float64)
    A57 = csc_matrix((nogen,nogen), dtype=np.float64)
    A58 = csc_matrix((nogen,nogen), dtype=np.float64)
    A59 = csc_matrix((nogen,nogen), dtype=np.float64)
    A5A = csc_matrix((nogen,nogen), dtype=np.float64)
    A5B = csc_matrix((nogen,nogen), dtype=np.float64)
    A5C = csc_matrix((nogen,nogen), dtype=np.float64)
    A5D = csc_matrix((nogen,nogen), dtype=np.float64)
    A5E = csc_matrix((nogen,nogen), dtype=np.float64)
    A5F = csc_matrix((nogen,nogen), dtype=np.float64)

    A61 = csc_matrix((nogen,nogen), dtype=np.float64)
    A62 = csc_matrix((nogen,nogen), dtype=np.float64)
    A63 = sp.dia_matrix((-KA*KF/(TA*TF), 0), shape=(len(KA),len(KA)), dtype=np.float64).tocsc()
    A64 = sp.dia_matrix((KA/TA, 0), shape=(len(KA),len(KA)), dtype=np.float64).tocsc()
    A65 = sp.dia_matrix((-KA/TA, 0), shape=(len(KA),len(KA)), dtype=np.float64).tocsc()
    A66 = sp.dia_matrix((-1/TA, 0), shape=(len(TA),len(TA)), dtype=np.float64).tocsc()
    A67 = csc_matrix((nogen,nogen), dtype=np.float64)
    A68 = csc_matrix((nogen,nogen), dtype=np.float64)
    A69 = csc_matrix((nogen,nogen), dtype=np.float64)
    A6A = csc_matrix((nogen,nogen), dtype=np.float64)
    A6B = csc_matrix((nogen,nogen), dtype=np.float64)
    A6C = csc_matrix((nogen,nogen), dtype=np.float64)
    A6D = csc_matrix((nogen,nogen), dtype=np.float64)
    A6E = csc_matrix((nogen,nogen), dtype=np.float64)
    A6F = csc_matrix((nogen,nogen), dtype=np.float64)

    A71 = csc_matrix((nogen,nogen), dtype=np.float64)
    A72 = csc_matrix((nogen,nogen), dtype=np.float64)
    A73 = csc_matrix((nogen,nogen), dtype=np.float64)
    A74 = csc_matrix((nogen,nogen), dtype=np.float64)
    A75 = csc_matrix((nogen,nogen), dtype=np.float64)
    A76 = csc_matrix((nogen,nogen), dtype=np.float64)
    A77 = csc_matrix((nogen,nogen), dtype=np.float64)
    A78 = sp.dia_matrix((wB, 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    A79 = csc_matrix((nogen,nogen), dtype=np.float64)
    A7A = csc_matrix((nogen,nogen), dtype=np.float64)
    A7B = csc_matrix((nogen,nogen), dtype=np.float64)
    A7C = csc_matrix((nogen,nogen), dtype=np.float64)
    A7D = csc_matrix((nogen,nogen), dtype=np.float64)
    A7E = csc_matrix((nogen,nogen), dtype=np.float64)
    A7F = csc_matrix((nogen,nogen), dtype=np.float64)

    A81 = csc_matrix((nogen,nogen), dtype=np.float64)
    A82 = sp.dia_matrix((1./(2*H), 0), shape=(len(H),len(H)), dtype=np.float64).tocsc()
    A83 = csc_matrix((nogen,nogen), dtype=np.float64)
    A84 = csc_matrix((nogen,nogen), dtype=np.float64)
    A85 = csc_matrix((nogen,nogen), dtype=np.float64)
    A86 = csc_matrix((nogen,nogen), dtype=np.float64)
    A87 = sp.dia_matrix(((1.0/(2*H)*(-Xadspp*(sif/Xfl+sih/Xhl)*(Ra*V*np.sin(-theta+delta)+(Xadspp+Xl)*V*np.cos(-theta+delta))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))+Xaqspp*(sig/Xgl+sik/Xkl)*(Ra*V*np.cos(-theta+delta)-(Xaqspp+Xl)*V*np.sin(-theta+delta))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))-(Ra*V*np.cos(-theta+delta)-(Xaqspp+Xl)*V*np.sin(-theta+delta))*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))*(Xadspp-Xaqspp)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2-(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))*(Ra*V*np.sin(-theta+delta)+(Xadspp+Xl)*V*np.cos(-theta+delta))*(Xadspp-Xaqspp)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2)), 0), shape=(len(H),len(H)), dtype=np.float64).tocsc()
    A88 = sp.dia_matrix((1.0/(2*H)*-D, 0), shape=(len(H),len(H)), dtype=np.float64).tocsc()
    A89 = sp.dia_matrix(((1.0/(2*H)*(-Xadspp*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))/(Xfl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl)))-Xadspp**2*(sif/Xfl+sih/Xhl)*Ra/(Xfl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl)))-Xaqspp*(sig/Xgl+sik/Xkl)*(Xaqspp+Xl)*Xadspp/(Xfl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl)))+(Xaqspp+Xl)*Xadspp*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))*(Xadspp-Xaqspp)/(Xfl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2)-(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))*Ra*Xadspp*(Xadspp-Xaqspp)/((Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2*Xfl))), 0), shape=(len(H),len(H)), dtype=np.float64).tocsc()
    A8A = sp.dia_matrix(((1.0/(2*H)*(-Xadspp*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))/(Xhl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl)))-Xadspp**2*(sif/Xfl+sih/Xhl)*Ra/(Xhl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl)))-Xaqspp*(sig/Xgl+sik/Xkl)*(Xaqspp+Xl)*Xadspp/(Xhl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl)))+(Xaqspp+Xl)*Xadspp*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))*(Xadspp-Xaqspp)/(Xhl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2)-(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))*Ra*Xadspp*(Xadspp-Xaqspp)/((Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2*Xhl))), 0), shape=(len(H),len(H)), dtype=np.float64).tocsc()
    A8B = sp.dia_matrix(((1.0/(2*H)*(Xadspp*(sif/Xfl+sih/Xhl)*(Xadspp+Xl)*Xaqspp/(Xgl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl)))+Xaqspp*(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))/(Xgl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl)))-Xaqspp**2*(sig/Xgl+sik/Xkl)*Ra/(Xgl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl)))+Ra*Xaqspp*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))*(Xadspp-Xaqspp)/(Xgl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2)+(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))*(Xadspp+Xl)*Xaqspp*(Xadspp-Xaqspp)/((Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2*Xgl))), 0), shape=(len(H),len(H)), dtype=np.float64).tocsc()
    A8C = sp.dia_matrix(((1.0/(2*H)*(Xadspp*(sif/Xfl+sih/Xhl)*(Xadspp+Xl)*Xaqspp/(Xkl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl)))+Xaqspp*(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))/(Xkl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl)))-Xaqspp**2*(sig/Xgl+sik/Xkl)*Ra/(Xkl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl)))+Ra*Xaqspp*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))*(Xadspp-Xaqspp)/(Xkl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2)+(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))*(Xadspp+Xl)*Xaqspp*(Xadspp-Xaqspp)/((Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2*Xkl))), 0), shape=(len(H),len(H)), dtype=np.float64).tocsc()
    A8D = csc_matrix((nogen,nogen), dtype=np.float64)
    A8E = sp.dia_matrix((1.0/(2*H)*(-(sif/Xfl+sih/Xhl)*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))-Xadspp*(sif/Xfl+sih/Xhl)*(Ra*(sif/Xfl+sih/Xhl)-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))+Xadspp*(sif/Xfl+sih/Xhl)*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))*(Xaqspp+Xl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2-Xaqspp*(sig/Xgl+sik/Xkl)*(Xaqspp+Xl)*(sif/Xfl+sih/Xhl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))-Xaqspp*(sig/Xgl+sik/Xkl)*(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))*(Xaqspp+Xl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2+(Xaqspp+Xl)*(sif/Xfl+sih/Xhl)*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))*(Xadspp-Xaqspp)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2+2.0*(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))*(Xadspp-Xaqspp)*(Xaqspp+Xl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**3-(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))*(Ra*(sif/Xfl+sih/Xhl)-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))*(Xadspp-Xaqspp)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2-(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2), 0), shape=(len(H),len(H)), dtype=np.float64).tocsc()
    A8F = sp.dia_matrix((1.0/(2*H)*(-Xadspp*(sif/Xfl+sih/Xhl)*(Xadspp+Xl)*(-sig/Xgl-sik/Xkl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))+Xadspp*(sif/Xfl+sih/Xhl)*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))*(Xadspp+Xl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2+(sig/Xgl+sik/Xkl)*(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))+Xaqspp*(sig/Xgl+sik/Xkl)*(Ra*(-sig/Xgl-sik/Xkl)-Xadspp*(sif/Xfl+sih/Xhl)+V*np.cos(-theta+delta))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))-Xaqspp*(sig/Xgl+sik/Xkl)*(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))*(Xadspp+Xl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2-(Ra*(-sig/Xgl-sik/Xkl)-Xadspp*(sif/Xfl+sih/Xhl)+V*np.cos(-theta+delta))*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))*(Xadspp-Xaqspp)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2+2.0*(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))*(Xadspp-Xaqspp)*(Xadspp+Xl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**3-(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))*(Xadspp+Xl)*(-sig/Xgl-sik/Xkl)*(Xadspp-Xaqspp)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2+(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2), 0), shape=(len(H),len(H)), dtype=np.float64).tocsc()

    A91 = csc_matrix((nogen,nogen), dtype=np.float64)
    A92 = csc_matrix((nogen,nogen), dtype=np.float64)
    A93 = sp.dia_matrix((wB*Rf/Xadu, 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    A94 = csc_matrix((nogen,nogen), dtype=np.float64)
    A95 = csc_matrix((nogen,nogen), dtype=np.float64)
    A96 = csc_matrix((nogen,nogen), dtype=np.float64)
    A97 = sp.dia_matrix((wB*Rf*Xadspp*(Ra*V*np.cos(-theta+delta)-(Xaqspp+Xl)*V*np.sin(-theta+delta))/(Xfl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    A98 = csc_matrix((nogen,nogen), dtype=np.float64)
    A99 = sp.dia_matrix((-wB*Rf/Xfl+wB*Rf*(-Xadspp**2*(Xaqspp+Xl)/(Xfl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl)))+Xadspp/Xfl)/Xfl, 0), shape=(len(H),len(H)), dtype=np.float64).tocsc()
    A9A = sp.dia_matrix((wB*Rf*(-Xadspp**2*(Xaqspp+Xl)/(Xhl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl)))+Xadspp/Xhl)/Xfl, 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    A9B = sp.dia_matrix((-wB*Rf*Xadspp*Ra*Xaqspp/(Xfl*Xgl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    A9C = sp.dia_matrix((-wB*Rf*Xadspp*Ra*Xaqspp/(Xfl*Xkl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    A9D = csc_matrix((nogen,nogen), dtype=np.float64)
    A9E = sp.dia_matrix((wB*Rf*((Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))-Xadspp*(Xaqspp+Xl)*(sif/Xfl+sih/Xhl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))-Xadspp*(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))*(Xaqspp+Xl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2+sif/Xfl+sih/Xhl)/Xfl, 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    A9F = sp.dia_matrix((wB*Rf*(Xadspp*(Ra*(-sig/Xgl-sik/Xkl)-Xadspp*(sif/Xfl+sih/Xhl)+V*np.cos(-theta+delta))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))-Xadspp*(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))*(Xadspp+Xl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2)/Xfl, 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()

    AA1 = csc_matrix((nogen,nogen), dtype=np.float64)
    AA2 = csc_matrix((nogen,nogen), dtype=np.float64)
    AA3 = csc_matrix((nogen,nogen), dtype=np.float64)
    AA4 = csc_matrix((nogen,nogen), dtype=np.float64)
    AA5 = csc_matrix((nogen,nogen), dtype=np.float64)
    AA6 = csc_matrix((nogen,nogen), dtype=np.float64)
    AA7 = sp.dia_matrix((wB*Rh*Xadspp*(Ra*V*np.cos(-theta+delta)-(Xaqspp+Xl)*V*np.sin(-theta+delta))/(Xhl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    AA8 = csc_matrix((nogen,nogen), dtype=np.float64)
    AA9 = sp.dia_matrix((wB*Rh*(-Xadspp**2*(Xaqspp+Xl)/(Xfl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl)))+Xadspp/Xfl)/Xhl, 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    AAA = sp.dia_matrix((-wB*Rh/Xhl+wB*Rh*(-Xadspp**2*(Xaqspp+Xl)/(Xhl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl)))+Xadspp/Xhl)/Xhl, 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    AAB = sp.dia_matrix((-wB*Rh*Xadspp*Ra*Xaqspp/(Xhl*Xgl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    AAC = sp.dia_matrix((-wB*Rh*Xadspp*Ra*Xaqspp/(Xhl*Xkl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    AAD = csc_matrix((nogen,nogen), dtype=np.float64)
    AAE = sp.dia_matrix((wB*Rh*((Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))-Xadspp*(Xaqspp+Xl)*(sif/Xfl+sih/Xhl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))-Xadspp*(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))*(Xaqspp+Xl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2+sif/Xfl+sih/Xhl)/Xhl, 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    AAF = sp.dia_matrix((wB*Rh*(Xadspp*(Ra*(-sig/Xgl-sik/Xkl)-Xadspp*(sif/Xfl+sih/Xhl)+V*np.cos(-theta+delta))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))-Xadspp*(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))*(Xadspp+Xl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2)/Xhl, 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()

    AB1 = csc_matrix((nogen,nogen), dtype=np.float64)
    AB2 = csc_matrix((nogen,nogen), dtype=np.float64)
    AB3 = csc_matrix((nogen,nogen), dtype=np.float64)
    AB4 = csc_matrix((nogen,nogen), dtype=np.float64)
    AB5 = csc_matrix((nogen,nogen), dtype=np.float64)
    AB6 = csc_matrix((nogen,nogen), dtype=np.float64)
    AB7 = sp.dia_matrix((wB*Rg*Xaqspp*(Ra*V*np.sin(-theta+delta)+(Xadspp+Xl)*V*np.cos(-theta+delta))/(Xgl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    AB8 = csc_matrix((nogen,nogen), dtype=np.float64)
    AB9 = sp.dia_matrix((wB*Rg*Xaqspp*Ra*Xadspp/(Xgl*Xfl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    ABA = sp.dia_matrix((wB*Rg*Xaqspp*Ra*Xadspp/(Xgl*Xhl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    ABB = sp.dia_matrix((-wB*Rg/Xgl+wB*Rg*(-Xaqspp**2*(Xadspp+Xl)/(Xgl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl)))+Xaqspp/Xgl)/Xgl, 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    ABC = sp.dia_matrix((wB*Rg*(-Xaqspp**2*(Xadspp+Xl)/(Xkl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl)))+Xaqspp/Xkl)/Xgl, 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    ABD = csc_matrix((nogen,nogen), dtype=np.float64)
    ABE = sp.dia_matrix((wB*Rg*(Xaqspp*(Ra*(sif/Xfl+sih/Xhl)-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))-Xaqspp*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))*(Xaqspp+Xl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2)/Xgl, 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    ABF = sp.dia_matrix((wB*Rg*((Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))+Xaqspp*(Xadspp+Xl)*(-sig/Xgl-sik/Xkl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))-Xaqspp*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))*(Xadspp+Xl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2+sig/Xgl+sik/Xkl)/Xgl, 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()

    AC1 = csc_matrix((nogen,nogen), dtype=np.float64)
    AC2 = csc_matrix((nogen,nogen), dtype=np.float64)
    AC3 = csc_matrix((nogen,nogen), dtype=np.float64)
    AC4 = csc_matrix((nogen,nogen), dtype=np.float64)
    AC5 = csc_matrix((nogen,nogen), dtype=np.float64)
    AC6 = csc_matrix((nogen,nogen), dtype=np.float64)
    AC7 = sp.dia_matrix((wB*Rk*Xaqspp*(Ra*V*np.sin(-theta+delta)+(Xadspp+Xl)*V*np.cos(-theta+delta))/(Xkl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    AC8 = csc_matrix((nogen,nogen), dtype=np.float64)
    AC9 = sp.dia_matrix((wB*Rk*Xaqspp*Ra*Xadspp/(Xkl*Xfl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    ACA = sp.dia_matrix((wB*Rk*Xaqspp*Ra*Xadspp/(Xkl*Xhl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    ACB = sp.dia_matrix((wB*Rk*(-Xaqspp**2*(Xadspp+Xl)/(Xgl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl)))+Xaqspp/Xgl)/Xkl, 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    ACC = sp.dia_matrix((-wB*Rk/Xkl+wB*Rk*(-Xaqspp**2*(Xadspp+Xl)/(Xkl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl)))+Xaqspp/Xkl)/Xkl, 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    ACD = csc_matrix((nogen,nogen), dtype=np.float64)
    ACE = sp.dia_matrix((wB*Rk*(Xaqspp*(Ra*(sif/Xfl+sih/Xhl)-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))-Xaqspp*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))*(Xaqspp+Xl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2)/Xkl, 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    ACF = sp.dia_matrix((wB*Rk*((Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))+Xaqspp*(Xadspp+Xl)*(-sig/Xgl-sik/Xkl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))-Xaqspp*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))*(Xadspp+Xl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2+sig/Xgl+sik/Xkl)/Xkl, 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()

    AD1 = csc_matrix((nogen,nogen), dtype=np.float64)
    AD2 = csc_matrix((nogen,nogen), dtype=np.float64)
    AD3 = csc_matrix((nogen,nogen), dtype=np.float64)
    AD4 = csc_matrix((nogen,nogen), dtype=np.float64)
    AD5 = csc_matrix((nogen,nogen), dtype=np.float64)
    AD6 = csc_matrix((nogen,nogen), dtype=np.float64)
    AD7 = sp.dia_matrix((-(Xaqspp-Xadspp)*(Ra*V*np.sin(-theta+delta)+(Xadspp+Xl)*V*np.cos(-theta+delta))/(Tc*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(Xaqspp),len(Xaqspp)), dtype=np.float64).tocsc()
    AD8 = csc_matrix((nogen,nogen), dtype=np.float64)
    AD9 = sp.dia_matrix((-(Xaqspp-Xadspp)*Ra*Xadspp/(Tc*Xfl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(Xaqspp),len(Xaqspp)), dtype=np.float64).tocsc()
    ADA = sp.dia_matrix((-(Xaqspp-Xadspp)*Ra*Xadspp/(Tc*Xhl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(Xaqspp),len(Xaqspp)), dtype=np.float64).tocsc()
    ADB = sp.dia_matrix(((Xaqspp-Xadspp)*(Xadspp+Xl)*Xaqspp/(Tc*Xgl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(Xaqspp),len(Xaqspp)), dtype=np.float64).tocsc()
    ADC = sp.dia_matrix(((Xaqspp-Xadspp)*(Xadspp+Xl)*Xaqspp/(Tc*Xkl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(Xaqspp),len(Xaqspp)), dtype=np.float64).tocsc()
    ADD = sp.dia_matrix((-1./Tc, 0), shape=(len(Tc),len(Tc)), dtype=np.float64).tocsc()
    ADE = sp.dia_matrix((((Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))-(Xaqspp-Xadspp)*(Ra*(sif/Xfl+sih/Xhl)-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))+(Xaqspp-Xadspp)*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))*(Xaqspp+Xl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2)/Tc, 0), shape=(len(Tc),len(Tc)), dtype=np.float64).tocsc()
    ADF = sp.dia_matrix(((-(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))-(Xaqspp-Xadspp)*(Xadspp+Xl)*(-sig/Xgl-sik/Xkl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))+(Xaqspp-Xadspp)*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))*(Xadspp+Xl)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2)/Tc, 0), shape=(len(Tc),len(Tc)), dtype=np.float64).tocsc()

    Txd_tem = ((-1./Txd)*(np.ones([nogen,1],dtype=int))).reshape(-1)
    AE1 = csc_matrix((nogen,nogen), dtype=np.float64)
    AE2 = csc_matrix((nogen,nogen), dtype=np.float64)
    AE3 = csc_matrix((nogen,nogen), dtype=np.float64)
    AE4 = csc_matrix((nogen,nogen), dtype=np.float64)
    AE5 = csc_matrix((nogen,nogen), dtype=np.float64)
    AE6 = csc_matrix((nogen,nogen), dtype=np.float64)
    AE7 = csc_matrix((nogen,nogen), dtype=np.float64)
    AE8 = csc_matrix((nogen,nogen), dtype=np.float64)
    AE9 = csc_matrix((nogen,nogen), dtype=np.float64)
    AEA = csc_matrix((nogen,nogen), dtype=np.float64)
    AEB = csc_matrix((nogen,nogen), dtype=np.float64)
    AEC = csc_matrix((nogen,nogen), dtype=np.float64)
    AED = csc_matrix((nogen,nogen), dtype=np.float64)
    AEE = sp.dia_matrix((Txd_tem, 0), shape=(nogen,nogen), dtype=np.float64).tocsc()
    AEF = csc_matrix((nogen,nogen), dtype=np.float64)

    Txq_tem = ((-1./Txq)*(np.ones([nogen,1],dtype=int))).reshape(-1)
    AF1 = csc_matrix((nogen,nogen), dtype=np.float64)
    AF2 = csc_matrix((nogen,nogen), dtype=np.float64)
    AF3 = csc_matrix((nogen,nogen), dtype=np.float64)
    AF4 = csc_matrix((nogen,nogen), dtype=np.float64)
    AF5 = csc_matrix((nogen,nogen), dtype=np.float64)
    AF6 = csc_matrix((nogen,nogen), dtype=np.float64)
    AF7 = csc_matrix((nogen,nogen), dtype=np.float64)
    AF8 = csc_matrix((nogen,nogen), dtype=np.float64)
    AF9 = csc_matrix((nogen,nogen), dtype=np.float64)
    AFA = csc_matrix((nogen,nogen), dtype=np.float64)
    AFB = csc_matrix((nogen,nogen), dtype=np.float64)
    AFC = csc_matrix((nogen,nogen), dtype=np.float64)
    AFD = csc_matrix((nogen,nogen), dtype=np.float64)
    AFE = csc_matrix((nogen,nogen), dtype=np.float64)
    AFF = sp.dia_matrix((Txq_tem, 0), shape=(nogen,nogen), dtype=np.float64).tocsc()

    A_111 = sparse.hstack([A11.tocoo(), A12.tocoo(), A13.tocoo(), A14.tocoo(), A15.tocoo(), A16.tocoo(), A17.tocoo(), A18.tocoo(), A19.tocoo(), A1A.tocoo(), A1B.tocoo(), A1C.tocoo(), A1D.tocoo(), A1E.tocoo(), A1F.tocoo()])
    A_112 = sparse.hstack([A21.tocoo(), A22.tocoo(), A23.tocoo(), A24.tocoo(), A25.tocoo(), A26.tocoo(), A27.tocoo(), A28.tocoo(), A29.tocoo(), A2A.tocoo(), A2B.tocoo(), A2C.tocoo(), A2D.tocoo(), A2E.tocoo(), A2F.tocoo()])
    A_113 = sparse.hstack([A31.tocoo(), A32.tocoo(), A33.tocoo(), A34.tocoo(), A35.tocoo(), A36.tocoo(), A37.tocoo(), A38.tocoo(), A39.tocoo(), A3A.tocoo(), A3B.tocoo(), A3C.tocoo(), A3D.tocoo(), A3E.tocoo(), A3F.tocoo()])
    A_114 = sparse.hstack([A41.tocoo(), A42.tocoo(), A43.tocoo(), A44.tocoo(), A45.tocoo(), A46.tocoo(), A47.tocoo(), A48.tocoo(), A49.tocoo(), A4A.tocoo(), A4B.tocoo(), A4C.tocoo(), A4D.tocoo(), A4E.tocoo(), A4F.tocoo()])
    A_115 = sparse.hstack([A51.tocoo(), A52.tocoo(), A53.tocoo(), A54.tocoo(), A55.tocoo(), A56.tocoo(), A57.tocoo(), A58.tocoo(), A59.tocoo(), A5A.tocoo(), A5B.tocoo(), A5C.tocoo(), A5D.tocoo(), A5E.tocoo(), A5F.tocoo()])
    A_116 = sparse.hstack([A61.tocoo(), A62.tocoo(), A63.tocoo(), A64.tocoo(), A65.tocoo(), A66.tocoo(), A67.tocoo(), A68.tocoo(), A69.tocoo(), A6A.tocoo(), A6B.tocoo(), A6C.tocoo(), A6D.tocoo(), A6E.tocoo(), A6F.tocoo()])
    A_117 = sparse.hstack([A71.tocoo(), A72.tocoo(), A73.tocoo(), A74.tocoo(), A75.tocoo(), A76.tocoo(), A77.tocoo(), A78.tocoo(), A79.tocoo(), A7A.tocoo(), A7B.tocoo(), A7C.tocoo(), A7D.tocoo(), A7E.tocoo(), A7F.tocoo()])
    A_118 = sparse.hstack([A81.tocoo(), A82.tocoo(), A83.tocoo(), A84.tocoo(), A85.tocoo(), A86.tocoo(), A87.tocoo(), A88.tocoo(), A89.tocoo(), A8A.tocoo(), A8B.tocoo(), A8C.tocoo(), A8D.tocoo(), A8E.tocoo(), A8F.tocoo()])
    A_119 = sparse.hstack([A91.tocoo(), A92.tocoo(), A93.tocoo(), A94.tocoo(), A95.tocoo(), A96.tocoo(), A97.tocoo(), A98.tocoo(), A99.tocoo(), A9A.tocoo(), A9B.tocoo(), A9C.tocoo(), A9D.tocoo(), A9E.tocoo(), A9F.tocoo()])
    A_11A = sparse.hstack([AA1.tocoo(), AA2.tocoo(), AA3.tocoo(), AA4.tocoo(), AA5.tocoo(), AA6.tocoo(), AA7.tocoo(), AA8.tocoo(), AA9.tocoo(), AAA.tocoo(), AAB.tocoo(), AAC.tocoo(), AAD.tocoo(), AAE.tocoo(), AAF.tocoo()])
    A_11B = sparse.hstack([AB1.tocoo(), AB2.tocoo(), AB3.tocoo(), AB4.tocoo(), AB5.tocoo(), AB6.tocoo(), AB7.tocoo(), AB8.tocoo(), AB9.tocoo(), ABA.tocoo(), ABB.tocoo(), ABC.tocoo(), ABD.tocoo(), ABE.tocoo(), ABF.tocoo()])
    A_11C = sparse.hstack([AC1.tocoo(), AC2.tocoo(), AC3.tocoo(), AC4.tocoo(), AC5.tocoo(), AC6.tocoo(), AC7.tocoo(), AC8.tocoo(), AC9.tocoo(), ACA.tocoo(), ACB.tocoo(), ACC.tocoo(), ACD.tocoo(), ACE.tocoo(), ACF.tocoo()])
    A_11D = sparse.hstack([AD1.tocoo(), AD2.tocoo(), AD3.tocoo(), AD4.tocoo(), AD5.tocoo(), AD6.tocoo(), AD7.tocoo(), AD8.tocoo(), AD9.tocoo(), ADA.tocoo(), ADB.tocoo(), ADC.tocoo(), ADD.tocoo(), ADE.tocoo(), ADF.tocoo()])
    A_11E = sparse.hstack([AE1.tocoo(), AE2.tocoo(), AE3.tocoo(), AE4.tocoo(), AE5.tocoo(), AE6.tocoo(), AE7.tocoo(), AE8.tocoo(), AE9.tocoo(), AEA.tocoo(), AEB.tocoo(), AEC.tocoo(), AED.tocoo(), AEE.tocoo(), AEF.tocoo()])
    A_11F = sparse.hstack([AF1.tocoo(), AF2.tocoo(), AF3.tocoo(), AF4.tocoo(), AF5.tocoo(), AF6.tocoo(), AF7.tocoo(), AF8.tocoo(), AF9.tocoo(), AFA.tocoo(), AFB.tocoo(), AFC.tocoo(), AFD.tocoo(), AFE.tocoo(), AFF.tocoo()])
    A_11 = sparse.vstack([A_111.tocoo(), A_112.tocoo(), A_113.tocoo(), A_114.tocoo(), A_115.tocoo(), A_116.tocoo(), A_117.tocoo(), A_118.tocoo(), A_119.tocoo(), A_11A.tocoo(), A_11B.tocoo(), A_11C.tocoo(), A_11D.tocoo(), A_11E.tocoo(), A_11F.tocoo()])

    A_121 = csc_matrix((nogen,2*noload), dtype=np.float64)
    A_122 = csc_matrix((nogen,2*noload), dtype=np.float64)
    A_123 = csc_matrix((nogen,2*noload), dtype=np.float64)
    A_124 = csc_matrix((nogen,2*noload), dtype=np.float64)
    A_125 = csc_matrix((nogen,2*noload), dtype=np.float64)
    A_126 = csc_matrix((nogen,2*noload), dtype=np.float64)
    A_127 = csc_matrix((nogen,2*noload), dtype=np.float64)
    A_128 = csc_matrix((nogen,2*noload), dtype=np.float64)
    A_129 = csc_matrix((nogen,2*noload), dtype=np.float64)
    A_12A = csc_matrix((nogen,2*noload), dtype=np.float64)
    A_12B = csc_matrix((nogen,2*noload), dtype=np.float64)
    A_12C = csc_matrix((nogen,2*noload), dtype=np.float64)
    A_12D = csc_matrix((nogen,2*noload), dtype=np.float64)
    A_12E = csc_matrix((nogen,2*noload), dtype=np.float64)
    A_12F = csc_matrix((nogen,2*noload), dtype=np.float64)
    A_12 = sparse.vstack([A_121.tocoo(), A_122.tocoo(), A_123.tocoo(), A_124.tocoo(), A_125.tocoo(), A_126.tocoo(), A_127.tocoo(), A_128.tocoo(), A_129.tocoo(), A_12A.tocoo(), A_12B.tocoo(), A_12C.tocoo(), A_12D.tocoo(), A_12E.tocoo(), A_12F.tocoo()])

    A_211 = csc_matrix((2*noload,nogen), dtype=np.float64)
    A_212 = csc_matrix((2*noload,nogen), dtype=np.float64)
    A_213 = csc_matrix((2*noload,nogen), dtype=np.float64)
    A_214 = csc_matrix((2*noload,nogen), dtype=np.float64)
    A_215 = csc_matrix((2*noload,nogen), dtype=np.float64)
    A_216 = csc_matrix((2*noload,nogen), dtype=np.float64)
    A_217 = csc_matrix((2*noload,nogen), dtype=np.float64)
    A_218 = csc_matrix((2*noload,nogen), dtype=np.float64)
    A_219 = csc_matrix((2*noload,nogen), dtype=np.float64)
    A_21A = csc_matrix((2*noload,nogen), dtype=np.float64)
    A_21B = csc_matrix((2*noload,nogen), dtype=np.float64)
    A_21C = csc_matrix((2*noload,nogen), dtype=np.float64)
    A_21D = csc_matrix((2*noload,nogen), dtype=np.float64)
    A_21E = csc_matrix((2*noload,nogen), dtype=np.float64)
    A_21F = csc_matrix((2*noload,nogen), dtype=np.float64)
    A_21 = sparse.hstack([A_211.tocoo(), A_212.tocoo(), A_213.tocoo(), A_214.tocoo(), A_215.tocoo(), A_216.tocoo(), A_217.tocoo(), A_218.tocoo(), A_219.tocoo(), A_21A.tocoo(), A_21B.tocoo(), A_21C.tocoo(), A_21D.tocoo(), A_21E.tocoo(), A_21F.tocoo()])

    TLr_tem = ((-1./TLr)*(np.ones([noload,1],dtype=int))).reshape(-1)
    TLi_tem = ((-1./TLi)*(np.ones([noload,1],dtype=int))).reshape(-1)
    A_22_11 = sp.dia_matrix((TLr_tem, 0), shape=(noload,noload), dtype=np.float64).tocsc()
    A_22_12 = csc_matrix((noload,noload), dtype=np.float64)
    A_22_21 = csc_matrix((noload,noload), dtype=np.float64)
    A_22_22 = sp.dia_matrix((TLi_tem, 0), shape=(noload,noload), dtype=np.float64).tocsc()
    A_221 = sparse.hstack([A_22_11.tocoo(), A_22_12.tocoo()])
    A_222 = sparse.hstack([A_22_21.tocoo(), A_22_22.tocoo()])
    A_22 = sparse.vstack([A_221.tocoo(), A_222.tocoo()])

    A1 = sparse.hstack([A_11.tocoo(), A_12.tocoo()])
    A2 = sparse.hstack([A_21.tocoo(), A_22.tocoo()])
    A = sparse.vstack([A1.tocoo(), A2.tocoo()])

    B11 = csc_matrix((nogen,nogen), dtype=np.float64)
    B12 = csc_matrix((nogen,nogen), dtype=np.float64)
    B21 = csc_matrix((nogen,nogen), dtype=np.float64)
    B22 = csc_matrix((nogen,nogen), dtype=np.float64)
    B31 = csc_matrix((nogen,nogen), dtype=np.float64)
    B32 = csc_matrix((nogen,nogen), dtype=np.float64)
    B41 = csc_matrix((nogen,nogen), dtype=np.float64)
    B42 = csc_matrix((nogen,nogen), dtype=np.float64)
    B51 = csc_matrix((nogen,nogen), dtype=np.float64)
    B52 = sp.dia_matrix((1./TR, 0), shape=(len(TR),len(TR)), dtype=np.float64).tocsc()
    B61 = csc_matrix((nogen,nogen), dtype=np.float64)
    B62 = csc_matrix((nogen,nogen), dtype=np.float64)
    B71 = csc_matrix((nogen,nogen), dtype=np.float64)
    B72 = csc_matrix((nogen,nogen), dtype=np.float64)
    B81 = sp.dia_matrix((1./(2*H)*(-Xadspp*(sif/Xfl+sih/Xhl)*(-Ra*V*np.sin(-theta+delta)-(Xadspp+Xl)*V*np.cos(-theta+delta))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))+Xaqspp*(sig/Xgl+sik/Xkl)*(-Ra*V*np.cos(-theta+delta)+(Xaqspp+Xl)*V*np.sin(-theta+delta))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))-(-Ra*V*np.cos(-theta+delta)+(Xaqspp+Xl)*V*np.sin(-theta+delta))*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))*(Xadspp-Xaqspp)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2-(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))*(-Ra*V*np.sin(-theta+delta)-(Xadspp+Xl)*V*np.cos(-theta+delta))*(Xadspp-Xaqspp)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2), 0), shape=(len(H),len(H)), dtype=np.float64).tocsc()
    B82 = sp.dia_matrix((1./(2*H)*(-Xadspp*(sif/Xfl+sih/Xhl)*(-Ra*np.cos(-theta+delta)+(Xadspp+Xl)*np.sin(-theta+delta))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))+Xaqspp*(sig/Xgl+sik/Xkl)*(Ra*np.sin(-theta+delta)+(Xaqspp+Xl)*np.cos(-theta+delta))/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))-(Ra*np.sin(-theta+delta)+(Xaqspp+Xl)*np.cos(-theta+delta))*(Ra*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta))+(Xadspp+Xl)*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta)))*(Xadspp-Xaqspp)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2-(Ra*(-Xaqspp*(sig/Xgl+sik/Xkl)+V*np.sin(-theta+delta))-(Xaqspp+Xl)*(Xadspp*(sif/Xfl+sih/Xhl)-V*np.cos(-theta+delta)))*(-Ra*np.cos(-theta+delta)+(Xadspp+Xl)*np.sin(-theta+delta))*(Xadspp-Xaqspp)/(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))**2), 0), shape=(len(H),len(H)), dtype=np.float64).tocsc()
    B91 = sp.dia_matrix((wB*Rf*Xadspp*(-Ra*V*np.cos(-theta+delta)+(Xaqspp+Xl)*V*np.sin(-theta+delta))/(Xfl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    B92 = sp.dia_matrix((wB*Rf*Xadspp*(Ra*np.sin(-theta+delta)+(Xaqspp+Xl)*np.cos(-theta+delta))/(Xfl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    BA1 = sp.dia_matrix((wB*Rh*Xadspp*(-Ra*V*np.cos(-theta+delta)+(Xaqspp+Xl)*V*np.sin(-theta+delta))/(Xhl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    BA2 = sp.dia_matrix((wB*Rh*Xadspp*(Ra*np.sin(-theta+delta)+(Xaqspp+Xl)*np.cos(-theta+delta))/(Xhl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    BB1 = sp.dia_matrix((wB*Rg*Xaqspp*(-Ra*V*np.sin(-theta+delta)-(Xadspp+Xl)*V*np.cos(-theta+delta))/(Xgl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    BB2 = sp.dia_matrix((wB*Rg*Xaqspp*(-Ra*np.cos(-theta+delta)+(Xadspp+Xl)*np.sin(-theta+delta))/(Xgl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    BC1 = sp.dia_matrix((wB*Rk*Xaqspp*(-Ra*V*np.sin(-theta+delta)-(Xadspp+Xl)*V*np.cos(-theta+delta))/(Xkl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    BC2 = sp.dia_matrix((wB*Rk*Xaqspp*(-Ra*np.cos(-theta+delta)+(Xadspp+Xl)*np.sin(-theta+delta))/(Xkl*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(wB),len(wB)), dtype=np.float64).tocsc()
    BD1 = sp.dia_matrix((-(Xaqspp-Xadspp)*(-Ra*V*np.sin(-theta+delta)-(Xadspp+Xl)*V*np.cos(-theta+delta))/(Tc*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(Xaqspp),len(Xaqspp)), dtype=np.float64).tocsc()
    BD2 = sp.dia_matrix((-(Xaqspp-Xadspp)*(-Ra*np.cos(-theta+delta)+(Xadspp+Xl)*np.sin(-theta+delta))/(Tc*(Ra**2+(Xadspp+Xl)*(Xaqspp+Xl))), 0), shape=(len(Xaqspp),len(Xaqspp)), dtype=np.float64).tocsc()
    BE1 = csc_matrix((nogen,nogen), dtype=np.float64)
    BE2 = csc_matrix((nogen,nogen), dtype=np.float64)
    BF1 = csc_matrix((nogen,nogen), dtype=np.float64)
    BF2 = csc_matrix((nogen,nogen), dtype=np.float64)

    B_111 = sparse.hstack([B11.tocoo(), B12.tocoo()])
    B_112 = sparse.hstack([B21.tocoo(), B22.tocoo()])
    B_113 = sparse.hstack([B31.tocoo(), B32.tocoo()])
    B_114 = sparse.hstack([B41.tocoo(), B42.tocoo()])
    B_115 = sparse.hstack([B51.tocoo(), B52.tocoo()])
    B_116 = sparse.hstack([B61.tocoo(), B62.tocoo()])
    B_117 = sparse.hstack([B71.tocoo(), B72.tocoo()])
    B_118 = sparse.hstack([B81.tocoo(), B82.tocoo()])
    B_119 = sparse.hstack([B91.tocoo(), B92.tocoo()])
    B_11A = sparse.hstack([BA1.tocoo(), BA2.tocoo()])
    B_11B = sparse.hstack([BB1.tocoo(), BB2.tocoo()])
    B_11C = sparse.hstack([BC1.tocoo(), BC2.tocoo()])
    B_11D = sparse.hstack([BD1.tocoo(), BD2.tocoo()])
    B_11E = sparse.hstack([BE1.tocoo(), BE2.tocoo()])
    B_11F = sparse.hstack([BF1.tocoo(), BF2.tocoo()])
    B_11 = sparse.vstack([B_111.tocoo(), B_112.tocoo(), B_113.tocoo(), B_114.tocoo(), B_115.tocoo(), B_116.tocoo(), B_117.tocoo(), B_118.tocoo(), B_119.tocoo(), B_11A.tocoo(), B_11B.tocoo(), B_11C.tocoo(), B_11D.tocoo(), B_11E.tocoo(), B_11F.tocoo()])

    B_21 = csc_matrix((2*noload,2*nogen), dtype=np.float64)
    B = sparse.vstack([B_11.tocoo(), B_21.tocoo()])

    d_delta = np.zeros((nogen,1), dtype=np.float64)

    for i in range(0,nogen):
        d_delta[i] = delta[i] - delta[d_ref]

    return A, B, states, inputs, d_delta


