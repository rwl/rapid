import numpy as np
from initialization.PowerModel import DistLoad_Dynamic, Sbase,  idxDistBus, First_Iteration 
from initialization.PowerModel import idxGov, idxTurb, idxExc, idxGen, idxLoad, idxV, idxI, idxTe, idxs, idxEfd, idxV2, idxV1, idxVR, idxdelta, idxsm, idxsif, idxsih, idxsig, idxsik, idxEdc, idxXadspp, idxXaqspp, idxILr, idxILi, idxGenBus, idxGenBus_new, idxLoadBus, \
                       Xfl, Xhl, Xkl, Xgl, Xl, Ra, Asd, Bsd, siTd, Asq, Bsq, siTq, Xadu, Xaqu, Rf, Rh, Rk, Rg, H, D, Tc, wB, Txd, Txq, \
                       RD, TSV, Psvmax, Psvmin, TCH, Pc, \
                       KA, TA, KE, TE, KF, TF, AE, BE, VRmax, VRmin, TR, Vref, \
                       VL0, PL0, QL0, a1, a2, a3, b1, b2, b3, TLr, TLi, \
                       ui, xi, A, B, states, inputs, d_delta, d_ref, d_thres, \
                       nobus, nogen, GENMVA, GenStatus, \
                       get_Ybus

def Eq_StatorAlgebraic(X, Vbus):
##        stator algebraic equations solver
        delta = X[idxdelta]
        sif = X[idxsif]
        sih = X[idxsih]
        sig = X[idxsig]
        sik = X[idxsik]
        Xadspp = X[idxXadspp]
        Xaqspp = X[idxXaqspp]
        Eqpp = Xadspp*(sif/Xfl+sih/Xhl)
        Edpp = -Xaqspp*(sig/Xgl+sik/Xkl)
        Xdspp = Xadspp+Xl
        Xqspp = Xaqspp+Xl
        Vgen = Vbus[idxGenBus]*np.exp(-1j*delta)
        vq = Vgen.real
        vd = Vgen.imag
        iq =1./(Ra**2+Xdspp*Xqspp)*(Ra*(Eqpp-vq)+Xdspp*(Edpp-vd))
        id = 1./(Ra**2+Xdspp*Xqspp)*(-Xqspp*(Eqpp-vq)+Ra*(Edpp-vd))
        Te = Eqpp*iq+Edpp*id+id*iq*(Xadspp-Xaqspp)
        siad = Xadspp*id+Eqpp
        siaq = Xaqspp*iq-Edpp
        siat = np.absolute(Vgen+(Ra+1j*Xl)*(iq+1j*id))
        siId = Asd*np.exp(Bsd*(siat-siTd))
        siIq = Asq*np.exp(Bsq*(siat-siTq))
        Ksd = siat/(siat+siId)
        Ksq = siat/(siat+siIq)
        Xads = Ksd*Xadu
        Xaqs = Ksq*Xaqu
        Fd = 1./(1./Xads+1./Xfl+1./Xhl)
        Fq = 1./(1./Xaqs+1./Xkl+1./Xgl)
        return Te, iq, id, siaq, siad, Fq, Fd

def Eq_LoadAlgebraic(Vbus):
##        load algebraic equations solver
        Vload = Vbus[idxLoadBus]
        VL = np.absolute(Vload)
        YL = (PL0-1j*QL0)/VL0**2
        PL = a1*PL0+a2*(PL0/VL0)*VL+a3*(PL0/VL0**2)*VL**2
        QL = b1*QL0+b2*(QL0/VL0)*VL+b3*(QL0/VL0**2)*VL**2
        SL = PL+1j*QL
        iL = Vload*YL-np.conj(SL/Vload)
        Fr = iL.real
        Fi = iL.imag
        return Fr, Fi

def NWAlgebraic(X, I_Dist_Dynamic):
##        network algebraic equations solver
        SYSMVA = 100.
        Ibus = np.zeros((nobus,1), dtype=complex)
        delta = X[idxdelta]
        sif = X[idxsif]
        sih = X[idxsih]
        sig = X[idxsig]
        sik = X[idxsik]
        Edc = X[idxEdc]
        Xadspp = X[idxXadspp]
        Xaqspp = X[idxXaqspp]
        ILr = X[idxILr]
        ILi = X[idxILi]
        IL = ILr+1j*ILi
        Ibus[idxLoadBus] = -IL

        if len(idxDistBus) > 0: 
                Ibus[idxDistBus] = Ibus[idxDistBus] + I_Dist_Dynamic  

        Eqpp = Xadspp*(sif/Xfl+sih/Xhl)
        Edpp = -Xaqspp*(sig/Xgl+sik/Xkl)
        Xdspp = Xadspp+Xl
        I1 = (Eqpp+1j*(Edpp+Edc))/(Ra+1j*Xdspp)*np.exp(1j*delta)*GENMVA/SYSMVA
        I1 = I1*GenStatus

        I1r = I1.real[:,0]
        I1i = I1.imag[:,0]
        I1r_new = np.bincount(idxGenBus, weights=I1r)
        I1i_new = np.bincount(idxGenBus, weights=I1i)
        I1_new = I1r_new+1j*I1i_new
        Ibus[idxGenBus_new] = np.reshape(I1_new[idxGenBus_new], (-1, 1))

        sbus = get_Ybus()
        Vbus = sbus.solve(Ibus)

        return Vbus , Ibus

def Eq_SteamGov(X,sm):
##        governor differential euqation
        Psv = X[idxGov]
        Psvdot = 1./TSV*(-Psv+Pc-1./RD*sm)
        return Psvdot

def Eq_SteamTurb(X,Psv):
##        turbine differential euqation
        Tm = X[idxTurb]
        Tmdot = 1./TCH*(-Tm+Psv)
        return Tmdot

def Eq_Exc(X,Vbus):
##        exciter differential euqations
        Efd = X[idxEfd]
        V2 = X[idxV2]
        V1 = X[idxV1]
        VR = X[idxVR]
        #for i in range(0,nogen):
        #        if VR[i]>=VRmax[i]:
        #                VR[i] = VRmax[i]
        #        elif VR[i]<=VRmin[i]:
        #                VR[i] = VRmin[i]
        VR=np.maximum(np.minimum(VR, VRmax), VRmin)
        Efddot = 1./TE*(-(KE+AE*np.exp(BE*Efd))*Efd+VR)
        V2dot = 1./TF*(-V2+KF/TF*Efd)
        Vt = np.absolute(Vbus[idxGenBus])
        V1dot = np.zeros((nogen,1))
        #for i in range(0,nogen):
        #        if TR[i]!=0:
        #                V1dot[i] = 1./TR[i]*(-V1[i]+Vt[i])
        #        else:
        #                V1dot[i] = 0
        #                V1[i] = Vt[i]
        idxTRN0 = TR[:,0]!=0
        V1dot[idxTRN0] = 1./TR[idxTRN0]*(-V1[idxTRN0]+Vt[idxTRN0])
        idxTR0=TR[:,0] == 0
        V1[idxTR0] = Vt[idxTR0]
        Verr = Vref-V1
        VF = KF/TF*Efd-V2
        FR = 1./TA*(-VR+KA*(Verr-VF))
        VRdot = FR
        for i in range(0,nogen):
                if VR[i]==VRmax[i] and FR[i] > 0 or VR[i]==VRmin[i] and FR[i] < 0:
                        VRdot[i] = 0

        Excdot=np.r_[Efddot,V2dot,V1dot,VRdot]
        return Excdot

def Eq_Gen(X,Exc,Tm,Te,iq,id,siaq,siad,Fq,Fd):
##        generator differential euqations
        sm = X[idxsm]
        sif = X[idxsif]
        sih = X[idxsih]
        sig = X[idxsig]
        sik = X[idxsik]
        Edc = X[idxEdc]
        Xadspp = X[idxXadspp]
        Xaqspp = X[idxXaqspp]
        Efd = Exc[0:nogen,:]
        Xdspp = Xadspp+Xl
        Xqspp = Xaqspp+Xl
        deltadot = wB*sm
        smdot = (1./(2*H))*(-D*sm+(Tm-Te))
        sifdot = -wB*Rf/Xfl*sif + wB*Rf/Xfl*siad + wB*Rf/Xadu*Efd
        sihdot = -wB*Rh/Xhl*sih+wB*Rh/Xhl*siad
        sigdot = -wB*Rg/Xgl*sig+wB*Rg/Xgl*siaq
        sikdot = -wB*Rk/Xkl*sik+wB*Rk/Xkl*siaq
        Edc_dot = (1./Tc)*(-Edc-(Xqspp-Xdspp)*iq)
        Xadspp = (1./Txd)*(-Xadspp+Fd)
        Xaqspp = (1./Txq)*(-Xaqspp+Fq)
        Gendot=np.r_[deltadot,smdot,sifdot,sihdot,sigdot,sikdot,Edc_dot,Xadspp,Xaqspp]

        return Gendot

def Eq_Load(X,Fr,Fi):
##        load differential euqations
        ILr = X[idxILr]
        ILi = X[idxILi]
        ILrdot = (1./TLr)*(-ILr+Fr)
        ILidot = (1./TLi)*(-ILi+Fi)
        Loaddot = np.r_[ILrdot,ILidot]
        return Loaddot

def fnRK4(t0,t1,nsteps,X0,Vbus0):
##                start the RK4 method
        dt = (t1-t0)/float(nsteps)
        for i in range(0, nsteps):
                Psv0 = X0[idxGov]
                Tm0 = X0[idxTurb]
                Exc0 = X0[idxExc]
                Gen0 = X0[idxGen]
                Load0 = X0[idxLoad]

##                step k1 calculation        
                Te, iq, id, siaq, siad, Fq, Fd = Eq_StatorAlgebraic(X0, Vbus0)
                Fr, Fi = Eq_LoadAlgebraic(Vbus0)
                if len(idxDistBus) > 0: 
                        I_Dist_Dynamic = DistLoad_Dynamic(idxDistBus, np.abs(Vbus0), Sbase, dt) 
                else:
                        I_Dist_Dynamic = 0 

                Vbus, Ibus = NWAlgebraic(X0, I_Dist_Dynamic) 
                sm0 = X0[idxsm]
                k1_Psv = Eq_SteamGov(X0,sm0)
                Psv = Psv0+dt/2.*k1_Psv
                #for i in range(0,nogen):
                #        if Psv[i]>=Psvmax[i]:
                #                Psv[i] = Psvmax[i]
                #        elif Psv[i]<=Psvmin[i]:
                #                Psv[i] = Psvmin[i]
                Psv = np.maximum(np.minimum(Psv, Psvmax), Psvmin)
                k1_Tm = Eq_SteamTurb(X0,Psv)
                Tm = Tm0+dt/2.*k1_Tm
                k1_Exc = Eq_Exc(X0,Vbus)
                Exc = Exc0 + dt/2.*k1_Exc
                k1_Gen = Eq_Gen(X0,Exc,Tm,Te,iq,id,siaq,siad,Fq,Fd)
                Gen = Gen0 + dt/2.*k1_Gen
                k1_Load = Eq_Load(X0,Fr,Fi)
                Load = Load0+dt/2.*k1_Load
                X = np.r_[Psv,Tm,Exc,Gen,Load]
##                step k2 calculation
                Te, iq, id, siaq, siad, Fq, Fd = Eq_StatorAlgebraic(X, Vbus)
                Fr, Fi = Eq_LoadAlgebraic(Vbus)
                if len(idxDistBus) > 0: 
                        I_Dist_Dynamic = DistLoad_Dynamic(idxDistBus, np.abs(Vbus), Sbase,dt)  

                Vbus, Ibus = NWAlgebraic(X,I_Dist_Dynamic) 
                sm = X[idxsm]
                k2_Psv = Eq_SteamGov(X,sm)
                Psv = Psv0+dt/2.*k2_Psv
                #for i in range(0,nogen):
                #        if Psv[i]>=Psvmax[i]:
                #                Psv[i] = Psvmax[i]
                #        elif Psv[i]<=Psvmin[i]:
                #                Psv[i] = Psvmin[i]
                Psv = np.maximum(np.minimum(Psv, Psvmax), Psvmin)
                k2_Tm = Eq_SteamTurb(X,Psv)
                Tm = Tm0 + dt/2.*k2_Tm
                k2_Exc = Eq_Exc(X,Vbus)
                Exc = Exc0+dt/2.*k2_Exc
                k2_Gen = Eq_Gen(X,Exc,Tm,Te,iq,id,siaq,siad,Fq,Fd)
                Gen = Gen0+dt/2.*k2_Gen
                k2_Load = Eq_Load(X,Fr,Fi)
                Load = Load0+dt/2.*k2_Load
                X = np.r_[Psv,Tm,Exc,Gen,Load]
##                step k3 calculation   
                Te, iq, id, siaq, siad, Fq, Fd = Eq_StatorAlgebraic(X, Vbus)
                Fr, Fi = Eq_LoadAlgebraic(Vbus)
                if len(idxDistBus) > 0: 
                        I_Dist_Dynamic = DistLoad_Dynamic(idxDistBus, np.abs(Vbus), Sbase,dt)  

                Vbus, Ibus = NWAlgebraic(X,I_Dist_Dynamic) 
                sm = X[idxsm]
                k3_Psv = Eq_SteamGov(X,sm)
                Psv = Psv0+dt*k3_Psv
                #for i in range(0,nogen):
                #        if Psv[i]>=Psvmax[i]:
                #                Psv[i] = Psvmax[i]
                #        elif Psv[i]<=Psvmin[i]:
                #                Psv[i] = Psvmin[i]
                Psv=np.maximum(np.minimum(Psv, Psvmax), Psvmin)
                k3_Tm = Eq_SteamTurb(X,Psv)
                Tm = Tm0+dt*k3_Tm
                k3_Exc = Eq_Exc(X,Vbus)
                Exc = Exc0+dt*k3_Exc
                k3_Gen = Eq_Gen(X,Exc,Tm,Te,iq,id,siaq,siad,Fq,Fd)
                Gen = Gen0+dt*k3_Gen
                k3_Load = Eq_Load(X,Fr,Fi)
                Load = Load0+dt*k3_Load
                X = np.r_[Psv,Tm,Exc,Gen,Load]
##                step k4 calculation   
                Te, iq, id, siaq, siad, Fq, Fd = Eq_StatorAlgebraic(X, Vbus)
                Fr, Fi = Eq_LoadAlgebraic(Vbus)
                if len(idxDistBus) > 0: 
                        I_Dist_Dynamic = DistLoad_Dynamic(idxDistBus, np.abs(Vbus), Sbase,dt)  

                Vbus, Ibus = NWAlgebraic(X,I_Dist_Dynamic) 
                sm = X[idxsm]
                k4_Psv = Eq_SteamGov(X,sm)
                Psv = Psv0+dt*(k1_Psv+2.*k2_Psv+2.*k3_Psv+k4_Psv)/6.
                #for i in range(0,nogen):
                #        if Psv[i]>=Psvmax[i]:
                #                Psv[i] = Psvmax[i]
                #        elif Psv[i]<=Psvmin[i]:
                #                Psv[i] = Psvmin[i]
                Psv=np.maximum(np.minimum(Psv, Psvmax), Psvmin)
                k4_Tm = Eq_SteamTurb(X,Psv)
                Tm = Tm0+dt*(k1_Tm+2.*k2_Tm+2.*k3_Tm+k4_Tm)/6.
                k4_Exc = Eq_Exc(X,Vbus)
                Exc = Exc0+dt*(k1_Exc+2.*k2_Exc+2.*k3_Exc+k4_Exc)/6.
                k4_Gen = Eq_Gen(X,Exc,Tm,Te,iq,id,siaq,siad,Fq,Fd)
                Gen = Gen0+dt*(k1_Gen+2.*k2_Gen+2.*k3_Gen+k4_Gen)/6.
                k4_Load = Eq_Load(X,Fr,Fi)
                Load = Load0+dt*(k1_Load+2.*k2_Load+2.*k3_Load+k4_Load)/6.
                X = np.r_[Psv,Tm,Exc,Gen,Load]
##                solution of algebraic equations
                if len(idxDistBus) > 0: 
                        I_Dist_Dynamic = DistLoad_Dynamic(idxDistBus, np.abs(Vbus), Sbase,dt)  

                Vbus, Ibus = NWAlgebraic(X,I_Dist_Dynamic) 
                X0 = X
                Vbus0 = Vbus

        return X, Vbus, Ibus

def fnTrap(t0,t1,nsteps,X0,Vbus0):
##                start the TRAP method
        dt = (t1-t0)/float(nsteps)
        for i in range(0, nsteps):
                Psv0 = X0[idxGov]
                Tm0 = X0[idxTurb]
                Exc0 = X0[idxExc]
                Gen0 = X0[idxGen]
                Load0 = X0[idxLoad]
                sm0 = X0[idxsm]
                Te, iq, id, siaq, siad, Fq, Fd = Eq_StatorAlgebraic(X0, Vbus0)
                Fr, Fi = Eq_LoadAlgebraic(Vbus0)
                if len(idxDistBus) > 0:
                        I_Dist_Dynamic = DistLoad_Dynamic(idxDistBus, np.abs(Vbus0), Sbase,dt)
                else:
                        I_Dist_Dynamic = 0

                Vbus, Ibus = NWAlgebraic(X0,I_Dist_Dynamic)
##                prediction 
                k1_Psv = Eq_SteamGov(X0,sm0)
                k1_Tm = Eq_SteamTurb(X0,Psv0)
                k1_Exc = Eq_Exc(X0,Vbus)
                k1_Gen = Eq_Gen(X0,Exc0,Tm0,Te,iq,id,siaq,siad,Fq,Fd)
                k1_Load = Eq_Load(X0,Fr,Fi)
                k1_Psv1 = Psv0+dt/2.*k1_Psv
                k1_Tm1 = Tm0+dt/2.*k1_Tm
                k1_Exc1 = Exc0+dt/2.*k1_Exc
                k1_Gen1 = Gen0+dt/2.*k1_Gen
                k1_Load1 = Load0+dt/2.*k1_Load
                k1_X1 = np.r_[k1_Psv1,k1_Tm1,k1_Exc1,k1_Gen1,k1_Load1]
                k1_Psv = Eq_SteamGov(k1_X1,sm0)
                k1_Tm = Eq_SteamTurb(k1_X1,Psv0)
                k1_Exc = Eq_Exc(k1_X1,Vbus)
                k1_Gen = Eq_Gen(k1_X1,Exc0,Tm0,Te,iq,id,siaq,siad,Fq,Fd)
                k1_Load = Eq_Load(k1_X1,Fr,Fi)
                Psv1 = Psv0+dt*k1_Psv
                Tm1 = Tm0+dt*k1_Tm
                Exc1 = Exc0+dt*k1_Exc
                Gen1 = Gen0+dt*k1_Gen
                Load1 = Load0+dt*k1_Load
                X1 = np.r_[Psv1,Tm1,Exc1,Gen1,Load1]
##                correction        
                k2_Psv = Eq_SteamGov(X1,sm0)
                k2_Tm = Eq_SteamTurb(X1,Psv0)
                k2_Exc = Eq_Exc(X1,Vbus)
                k2_Gen = Eq_Gen(X1,Exc0,Tm0,Te,iq,id,siaq,siad,Fq,Fd)
                k2_Load = Eq_Load(X1,Fr,Fi)
                Psv = Psv0+dt*(k1_Psv+k2_Psv)/2.
                #for i in range(0,nogen):
                #        if Psv[i]>=Psvmax[i]:
                #                Psv[i] = Psvmax[i]
                #        elif Psv[i]<=Psvmin[i]:
                #                Psv[i] = Psvmin[i]
                Psv=np.maximum(np.minimum(Psv, Psvmax), Psvmin)
                Tm = Tm0+dt*(k1_Tm+k2_Tm)/2.
                Exc = Exc0+dt*(k1_Exc+k2_Exc)/2.
                Gen = Gen0+dt*(k1_Gen+k2_Gen)/2.
                Load = Load0+dt*(k1_Load+k2_Load)/2.
                X = np.r_[Psv,Tm,Exc,Gen,Load]
##                solution of algebraic equations
                if len(idxDistBus) > 0: 
                        I_Dist_Dynamic = DistLoad_Dynamic(idxDistBus, np.abs(Vbus), Sbase,dt)  

                Vbus, Ibus = NWAlgebraic(X,I_Dist_Dynamic) 
                X0 = X
                Vbus0 = Vbus

        return X, Vbus, Ibus

def Eq_ADM(X0,Vbus0,Nt):
##                derive ADM terms
        Psv = X0[idxGov]
        Tm = X0[idxTurb]
        Efd = X0[idxEfd]
        V2 = X0[idxV2]
        V1 = X0[idxV1]
        VR = X0[idxVR]
        delta = X0[idxdelta]
        sm = X0[idxsm]
        sif = X0[idxsif]
        sih = X0[idxsih]
        sig = X0[idxsig]
        sik = X0[idxsik]
        Edc = X0[idxEdc]
        Xadspp = X0[idxXadspp]
        Xaqspp = X0[idxXaqspp]
        ILr = X0[idxILr]
        ILi = X0[idxILi]
        ExcSize = Efd.size

        zz0 = np.r_[Psv, Tm, Efd, V2, V1, VR, delta, sm, sif, sih, sig, sik, Edc, Xadspp, Xaqspp, ILr, ILi]

        # algebraic equations
        Te, iq, id, siaq, siad, Fq, Fd = Eq_StatorAlgebraic(X0, Vbus0)
        Fr, Fi = Eq_LoadAlgebraic(Vbus0)
        if len(idxDistBus) > 0:
                I_Dist_Dynamic = DistLoad_Dynamic(idxDistBus, np.abs(Vbus0), Sbase, dt) 
        else:
                I_Dist_Dynamic = 0

        Vbus, Ibus = NWAlgebraic(X0, I_Dist_Dynamic)

        ## Calculate ADM terms - 1
        ## turbine-governor
        Psv1 = (1./TSV) * (-Psv + Pc - (1./RD)*(sm) )
        Tm1 = (1./TCH) * (-Tm + Psv)

        ## excitation
        VR = np.maximum(np.minimum(VR, VRmax), VRmin)
        A1 = AE * np.exp(BE * Efd)
        Efd1 = (1./TE) * (-1*KE*Efd - A1*Efd + VR)
        V21 = (1./TF) * (-V2 + (KF/TF)*Efd)
        Vt = np.absolute(Vbus[idxGenBus])
        V11 = np.zeros((nogen, 1))
        idxTRN0 = TR[:, 0] != 0
        idxTR0 = TR[:,0]==0
        V1[idxTR0] = Vt[idxTR0]
        V11 = (1./TR) * (-V1 + Vt)
        VR1 = (1./TA)*(-VR + KA*Vref -1*KA*V1 -1*KA*(KF/TF)*Efd + KA*V2)
        for i in range(0, nogen):
                if VR[i] == VRmax[i] and VR1[i] > 0 or VR[i] == VRmin[i] and VR1[i] < 0:
                        VR1[i] = 0

        ## generator
        Xdspp = Xadspp + Xl
        Xqspp = Xaqspp + Xl
        delta1 = wB * sm
        sm1 = (1./(2*H)) * (-1*D*sm + (Tm - Te))
        sif1 = -wB*Rf/Xfl*(sif) + wB*Rf/Xfl*(siad) + wB*Rf/Xadu*(Efd)
        sih1 = -wB*Rh/Xhl*(sih) + wB*Rh/Xhl*(siad)
        sig1 = -wB*Rg/Xgl*(sig) + wB*Rg/Xgl*(siaq)
        sik1 = -wB*Rk/Xkl*(sik) + wB*Rk/Xkl*(siaq)
        Edc1 = (1./Tc) * (-Edc - (Xqspp-Xdspp)*iq)
        Xadspp1 = (1./Txd) * (-Xadspp + Fd)
        Xaqspp1 = (1./Txq) * (-Xaqspp + Fq)

        ## load
        ILr1 = (1./TLr)*(-ILr + Fr)
        ILi1 = (1./TLi)*(-ILi + Fi)

        zz1 = np.r_[Psv1, Tm1, Efd1, V21, V11, VR1, delta1, sm1, sif1, sih1, sig1, sik1, Edc1, Xadspp1, Xaqspp1, ILr1, ILi1]

        ## Calculate ADM terms - 2
        ## turbine-governor
        Psv2 = (1/2) *(1./TSV) * (-Psv1 - (1./RD)*(sm1))
        Tm2 = (1/2)*(1./TCH) * (-Tm1 + Psv1)

        ## excitation
        VR1 = np.maximum(np.minimum(VR1, VRmax), VRmin)
        A2 = (A1 + BE*A1*Efd)
        Efd2 = (1./TE)*(-1*KE*Efd1*(1/2) - (1/2)*A2*Efd1 + VR1*(1/2))
        V22 = (1./TF)*(-V21*(1/2) + (1/2)*(KF/TF)*Efd1)
        V12 = (1./TR) * (-V11*(1/2))
        VR2 = (1./TA)*(-VR1*(1/2) -1*KA*V11*(1/2) -1*KA*(KF/TF)*Efd1*(1/2) + KA*V21*(1/2))
        for i in range(0, nogen):
                if VR1[i] == 0:
                        VR2[i] = 0

        ## generator
        Xdspp1 = Xadspp1
        Xqspp1 = Xaqspp1
        delta2 = wB*sm1*(1/2)
        sm2 = (1/2)*(1./(2*H)) * (-D*sm1 + Tm1)
        sif2 = -(wB*Rf)/Xfl * (sif1)*(1/2) + (wB*Rf)/Xadu*(Efd1)*(1/2)
        sih2 = -(wB*Rh)/Xhl * (sih1)*(1/2)
        sig2 = -(wB*Rg)/Xgl * (sig1)*(1/2)
        sik2 = -(wB*Rk)/Xkl * (sik1)*(1/2)
        Edc2 = (1/2)*(1./Tc) * (-Edc1 - (Xqspp1 - Xdspp1)*iq)
        Xadspp2 = (1/2)*(1./Txd) * (-Xadspp1)
        Xaqspp2 = (1/2)*(1./Txq) * (-Xaqspp1)

        ## load
        ILr2 = (1/2)*(1./TLr) * (-ILr1)
        ILi2 = (1/2)*(1./TLi) * (-ILi1)

        zz2 = np.r_[Psv2, Tm2, Efd2, V22, V12, VR2, delta2, sm2, sif2, sih2, sig2, sik2, Edc2, Xadspp2, Xaqspp2, ILr2, ILi2]

        ## Calculate ADM terms - 3
        ## turbine-governor
        Psv3 = (1/3)*(1./TSV)*(-Psv2 - (1./RD)*sm2)
        Tm3 = (1/3)*(1./TCH)*(-Tm2 + Psv2)

        ## excitation
        VR2 = np.maximum(np.minimum(VR2, VRmax), VRmin)
        A3 = (1/2)*(2*BE*A1 + BE**2*A1*Efd)*Efd1**2 + (A1 + BE*A1*Efd)*Efd2
        Efd3 = 1. / TE * (-1*KE*Efd2*(1/3) - (1/3)*A3 + VR2*(1/3))
        V23 = 1. / TF * (-V22*(1/3) + (1/3)*(KF/TF)*Efd2)
        V13 = (1./TR)*(-V12*(1/3))
        VR3 = (1./TA)*(-VR2*(1/3) -1*KA*V12*(1/3) -1*KA*(KF/TF)*Efd2*(1/3) + KA*V22*(1/3))
        for i in range(0, nogen):
                if VR2[i] == 0:
                        VR3[i] = 0

        ## generator
        Xdspp2 = Xadspp2
        Xqspp2 = Xaqspp2
        delta3 = wB*sm2*(1/3)
        sm3 = (1/3)*(1./(2*H)) * (-D*sm2 + Tm2)
        sif3 = -(wB*Rf)/Xfl*(sif2)*(1/3) + (wB*Rf)/Xadu*(Efd2)*(1/3)
        sih3 = -(wB*Rh)/Xhl*(sih2)*(1/3)
        sig3 = -(wB*Rg)/Xgl*(sig2)*(1/3)
        sik3 = -(wB*Rk)/Xkl*(sik2)*(1/3)
        Edc3 = (1/3)*(1./Tc)*(-Edc2 - (Xqspp2 - Xdspp2)*iq)
        Xadspp3 = (1/3)*(1./Txd)*(-Xadspp2)
        Xaqspp3 = (1/3)*(1./Txq)*(-Xaqspp2)

        ## load
        ILr3 = (1/3)*(1./TLr)*(-ILr2)
        ILi3 = (1/3)*(1./TLi)*(-ILi2)

        zz3 = np.r_[Psv3, Tm3, Efd3, V23, V13, VR3, delta3, sm3, sif3, sih3, sig3, sik3, Edc3, Xadspp3, Xaqspp3, ILr3, ILi3]

        ## Calculate ADM terms - 4
        ## turbine-governor
        Psv4 = (1/4)*(1./TSV)*(-Psv3 - (1./RD)*(sm3))
        Tm4 = (1/4)*(1./TCH)*(-Tm3 + Psv3)

        ## excitation
        VR3 = np.maximum(np.minimum(VR3, VRmax), VRmin)
        A4 = (A1 + BE*A1*Efd)*Efd3 + Efd1*Efd2*(2*BE*A1 + BE**2*A1*Efd) + (1/6)*(3*BE**2*A1 + BE**3*A1*Efd)*Efd1**3
        Efd4 = 1./TE * (-1*KE*Efd3*(1/4) - (1/4)*A4 + VR3*(1/4))
        V24 = 1./TF * (-V23*(1/4) + (1/4)*(KF/TF)*Efd3)
        V14 = 1./TR * (-V13*(1/4))
        VR4 = 1./TA*(-VR3*(1/4) -1*KA*V13*(1/4) -1*KA*(KF/TF)*Efd3*(1/4) + KA*V23*(1/4))
        for i in range(0, nogen):
                if VR3[i] == 0:
                        VR4[i] = 0

        ## generator
        Xdspp3 = Xadspp3
        Xqspp3 = Xaqspp3
        delta4 = wB*sm3*(1/4)
        sm4 = (1/4)*(1./(2*H)) * (-D*sm3 + Tm3)
        sif4 = -(wB*Rf)/Xfl*(sif3)*(1/4) + (wB*Rf)/Xadu*(Efd3)*(1/4)
        sih4 = -(wB*Rh)/Xhl*(sih3)*(1/4)
        sig4 = -(wB*Rg)/Xgl*(sig3)*(1/4)
        sik4 = -(wB*Rk)/Xkl*(sik3)*(1/4)
        Edc4 = (1/4)*(1./Tc)*(-Edc3 - (Xqspp3 - Xdspp3)*iq)
        Xadspp4 = (1/4)*(1./Txd)*(-Xadspp3)
        Xaqspp4 = (1/4)*(1./Txq)*(-Xaqspp3)

        ## load
        ILr4 = (1/4)*(1./TLr)*(-ILr3)
        ILi4 = (1/4)*(1./TLi)*(-ILi3)

        zz4 = np.r_[Psv4, Tm4, Efd4, V24, V14, VR4, delta4, sm4, sif4, sih4, sig4, sik4, Edc4, Xadspp4, Xaqspp4, ILr4, ILi4]

        ## Calculate ADM terms - 5
        ## turbine-governor
        Psv5 = (1/5)*(1./TSV)*(-Psv4 - (1./RD)*(sm4))
        Tm5 = (1/5)*(1./TCH)*(-Tm4 + Psv4)

        ## excitation
        Dum_Exc = np.zeros([4*ExcSize, 1], dtype=np.float64)

        ## generator
        Xdspp4 = Xadspp4
        Xqspp4 = Xaqspp4
        delta5 = wB*sm4*(1/5)
        sm5 = (1/5)*(1./(2*H)) * (-D*sm4 + Tm4)
        sif5 = -(wB*Rf)/Xfl*(sif4)*(1/5) + (wB*Rf)/Xadu*(Efd4)*(1/5)
        sih5 = -(wB*Rh)/Xhl*(sih4)*(1/5)
        sig5 = -(wB*Rg)/Xgl*(sig4)*(1/5)
        sik5 = -(wB*Rk)/Xkl*(sik4)*(1/5)
        Edc5 = (1/5)*(1./Tc)*(-Edc4 - (Xqspp4 - Xdspp4)*iq)
        Xadspp5 = (1/5)*(1./Txd)*(-Xadspp4)
        Xaqspp5 = (1/5)*(1./Txq)*(-Xaqspp4)

        ## load
        ILr5 = (1/5)*(1./TLr)*(-ILr4)
        ILi5 = (1/5)*(1./TLi)*(-ILi4)

        zz5 = np.r_[Psv5, Tm5, Dum_Exc, delta5, sm5, sif5, sih5, sig5, sik5, Edc5, Xadspp5, Xaqspp5, ILr5, ILi5]

        ##
        # Determine the number of terms
        nVector = zz1.size
        zz = np.zeros([nVector, Nt], dtype=np.float64)
        zz[:,0:1] = zz0
        zz[:,1:2] = zz1
        zz[:,2:3] = zz2
        zz[:,3:4] = zz3
        zz[:,4:5] = zz4
        zz[:,5:6] = zz5

        # To derive ADM terms more than 5 without the excitation
        for i in range(6, Nt):
                #print(zz[idxGov,i-1:i])  # take i-1 value but include i to not drop the dimension
                zz[idxGov,i:i+1] = (1/i) * (1./TSV)*(-zz[idxGov,i-1:i] - (1./RD)*zz[idxsm,i-1:i])
                zz[idxTurb,i:i+1] = (1/i) * (1./TCH)*(-zz[idxTurb,i-1:i] + zz[idxGov,i-1:i])
                zz[idxExc,i:i+1] = np.zeros([4*ExcSize, 1], dtype=np.float64)

                Xdspp = zz[idxXadspp,i-1:i]
                Xqspp = zz[idxXaqspp,i-1:i]

                zz[idxdelta,i:i+1] = wB*zz[idxsm,i-1:i]*(1/i)
                zz[idxsm,i:i+1] = (1/i)*(1./(2*H))*(-D*zz[idxsm,i-1:i] + zz[idxTurb,i-1:i])
                zz[idxsif,i:i+1] = -(wB*Rf)/Xfl*zz[idxsif,i-1:i]*(1/i)
                zz[idxsih,i:i+1] = -(wB*Rh)/Xhl*zz[idxsih,i-1:i]*(1/i)
                zz[idxsig,i:i+1] = -(wB*Rg)/Xgl*zz[idxsig,i-1:i]*(1/i)
                zz[idxsik,i:i+1] = -(wB*Rk)/Xkl*zz[idxsik,i-1:i]*(1/i)
                zz[idxEdc,i:i+1] = (1/i)*(1./Tc)*(-zz[idxEdc,i-1:i] - (Xqspp - Xdspp)*iq)
                zz[idxXadspp,i:i+1] = (1/i)*(1./Txd)*(-zz[idxXadspp,i-1:i])
                zz[idxXaqspp,i:i+1] = (1/i)*(1./Txq)*(-zz[idxXaqspp,i-1:i])

                ## load
                zz[idxILr,i:i+1] = (1/i) * (1./TLr) * (-zz[idxILr,i-1:i])
                zz[idxILi,i:i+1] = (1/i) * (1./TLi) * (-zz[idxILi,i-1:i])
        #

        return zz

def fnADM(t0,t1,nsteps,X0,Vbus0):
##                start the ADM method
        dt = (t1 - t0) / float(nsteps)

        for i in range(0, nsteps):
                Nt = 11
                zz = Eq_ADM(X0, Vbus0, Nt)

                X = zz[:,0:1]
                for i in range(1, Nt):
                        X = X + zz[:,i:i+1]*dt**i

#                X = zz0 + zz1*dt + zz2*dt**2 + zz3*dt**3 + zz4*dt**4 + zz5*dt**
                Psv = X[idxGov]
                VR = X[idxVR]
                Psv = np.maximum(np.minimum(Psv, Psvmax), Psvmin)
                VR = np.maximum(np.minimum(VR, VRmax), VRmin)
                X[idxGov] = Psv
                X[idxVR] = VR

                ##  solution of algebraic equations
                if len(idxDistBus) > 0: 
                        I_Dist_Dynamic = DistLoad_Dynamic(idxDistBus, np.abs(Vbus0), Sbase,dt)
                else:
                        I_Dist_Dynamic = 0

                Vbus, Ibus = NWAlgebraic(X,I_Dist_Dynamic) 
                X0 = X
                Vbus0 = Vbus

        return X, Vbus, Ibus

def Eq_HAM(X0,Vbus0,h,dt):
##                derive HAM terms
        Psv = X0[idxGov]
        Tm = X0[idxTurb]
        Efd = X0[idxEfd]
        V2 = X0[idxV2]
        V1 = X0[idxV1]
        VR = X0[idxVR]
        delta = X0[idxdelta]
        sm = X0[idxsm]
        sif = X0[idxsif]
        sih = X0[idxsih]
        sig = X0[idxsig]
        sik = X0[idxsik]
        Edc = X0[idxEdc]
        Xadspp = X0[idxXadspp]
        Xaqspp = X0[idxXaqspp]
        ILr = X0[idxILr]
        ILi = X0[idxILi]

        zz0 = np.r_[Psv, Tm, Efd, V2, V1, VR, delta, sm, sif, sih, sig, sik, Edc, Xadspp, Xaqspp, ILr, ILi]

        # algebraic equations
        Te, iq, id, siaq, siad, Fq, Fd = Eq_StatorAlgebraic(X0, Vbus0)
        Fr, Fi = Eq_LoadAlgebraic(Vbus0)

        if len(idxDistBus) > 0:
                I_Dist_Dynamic = DistLoad_Dynamic(idxDistBus, np.abs(Vbus0), Sbase, dt)
        else:
                I_Dist_Dynamic = 0

        Vbus, Ibus = NWAlgebraic(X0, I_Dist_Dynamic)

        ## Calculate HAM terms - 1
        ## turbine-governor
        Psv1 = h*(1./TSV) * (Psv - Pc + (1./RD)*(sm) )*dt
        Tm1 = h*(1./TCH) * (Tm - Psv)*dt

        ## excitation
        VR = np.maximum(np.minimum(VR, VRmax), VRmin)
        A1 = AE * np.exp(BE * Efd)*Efd

        Efd1 = h*(1./TE) * (1*KE*Efd + A1 - VR)*dt
        V21 = h*(1./TF) * (V2 - (KF/TF)*Efd)*dt
        Vt = np.absolute(Vbus[idxGenBus])
        V11 = np.zeros((nogen, 1))
        idxTRN0 = TR[:, 0] != 0
        V11[idxTRN0] = h*(1./TR[idxTRN0]) * (V1[idxTRN0] - Vt[idxTRN0])*dt
        idxTR0 = TR[:,0]==0
        V1[idxTR0] = Vt[idxTR0]
        VR1 = h*(1./TA)*(VR - KA*Vref + 1*KA*V1 + 1*KA*(KF/TF)*Efd - KA*V2)*dt
        for i in range(0, nogen):
                if VR[i] == VRmax[i] and VR1[i] > 0 or VR[i] == VRmin[i] and VR1[i] < 0:
                        VR1[i] = 0

        ## generator
        Xdspp = Xadspp + Xl
        Xqspp = Xaqspp + Xl
        delta1 = -h*wB * sm  * dt
        sm1 = h *(1./(2*H)) * (1*D*sm - (Tm - Te)) * dt
        sif1 = h*(wB*Rf/Xfl*(sif) - wB*Rf/Xfl*(siad) - wB*Rf/Xadu*(Efd))*dt
        sih1 = h*(wB*Rh/Xhl*(sih) - wB*Rh/Xhl*(siad))*dt
        sig1 = h*(wB*Rg/Xgl*(sig) - wB*Rg/Xgl*(siaq))*dt
        sik1 = h*(wB*Rk/Xkl*(sik) - wB*Rk/Xkl*(siaq))*dt
        Edc1 = h*(1./Tc) * (Edc + (Xqspp-Xdspp)*iq)*dt
        Xadspp1 = h*(1./Txd) * (Xadspp - Fd)*dt
        Xaqspp1 = h*(1./Txq) * (Xaqspp - Fq)*dt

        ## load
        ILr1 = h*(1./TLr)*(ILr - Fr)*dt
        ILi1 = h*(1./TLi)*(ILi - Fi)*dt

        zz1 = np.r_[Psv1, Tm1, Efd1, V21, V11, VR1, delta1, sm1, sif1, sih1, sig1, sik1, Edc1, Xadspp1, Xaqspp1, ILr1, ILi1]

        ## Calculate HAM terms - 2
        ## turbine-governor
        Psv2_1 = (1+h)*Psv1
        Psv2_2 = h*(1./TSV) * ((1/2)*Psv1 + (1./RD)*(sm1)*(1/2))*dt
        Psv2 = Psv2_1 + Psv2_2

        Tm2_1 = (1+h)*Tm1
        Tm2_2 = h*(1./TCH) * ((1/2)*Tm1 - (1/2)*Psv1)*dt
        Tm2 = Tm2_1 + Tm2_2

        ## excitation
        A2 = (AE*np.exp(BE*Efd) + BE*AE*np.exp(BE*Efd)*Efd)*Efd1

        Efd2_1 = (1+h)*Efd1
        Efd2_2 = h*(1./TE)*(1*KE*Efd1*(1/2) + (1/2)*A2 - VR1*(1/2))*dt
        Efd2 = Efd2_1 + Efd2_2

        V22_1 = (1+h)*V21
        V22_2 = h*(1./TF)*(V21*(1/2) - (1/2)*(KF/TF)*Efd1)*dt
        V22 = V22_1 + V22_2

        V12_2 = np.zeros((nogen, 1))
        V12_1 = (1+h)*V11
        V12_2[idxTRN0] = h*(1./TR[idxTRN0])*(V11[idxTRN0]*(1/2))*dt
        V12 = V12_1 + V12_2

 #       V11[idxTR0] = Vt[idxTR0]
        VR2_1 = (1+h)*VR1
        VR2_2 = h*((1./TA)*(VR1*(1/2) +1*KA*V11*(1/2) +1*KA*(KF/TF)*Efd1*(1/2) - KA*V21*(1/2)))*dt
        VR2 = VR2_1 + VR2_2
        for i in range(0, nogen):
                if VR1[i] == 0:
                        VR2[i] = 0

        ## generator
        Xdspp1 = Xadspp1 + Xl
        Xqspp1 = Xaqspp1 + Xl

        delta2_1 = (1+h)*delta1
        delta2_2 = - h*wB*sm1*(1/2)*dt
        delta2 = delta2_1 + delta2_2

        sm2_1 = (1+h)*sm1
        sm2_2 = h*(1./(2 * H)) * (D*sm1*(1/2) - Tm1*(1/2) )*dt
        sm2 = sm2_1 + sm2_2

        sif2_1 = (1+h)*sif1
        sif2_2 = h*((wB*Rf)/Xfl * (sif1)*(1/2) - (wB*Rf)/Xadu*(Efd1)*(1/2))*dt
        sif2 = sif2_1 + sif2_2

        sih2_1 = (1+h)*sih1
        sih2_2 = h*((wB*Rh)/Xhl * (sih1)*(1/2))*dt
        sih2 = sih2_1 + sih2_2

        sig2_1 = (1+h)*sig1
        sig2_2 = h*((wB*Rg)/Xgl * (sig1)*(1/2))*dt
        sig2 = sig2_1 + sig2_2

        sik2_1 = (1+h)*sik1
        sik2_2 = h*((wB*Rk)/Xkl * (sik1)*(1/2))*dt
        sik2 = sik2_1 + sik2_2

        Edc2_1 = (1+h)*Edc1
        Edc2_2 = h*((1./Tc) * (Edc1*(1/2) + (Xqspp1 - Xdspp1)*(1/2)*iq))*dt
        Edc2 = Edc2_1 + Edc2_2

        Xadspp2_1 = (1+h)*Xadspp1
        Xadspp2_2 = h*((1./Txd) * (Xadspp1*(1/2)))*dt
        Xadspp2 = Xadspp2_1 + Xadspp2_2

        Xaqspp2_1 = (1+h)*Xaqspp1
        Xaqspp2_2 = h*((1./Txq) * (Xaqspp1*(1/2)))*dt
        Xaqspp2 = Xaqspp2_1 + Xaqspp2_2

        ## load
        ILr2_1 = (1+h)*ILr1
        ILr2_2 = h*((1./TLr) * (ILr1*(1/2)))*dt
        ILr2 = ILr2_1 + ILr2_2

        ILi2_1 = (1+h)*ILi1
        ILi2_2 = h*((1./TLi) * (ILi1*(1/2)))*dt
        ILi2 = ILi2_1 + ILi2_2

        zz2 = np.r_[Psv2, Tm2, Efd2, V22, V12, VR2, delta2, sm2, sif2, sih2, sig2, sik2, Edc2, Xadspp2, Xaqspp2, ILr2, ILi2]

        ## Calculate HAM terms - 3
        ## turbine-governor
        Psv3 = (1+h)*Psv2 +  h*(1./TSV)*(1/2)*Psv2_1*dt + h*(1./TSV)*(1/3)*Psv2_2*dt + h*(1./TSV)*(1./RD)*(1/2)*sm2_1*dt + h*(1./TSV)*(1./RD)*(1/3)*sm2_2*dt
        Tm3 = (1+h)*Tm2 + h*(1./TCH)*(1/2)*Tm2_1*dt + h*(1./TCH)*(1/3)*Tm2_2*dt - h*(1./TCH)*(1/2)*Psv2_1*dt - h*(1./TCH)*(1/3)*Psv2_2*dt

        ## excitation
#        A3 = (1/2)*(2*BE*AE*np.exp(BE*Efd) + BE**2*AE*np.exp(BE*Efd)*Efd)*Efd1**2 + (AE*np.exp(BE*Efd) + BE*AE*np.exp(BE*Efd)*Efd)*Efd2
        Efd3 = (1+h)*Efd2 + h*(1./TE)*(1/2)*KE*Efd2_1*dt  + h*(1./TE)*(1/3)*KE*Efd2_2*dt + h*(1./TE)*(1/2)*( (1/3)*2*BE*AE*np.exp(BE*Efd)*Efd1**2*dt + (1/3)*BE**2*AE*np.exp(BE*Efd)*Efd*Efd1**2*dt + (1/2)*2*AE*np.exp(BE*Efd)*Efd2_1*dt + (1/3)*2*AE*np.exp(BE*Efd)*Efd2_2*dt + (1/2)*2*BE*AE*np.exp(BE*Efd)*Efd*Efd2_1*dt + (1/3)*2*BE*AE*np.exp(BE*Efd)*Efd*Efd2_2*dt ) \
               - h*(1./TE)*(1/2)*VR2_1*dt - h*(1./TE)*(1/3)*VR2_2*dt
        V23 =  (1+h)*V22  + h*(1./TF)*(1/2)*V22_1*dt      + h*(1./TF)*(1/3)*V22_2*dt - h*(1./TF)*(KF/TF)*(1/2)*Efd2_1*dt - h*(1./TF)*(KF/TF)*(1/3)*Efd2_2*dt
        V13 = np.zeros((nogen, 1))
        V13[idxTRN0] = (1+h)*V12[idxTRN0] + h*(1./TR[idxTRN0])*(1/2)*V12_1[idxTRN0]*dt + h*(1./TR[idxTRN0])*(1/3)*V12_2[idxTRN0]*dt
        VR3 = (1+h)*VR2 + h*(1./TA)*(1/2)*VR2_1*dt + h*(1./TA)*(1/3)*VR2_2*dt + h*(1./TA)*KA*(1/2)*V12_1*dt + h*(1./TA)*KA*(1/3)*V12_2*dt + h*(1./TA)*KA*(KF/TF)*(1/2)*Efd2_1*dt + h*(1./TA)*KA*(KF/TF)*(1/3)*Efd2_2*dt - h*(1./TA)*KA*(1/2)*V22_1*dt - h*(1./TA)*KA*(1/3)*V22_2*dt
        for i in range(0, nogen):
                if VR2[i] == 0:
                        VR3[i] = 0

        ## generator
        Xdspp2 = Xadspp2 + Xl
        Xqspp2 = Xaqspp2 + Xl

        delta3 = (1+h)*delta2 - h*wB*(1/2)*sm2_1*dt - h*wB*(1/3)*sm2_2*dt
        sm3 = (1+h)*sm2 + h*(1./(2*H))*(1/2)*(D*sm2_1*dt) + h*(1./(2*H))*(1/3)*(D*sm2_2*dt) - h*(1./(2*H))*(1/2)*(Tm2_1*dt)  - h*(1./(2*H))*(1/3)*(Tm2_2*dt)
        sif3 = (1+h)*sif2 + h*(wB*Rf)/Xfl*(1/2)*(sif2_1)*dt + h*(wB*Rf)/Xfl*(1/3)*(sif2_2)*dt + h*(wB*Rf)/Xadu*(1/2)*(Efd2_1)*dt + h*(wB*Rf)/Xadu*(1/3)*(Efd2_2)*dt
        sih3 = (1+h)*sih2 + h*(wB*Rh)/Xhl*(1/2)*sih2_1*dt + h*(wB*Rh)/Xhl*(1/3)*sih2_2*dt
        sig3 = (1+h)*sig2 + h*(wB*Rg)/Xgl*(1/2)*sig2_1*dt + h*(wB*Rg)/Xgl*(1/3)*sig2_2*dt
        sik3 = (1+h)*sik2 + h*(wB*Rk)/Xkl*(1/2)*sik2_1*dt + h*(wB*Rk)/Xkl*(1/3)*sik2_2*dt
        Edc3 = (1+h)*Edc2 + h*(1./Tc)*(1/2)*Edc2_1*dt + h*(1./Tc)*(1/3)*Edc2_2*dt + h*(1./Tc)*iq*(1/2)*(Xaqspp2_1 - Xadspp2_1)*dt + h*(1./Tc)*iq*(1/3)*(Xaqspp2_2 - Xadspp2_2)*dt
        Xadspp3 = (1+h)*Xadspp2 + h*(1./Txd)*(1/2)*Xadspp2_1*dt + h*(1./Txd)*(1/3)*Xadspp2_2*dt
        Xaqspp3 = (1+h)*Xaqspp2 + h*(1./Txd)*(1/2)*Xaqspp2_1*dt + h*(1./Txd)*(1/3)*Xaqspp2_2*dt

        ## load
        ILr3 = (1+h)*ILr2 + h*(1./TLr)*(1/2)*(ILr2_1)*dt + h*(1./TLr)*(1/3)*(ILr2_2)*dt
        ILi3 = (1+h)*ILi2 + h*(1./TLi)*(1/2)*(ILi2_1)*dt + h*(1./TLi)*(1/3)*(ILi2_2)*dt

        zz3 = np.r_[Psv3, Tm3, Efd3, V23, V13, VR3, delta3, sm3, sif3, sih3, sig3, sik3, Edc3, Xadspp3, Xaqspp3, ILr3, ILi3]

        return zz0, zz1, zz2, zz3

def fnHAM(t0,t1,nsteps,X0,Vbus0):
##                start the HAM method
        dt = (t1 - t0) / float(nsteps)
        for i in range(0, nsteps):
                h = -1;    # one can control h here
                zz0, zz1, zz2, zz3 = Eq_HAM(X0, Vbus0, h, dt)
                X = zz0 + zz1 + zz2 + zz3

                Psv = X[idxGov]
                VR = X[idxVR]
                Psv = np.maximum(np.minimum(Psv, Psvmax), Psvmin)
                VR = np.maximum(np.minimum(VR, VRmax), VRmin)
                X[idxGov] = Psv
                X[idxVR] = VR

                ##  solution of algebraic equations
                if len(idxDistBus) > 0: 
                        I_Dist_Dynamic = DistLoad_Dynamic(idxDistBus, np.abs(Vbus0), Sbase,dt)
                else:
                        I_Dist_Dynamic = 0

                Vbus, Ibus = NWAlgebraic(X,I_Dist_Dynamic) 
                X0 = X
                Vbus0 = Vbus

        return X, Vbus, Ibus


def Eq_Algebraic_adap(X, Vbus):
        delta = X[idxdelta]
        sif = X[idxsif]
        sih = X[idxsih]
        sig = X[idxsig]
        sik = X[idxsik]
        Xadspp = X[idxXadspp]
        Xaqspp = X[idxXaqspp]        
        Eqpp = Xadspp*(sif/Xfl+sih/Xhl)
        Edpp = -Xaqspp*(sig/Xgl+sik/Xkl)
        Xdspp = Xadspp+Xl
        Xqspp = Xaqspp+Xl
        Vgen = Vbus[idxGenBus]*np.exp(-1j*delta)
        vq = Vgen.real
        vd = Vgen.imag
        V = np.absolute(Vbus[idxGenBus])
        theta = np.angle(Vbus[idxGenBus])
        du = np.r_[theta,V]-ui
        iq =1./(Ra**2+Xdspp*Xqspp)*(Ra*(Eqpp-vq)+Xdspp*(Edpp-vd))
        id = 1./(Ra**2+Xdspp*Xqspp)*(-Xqspp*(Eqpp-vq)+Ra*(Edpp-vd))
        Te = Eqpp*iq+Edpp*id+id*iq*(Xadspp-Xaqspp)
        siad = Xadspp*id+Eqpp
        siaq = Xaqspp*iq-Edpp
        Vdc = (Xqspp-Xdspp)*iq
        Vload = Vbus[idxLoadBus]
        VL = np.absolute(Vload)
        YL = (PL0-1j*QL0)/VL0**2
        PL = a1*PL0+a2*(PL0/VL0)*VL+a3*(PL0/VL0**2)*VL**2
        QL = b1*QL0+b2*(QL0/VL0)*VL+b3*(QL0/VL0**2)*VL**2
        SL = PL+1j*QL
        iL = Vload*YL-np.conj(SL/Vload)
        Fr = iL.real
        Fi = iL.imag
        return du, Te, siaq, siad, Vdc, Fr, Fi

def Eq_adap(dx,du,Te,siaq,siad,Vdc,Fr,Fi,d_flag):
        if d_flag:
                X = dx + xi
                Tm = X[idxTurb]
                Efd = X[idxEfd]
                VR = X[idxVR]
                sm = X[idxsm]
                sif = X[idxsif]
                sih = X[idxsih]
                sig = X[idxsig]
                sik = X[idxsik]
                Edc = X[idxEdc]
                ILr = X[idxILr]
                ILi = X[idxILi]
                dxdot = A.dot(dx)+B.dot(du)
                dxdot[idxEfd] = 1./TE*(-(KE+AE*np.exp(BE*Efd))*Efd+VR)
                dxdot[idxsm] = (1./(2.*H))*(-D*sm+(Tm-Te))
                dxdot[idxsif] = -wB*Rf/Xfl*sif + wB*Rf/Xfl*siad + wB*Rf/Xadu*Efd
                dxdot[idxsih] = -wB*Rh/Xhl*sih+wB*Rh/Xhl*siad
                dxdot[idxsig] = -wB*Rg/Xgl*sig+wB*Rg/Xgl*siaq
                dxdot[idxsik] = -wB*Rk/Xkl*sik+wB*Rk/Xkl*siaq
                dxdot[idxEdc] = (1./Tc)*(-Edc-Vdc)
                dxdot[idxILr] = (1./TLr)*(-ILr+Fr)
                dxdot[idxILi] = (1./TLi)*(-ILi+Fi)
        else:
                dxdot = A.dot(dx)+B.dot(du)
        return dxdot

def fnTrap_adap(t0,t1,nsteps,X0,Vbus0):
        dt = (t1-t0)/float(nsteps)
        for i in range(0, nsteps):
                delta0 = X0[idxdelta]
                du, Te, siaq, siad, Vdc, Fr, Fi = Eq_Algebraic_adap(X0, Vbus0)
##                rotor angle deviation check
                d_flag = 1
                if max(abs(delta0-d_delta-delta0[d_ref]))<d_thres:
                        d_flag = 0
##                prediction
                dx0 = X0-xi
                k1_dx = Eq_adap(dx0,du,Te,siaq,siad,Vdc,Fr,Fi,d_flag)
                k1_dx1 = dx0+dt/2.*k1_dx
                k1_dx = Eq_adap(k1_dx1,du,Te,siaq,siad,Vdc,Fr,Fi,d_flag)
                dx1 = dx0+dt*k1_dx
##                correction
                k2_dx = Eq_adap(dx1,du,Te,siaq,siad,Vdc,Fr,Fi,d_flag)
                dx = dx0+dt*(k1_dx+k2_dx)/2.
                X = dx+xi
                Psv = X[idxGov]
                Psv = np.maximum(np.minimum(Psv, Psvmax), Psvmin)
                X[idxGov] = Psv
##                solution of algebraic equations
                if len(idxDistBus) > 0: 
                        I_Dist_Dynamic = DistLoad_Dynamic(idxDistBus, np.abs(Vbus0), Sbase,dt)
                else:
                        I_Dist_Dynamic = 0

                Vbus, Ibus = NWAlgebraic(X,I_Dist_Dynamic) 

                X0 = X
                Vbus0 = Vbus
        return X, Vbus, Ibus

if __name__ == "__main__":
        X, Vbus = fnTrap_adap(0.0, 0.002,1,X0,Vbus0)
##        X, Vbus = fnTrap(0.,0.02,1,X0,Vbus0)
##        X, Vbus = fnTrap_adap(0.,0.02,1,X0,Vbus0)
        print(X)

