'''
Solves the power flow using a full Newton's method.
Solves for bus voltages given the full system admittance matrix (for
all buses), the complex bus power injection vector (for all buses),
the initial vector of complex bus voltages, and column vectors with
the lists of bus indices for the swing bus, PV buses, and PQ buses,
respectively. The bus voltage vector contains the set point for
generator (including ref bus) buses, and the reference angle of the
swing bus, as well as an initial guess for remaining magnitudes and
angles.
@see: L{runpf}
@author: Ray Zimmerman (PSERC Cornell)
@author: Richard Lincoln
Modified by University of Kassel (Florian Schaefer) to use numba
Modified by Oak Ridge National Laboratory (Byungkwon Park) to be used in the parareal algoritm
'''

import numpy as np
from scipy.sparse.linalg import spsolve
from initialization.create_jacobian import create_jacobian_matrix, get_fastest_jacobian_function

def newtonpf(baseMVA, Ybus, Sbus, V0, ref, pv, pq, max_it):
    # initialize
    #max_it = 20
    tol = 1e-8
    i = 0
    V = V0
    Va = np.angle(V)
    Vm = abs(V)

    # set up indexing for updating V
    pvpq = np.r_[pv, pq]
    # generate lookup pvpq -> index pvpq (used in createJ)
    pvpq_lookup = np.zeros(max(Ybus.indices) + 1, dtype=int)
    pvpq_lookup[pvpq] = np.arange(len(pvpq))

    # get jacobian function
    createJ = get_fastest_jacobian_function(pvpq, pq)

    npv = len(pv)
    npq = len(pq)
    j1 = 0
    j2 = npv       # j1:j2 - V angle of pv buses
    j3 = j2
    j4 = j2 + npq  # j3:j4 - V angle of pq buses
    j5 = j4
    j6 = j4 + npq  # j5:j6 - V mag of pq buses

    # evaluate F(x0)
    F = _evaluate_Fx(Ybus, V, Sbus, pv, pq)
    converged = _check_for_convergence(F, tol)

    Ybus = Ybus.tocsr()
    J = None

    # do Newton iterations
    while (not converged and i < max_it):
        # update iteration counter
        i = i + 1

        J = create_jacobian_matrix(Ybus, V, pvpq, pq, createJ, pvpq_lookup, npv, npq, numba=1)

        dx = -1*spsolve(J, F, permc_spec=None, use_umfpack=True)

        # update voltage
        if npv:
            Va[pv] = Va[pv] + dx[j1:j2]
        if npq:
            Va[pq] = Va[pq] + dx[j3:j4]
            Vm[pq] = Vm[pq] + dx[j5:j6]

        V = Vm*np.exp(1j*Va)
        Vm = abs(V)       # update Vm and Va again in case
        Va = np.angle(V)  # we wrapped around with a negative Vm

        F = _evaluate_Fx(Ybus, V, Sbus, pv, pq)
        converged = _check_for_convergence(F, tol)

    return V, converged, i


def _evaluate_Fx(Ybus, V, Sbus, pv, pq):
    mis = V * np.conj(Ybus*V) - Sbus
    F = np.r_[mis[pv].real,
              mis[pq].real,
              mis[pq].imag]
    return F

def _check_for_convergence(F, tol):
    # calculate infinity norm
    return np.linalg.norm(F, np.Inf) < tol


