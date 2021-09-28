'''
Builds the bus admittance matrix and branch admittance matrices.
Returns the full bus admittance matrix (i.e. for all buses) and the
matrices C{Yf} and C{Yt} which, when multiplied by a complex voltage
vector, yield the vector currents injected into each line from the
"from" and "to" buses respectively of each line. Does appropriate
conversions to p.u.
@see: L{makeSbus}
@author: Ray Zimmerman (PSERC Cornell)
@author: Richard Lincoln
Modified by Oak Ridge National Laboratory (Byungkwon Park) to be used in the parareal algorithm
'''

import numpy as np
from scipy.sparse import csr_matrix, csc_matrix, lil_matrix, identity
from numpy import ones, conj, nonzero, any, exp, pi, hstack, real, vstack

def makeYbus(baseMVA, bus, branch):
    ## constants
    nb = bus.shape[0]     ## number of buses
    nl = branch.shape[0]  ## number of lines

    ## for each branch, compute the elements of the branch admittance matrix where
    ##
    ##      | If |   | Yff  Yft |   | Vf |
    ##      |    | = |          | * |    |
    ##      | It |   | Ytf  Ytt |   | Vt |
    ##
    Ytt, Yff, Yft, Ytf = branch_vectors(branch, nl)
    ## compute shunt admittance
    ## if Psh is the real power consumed by the shunt at V = 1.0 p.u.
    ## and Qsh is the reactive power injected by the shunt at V = 1.0 p.u.
    ## then Psh - j Qsh = V * conj(Ysh * V) = conj(Ysh) = Gs - j Bs,
    ## i.e. Ysh = Psh + j Qsh, so ...
    ## vector of shunt admittances
    Ysh = (bus[:, 4].toarray() + 1j * bus[:, 5].toarray()) / baseMVA

    ## build connection matrices
    f = real(branch[:, 0].toarray().reshape(-1)-1).astype(int)  ## list of "from" buses
    t = real(branch[:, 1].toarray().reshape(-1)-1).astype(int)  ## list of "to" buses
    ## connection matrix for line & from buses
    Cf = csr_matrix((ones(nl), (np.arange(nl), f)), (nl, nb))
    ## connection matrix for line & to buses
    Ct = csr_matrix((ones(nl), (range(nl), t)), (nl, nb))

    ## build Yf and Yt such that Yf * V is the vector of complex branch currents injected
    ## at each branch's "from" bus, and Yt is the same for the "to" bus end
    i = hstack([range(nl), range(nl)])  ## double set of row indices
    Yf = csr_matrix( (hstack([Yff.reshape(-1), Yft.reshape(-1)]), (i, hstack([f, t])) ), (nl, nb))
    Yt = csr_matrix( (hstack([Ytf.reshape(-1), Ytt.reshape(-1)]), (i, hstack([f, t])) ), (nl, nb))
    # Yf = spdiags(Yff, 0, nl, nl) * Cf + spdiags(Yft, 0, nl, nl) * Ct
    # Yt = spdiags(Ytf, 0, nl, nl) * Cf + spdiags(Ytt, 0, nl, nl) * Ct

    ## build Ybus
    Ybus = Cf.T * Yf + Ct.T*Yt + csr_matrix((Ysh.reshape(-1), (range(nb), range(nb))), (nb, nb))

    return Ybus, Yf, Yt


def branch_vectors(branch, nl):
    stat = branch[:, 10]   ## ones at in-service branches
    Ysf = stat.toarray() / (branch[:, 2].toarray() + 1j*branch[:, 3].toarray())  ## series admittance
    Yst = Ysf
    Bc = stat.toarray()*branch[:, 4].toarray() ## line charging susceptance
    tap = ones(nl)          ## default tap ratio = 1
    i = nonzero(real(branch[:, 8]))[0]  ## indices of non-zero tap ratios
    tap[i] = real(branch[i, 8].toarray().reshape(-1))         ## assign non-zero tap ratios
    tap = tap.reshape(-1,1)*exp(1j*pi/180*(branch[:, 9].toarray()))  ## add phase shifters

    Ytt = (Yst + 1j*Bc/2)
    Yff = (Ysf + 1j*Bc/2)/(tap*conj(tap))
    Yft = - Ysf/conj(tap)
    Ytf = - Yst/tap

    return Ytt, Yff, Yft, Ytf


