
import time
from datetime import datetime
txst = time.perf_counter()
print("START  TIME = ", datetime.now().strftime("%m/%d/%Y %H:%M:%S"))

from pyscf import scf, symm
import numpy as np
mfchk = "../mf-1/mf.chk"
mol, mfx = scf.chkfile.load_scf(mfchk)
mo_coeff = mfx["mo_coeff"]
mo_energy = mfx["mo_energy"]
mo_occ = mfx["mo_occ"]

if mo_coeff[0].ndim == 2:
    ma, mb = mo_coeff

    nalpha = (mol.nelectron + mol.spin) // 2
    nbeta = (mol.nelectron - mol.spin) // 2
    pTa = np.dot(ma[:, :nalpha], ma[:, :nalpha].T)
    pTb = np.dot(mb[:, :nbeta], mb[:, :nbeta].T)
    pav = 0.5 * (pTa + pTb)

    fa = ma @ np.diag(mo_energy[0]) @ ma.T
    fb = mb @ np.diag(mo_energy[1]) @ mb.T
    fav = 0.5 * (fa + fb)
else:
    pav = 0.5 * (mo_coeff @ np.diag(mo_occ) @ mo_coeff.T)
    fav = mo_coeff @ np.diag(mo_energy) @ mo_coeff.T

ova = mol.intor_symmetric('cint1e_ovlp_sph')
nactorb = None
nactelec = None
split_low = 0.0
split_high = 0.0
alpha = False
beta = False
uno = False
average_occ = False
loc_with_pg = False
lde = '../cc-5'
cas_list = [58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120]
split_low = 0.1
split_high = 1.9
uno = True
do_loc = True

import numpy as np
import scipy.linalg

# Active space procedure and Boys/PM-Localization
# Original author : Zhendong Li

def sqrtm(s):
    e, v = np.linalg.eigh(s)
    return np.dot(v * np.sqrt(e), v.T.conj())

def lowdin(s):
    e, v = np.linalg.eigh(s)
    return np.dot(v / np.sqrt(e), v.T.conj())

def loc(mol, mocoeff, tol=1E-6, maxcycle=1000, iop=0):
    part = {}
    for iatom in range(mol.natm):
        part[iatom] = []
    ncgto = 0
    for binfo in mol._bas:
        atom_id = binfo[0]
        lang = binfo[1]
        ncntr = binfo[3]
        nbas = ncntr * (2 * lang + 1)
        part[atom_id] += range(ncgto, ncgto + nbas)
        ncgto += nbas
    partition = []
    for iatom in range(mol.natm):
        partition.append(part[iatom])
    ova = mol.intor_symmetric("cint1e_ovlp_sph")
    print()
    print('[pm_loc_kernel]')
    print(' mocoeff.shape=',mocoeff.shape)
    print(' tol=',tol)
    print(' maxcycle=',maxcycle)
    print(' partition=',len(partition),'\n',partition)
    k = mocoeff.shape[0]
    n = mocoeff.shape[1]
    natom = len(partition)
 
    def genPaij(mol,mocoeff,ova,partition,iop):
        c = mocoeff.copy()
        # Mulliken matrix
        if iop == 0:
            cts = c.T.dot(ova)
            natom = len(partition)
            pija = np.zeros((natom,n,n))
            for iatom in range(natom):
                idx = partition[iatom]
                tmp = np.dot(cts[:,idx],c[idx,:])
                pija[iatom] = 0.5*(tmp+tmp.T)
        # Lowdin
        elif iop == 1:
            s12 = sqrtm(ova)
            s12c = s12.T.dot(c)
            natom = len(partition)
            pija = np.zeros((natom,n,n))
            for iatom in range(natom):
                idx = partition[iatom]
                pija[iatom] = np.dot(s12c[idx,:].T,s12c[idx,:])
        # Boys
        elif iop == 2:
            rmat = mol.intor_symmetric('cint1e_r_sph',3)
            pija = np.zeros((3,n,n))
            for icart in range(3):
                pija[icart] = c.T @ rmat[icart] @ c
        # P[i,j,a]
        pija = pija.transpose(1,2,0).copy()
        return pija
 
    u = np.identity(n)
    pija = genPaij(mol,mocoeff,ova,partition,iop)
 
    # Start
    def funval(pija):
        return np.einsum('iia,iia',pija,pija)
 
    fun = funval(pija)
    print(' initial funval = ',fun)
    for icycle in range(maxcycle):
        delta = 0.0
        # i>j
        ijdx = []
        for i in range(n-1):
            for j in range(i+1,n):
                bij = abs(np.sum(pija[i,j]*(pija[i,i]-pija[j,j])))
                ijdx.append((i,j,bij))
        ijdx = sorted(ijdx,key=lambda x:x[2], reverse=True)
        for i,j,bij in ijdx:
            # determine angle
            vij = pija[i,i]-pija[j,j]
            aij = np.dot(pija[i,j],pija[i,j]) - 0.25*np.dot(vij,vij)
            bij = np.dot(pija[i,j],vij)
            if abs(aij)<1.e-10 and abs(bij)<1.e-10: continue
            p1 = np.sqrt(aij**2+bij**2)
            cos4a = -aij/p1
            sin4a = bij/p1
            cos2a = np.sqrt((1+cos4a)*0.5)
            sin2a = np.sqrt((1-cos4a)*0.5)
            cosa  = np.sqrt((1+cos2a)*0.5)
            sina  = np.sqrt((1-cos2a)*0.5)
            # Why? Because we require alpha in [0,pi/2]
            if sin4a < 0.0:
                cos2a = -cos2a
                sina, cosa = cosa, sina
            # stationary condition
            if abs(cosa-1.0)<1.e-10: continue
            if abs(sina-1.0)<1.e-10: continue
            # incremental value
            delta += p1*(1-cos4a)
            # Transformation
            # Urot
            ui = u[:,i]*cosa+u[:,j]*sina
            uj = -u[:,i]*sina+u[:,j]*cosa
            u[:,i] = ui.copy()
            u[:,j] = uj.copy()
            # Bra-transformation of Integrals
            tmp_ip = pija[i,:,:]*cosa+pija[j,:,:]*sina
            tmp_jp = -pija[i,:,:]*sina+pija[j,:,:]*cosa
            pija[i,:,:] = tmp_ip.copy()
            pija[j,:,:] = tmp_jp.copy()
            # Ket-transformation of Integrals
            tmp_ip = pija[:,i,:]*cosa+pija[:,j,:]*sina
            tmp_jp = -pija[:,i,:]*sina+pija[:,j,:]*cosa
            pija[:,i,:] = tmp_ip.copy()
            pija[:,j,:] = tmp_jp.copy()
        fun = fun+delta
        print('icycle=', icycle, 'delta=', delta, 'fun=', fun)
        if delta < tol:
            break

    # Check
    ierr = 0
    if delta < tol: 
        print('CONG: PMloc converged!')
    else:
        ierr = 1
        print('WARNING: PMloc not converged')
    return ierr, u

def loc_pg(mol, mocoeff, orb_sym):
    assert mocoeff.shape[1] == len(orb_sym)
    ierr = 0
    ru = np.zeros((len(orb_sym), len(orb_sym)), dtype=mocoeff.dtype)
    for isym in set(orb_sym):
        mask = np.array(orb_sym) == isym
        jerr, u = loc(mol, mocoeff[:, mask])
        ru[np.outer(mask, mask)] = u.flatten()
        ierr = ierr | jerr
    return ierr, ru


from pyscf import tools
import numpy as np
import os

for fname in ["mo_coeff.npy", "lo_coeff.npy", "nat_coeff.npy"]:
    if os.path.isfile(lde + "/" + fname):
        print("use: " + lde + "/" + fname)
        coeff = np.load(lde + "/" + fname)
        break

for fname in ["mf_occ.npy", "lo_occ.npy", "nat_occ.npy"]:
    if os.path.isfile(lde + "/" + fname):
        print("use: " + lde + "/" + fname)
        mo_occ = np.load(lde + "/" + fname)
        break

if alpha:
    coeff = coeff[0]
    mo_occ = mo_occ[0]
elif beta:
    coeff = coeff[1]
    mo_occ = mo_occ[1]
elif uno:
    # 1. Read UHF-alpha/beta orbitals from chkfile
    ma, mb = coeff
    norb = ma.shape[1]
    nalpha = (mol.nelectron + mol.spin) // 2
    nbeta  = (mol.nelectron - mol.spin) // 2
    print('Nalpha = %d, Nbeta %d, Sz = %d, Norb = %d' % (nalpha, nbeta, mol.spin, norb))

    # 2. Sanity check, using orthogonality

    ova = mol.intor_symmetric("cint1e_ovlp_sph")
    diff = ma.T @ ova @ ma - np.identity(norb)
    assert np.linalg.norm(diff) < 1E-7
    diff = mb.T @ ova @ mb - np.identity(norb)
    assert np.linalg.norm(diff) < 1E-7

    print('alpha occ = ', mo_occ[0])
    print('beta  occ = ', mo_occ[1])

    pTa = ma @ np.diag(mo_occ[0]) @ ma.T
    pTb = mb @ np.diag(mo_occ[1]) @ mb.T
    pT = 0.5 * (pTa + pTb)

    # Lowdin basis
    s12 = sqrtm(ova)
    s12inv = lowdin(ova)
    pT = s12 @ pT @ s12
    print('Idemponency of DM: %s' % np.linalg.norm(pT.dot(pT) - pT))

    # 'natural' occupations and orbitals
    mo_occ, coeff = np.linalg.eigh(pT)
    mo_occ = 2 * mo_occ
    mo_occ[abs(mo_occ) < 1E-14] = 0.0

    # Rotate back to AO representation and check orthogonality
    coeff = np.dot(s12inv, coeff)
    diff = coeff.T @ ova @ coeff - np.identity(norb)
    assert np.linalg.norm(diff) < 1E-7

    index = np.argsort(-mo_occ)
    mo_occ  = mo_occ[index]
    coeff = coeff[:, index]

    np.save("lo_coeff.npy", coeff)
    np.save("lo_occ.npy", mo_occ)

    if average_occ:
        for fname in ["cc_mo_coeff.npy"]:
            if os.path.isfile(lde + "/" + fname):
                print("use: " + lde + "/" + fname)
                mo_coeff = np.load(lde + "/" + fname)
                break

        for fname in ["cc_dmmo.npy"]:
            if os.path.isfile(lde + "/" + fname):
                print("use: " + lde + "/" + fname)
                dmmo = np.load(lde + "/" + fname)
                break
        
        dmao = np.einsum('...pi,...ij,...qj->...pq', mo_coeff, dmmo, mo_coeff, optimize=True)
        if dmao.ndim == 3:
            dmao = np.einsum('spq->pq', dmao, optimize=True)

        coeff_inv = np.linalg.pinv(coeff)
        dmmo = np.einsum('ip,pq,jq->ij', coeff_inv, dmao, coeff_inv, optimize=True)

        print('AVERAGE TRACE = %8.5f' % np.trace(dmmo))

        nat_occ, u = np.linalg.eigh(dmmo)
        print('AVERAGE NAT OCC = ', ''.join(['%8.5f,' % x for x in nat_occ[::-1]]))



def psort(ova, fav, pT, coeff, orb_sym=None):
    pTnew = 2.0 * (coeff.T @ ova @ pT @ ova @ coeff)
    nocc  = np.diag(pTnew)
    index = np.argsort(-nocc)
    ncoeff = coeff[:, index]
    nocc = nocc[index]
    enorb = np.diag(coeff.T @ ova @ fav @ ova @ coeff)
    enorb = enorb[index]
    if orb_sym is not None:
        orb_sym = orb_sym[index]
    return ncoeff, nocc, enorb, orb_sym

if cas_list is None:
    assert nactorb is not None
    assert nactelec is not None
    ncore = (mol.nelectron - nactelec) // 2
    cas_list = list(range(ncore, ncore + nactorb))

print('cas list = ', cas_list)

orb_sym = None
if loc_with_pg:
    orb_sym = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, coeff, tol=1e-2)
    orb_sym = np.array([symm.irrep_name2id(mol.groupname, ir) for ir in orb_sym])

if split_low == 0.0 and split_high == 0.0:

    print('simple localization')

    actmo = coeff[:, np.array(cas_list, dtype=int)]
    if do_loc:
        if not loc_with_pg:
            ierr, ua = loc(mol, actmo)
            actmo = actmo.dot(ua)
        else:
            act_orb_sym = orb_sym[np.array(cas_list, dtype=int)]
            ierr, ua = loc_pg(mol, actmo, act_orb_sym)
            actmo = actmo.dot(ua)
    if not loc_with_pg:
        actmo, actocc, e_o, _ = psort(ova, fav, pav, actmo, None)
    else:
        act_orb_sym = orb_sym[np.array(cas_list, dtype=int)]
        actmo, actocc, e_o, act_orb_sym = psort(ova, fav, pav, actmo, act_orb_sym)

else:

    print('split localization at', split_low, '~', split_high)
    assert do_loc
    assert split_high >= split_low
    actmo = coeff[:, np.array(cas_list, dtype=int)]
    if loc_with_pg:
        act_orb_sym = orb_sym[np.array(cas_list, dtype=int)]
    actocc = mo_occ[np.array(cas_list, dtype=int)]
    print('active occ = ', np.sum(actocc, axis=-1), actocc)
    lidx = actocc <= split_low
    midx = (actocc > split_low) & (actocc <= split_high)
    hidx = actocc > split_high

    if len(actmo[:, lidx]) != 0:
        print('low orbs = ', np.array(list(range(len(lidx))))[lidx])
        if not loc_with_pg:
            ierr, ua = loc(mol, actmo[:, lidx])
            actmo[:, lidx] = actmo[:, lidx].dot(ua)
            actmo[:, lidx], actocc[lidx], _, _ = psort(ova, fav, pav, actmo[:, lidx], None)
        else:
            ierr, ua = loc_pg(mol, actmo[:, lidx], act_orb_sym[lidx])
            actmo[:, lidx] = actmo[:, lidx].dot(ua)
            actmo[:, lidx], actocc[lidx], _, act_orb_sym[lidx] = psort(ova, fav, pav, actmo[:, lidx], act_orb_sym[lidx])

    if len(actmo[:, midx]) != 0:
        print('mid orbs = ', np.array(list(range(len(midx))))[midx])
        if not loc_with_pg:
            ierr, ua = loc(mol, actmo[:, midx])
            actmo[:, midx] = actmo[:, midx].dot(ua)
            actmo[:, midx], actocc[midx], _, _ = psort(ova, fav, pav, actmo[:, midx])
        else:
            ierr, ua = loc_pg(mol, actmo[:, midx], act_orb_sym[midx])
            actmo[:, midx] = actmo[:, midx].dot(ua)
            actmo[:, midx], actocc[midx], _, act_orb_sym[midx] = psort(ova, fav, pav, actmo[:, midx], act_orb_sym[midx])


    if len(actmo[:, hidx]) != 0:
        print('high orbs = ', np.array(list(range(len(hidx))))[hidx])
        if not loc_with_pg:
            ierr, ua = loc(mol, actmo[:, hidx])
            actmo[:, hidx] = actmo[:, hidx].dot(ua)
            actmo[:, hidx], actocc[hidx], _, _ = psort(ova, fav, pav, actmo[:, hidx])
        else:
            ierr, ua = loc_pg(mol, actmo[:, hidx], act_orb_sym[hidx])
            actmo[:, hidx] = actmo[:, hidx].dot(ua)
            actmo[:, hidx], actocc[hidx], _, act_orb_sym[hidx] = psort(ova, fav, pav, actmo[:, hidx], act_orb_sym[hidx])

coeff[:, np.array(sorted(cas_list), dtype=int)] = actmo
mo_occ[np.array(sorted(cas_list), dtype=int)] = actocc
if loc_with_pg:
    orb_sym[np.array(sorted(cas_list), dtype=int)] = act_orb_sym

# sort_mo from pyscf.mcscf.addons

cas_list = np.array(sorted(cas_list), dtype=int)
mask = np.ones(coeff.shape[1], dtype=bool)
mask[cas_list] = False
idx = np.where(mask)[0]
nactorb = len(cas_list)
nactelec = int(np.round(sum(mo_occ[cas_list])) + 0.1)
assert (mol.nelectron - nactelec) % 2 == 0
ncore = (mol.nelectron - nactelec) // 2
print("NACTORB = %d NACTELEC = %d NCORE = %d" % (nactorb, nactelec, ncore))
coeff = np.hstack((coeff[:, idx[:ncore]], coeff[:, cas_list], coeff[:, idx[ncore:]]))
print('lo occ =', mo_occ[cas_list])
mo_occ = np.hstack((mo_occ[idx[:ncore]], mo_occ[cas_list], mo_occ[idx[ncore:]]))
if loc_with_pg:
    print('loc orb_sym =', orb_sym[cas_list])
    orb_sym = np.hstack((orb_sym[idx[:ncore]], orb_sym[cas_list], orb_sym[idx[ncore:]]))

np.save("lo_coeff.npy", coeff)
np.save("lo_occ.npy", mo_occ)
np.save("active_space.npy", (nactorb, nactelec))
if loc_with_pg:
    np.save("lo_orb_sym.npy", orb_sym)

txed = time.perf_counter()
print("FINISH TIME = ", datetime.now().strftime("%m/%d/%Y %H:%M:%S"))
print("TOTAL TIME  = %20.3f" % (txed - txst))
