
import time
from datetime import datetime
txst = time.perf_counter()
print("START  TIME = ", datetime.now().strftime("%m/%d/%Y %H:%M:%S"))

from pyscf import scf
import numpy as np
mfchk = "../mf-1/mf.chk"
mol, mfx = scf.chkfile.load_scf(mfchk)
x2c = True
d3 = False
nactorb = None
nactelec = None
semi_canonical = False
cascc_mf_relax = False
lde = '../select-17'
spin = 0

from pyscf import mcscf
import os

mcscf.casci.FRAC_OCC_THRESHOLD = 1E-6

for fname in ["mo_coeff.npy", "lo_coeff.npy", "nat_coeff.npy"]:
    if os.path.isfile(lde + "/" + fname):
        print("use: " + lde + "/" + fname)
        coeff = np.load(lde + "/" + fname)
        break

if nactelec is None or nactelec is None:
    print("use: " + lde + "/active_space.npy")
    nactorb, nactelec = np.load(lde + "/active_space.npy")

print("act: orb = %d elec = %d spin = %d" % (nactorb, nactelec, spin))

nacta = (nactelec + spin) // 2
nactb = (nactelec - spin) // 2

if coeff.ndim == 3:
    print('use UHF')
    mf = scf.UHF(mol)
else:
    print('use RHF')
    mf = scf.RHF(mol)
if x2c:
    mf = scf.sfx2c(mf)
if d3:
    from pyscf import dftd3
    mf = dftd3.dftd3(mf)

mf.chkfile = "mf.chk"
mf.mo_coeff = coeff

mc = mcscf.CASCI(mf, nactorb, (nacta, nactb))
mc.conv_tol = 1E-8
mc.max_cycle_macro = 50
mc.canonicalization = False
mc.natorb = False

mcfs = [mc.fcisolver]

import numpy as np
from pyscf import cc, gto
from libdmet.basis_transform import make_basis

scf_dmao = np.load("../mf-1/mf_dmao.npy")
scf_dmlo = make_basis.transform_rdm1_to_lo_mol(scf_dmao, coeff, mf.get_ovlp())
dmcas = scf_dmlo[:, mc.ncore:mc.ncore + mc.ncas, mc.ncore:mc.ncore + mc.ncas]

print('idemponency of dmcas[0]: %s' % np.linalg.norm(dmcas[0].dot(dmcas[0]) - dmcas[0]))
print('idemponency of dmcas[1]: %s' % np.linalg.norm(dmcas[1].dot(dmcas[1]) - dmcas[1]))
print('trace of dmcas[0]: %s' % np.trace(dmcas[0]))
print('trace of dmcas[1]: %s' % np.trace(dmcas[1]))

class UCCSolver:
    def __init__(self, ccsd_t=False, dmcas=None):
        self.ccsd_t = ccsd_t
        self.dmcas = dmcas

    def kernel(self, h1, h2, norb, nelec, ci0=None, ecore=0, **kwargs):
        mol = gto.M(verbose=4)
        mol.nelectron = sum(nelec)
        mol.spin = nelec[0] - nelec[1]
        mf = mol.UHF()
        mf._eri = h2
        mf.get_hcore = lambda *args: h1
        mf.get_ovlp = lambda *args: np.identity(norb)
        mf.max_cycle = 200
        mf.kernel(dm0=self.dmcas)

        if semi_canonical:
            semi_canon(mf)

        self.cc = cc.UCCSD(mf)
        self.cc.level_shift = 0.0
        self.cc.max_cycle = 200
        self.cc.incore_complete = True
        self.cc.run()
        if self.ccsd_t:
            e_ccsd_t = self.cc.e_tot + self.cc.ccsd_t()
        else:
            e_ccsd_t = self.cc.e_tot
        return e_ccsd_t + ecore, dict(t1=self.cc.t1, t2=self.cc.t2)

    def make_rdm1(self, t12, norb, nelec):
        dms = self.cc.make_rdm1(**t12)
        if isinstance(dms, tuple):
            return dms[0] + dms[1]
        else:
            return dms

mc.fcisolver = UCCSolver(ccsd_t=True, dmcas=dmcas)
mcfs = [mc.fcisolver]

for mcf in mcfs:
    mcf.conv_tol = 1E-10
mc.kernel()
for _ in range(2):
    if not mc.converged:
        mc.kernel()
dmao = mc.make_rdm1()

import numpy as np
coeff_inv = np.linalg.pinv(mc.mo_coeff)
dmmo = np.einsum('...ip,...pq,...jq->...ij', coeff_inv, dmao, coeff_inv, optimize=True)
if dmmo[0].ndim == 2:
    mc_occ = np.diag(dmmo[0]) + np.diag(dmmo[1])
else:
    mc_occ = np.diag(dmmo)

np.save("mc_occ.npy", mc_occ)
np.save("mc_mo_coeff.npy", mc.mo_coeff)
np.save("mc_e_tot.npy", mc.e_tot)
np.save("mc_dmao.npy", dmao)
np.save("mc_dmmo.npy", dmmo)

nat_occ, u = np.linalg.eigh(dmmo)
nat_coeff = np.einsum('...pi,...ij->...pj', mc.mo_coeff, u, optimize=True)
np.save("nat_coeff.npy", nat_coeff[..., ::-1])
np.save("nat_occ.npy", nat_occ[..., ::-1])
np.save("active_space.npy", (nactorb, nactelec))
np.save("spin.npy", (spin, ))

txed = time.perf_counter()
print("FINISH TIME = ", datetime.now().strftime("%m/%d/%Y %H:%M:%S"))
print("TOTAL TIME  = %20.3f" % (txed - txst))
