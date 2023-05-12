
import time
from datetime import datetime
txst = time.perf_counter()
print("START  TIME = ", datetime.now().strftime("%m/%d/%Y %H:%M:%S"))
spin = 0

from pyscf import scf, lib
import numpy as np
mfchk = "../mf-1/mf.chk"
mol, mfx = scf.chkfile.load_scf(mfchk)
if spin is not None:
    mol.spin = spin
    mol.build()
mf = scf.sfx2c(scf.UHF(mol))
mf.chkfile = "mf.chk"
mf.mo_coeff = mfx["mo_coeff"]
mf.mo_energy = mfx["mo_energy"]
mf.mo_occ = mfx["mo_occ"]

non_canonical = False
do_spin_square = False
mf.max_memory = 120000
non_canonical = True

from sys import argv
import os
is_restart = len(argv) >= 2 and argv[1] == "1"

if not is_restart:
    for fname in ['/mpdiis.h5']:
        if os.path.isfile(lib.param.TMPDIR + fname):
            fid = 1
            while os.path.isfile(lib.param.TMPDIR + fname + '.%d' % fid):
                fid += 1
            os.rename(lib.param.TMPDIR + fname,
                lib.param.TMPDIR + fname + '.%d' % fid)

print('mf occ', np.sum(mf.mo_occ, axis=-1), mf.mo_occ)

from pyscf import mp
mc = mp.MP2(mf)
if non_canonical:
    mf.converged = False
mc.diis_file = lib.param.TMPDIR + '/mpdiis.h5'
mc.max_cycle = 1000

if is_restart and os.path.isfile(lib.param.TMPDIR + '/mpdiis.h5'):
    print("restart mp from ", lib.param.TMPDIR + '/mpdiis.h5')
    mc.restore_from_diis_(lib.param.TMPDIR + '/mpdiis.h5')
    t1, t2 = mc.t1, mc.t2
    mc.kernel(t1, t2)
else:
    mc.kernel()
e_mp2 = mc.e_tot
print('EMP2     = ', e_mp2)
print("PART TIME (MP2) = %20.3f" % (time.perf_counter() - txst))

if do_spin_square:
    S2 = mc.spin_square()[0]
    print('MP2 <S^2> = ', S2)
    print("PART TIME (MP2 S2) = %20.3f" % (time.perf_counter() - txst))

import numpy as np
dm = mc.make_rdm1()
if dm[0].ndim == 2:
    mc_occ = np.diag(dm[0]) + np.diag(dm[1])
else:
    mc_occ = np.diag(dm)
print("PART TIME (1PDM)  = %20.3f" % (time.perf_counter() - txst))


import numpy as np
np.save("cc_occ.npy", mc_occ)
np.save("cc_mo_coeff.npy", mc.mo_coeff)
np.save("cc_e_tot.npy", mc.e_tot)
np.save("cc_dmmo.npy", dm)

# dmao = np.einsum('xpi,xij,xqj->xpq', mc.mo_coeff, dm, mc.mo_coeff)
# coeff_inv = np.linalg.pinv(mc.mo_coeff)
# dmmo = np.einsum('xip,xpq,xjq->xij', coeff_inv, dmao, coeff_inv)

nat_occ, u = np.linalg.eigh(dm)
nat_coeff = np.einsum('...pi,...ij->...pj', mc.mo_coeff, u, optimize=True)
np.save("nat_coeff.npy", nat_coeff[..., ::-1])
np.save("nat_occ.npy", nat_occ[..., ::-1])

print('nat occ', np.sum(nat_occ, axis=-1), nat_occ)

txed = time.perf_counter()
print("FINISH TIME = ", datetime.now().strftime("%m/%d/%Y %H:%M:%S"))
print("TOTAL TIME  = %20.3f" % (txed - txst))
