
import time
from datetime import datetime
txst = time.perf_counter()
print("START  TIME = ", datetime.now().strftime("%m/%d/%Y %H:%M:%S"))

from pyscf import gto, scf
mol = gto.Mole()
mol.verbose = 4
mol.atom = '''
C    -1.4481492   -1.2918793   -4.0685921 
H    -1.5592125   -0.1859823   -4.0893224 
H    -2.3309980   -1.7079134   -3.5297987 
H    -1.4760151   -1.6885040   -5.1175674 
S     0.1243831   -1.7171690   -3.2045707 
C     1.0341701    1.7983129   -3.6269368 
H     1.9627435    1.9656559   -3.0387507 
H     0.9926332    2.5651117   -4.4445661 
H     1.1042400    0.7824705   -4.0723829 
S    -0.4446417    1.8649325   -2.5294114 
C    -1.1427404    2.6271747    1.3777265 
H    -0.7051544    2.6282955    0.3581567 
H    -1.2853067    3.6896251    1.7066593 
H    -2.1508112    2.1651199    1.2944240 
S    -0.0778785    1.6928495    2.5587062 
C    -0.7929520   -1.4324890    4.2467679 
H    -0.9625390   -0.3357608    4.1885620 
H    -0.5400427   -1.7119510    5.3045419 
H    -1.7588854   -1.9373028    4.0010839 
S     0.5386887   -1.8989304    3.0569739 
Fe   -0.0624496   -0.2245743   -1.3014183 
Fe    0.1805793   -0.3587873    1.2612831 
S     2.1356417    0.1952062   -0.2457191 
C    -1.2915962   -0.9617811    0.0689167 
H     1.9086364    1.5088483    0.0713138 
H    -2.2983991   -0.4932717    0.2239739 
H    -1.4219432   -2.0743068    0.0479475

'''
mol.basis = "ccpvqz-dk"
mol.spin = 0
mol.charge = -3
mol.max_memory = 82000

mol.build()
print("NAO   = ", mol.nao)
print("NELEC = ", mol.nelec)
dm = None

spin_bak = mol.spin
mol.spin = 0
mol.build()

import numpy as np

rmf = scf.RHF(mol)
dm0 = rmf.get_init_guess(key="atom")
dm0 = np.array([dm0, dm0]) * 0.5
idx0 = mol.search_ao_label("0 Fe 3d.*")
idx1 = mol.search_ao_label("1 Fe 3d.*")

from pyscf import lo
from libdmet.basis_transform import make_basis

ld = lo.orth_ao(mol, 'lowdin', pre_orth_ao='SCF')
dl0 = make_basis.transform_rdm1_to_lo_mol(dm0, ld, rmf.get_ovlp())

dspin = float(5) / 2
afm = True

if afm:
    dl0[0][idx0, idx0] += dspin / len(idx0)
    dl0[0][idx1, idx1] -= dspin / len(idx1)
    dl0[1][idx0, idx0] -= dspin / len(idx0)
    dl0[1][idx1, idx1] += dspin / len(idx1)
else:
    dl0[0][idx0, idx0] += dspin / len(idx0)
    dl0[0][idx1, idx1] += dspin / len(idx1)
    dl0[1][idx0, idx0] -= dspin / len(idx0)
    dl0[1][idx1, idx1] -= dspin / len(idx1)

dm = make_basis.transform_rdm1_to_ao_mol(dl0, ld)

mol.spin = spin_bak
mol.build()

mf = scf.sfx2c(scf.UHF(mol))
mf.chkfile = 'mf.chk'
mf.conv_tol = 1E-12
mf.max_cycle = 1000
mf = mf.newton()

mf.kernel(dm0=dm)
dm = mf.make_rdm1()

import numpy as np
np.save("mf_occ.npy", mf.mo_occ)
np.save("mo_coeff.npy", mf.mo_coeff)
np.save("mo_energy.npy", mf.mo_energy)
np.save("e_tot.npy", mf.e_tot)
np.save("mf_dmao.npy", dm)

txed = time.perf_counter()
print("FINISH TIME = ", datetime.now().strftime("%m/%d/%Y %H:%M:%S"))
print("TOTAL TIME  = %20.3f" % (txed - txst))
