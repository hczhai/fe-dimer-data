
import time
from datetime import datetime
txst = time.perf_counter()
print("START  TIME = ", datetime.now().strftime("%m/%d/%Y %H:%M:%S"))

from pyscf import gto, scf
mol = gto.Mole()
mol.verbose = 4
mol.atom = '''
C    -0.6004690   -1.8541396   -4.3931123 
H    -0.9885128   -0.8116508   -4.4441158 
H    -1.2233186   -2.5109691   -5.0545639 
H     0.4332653   -1.8614472   -4.8208224 
S    -0.6420813   -2.4169491   -2.6407698 
C     0.4911160    1.7119732   -3.4952934 
H     1.1454118    1.7908560   -2.6017265 
H     0.4940557    2.6834804   -4.0526033 
H     0.9297739    0.9375322   -4.1667404 
S    -1.2045031    1.2624527   -2.9490191 
C    -0.5170233    2.7735147    1.7890395 
H     0.1897341    2.5579946    0.9587402 
H    -0.1817476    3.6864329    2.3434608 
H    -1.5092436    2.9968065    1.3360800 
S    -0.5859628    1.3069271    2.8986914 
C    -0.5055631   -1.6440394    4.3765240 
H    -1.3125977   -0.8903412    4.5041769 
H     0.4429013   -1.1359999    4.6628636 
H    -0.6817372   -2.4716890    5.1124802 
S    -0.4582362   -2.2773941    2.6480953 
Fe   -0.5822491   -0.4261876   -1.2438221 
Fe   -0.4040092   -0.3452652    1.1749566 
S     1.3672017    0.2724370   -0.2315191 
C    -2.0078868    0.0019752    0.0519971 
H    -0.5166457   -1.6358452   -0.0228939 
H    -2.3493298    1.0582445   -0.0027923 
H    -2.8433420   -0.7077099    0.2366894 

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
