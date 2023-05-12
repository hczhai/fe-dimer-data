
import time
from datetime import datetime
txst = time.perf_counter()
print("START  TIME = ", datetime.now().strftime("%m/%d/%Y %H:%M:%S"))

from pyscf import gto, scf
mol = gto.Mole()
mol.verbose = 4
mol.atom = '''
C    -1.0654597   -1.6165964   -3.6938910 
H    -1.4383635   -0.5689779   -3.6741341 
H    -1.7807115   -2.2443224   -3.1121367 
H    -1.0773907   -1.9924332   -4.7494055 
S     0.6243829   -1.6537410   -2.9602556 
C     0.4156598    1.7876147   -4.0556776 
H     1.3488982    2.3835358   -3.9218054 
H    -0.1330190    2.2251636   -4.9293434 
H     0.7135941    0.7422418   -4.2919146 
S    -0.6070437    1.7896179   -2.5222937 
C    -0.8232495    2.8573739    1.4681234 
H    -0.0877151    2.8066116    0.6407689 
H    -0.9439809    3.9246222    1.7897946 
H    -1.7975667    2.5226444    1.0498580 
S    -0.2735647    1.7807874    2.8579066 
C    -1.0643237   -1.5242334    4.0337374 
H    -1.0304628   -1.8368555    5.1103107 
H    -1.8734238   -2.1231958    3.5494502 
H    -1.3541746   -0.4526895    3.9781077 
S     0.5588840   -1.7538671    3.1864833 
Fe    0.4943827    0.0464709   -1.2052434 
Fe    0.5810012   -0.0475200    1.4684925 
S     2.2702018    0.5987571    0.0818694 
C    -0.9804810   -0.8504359    0.1595883 
H    -0.9278775   -1.9603768    0.2106940 
H    -1.5301625   -0.5677658   -0.7796813 
H    -1.6330329   -0.4914309    0.9945982 

'''
mol.basis = "ccpvtz-dk"
mol.spin = 0
mol.charge = -3
mol.max_memory = 82000

mol.build()
print("NAO   = ", mol.nao)
print("NELEC = ", mol.nelec)
dm = None

from pyblock2._pyscf import scf as b2scf
dm = b2scf.get_metal_init_guess(mol, orb="Fe 3d", atom_idxs=[20, 21], coupling="+-", atomic_spin=4)


print("PG = ", mol.groupname)

mf = scf.sfx2c(scf.UKS(mol))
mf.chkfile = 'mf.chk'
mf.conv_tol = 1E-12
mf.xc = 'b3lyp'
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

from pyblock2._pyscf import scf as b2scf
b2scf.mulliken_pop_dmao(mol, dm)


txed = time.perf_counter()
print("FINISH TIME = ", datetime.now().strftime("%m/%d/%Y %H:%M:%S"))
print("TOTAL TIME  = %20.3f" % (txed - txst))
