
import time
from datetime import datetime
txst = time.perf_counter()
print("START  TIME = ", datetime.now().strftime("%m/%d/%Y %H:%M:%S"))

from pyscf import gto, scf
mol = gto.Mole()
mol.verbose = 4
mol.atom = '''
C    -1.7030912   -0.2687576   -3.8857510 
H    -1.4276567    0.6161659   -3.2739484 
H    -2.7282115   -0.5817080   -3.5791673 
H    -1.7544034    0.0472283   -4.9585791 
S    -0.4654014   -1.5971491   -3.6068456 
C     1.3427752    1.2911521   -3.6530591 
H     2.4027725    1.5361646   -3.9146487 
H     0.7051324    2.0932763   -4.0999022 
H     1.0611069    0.3183708   -4.1168671 
S     1.1257526    1.1531598   -1.8302475 
C    -0.5460690    2.6148722    1.2254132 
H     0.3817785    2.2964072    1.7410684 
H    -0.7581861    3.6865849    1.4682270 
H    -0.3516724    2.5166490    0.1383690 
S    -1.9524226    1.5437145    1.7437219 
C    -1.2593992   -1.2001762    4.2059136 
H    -0.1794134   -1.4477025    4.1290192 
H    -1.6496014   -1.5773749    5.1869269 
H    -1.3493637   -0.0933792    4.1760396 
S    -2.1872294   -1.9182902    2.7849438 
Fe    0.2293122   -1.0722584   -1.1389980 
Fe   -1.0895429   -0.6236298    1.1100321 
S     1.1201274   -0.6801252    1.0718109 
C    -1.6833719   -1.3940353   -0.5791815 
H     1.1279016   -2.5680078   -1.3656985 
H    -2.5031122   -0.8682227   -1.1171969 
H    -1.7445096   -2.4939284   -0.7433938 

'''
mol.basis = "ccpvqz-dk"
mol.spin = 0
mol.charge = -3
mol.max_memory = 82000
from pyscf import dftd3

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
mf.xc = 'pbe'
mf.max_cycle = 1000
mf = dftd3.dftd3(mf)
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
