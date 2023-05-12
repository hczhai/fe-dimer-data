
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
lde = '../select-19'
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

from pyscf import dmrgscf, lib
import os

dmrgscf.settings.BLOCKEXE = os.popen("which block2main").read().strip()
dmrgscf.settings.MPIPREFIX = "" if "PYSCF_MPIPREFIX" not in os.environ else os.environ["PYSCF_MPIPREFIX"]

mc.fcisolver = dmrgscf.DMRGCI(mol, maxM=5000, tol=1E-10)
mc.fcisolver.runtimeDir = lib.param.TMPDIR
mc.fcisolver.scratchDirectory = lib.param.TMPDIR
mc.fcisolver.threads = int(os.environ["OMP_NUM_THREADS"])
mc.fcisolver.memory = int(mol.max_memory / 1000)

mcfs = [mc.fcisolver]
for mcf in mcfs:
    mcf.scheduleSweeps = [0, 4, 8, 12, 16, 20, 24, 28, 32, 36]
    mcf.scheduleMaxMs = [400, 800, 1200, 1600, 2000, 2500, 3000, 4000, 5000, 5000]
    mcf.scheduleTols = [0.0001, 0.0001, 1e-05, 1e-05, 1e-06, 1e-06, 1e-06, 1e-06, 1e-06, 1e-07]
    mcf.scheduleNoises = [0.0001, 0.0001, 0.0001, 5e-05, 5e-05, 5e-05, 5e-05, 5e-05, 5e-05, 0.0]
    mcf.maxIter = 48
    mcf.twodot_to_onedot = 40
    mcf.tol = 1E-8
    mcf.twopdm = False
    mcf.block_extra_keyword = ['onepdm']

for mcf in mcfs:
    mcf.conv_tol = 1E-10
from pyscf import dmrgscf
dmrgscf.dryrun(mc)
for mcf in mcfs:
    mcf.scheduleSweeps = [0, 4, 8, 12, 16, 20, 24, 28]
    mcf.scheduleMaxMs = [5000, 4500, 4000, 3500, 3000, 2500, 2000, 1500]
    mcf.scheduleTols = [1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07]
    mcf.scheduleNoises = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    mcf.maxIter = 32
    mcf.twodot_to_onedot = 0
    mcf.tol = 0.0
    mcf.twopdm = False
    mcf.block_extra_keyword = ['fullrestart', 'twodot', 'extrapolation']
    mcf.configFile = "dmrg-rev.conf"

for mcf in mcfs:
    mcf.conv_tol = 1E-10
from pyscf import dmrgscf
dmrgscf.dryrun(mc)
for mcf in mcfs:
    mcf.scheduleSweeps = [0]
    mcf.scheduleMaxMs = [5000]
    mcf.scheduleTols = [5E-6]
    mcf.scheduleNoises = [0.0]
    mcf.maxIter = 1
    mcf.twodot_to_onedot = 0
    mcf.tol = 0.0
    mcf.twopdm = False
    mcf.block_extra_keyword = ['trans_mps_to_singlet_embedding']
    mcf.block_extra_keyword += ['restart_copy_mps SEKET']
    mcf.block_extra_keyword += ['restart_sample 0.1']
    mcf.configFile = "dmrg-csf.conf"

for mcf in mcfs:
    mcf.conv_tol = 1E-10
from pyscf import dmrgscf
dmrgscf.dryrun(mc)

txed = time.perf_counter()
print("FINISH TIME = ", datetime.now().strftime("%m/%d/%Y %H:%M:%S"))
print("TOTAL TIME  = %20.3f" % (txed - txst))
