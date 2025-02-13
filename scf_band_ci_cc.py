#!/usr/bin/env python

'''
SCF module currently does not apply SO-ECP automatically. SO-ECP contributions
can be added to GHF/GKS core Hamiltonian by overwriding the method get_hcore.
Since pyscf-2.0 setting attribte with_soc in GHF object can include the
ECP-SOC integrals in core Hamiltonian.

See also examples/gto/20-soc_ecp.py
'''

import pyscf
from pyscf.pbc import gto, cc

cell = gto.Cell()
cell.build(
    a = [[0.00000000000, 1.77549387821, 1.77549387821],
         [1.77549387821, 0.00000000000, 1.77549387821],
         [1.77549387821, 1.77549387821, 0.00000000000]],
    atom = 'C 0 0 0',  # in Angstrom
#    basis = 'gth-cc-qzvp', #'631g',
    basis = '321g', #'631g',
    ke_cutoff = 550,
    mesh = [10]*3,
    space_group_symmetry = True,
    symmetry = True,
    verbose = 5
)


kp = 2
kpts = cell.make_kpts([kp,kp,kp],wrap_around=True, with_gamma_point=True, scaled_center=None, space_group_symmetry=True, time_reversal_symmetry=True)
print('kp = ', kp) 

##################
# SCF
##################

cscf_set = pyscf.pbc.scf.KRHF(cell,kpts=kpts)
# https://pyscf.org/user/pbc/scf.html#smearing
cscf_set = pyscf.pbc.scf.addons.smearing_(cscf_set, sigma=0.0001, method='fermi')

cscf = cscf_set.run(
    init_guess = 'minao', # chk to read from chkfile
    chkfile    = 'scf_KRHF.chk',
    max_memory = 250000,
    verbose    = 4,
    diis_start_cycle = 1,
    diis_space = 8,
    diis_damp = 0.5,
    conv_tol = 1e-6,
    conv_tol_grad = None,
    max_cycle = 150,
    direct_scf = True
)

# Code to print eigenvalues and occupancy for each k-point in blocks with k-point vector components
mo_energy_kpts = cscf.mo_energy
mo_occ_kpts = cscf.mo_occ

for k, (mo_energy, mo_occ) in enumerate(zip(mo_energy_kpts, mo_occ_kpts)):
    kpt_vector = kpts.kpts[k]  # Correct way to access k-point vector
    print(f"k-point {k+1} ({kpt_vector[0]:.4f}, {kpt_vector[1]:.4f}, {kpt_vector[2]:.4f}):")
    print(f"{'Eigenvalues':<20} {'Occupancy':<20}")
    for energy, occ in zip(mo_energy, mo_occ):
        print(f"{energy:<20} {occ:<20}")
    print()
    


#cisolver = pyscf.pbc.ci.KCIS(cscf)
#cisolver.run(
#    threads = 30,
#    verbose = 4,
#    with_soc = True,
#)

# Disable smearing and recalculate occupation numbers
cscf.smearing_method = False
cscf.mo_occ = cscf.get_occ()

ccsolver = pyscf.pbc.cc.KRCCSD(cscf)
ccsolver.kpts = kpts.kpts  # Ensure kpts are correctly set
ccsolver.run(
    diis_start_cycle = 3,
    diis_space = 8,
    verbose = 4,
    with_soc = True,
    max_cycle = 2
)

##################################
#  Band structure calculation
##################################

import pyscf.pbc.tools.pyscf_ase as pyscf_ase
from ase.build import bulk
from ase.dft.kpoints import sc_special_points as special_points, get_bandpath

#points = special_points['fcc']
#G = points['G']
#X = points['X']
#W = points['W']
#K = points['K']
#L = points['L']
#band_kpts, kpath, sp_points = get_bandpath([L, G, X, W, K, G], cell.a, npoints=50)

G=[0.0000,   0.0000,   0.0000] 
X=[0.5000,   0.0000,   0.5000]
W=[0.5000,   0.2500,   0.7500]
K=[0.3750,   0.3750,   0.7500]
G=[0.0000,   0.0000,   0.0000]
L=[0.5000,   0.5000,   0.5000]
U=[0.6250,   0.2500,   0.6250]
W=[0.5000,   0.2500,   0.7500]
L=[0.5000,   0.5000,   0.5000]
K=[0.3750,   0.3750,   0.7500]
U=[0.6250,   0.2500,   0.6250]
X=[0.5000,   0.0000,   0.5000]


band_kpts, kpath, sp_points = get_bandpath([G, X, W, K, G, L, U, W, L, K, U, X], cell.a, npoints=50)

band_kpts = cell.get_abs_kpts(band_kpts)

#e_kn = cscf.get_bands(band_kpts)[0]
# Replace line 107 with the following code
e_kn = []
for kpt in band_kpts:
    # Update k-points in ccsolver
    ccsolver.kpts = [kpt]
    # Run the coupled-cluster solver for the specific k-point
    ccsolver.run()
    # Get the band energies for this k-point
    e_kn.append(ccsolver.mo_energy[0])

# Continue with the rest of the code
Ef = cscf.get_fermi()
print("Fermi level: ", Ef)
e_kn_shift = [en - Ef for en in e_kn]

au2ev = 27.21138624598

# Band structure plot
nbands = cell.nao_nr()
# Write bandstructure to file band.dat

outf = "band.dat"
out = open(outf,"w")

print("Writing band structure to file ", outf)

print("# Fermi level is set to 0; energy levels in eV", file=out)
print("# G, X, W, K, G, L, U, W, L, K, U, X", file=out)
print("#", *sp_points, file=out)

for i in range(nbands):
    for j in range(len(kpath)):
        print(kpath[j], e_kn_shift[j][i] * au2ev, file=out)
    print(file=out)

out.close()

Ef = cscf.get_fermi()
print("Fermi level: ", Ef)
e_kn_shift = [en - Ef for en in e_kn]

au2ev = 27.21138624598

#########################
# Band structure plot
#########################
nbands = cell.nao_nr()
# Write bandstructure to file band.dat

outf = "band.dat"
out = open(outf,"w")

print(" Writing band structure to file ", outf)

print("# Fermi level is set to 0; energy levels in eV", file=out)
print("# G, X, W, K, G, L, U, W, L, K, U, X", file=out)
print("#", *sp_points, file=out)

for i in range(nbands):
    for j in range(len(kpath)):
        print(kpath[j], e_kn_shift[j][i]*au2ev, file=out)
    print(file=out)

out.close()

