from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from ase.io import read, write
from ase import *

from pymatgen.core.structure import Composition, Lattice, Structure

import os, argparse, pymatgen, yaml, sys, json
import numpy as np
from pathlib import Path

def write_inputs(atoms, user_incar_settings={}, mags=None, Ucorr=False, 
                 dipolcorr=False, kmult=30):
    """
    Generate VASP inputs for running calculations of RuO2 doped slab w/ 
        and w/o adsorbates 
    atoms (ASE Atoms): Atoms representation of a RuO2 slab, doped slab, 
        or adsorbed slab
    user_incar_settings (dict): User defined VASP inputs, will override 
        any default settings
    mags (list): Magnetic moment of each atom
    Ucorr (bool): Whether or not to apply Materials Project Hubbard U 
        corrections to atoms
    dipolcorr (bool): Whether or not to apply dipole corrections, useful 
        for converging slabs with assymetric surfaces
    kmult (int): Kpoints are set to int(kmult/|a|) x int(kmult/|b|) x 1
    """
    

    # default settings for oxide
    vasp_params = dict(xc='PBE', gga='PE', lreal=False, encut=500, ediff=1e-4,
                       ediffg=-0.03, ispin=2, isif=2, lcharg=False,
                       lwave=False, ismear=0, sigma=0.2, isym=2, nsw=300, lorbit=11,
                       lvtot=False, ibrion=2, potim=0.5, nelm=300, algo='Fast')

    if Ucorr:
        # add Hubbard U corrections
        hubbard_u_dict = {'Co': 3.32, 'Cr': 3.7, 'Fe': 5.3, 'Mn': 3.9,
                          'Mo': 4.38, 'Ni': 6.2, 'V': 3.25, 'W': 6.2}
        if any([site.symbol in hubbard_u_dict.keys() for site in atoms]):
            vasp_params['ldau'] = True
            ldau_luj = {site.symbol: {} for site in atoms}
            for el in ldau_luj.keys():
                ldau_luj[el]['U'] = hubbard_u_dict[el] if el \
                in hubbard_u_dict.keys() else 0
                ldau_luj[el]['J'] = 0
                ldau_luj[el]['L'] = 2 if el in hubbard_u_dict.keys() else 0
            vasp_params['ldau_luj'] = ldau_luj
            vasp_params['ldauprint'] = 0
            vasp_params['ldautype'] = 2
        

            # contains f-electrons
            if any(z > 56 for z in atoms.get_atomic_numbers()):
                vasp_params["lmaxmix"] = 6

            # contains d-electrons
            elif any(z > 20 for z in atoms.get_atomic_numbers()):
                vasp_params["lmaxmix"] = 4

    if dipolcorr:
        # Add dipole correction if adslab
        if 'adslab-' in os.getcwd():
            vasp_params['ldipol'] = True
            vasp_params['idipol'] = 3
            weights = [Composition(site.symbol).weight for site in atoms]
            # center of mass for the slab
            vasp_params['dipol'] = np.average(atoms.get_scaled_positions(),
                                              weights=weights, axis=0)

    # Calculate appropriate kpoints
    vasp_params['kpts'] = (np.ceil(kmult/atoms.cell.lengths()[0]),
                           np.ceil(kmult/atoms.cell.lengths()[1]), 1)

    # set magmoms
    if not mags:
        mags = {'Ce': 5, 'Ce3+': 1, 'Co': 0.6, 'Co3+': 0.6, 'Co4+': 1, 'Cr': 5, 
                'Dy3+': 5, 'Er3+': 3, 'Eu': 10, 'Eu2+': 7, 'Eu3+': 6, 'Fe': 5, 
                'Gd3+': 7, 'Ho3+': 4, 'La3+': 0.6, 'Lu3+': 0.6, 'Mn': 5, 
                'Mn3+': 4, 'Mn4+': 3, 'Mo': 5, 'Nd3+': 3, 'Ni': 5, 'Pm3+': 4, 
                'Pr3+': 2, 'Sm3+': 5, 'Tb3+': 6, 'Tm3+': 2, 'V': 5, 'W': 5, 'Yb3+': 1}   
        vasp_params['magmom'] = [mags[site.symbol] if site.symbol in \
                                 mags.keys() else 0.6 for site in atoms]
    else:
        vasp_params['magmom'] = mags
        
    vasp_params.update(user_incar_settings)

    # add selective dynamics
    atoms.set_constraint(FixAtoms([i for i, t in enumerate(atoms.get_tags()) if t == 0]))
    calc = Vasp(**vasp_params)
    calc.write_input(atoms)

