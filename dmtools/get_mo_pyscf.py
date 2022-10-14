#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 10:23:10 2021

@author: nico
"""
import numpy as np
import sys
#import pyscf
from pyscf import scf, gto


def get_pyscf_coeffs(xyz, basname):
    basfile = basname + ".nwchem"
    with open(basfile, "r") as f:
        bas = f.read()
    mol = gto.M(atom=xyz, basis=bas)
    scfres = scf.RHF(mol)
    scfres.conv_tol = 1e-12
    scfres.kernel()
    coeffs = scfres.mo_coeff
    np.savetxt("{}_coeffs_pyscf.txt".format(basname), coeffs)
    return coeffs

if __name__ == "__main__":
    _, xyz, basname = sys.argv
    get_pyscf_coeffs(xyz, basname)