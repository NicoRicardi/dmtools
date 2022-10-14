#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 15:25:34 2021

@author: nico
"""

# imports
import json as js
import sys
import shutil as sh
import os
from pyscf import gto
import numpy as np
from pyscf.gto import ANG_OF, NCTR_OF
import mendeleev as md

# functions
def mkdif(path):
    """
    Parameters
    ----------
    path: str
        the path to create if not present
    """
    if not os.path.isdir(path):
        os.mkdir(path)
# to change in case
repeat_elements = True  # whether more atoms of the same elements should be analysed. must be True for mixed basis!
rows_equal=True  # Whether atoms in the same period have the same order. must be False if mixed basis!
# getting started
root = os.getcwd()
args = sys.argv[1:]
if len(args) != 2:
    raise ValueError("Run this script as \"get_correspondence.py file.xyz basis_name\"")
xyz, basname = args
nwchem = basname + ".nwchem"
basfile = basname + ".bas"
mkdif(basname)
if not os.path.isfile(os.path.join(root, basname, xyz)):
    if os.path.isfile(os.path.join(root, xyz)):
        sh.copy(os.path.join(root, xyz), os.path.join(root, basname, xyz))
    else:
        raise FileNotFoundError("Cannot find your xyz file")
if not os.path.isfile(os.path.join(root, basname, nwchem)):
    if os.path.isfile(os.path.join(root, nwchem)):
        sh.move(os.path.join(root, nwchem), os.path.join(root, basname, nwchem))
    else:
        raise FileNotFoundError("Cannot find an nwchem file for your basis. Check your spelling.\n\
                                You can get nwchem basis files from https://www.basissetexchange.org/")       
if not os.path.isfile(os.path.join(root, basname, basfile)):
    if os.path.isfile(os.path.join(root, basfile)):
        sh.move(os.path.join(root, basfile), os.path.join(root, basname, basfile))
    else:
        raise FileNotFoundError("Cannot find a .bas file for your basis. Check your spelling.\n\
                                You can get nwchem basis files from https://www.basissetexchange.org/\
                                asking for format \"gaussian\", or in your $QCAUX")          
# getting coeffs 
os.chdir(basname)
qcoeff_file = "{}_coeffs_qchem.txt".format(basname)
pcoeff_file = "{}_coeffs_pyscf.txt".format(basname)
if not os.path.isfile(pcoeff_file):
    from get_mo_pyscf import get_pyscf_coeffs
    get_pyscf_coeffs(xyz, basname)
if not os.path.isfile(qcoeff_file):
    print("getting qchem coefficients")
    from get_mo_qchem import get_qchem_coeffs
    get_qchem_coeffs(xyz, basname)   
qcoeffs = abs(np.loadtxt(qcoeff_file)) # abs because phase can change
pcoeffs = abs(np.loadtxt(pcoeff_file).T)  # transpose so that it is [MO, basis_func]

# let's go!
with open(nwchem, "r") as f:
    bas = f.read()
mol = gto.M(atom=xyz, basis=bas)

atm_numbers = [i[0] for  i in mol._atm] 
set_atm_numbers = list(set(atm_numbers))
Nbas = {}  # {atom_weight: nbas_atom}
Nbas_tot = pcoeffs.shape[0]
for a in atm_numbers:
    idx = [n for n,i in enumerate(mol._atm) if i[0] == a][0] #get only the first atom of that type
    Nbas[a] = ((mol._bas[mol._bas[:,0]==idx][:,ANG_OF]*2+1) * mol._bas[mol._bas[:,0]==idx][:,NCTR_OF]).sum()  
    
cntr=0
js_o={}
MOs=min(qcoeffs.shape[0], pcoeffs.shape[0])
for c in atm_numbers:
    b = cntr  # beginning
    e = b + Nbas[c]  # end
    q_slc = qcoeffs[:, b:e]
    p_slc = pcoeffs[:, b:e]
    cntr += Nbas[c]
    if str(c) in js_o.keys():
        if not repeat_elements:
            continue
    d_arr = np.zeros([Nbas[c],Nbas[c]], dtype=float)  # difference array
    for i in range(MOs):
        for j in range(Nbas[c]):
            d_arr[j] += abs(p_slc[i,j] - q_slc[i,:])
    ID= "p{}".format(md.element(int(c)).period) if rows_equal else md.element(int(c)).symbol
    if ID in js_o.keys():
        assert (list(map(int,np.argmin(d_arr, axis=1))) == js_o[ID]), "BIG ISSUE!! two centers of atomic weight "+str(c)+" have different assignment!!"
    js_o[ID] = list(map(int,np.argmin(d_arr, axis=1)))

d=dict(source="qchem", destination="pyscf", basisname="cc-pVDZ(seg-opt)" ,
                   basisfile=nwchem, order=js_o)
with open(os.path.join(root,"{}.corr".format(basname)),"w") as f:
    js.dump(d,f)


