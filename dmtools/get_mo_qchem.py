#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 13:52:30 2021

@author: nico
"""

import CCJob as ccj
import CCJob.utils as ut
from CCJob.Composable_templates import Tdefaults, Tinp, Trem_kw, Tmolecule, Tbasis
#import slurm_stuff_yggdrasil as slrm
import slurm_stuff_baobab as slrm
import sys
import os
import CCParser as ccp
import numpy as np

def get_qchem_coeffs(xyz, basname, print_orbitals=50):
    root = os.getcwd()
    basfile = basname + ".bas"
    meta_HF  = dict(method="HF", status=None, basename="HF")
    meta_HF["path"] = os.path.join(root, "qchem")
    already_done = ut.status_ok(path=meta_HF["path"])
    if already_done == False:
        memory = 14000
        bs_kw = "gen"  
        bs_string = ut.read_file(basfile)
        if bs_string[-1] == "\n":
            bs_string = bs_string[:-1]
        rem_extras = ["thresh = 14", "basis_lin_dep_thresh = 5"]
        extra_basic = Tbasis.substitute(**{"basis_specs": bs_string})
        rem_kw = Trem_kw.substitute(Tdefaults["rem_kw"], **dict(memory=memory, 
                                    method="hf", basis=bs_kw, print_orbitals=print_orbitals,
                                    rem_extras="\n".join(rem_extras)))
        with open(xyz, "r") as f:
            lines = f.readlines()
        coords = "".join(lines[2:])
        molecule = Tmolecule.substitute(Tdefaults["molecule"], **dict(xyz=coords))
        specs_HF = ut.myupd(Tdefaults["inp"], rem_kw=rem_kw, molecule=molecule, extras=extra_basic)
        queue_HF = dict(** slrm.shabug_S)  
        ut.run_job(specs_HF, queue_HF, meta_HF, Tinp,  
                batch_mode=False)  # because we want to extract data and copy matrices
        
        ut.save_status(meta_HF)
    if not os.path.isfile(os.path.join(meta_HF["path"], "CCParser.json")):  # generally it is there because of save_status
        ccp.Parser(os.path.join(meta_HF["path"], meta_HF["basename"]+".out"),
                                to_json=True, to_console=False, json_file="CCParser.json",
                                large_fn="matrices.npz")
    d = ut.load_js(os.path.join(meta_HF["path"], "CCParser.json"))
    coeffs = d["C"][-1][0]
    if type(coeffs) == str:
        coeffs = np.load(os.path.join(meta_HF["path"], coeffs), allow_pickle=True)["C"][-1]
    np.savetxt("{}_coeffs_qchem.txt".format(basname), coeffs)

if __name__ == "__main__":
    _, xyz, basname = sys.argv
    get_qchem_coeffs(xyz, basname)
