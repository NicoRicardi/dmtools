#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 11:38:10 2019

@author: nico
"""

import numpy as np

class BasisSet():
    def __init__(self, bas_name, source="", target="", cart=False):
        self.basis = bas_name
        self.isCartesian = cart #also not 100% practical, should be per shell
        self.source = source
        self.target = target
        self.first_row = ['H', 'He']
        self.second_row = ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']
        self.row1 = None
        self.row2 = None
        self.str1 = None
        self.str2 = None
        self.str1_ord = None
        self.str2_ord= None
        # set order based in basis name, source and target
        self.distributor()
    def distributor(self):
        if self.basis == 'cc-pVDZ':
            self.cc_pvdz()
            self.Nrow1 = len(self.row1)
            self.Nrow2 = len(self.row2)
        elif self.basis == 'cc-pVTZ':
            self.cc_pvtz()
            self.Nrow1 = len(self.row1)
            self.Nrow2 = len(self.row2)
        else:
            raise NotImplementedError("""[!] Only [cc-pVDZ] and \
                                      [cc-pVTZ] supported for now.""") 
    
    def cc_pvdz(self):      # Basis function library [cc-pVDZ]
        self.row1 = [n for n in range(5)]
        # strings of AOs in Q-Chem format as reference
        self.str1 = ["1s", "2s", "3px", "3py", "3pz"]
        if self.isCartesian: 
            # Gaussian --> Q-Chem
            if self.source == "gaussian" and self.target == "qchem":
                # AO strings in source format
                self.str2 = ["1s", "2s", "3s", "4px", "4py", "4pz", "5px",
                             "5py", "5pz", "6dxx", "6dyy", "6dzz", "6dxy",
                             "6dxz", "6dyz"]
                # d-orbitals in Gaussian: dxx, dyy, dzz, dxy, dxz, dyz
                self.row2 = [n for n in range(9)] + [9,12,10,13,14,11]
            # Q-Chem --> Gaussian
            elif self.source == "qchem" and self.target == "gaussian":
                # AO strings in source format
                self.str2 = ["1s", "2s", "3s", "4px", "4py", "4pz", "5px",
                             "5py", "5pz", "6dxx", "6dxy", "6dyy", "6dxz",
                             "6dyz", "6dzz"]
                # d-orbitals in Q-Chem: dxx, dxy, dyy, dxz, dyz, dzz
                self.row2 = [n for n in range(9)] + [9,11,14,10,12,13]
            # Molcas --> Q-Chem
            elif self.source == "molcas" and self.target == "qchem":
                # AO strings in source format
                self.str2 = ["1s", "2s", "3s", "4px", "5px", "4py", "5py",
                             "4pz", "5pz", "6dxx", "6dxy", "6dxz", "6dyy",
                             "6dyz", "6dzz"]
                # d-orbitals in Molcas: dxx, dxy, dxz, dyy, dyz, dzz
                self.row2 = [0,1,2,3,5,7,4,6,8,9,10,12,11,13,14]
            # Q-Chem --> Molcas
            elif self.source == "qchem" and self.target == "molcas":
                 # AO strings in source format
                self.str2 = ["1s", "2s", "3s", "4px", "4py", "4pz", "5px",
                             "5py", "5pz", "6dxx", "6dxy", "6dyy", "6dxz",
                             "6dyz", "6dzz"]
                self.row2 = [0,1,2,3,6,4,7,5,8,9,10,12,11,13,14]
            else:
                raise NotImplementedError("""Molcas <-> Gaussian not impl.""")
            
        else: 
            # Q-Chem --> Molcas
            if self.source == "qchem" and self.target == "molcas":
                # p-orbitals in Q-Chem: 2px, 2py, 2pz, 3px, 3py, 3pz
                # d-orbitals correspond directly, e.g. 6d1 = 6d2-, etc.
                self.str2 = ["1s", "2s", "3s", "4px", "4py", "4pz", "5px",
                             "5py", "5pz", "6d2-", "6d1-", "6d0", "6d1+", "6d2+"]
                self.row2 = [0,1,2,3,6,4,7,5,8,9,10,11,12,13]
            # Molcas --> Q-Chem
            elif self.source == "molcas" and self.target == "qchem":
                # p-orbitals in Molcas: 2px, 3px, 2py, 3py, 2pz, 3pz
                self.str2 = ["1s", "2s", "3s", "4px", "5px", "4py", "5py", 
                             "4pz", "5pz", "6d2-", "6d1-", "6d0", "6d1+", "6d2+"]
                self.row2 = [0,1,2,3,5,7,4,6,8,9,10,11,12,13]
            # Q-Chem --> Gaussian
            elif self.source == "qchem" and self.target == "gaussian":
                self.str2 = ["1s", "2s", "3s", "4px", "4py", "4pz", "5px",
                             "5py", "5pz", "6d2-", "6d1-", "6d0", "6d1+", "6d2+"]
                self.row2 = [n for n in range(9)] + [11,12,10,13,9]
            # Gaussian --> Q-Chem
            elif self.source == "gaussian" and self.target == "qchem":
                self.str2 = ["1s", "2s", "3s", "4px", "4py", "4pz", "5px",
                             "5py", "5pz", "6d0", "6d1+", "6d1-", "6d2+", "6d2-"]
                self.row2 = [n for n in range(9)] + [13,11,9,10,12]
#            # Q-chem --> Horton
#            elif self.source == "qchem" and self.target == "horton":
#                self.str2 = ["1s", "2s", "3s", "4px", "4py", "4pz", "5px",
#                             "5py", "5pz", "6d2-", "6d1-", "6d0", "6d1+", "6d2+"]
#                self.row2 = [n for n in range(9)] + [13,11,9,10,12]                
#            elif self.source == "horton" and self.target == "qchem":
#                self.str2 = ["1s", "2s", "3s", "4px", "4py", "4pz", "5px",
#                             "5py", "5pz", "6d0", "6d1+", "6d1-", "6d2+", "6d2-"]
#                self.row2 = [n for n in range(9)] + [11,12,10,13,9]                
            else:
                raise ValueError("Invalid value for source/target format.")
            
            
    def cc_pvtz(self):      # Basis function library [cc-pVTZ]
        if self.isCartesian:
            raise NotImplementedError("Order of cartesian AOs not implemented yet!")
        else:
            # Q-Chem --> Molcas
            if self.source == "qchem" and self.target == "molcas":
                # AO strings in source format
                self.str1 = ["1s", "2s", "3s", "4px", "4py", "4pz", "5px", "5py",
                             "5pz", "6d1", "6d2", "6d3", "6d4", "6d5"]
                self.row1 = [0,1,2,3,6,4,7,5,8,9,10,11,12,13]
                self.str2 = ["1s", "2s", "3s", "4s", "5px", "5py", "5pz", "6px",
                             "6py", "6pz", "7px", "7py", "7pz", "8d1", "8d2",
                             "8d3", "8d4", "8d5", "9d1", "9d2", "9d3", "9d4",
                             "9d5", "10f1", "10f2", "10f3", "10f4", "10f5", "10f6",
                             "10f7"]
                self.row2 = [0,1,2,3,4,7,10,5,8,11,6,9,12,13,18,14,19,15,20,16,21,
                             17,22,23,24,25,26,27,28,29]
            # Molcas --> Q-Chem
            elif self.source == "molcas" and self.target == "qchem":
                # AO strings in source format
                self.str1 = ["1s", "2s", "3s", "4px", "5px", "4py", "5py","4pz", 
                             "5pz", "6d2-", "6d1-", "6d0", "6d1+", "6d2+"]
                self.row1 = [0,1,2,3,5,7,4,6,8,9,10,11,12,13]
                self.str2 = ["1s", "2s", "3s", "4s", "5px", "6px", "7px", "5py",
                             "6py", "7py", "5pz", "6pz", "7pz", "8d2-", "9d2-",
                             "8d1-", "9d1-", "8d0", "9d0", "8d1+", "9d1+", "8d2+",
                             "9d2+", "10f3-", "10f2-", "10f1-", "10f0", "10f1+",
                             "10f2+", "10f3+"]
                self.row2 = [0,1,2,3,4,7,10,5,8,11,6,9,12,13,15,17,19,21,14,16,18,
                             20,22,23,24,25,26,27,28,29]
            else:
                raise NotImplementedError("""Order of pure AOs in Gaussian \
                                          format notDMqa= implemented yet!""")
def get_SortArray(bas,atomlist):
    sort = []
    a = 0   # 1st row iterator
    b = 0   # 2nd row iterator
    for atom in atomlist:
        if atom in bas.first_row:
            for n in bas.row1:
                sort.append(n +a*bas.Nrow1 +b*bas.Nrow2)
            a = a+1
        if atom in bas.second_row:
            for m in bas.row2:
                sort.append(m +a*bas.Nrow1 +b*bas.Nrow2)
            b = b+1
    return sort
    
def transform(inparr, bas, atomlist, source="", target=""):
    """Perform symmetric transformation of a matrix.

    The symmetric transformation of a matrix :math:`\\mathbf{X}` uses a
    rearranged identity matrix :math:`\\mathbf{P}` such that the working
    equation is:
    :math:`\\mathbf{P}^{T}\\cdot \\mathbf{X} \\cdot\\mathbf{P}`
 
    Parameters:
    -----------
    inparr : np.ndarray
        Input matrix :math:`\\mathbf{X}`.
    source : string
        Source format ('qchem', 'molcas' or 'gaussian').
    target : string
        Target format ('qchem', 'molcas' or 'gaussian').
    
    Returns:
    --------
    Q : np.ndarray
        Transformed matrix according to target format.
    """
    # -------------------
    # Q = P.T * X * P
    # -------------------
    sort_array = get_SortArray(bas,atomlist)
    nAO=inparr.shape[0]
    idarr = np.identity(nAO)
    # do transformation
    P = idarr[:, sort_array] # black magic: rearrange columns of ID matrix
    M = np.dot(P.T, inparr)
    Q = np.dot(M,P)
    
    return Q