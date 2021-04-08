#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 11:24:26 2019

@author: nico
"""
import numpy as np
import itertools as ittl

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class otherError(Error):
    """Some error.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message="some error!"):
        self.message = message
        print(self.message)

def isfloat(value):
    """
    Parameters
    ----------
    value: float
        value/variable to check
    
    Returns
    -------
    bool
        whether it is a float or not
    """
    try:
        float(value)
        return True
    except ValueError:
        return False

class DM():
    """
    Note
    ----
    main class in this module. Stores alpha and beta density matrix
    """
    def __init__(self, coeffs_a, coeffs_b, comment=""):
        """
        Parameters
        ----------
        coeffs_a: array
            alpha DM coefficients
        coeffs_b: array
            beta DM coefficients
        comment: str
            an optional comment
        """
        self.coeffs = {"a" : coeffs_a, "b" : coeffs_b}
        self.closed_shell = (self.coeffs["a"] - self.coeffs["b"] < 0.000001).all()
        self.comment = comment 
        
    def get_dm_full(self):
        """
        Note
        ----
        for closed shell systems, provides the sum of alpha and beta DMs.
        Calculated only if not stored yet.
        
        Returns
        -------
        array
            the full DM
        """
        if not self.closed_shell:
            raise otherError("Not closed shell!")
        else:
            if not hasattr(self,"dm_full"):
                self.dm_full = np.multiply(2,self.coeffs["a"])
            return self.dm_full
        
    def copy(self):
        """
        Returns
        -------
        DM
            a copy of self
        """
        import copy as c
        return c.deepcopy(self)
    
    def __sub__(self, other):
        """
        Note
        ----
        Magic method for a-b
        
        Parameters
        ----------
        other: DM
            the other DM
        
        Returns
        -------
        array
            the difference array
        """
        return self.diff(other)
    
    def __div__(self, other):
        """
        Note
        ----
        Magic method for a/b
        
        Parameters
        ----------
        other: DM
            the other DM
        
        Returns
        -------
        float
            the maximum of the difference array
        """
        return self.maxdiff
    
    def __floordiv__(self, other):
        """
        Note
        ----
        Magic method for a//b
        
        Parameters
        ----------
        other: DM
            the other DM
        
        Returns
        -------
        float
            the maximum of the absolute difference array
        """
        return self.absmaxdiff
    
    def __or__(self, other):
        """
        Note
        ----
        Magic method for a|b
        
        Parameters
        ----------
        other: DM
            the other DM
        
        Returns
        -------
        float
            the difference of the traces
        """
        return self.tracediff(other)
    
    def __invert__(self):
        """
        Note
        ----
        Magic method for ~a
        
        Returns
        -------
        float
            the trace
        """
        return self.trace()
        
    def trace(self,spin="a"):
        """
        Parameters
        ----------
        spin: {"a", "b"}
            matrix to get the trace of
        
        Returns
        -------
        float
            the trace of the alpha/beta DM
        """
        return self.coeffs[spin].trace()
    
    def from_dmfile(fname):
        """
        Note
        ----
        Reads both "alpha_only" or "alpha,beta" files, header or not
        
        Parameters
        ----------
        fname: str
            filename or path
            
        Returns
        -------
        DM object
        """
        return read_dmfile(fname)
    
    def from_qcfile(fname):
        """
        Note
        ----
        Parses from qchem output
        
        Parameters
        ----------
        fname: str
            filename or path
            
        Returns
        -------
        DM object
        """
        print("This may take a while, please be patient")
        return DM_from_qcfile(fname)
    
    def from_fchk(fname):
        """
        Note
        ----
        DM from fchk file (uses SCF coefficients) #TODO use DM in upper triangular
        
        Parameters
        ----------
        fname: str
            filename or path
            
        Returns
        -------
        DM object
        """        
        return DM_from_fchk(fname)
    
    def from_HFcoeffs(coeffs,alpha_only_in=True):
        """
        Parameters
        ----------
        coeffs: array
            HFcoeffs
        alpha_only_in: bool
            whether these coeffs are alpha+beta or only alpha
            
        Returns
        -------
        DM
            DM object
        """
        return DM_from_HFcoeffs(coeffs,alpha_only_in=True,alpha_only_out=False)
    
    def to_file(self, fname,header=True, alpha_only=False, exponential=True, digits=10):
        """
        Note
        ----
        writes down to file. type of file can be chosen
        
        Parameters
        ----------
        fname: str
            filename or path
        header: bool
            whether there should be a header (2 Nbas Nbas)
        alpha_only: bool
            whether it should write only the alpha matrix or alpha,beta
        exponential: bool
            whether the values should be written in exponential form
        digits: int
            the number of digits desired
        """
        s = ["a","b"]
        N_bfs = self.coeffs["a"].shape[0]
        format_ = "e" if exponential else "f"
        with open(fname,"w") as out:
            spin = 1 if alpha_only else 2
            if  header:
                out.write(str(spin)+" {0:d} {1:d}\n".format(N_bfs, N_bfs))
            for k in range(spin):
                for i in range(self.coeffs[s[k]].shape[0]):
                    for j in range(self.coeffs[s[k]].shape[1]):
                        out.write('{: .{dgt}{frmt}} {}'.format(self.coeffs[s[k]][i,j],'\n', dgt=digits,frmt=format_))
    
    def to_numpy_txt(self, fname, alpha_only=False):
        """
        Note
        ----
        writes down to numpy txt. Either 2*alpha or only alpha
        
        Parameters
        ----------
        fname: str
            filename or path
        alpha_only: bool
            whether it should write only the alpha matrix or alpha,beta
        """
        if alpha_only:
            np.savetxt(fname,self.coeffs["a"])
        else:
            np.savetxt(fname,np.multiply(2,self.coeffs["a"]))
                        
    def diff(self, other, spin="a"):
        """
        Note
        ----
        Returns directly the array(diff)!
        
        Parameters
        ----------
        other: DM
            the dm to subtract
        spin: {"a", "b"}
        
        Returns
        -------
        array
            the difference array
        """
        return DM_diff(self, other, spin)
    def maxdiff(self, other, spin="a"):
        """
        Note
        ----
        Returns directly the array(diff)!
        
        Parameters
        ----------
        other: DM
            the dm to subtract
        spin: {"a", "b"}
        
        Returns
        -------
        float
            the maximum value of the difference array
        """
        return max_DM_diff(self, other, spin)
    def absmaxdiff(self, other, spin="a"):
        """
        Note
        ----
        Returns directly the array(diff)!
        
        Parameters
        ----------
        other: DM
            the dm to subtract
        spin: {"a", "b"}
        
        Returns
        -------
        float
            the maximum value of the absolute difference array
        """
        return absmax_DM_diff(self, other, spin)
    def tracediff(self, other, spin="a"):
        """
        Note
        ----
        Returns directly the array(diff)!
        
        Parameters
        ----------
        other: DM
            the dm to subtract
        spin: {"a", "b"}
        
        Returns
        -------
        float
            the difference of traces
        """
        return self.trace(spin) - other.trace(spin)
    
    def reorder(self, geom, order, basis):
        """
        Note
        ----
        Reorders the atoms in the density matrix
        
        Parameters
        ----------
        geom: geom_obj, str
            geometry object, xyz filename, or xyz str
        order: np.array
            the order to use to reorder (e.g. [1,0,3,2,5,4] to swap every two atoms)
        basis: str
            basis name or *.nwchem filepath
            
        Returns
        -------
        DM obj
            the reordered DM
        """
        if geom != str:
            geom = geom.__str__()
        sort_arr = get_atoms_sortarr(geom, order, basis)
        coeffs_a, coeffs_b = reorder(self.coeffs["a"], sort_arr), reorder(self.coeffs["b"], sort_arr)
        return DM(coeffs_a, coeffs_b, comment="reordered with arr {}".format(order))
    
class HFcoeffs:
    """
    HF coefficient matrix class
    """
    def __init__(self, coeffs):
        """
        Parameters
        ----------
        coeffs: array
            array of the HF coefficients
        """
        self.coeffs = coeffs
        
    def from_fchk(fname):
        """
        Parameters
        ----------
        fname: str
            fchk filename or path
            
        Returns
        -------
        HFcoeffs
            the HF coefficients object
        """        
        return coeffs_from_fchk(fname)
    
    def cdiff(self,other):
        """
        Parameters
        ----
        Reads both "alpha_only" or "alpha,beta" files, header or not
        ----------
        other: HFcoeffs
            HF coefficients to do the difference of
            
        Returns
        -------
        array
            the difference array
        """    
        return cdiff(self,other)
    
    def acdiff(self,other):
        """
        Parameters
        ----------
        other: HFcoeffs
            HF coefficients to do the difference of
            
        Returns
        -------
        array
            the absolut difference array
        """            
        return acdiff(self,other)
    def acadiff(self,other):
        """
        Parameters
        ----------
        other: HFcoeffs
            HF coefficients to do the difference of
            
        Returns
        -------
        array
            the absolut difference of absolute values array
        """     
        return acadiff(self,other)
    
def get_atoms_sortarr(geom, order, basis):
    """
    Note
    Returns the sorting array, to obtain a DM, where the atoms in "geom" have been rearranged according to "order".
    
    Parameters
    ----------
    geom: str
        filename or string of the xyz
    order: np.array
        the order to use to reorder (e.g. [1,0,3,2,5,4] to swap every two atoms)
    basis: str
        basis name or *.nwchem filepath
    
    Returns
    -------
    np.array
        the sorting array
    """
    from pyscf import gto
    mol = gto.M(atom=geom, basis=basis)
    tmp = gto.aoslice_by_atom(mol)
    sort_arr = []
    for o in order:
        slc = tmp[o][[2,3]]
        sort_arr.append(np.arange(slc[0], slc[1]))
    sort_arr = list(ittl.chain.from_iterable(sort_arr))
    return sort_arr

def reorder(arr, sort_arr):
    """
    Parameters
    ----------
    arr: np.array
        array of the density matrix
    sort_arr: np.array
        sorting array
    
    Returns
    -------
    np.array
        the DM array where the rows have been rearranged
    """
    I = np.identity(arr.shape[0])
    P = I[:, sort_arr] # black magic: rearrange columns of ID matrix
    M = np.dot(P.T, arr)
    Q = np.dot(M,P)
    return Q
    
    
    
def coeffs_from_fchk(fname):
    """
    Parameters
    ----------
    fname: str
        fchk filename or path
        
    Returns
    -------
    HFcoeffs
        HFcoeffs object
    """     
    with open(fname) as f:
        fl=f.readlines()
    MO_coeffs=[]
    for i in range(20):
        if "Number of alpha" in fl[i]:
            N_el=int(fl[i].split()[5])
        if "Number of beta" in fl[i]:
            N_el2=int(fl[i].split()[5])
    
    assert(N_el==N_el2), "Watch out: open shell!"  #checking it is closed shell
#    parse=False
    for i, line in enumerate(fl):
        if "Alpha MO" in line:
    #        parse=True
            N_bfs=int(np.sqrt(int(line.split()[5])))
            k = 0
            while True:
                tmp_list = fl[i+1+k].split()
                if not isfloat(tmp_list[0]):
                    break
                MO_coeffs.extend(list(map(float,tmp_list)))
                k += 1
    coeffs=np.array(MO_coeffs)
    coeffs.shape=(N_bfs,N_bfs)
    return HFcoeffs(coeffs)

def cdiff(a,b):
    """
    Parameters
    ----------
    a: HFcoeffs
        first HFcoeff object
    b: HFcoeffs
        second HFcoeff object    
        
    Returns
    -------
    array
        the difference array
    """  
    diff=a.coeffs-b.coeffs
    return diff

def acdiff(a,b):
    """
    Parameters
    ----------
    a: HFcoeffs
        first HFcoeff object
    b: HFcoeffs
        second HFcoeff object         cf2: cubefile
        the cubefile to compare to
        
    Returns
    ------
    array
        the absolute difference array
    """  
    diff=abs(a.coeffs-b.coeffs)
    return diff

def acadiff(a,b):
    """
    Parameters
    ----------
    a: HFcoeffs
        first HFcoeff object
    b: HFcoeffs
        second HFcoeff object    
        
    Returns
    -------
    array
        the absolute difference of absolute values array
    """  
    diff=abs(abs(a.coeffs)-abs(b.coeffs))
    return diff

    
def read_dmfile(fname):
    """
    Note
    ----
    Reads both "alpha_only" or "alpha,beta" files, header or not
    
    Parameters
    ----------
    fname: str
        fchk filename or path
        
    Returns
    -------
    DM 
        DM object
    """    
    with open(fname) as f:
        fl=f.readlines()
        tml=fl[0].split()
        if len(tml)==3:
            Nbas=int(tml[1])
            spinless = True if tml[0]=="1" else False
            fl=fl[1:]
        Float=np.sqrt(len(fl))
        Int=int(np.sqrt(len(fl)))
        if Float-Int < 0.000001:
            spinless=True
            Nbas=Int
        else:
            spinless=False
            Nbas=int(np.sqrt(0.5*len(fl)))
        if spinless==False:
            a,b = fl[:int(len(fl)/2)], fl[int(len(fl)/2):]
        else:
            a,b = fl, fl
        if len(tml) not in [1,3]:
            print("Something is wrong!")
        coeffs={"a":None,"b":None}
        l=[a,b]
        for n,s in enumerate(coeffs.keys()):
            tmp=l[n][0:int(Nbas)**2]
            tmp=[float(x) for x in tmp]
            tmp=np.array(tmp)
            coeffs[s]=tmp.reshape(Nbas, Nbas)
        dm=DM(coeffs["a"], coeffs["b"])
    return dm 
    
def DM_from_qcfile(fname):
    """
    Note
    ----
    Parses from qchem output
    
    Parameters    cf2: cubefile
        the cubefile to compare to
    ----------
    fname: str
        filename or path
        
    Returns
    -------
    DM object
    """
    import CCParser as cc
    data=cc.Parser(fname,software="qchem",to_console=False)
    dm=DM(data.results.P_alpha.get_last(),data.results.P_alpha.get_last())#only closed shell now!!!
    return dm

def DM_from_fchk(fname): #TODO parse UT DM
    """
    Note
    ----
    Parses from fchk. NB only SCF from coeffs!
    
    Parameters
    ----------
    fname: str
        filename or path
        
    Returns
    -------
    DM object
    """
    with open(fname) as f:
        fl=f.readlines()
    MO_coeffs=[]
    for i in range(20):
        if "Number of alpha" in fl[i]:
            N_el=int(fl[i].split()[5])
        if "Number of beta" in fl[i]:
            N_el2=int(fl[i].split()[5])
    assert(N_el==N_el2), "Watch out: open shell!"  #checking it is closed shell
    for i, line in enumerate(fl):
        if "Alpha MO" in line:
            N_bfs=int(np.sqrt(int(line.split()[5])))
            k = 0
            while True:
                tmp_list = fl[i+1+k].split()
                if not isfloat(tmp_list[0]):
                    break
                MO_coeffs.extend(tmp_list)
                k += 1
    coeffs=np.array(MO_coeffs)
    coeffs.shape=(N_bfs,N_bfs)
    P=np.array(np.zeros((N_bfs,N_bfs)))
    for i in range(N_bfs):
        for j in range(N_bfs):
            for k in range(N_el):
                P[i,j]+=float(coeffs[k,i])*float(coeffs[k,j])
    dm=DM(P,P)
    return dm

def DM_from_HFcoeffs(coeffs, N_el, alpha_only_in=True):
    """
    Note
    ----
    DM from coeffs and N_el
    
    Parameters
    ----------
    coeffs: array (or list)
        array of the MO coefficients
    N_el: int
        number of electrons
    alpha_only_in: bool
        whether the input is just alpha MO coeffs or alpha+beta
    Returns
    -------
    DM object
    """
    if type(coeffs)==list:
        coeffs=np.asarray(coeffs)
    shape=coeffs.shape
    if shape[0]==shape[1]:
        N_bfs=shape[0]
    else:
        N_bfs=np.sqrt(max(shape))
        coeffs=coeffs.reshape(N_bfs,N_bfs)
    if not alpha_only_in:
        coeffs=np.divide(coeffs,2)
        N_el=int(N_el/2)
    P=np.array(np.zeros((N_bfs,N_bfs)))
    for i in range(N_bfs):
        for j in range(N_bfs):
            for k in range(N_el):
                P[i,j]+=float(coeffs[k,i])*float(coeffs[k,j])
    dm=DM(P,P)
    return dm

def DM_diff(x,y,spin="a"):
    """
    Parameters
    ----------
    x: DM
        first DM object
    y: DM
        second DM object
    spin: {"a", "b"}
        what spin to calculate the diff for
    Returns
    -------
    array
        the difference array
    """
    diff=x.coeffs[spin]-y.coeffs[spin]
    return diff

def max_DM_diff(x,y,spin="a"):
    """
    Parameters
    ----------
    x: DM
        first DM object
    y: DM
        second DM object
    spin: {"a", "b"}
        what spin to calculate the diff for
    Returns
    -------
    float
        the maximum value of the difference array 
    """
    diff=x.coeffs[spin]-y.coeffs[spin]
    return diff.max()

def absmax_DM_diff(x,y,spin="a"):
    """
    Parameters
    ----------
    x: DM
        first DM object
    y: DM
        second DM object
    spin: {"a", "b"}
        what spin to calculate the diff for
    Returns
    -------
    float
        the maximum value of the absolute difference array 
    """
    diff=abs(x.coeffs[spin]-y.coeffs[spin])
    return diff.max()