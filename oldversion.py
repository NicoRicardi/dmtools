#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 11:24:26 2019

@author: nico
"""
import numpy as np


def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

class DM():
    def __init__(self,coeffs,alpha_only):
        self.coeffs=coeffs
        self.alpha_only=alpha_only
    def trace(self):
        return self.coeffs.trace()
    def from_dmfile(fname,in_alpha_only=False,out_alpha_only=True):
        return read_dmfile(fname,in_alpha_only)
    def from_qcfile(fname):
        print("This may take a while, please be patient")
        return DM_from_qcfile(fname)
    def from_fchk(fname):
        return DM_from_fchk(fname)
    def to_file(self,fname,header=True, twice=False,exponential=True, digits=10):
        N_bfs=self.coeffs.shape[0]
        format_ = "e" if exponential else "f"
        with open(fname,"w") as out:
            b = 2 if twice else 1
            if  header:
                out.write(str(b)+" {0:d} {1:d}\n".format(N_bfs, N_bfs))
            for k in range(b):
                for i in range(self.coeffs.shape[0]):
                    for j in range(self.coeffs.shape[1]):
                        print('{: .{dgt}{frmt}} {}'.format(self.coeffs[i,j],'\n', dgt=digits,frmt=format_))
    def diff(self,other):
        return DM_diff(self,other)
    def maxdiff(self,other):
        return max_DM_diff(self,other)
    def tracediff(self,other):
        return self.trace()-other.trace()
    def to_alphabeta(self):
        if self.alpha_only==False:
            print("It is already a 2-particle density matrix")
            return self
        elif self.alpha_only==True:
            print("doubling alpha density matrix")
            return DM(np.multiply(2,self.coeffs),alpha_only=False)
    
class HFcoeffs:
    def __init__(self, coeffs):
        self.coeffs = coeffs
    def from_fchk(fname):
        return coeffs_from_fchk(fname)
    def cdiff(self,other):
        return cdiff(self,other)
    def acdiff(self,other):
        return acdiff(self,other)
    def acadiff(self,other):
        return acadiff(self,other)
    
def coeffs_from_fchk(fname):
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
    diff=a.coeffs-b.coeffs
    return diff

def acdiff(a,b):
    diff=abs(a.coeffs-b.coeffs)
    return diff

def acadiff(a,b):
    diff=abs(abs(a.coeffs)-abs(b.coeffs))
    return diff

#def read_dmfile(file):
#    with open(file) as f:
#        Nbas=int(f.readline().split()[1])
#        coeffs=f.read().splitlines()[0:int(Nbas)**2]
#        coeffs=[float(x) for x in coeffs]
#        coeffs=np.array(coeffs)
#        coeffs=coeffs.reshape(Nbas,Nbas)
#        return DM(coeffs)
    
#def read_dmfile(file,in_alpha_only=False,out_alpha_only=False):
#    with open(file) as f:
#        fl=f.readlines()
#        if len(fl[0].split())==1:
#            Nbas=int(np.sqrt(len(fl)))
#        elif len(fl[0].split())==3:
#            Nbas=int(fl[0].split()[1])
#            fl=fl[1:]
#        else:
#            print("Something is wrong!")
#        if in_alpha_only==False:
#            fl=fl[:len(fl)/2]
#        coeffs=fl[0:int(Nbas)**2]
#        coeffs=[float(x) for x in coeffs]
#        coeffs=np.array(coeffs)
#        coeffs=coeffs.reshape(Nbas,Nbas)
#        dm=DM(coeffs,alpha_only=True)
#        if out_alpha_only==False:
#            dm.to_alphabeta()
#        return dm

def read_dmfile(file,in_alpha_only=False,out_alpha_only=True):
    import sys
    with open(file) as f:
        fl=f.readlines()
        tml=fl[0].split()
        if len(tml)==3:
            Nbas=int(tml[1])
            check=True if tml[0]=="1" else False
            if check!=in_alpha_only:
                print("It seems like the file is not how you expected it to be!!")
                sys.exit()
            fl=fl[1:]
        if in_alpha_only==False:
            fl=fl[:int(len(fl)/2)]
        if len(tml)==1:
            Nbas=int(np.sqrt(len(fl)))
        if len(tml)!=1 and len(tml)!=3:
            print("Something is wrong!")
        coeffs=fl[0:int(Nbas)**2]
        coeffs=[float(x) for x in coeffs]
        coeffs=np.array(coeffs)
        coeffs=coeffs.reshape(Nbas,Nbas)
        dm=DM(coeffs,alpha_only=True)
        if out_alpha_only==False:
            dm.to_alphabeta()
    return dm 
    
def DM_from_qcfile(fname,alpha_only=True):
    import CCParser as cc
    data=cc.Parser(fname,software="qchem",to_console=False)
    dm=DM(data.results.P_alpha.get_last(),alpha_only=True)
    if alpha_only==False:
        dm.to_alphabeta()
    return dm

def DM_from_fchk(fname,alpha_only=True):
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
    dm=DM(P,alpha_only=True)
    if alpha_only==False:
        dm.to_alphabeta()
    return dm

def DM_diff(a,b):
    diff=a.coeffs-b.coeffs
    return diff

def max_DM_diff(a,b):
    diff=a.coeffs-b.coeffs
    return diff.max()
