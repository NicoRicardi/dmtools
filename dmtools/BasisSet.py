#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 11:38:10 2019

@author: nico
"""

import numpy as np

class correspondence():
    def __init__(self, source, destination, basisname, basisfile, order):
        self.source = source
        self.destination = destination
        self.basisname = basisname
        self.basisfile = basisfile
        self.order = {i : np.array(order[i]) for i in order.keys()}
        
    def from_file(fname):
        return read_corresp_file(fname)
    
    def to_file(self, fname):
        write_corresp_file(self,fname)
        
    def get_reverse_order(self):
        if not hasattr(self, "reverse_order"):
            self.reverse_order=calc_reverse_order(self)
        return self.reverse_order
    
def read_corresp_file(fname):
    import json as js
    fname = fname+".corr" if "." not in fname else fname
    with open(fname,"r") as f:
        c = js.load(f)
    return correspondence(c["source"],c["destination"],c["basisname"],c["basisfile"],c["order"])

def write_corresp_file(corr_obj, fname):
    import json as js
    to_dump = dict(source = corr_obj.source, destination=corr_obj.destination, basisname=corr_obj.basisname ,
                   basisfile=corr_obj.basisfile, order=corr_obj.order.to_list(), rows_equal=corr_obj.rows_equal)
    fname = fname+".corr" if "." not in fname else fname
    with (fname,"w") as f:
        js.load(to_dump,f)
        
def calc_reverse_order(corr_obj):
    olist={i: corr_obj.order[i].tolist() for i in corr_obj.order.keys()}
    return {k: np.array([olist[k].index(i) for i in range(len(olist[k]))]) for k in olist.keys()}

def get_order(fname):
    import json as js
    fname = fname+".corr" if "." not in fname else fname
    with open(fname,"r") as f:
        c = js.load(f)
    return {i: np.array(c["order"][i]) for i in c["order"].keys()}

def get_rev_order(fname):
    import json as js
    fname = fname+".corr" if "." not in fname else fname
    with open(fname,"r") as f:
        olist = js.load(f)["order"]
    return {k: np.array([olist[k].index(i) for i in range(len(olist[k]))]) for k in olist.keys()}

def get_sort_arr(atomlist, order):
    import mendeleev as md
    arr = np.array([],dtype=int)
    cntr = 0
    for n,atom in enumerate(atomlist):
        if atom not in order.keys():
            try:
                ordr=order["p{}".format(md.element(atom).period)]
            except:
                print("you don't have the order for this period")
        else:
            ordr=order[atom]
        to_app=ordr+cntr
        arr=np.append(arr,to_app)
        cntr+=ordr.shape[0]
    return arr
        
def reorder(inp, sort_arr):
    I = np.identity(inp.shape[0])
    P = I[:, sort_arr] # black magic: rearrange columns of ID matrix
    M = np.dot(P.T, inp)
    Q = np.dot(M,P)
    return Q
            
        