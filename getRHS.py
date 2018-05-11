#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 23:16:37 2018

@author: emilyng
"""
import numpy as np
#from RHS_eqs import HI, HII, HM, HeI, HeII, HeIII, H2I, H2II, de


import sympy

HI, HII, HM, HeI, HeII, HeIII, H2I, H2II = sympy.sympify("HI, HII, HM, HeI, HeII, HeIII, H2I, H2II")
k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, k20, k21, k22, k23, k24, k25, k26, k27, k28, k29, k30, k31, k32 = sympy.sympify("k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, k20, k21, k22, k23, k24, k25, k26, k27, k28, k29, k30, k31, k32")    
  
r1 = (HI), (HII), k1
r2 = (HII), (HI), k2
r3 = (HeI), (HeII), k3
r4 = (HeII), (HeI), k4
r5 = (HeII), (HeIII), k5
r6 = (HeIII), (HeII), k6
r7 = (HI), (HM), k7
r8 = (HM + HI), (H2I), k8
r9 = (HI + HII), (H2II), k9
r10 = (H2II + HI), (H2I + HII), k10
r11 = (H2I + HII), (H2II + HI), k11
r12 = (H2I), (2*HI), k12
r13 = (H2I + HI), (3*HI), k13
r14 = (HM), (HI ), k14
r15 = (HM + HI), (2*HI), k15
r16 = (HM + HII), (2*HI), k16
r17 = (HM + HII), (H2II), k17
r18 = (H2II), (2*HI), k18
r19 = (H2II + HM), (HI + H2I), k19

all_reactions = [r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11 ,r12, r13, r14, r15, r16, r17, r18, r19]

def find_formation(species):
    rxns = []
    for r in all_reactions:
        if species in r[1].atoms():
            rxns.append(r)
    return rxns

def find_destruction(species):
    rxns = []
    for r in all_reactions:
        if species in r[0].atoms():
            rxns.append(r)
    return rxns

def get_rhs(species):
    dSdt = 0
    for lhs, rhs, coeff in find_formation(species):
        term = coeff
        for atom in list(lhs.atoms()):
            term *= atom
        dSdt += term
    for lhs, rhs, coeff in find_destruction(species):
        term = -coeff
        for atom in list(lhs.atoms()):
            term *= atom
        dSdt += term
    return dSdt
