#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 11 15:19:11 2018

@author: emilyng
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sint
import rate_coeff as rc
import timeit

'''RHS equations follow as: [destruction term] - [construction term]'''
def rhs(t, state):
    dnHIdt = state[5]*state[1]*rc.k11(state[9]) + state[5]*state[8]*rc.k12(state[9]) \
               - state[6]*state[0]*rc.k10(state[9]) + state[6]*state[7]*rc.k19(state[9]) \
               + state[6]*state[8]*rc.k18(state[9]) - state[0]*state[1]*rc.k9(state[9]) \
               - state[0]*state[7]*rc.k8(state[9]) - state[0]*state[8]*rc.k1(state[9]) \
                   + state[0]*state[8]*rc.k7(state[9]) + state[1]*state[7]*rc.k16(state[9]) \
                       + state[1]*state[8]*rc.k2(state[9]) + state[7]*state[8]*rc.k14(state[9])
    
    dnHIIdt = - state[5]*state[1]*rc.k11(state[9]) + state[6]*state[0]*rc.k10(state[9]) \
              - state[0]*state[1]*rc.k9(state[9]) + state[0]*state[8]*rc.k1(state[9]) \
              - state[1]*state[7]*rc.k16(state[9]) - state[1]*state[7]*rc.k17(state[9]) \
              - state[1]*state[8]*rc.k2(state[9])

    dnHeIdt = - state[2]*state[8]*rc.k3(state[9]) - state[3]*state[8]*rc.k4(state[9])

    dnHeIIdt = state[2]*state[8]*rc.k3(state[9]) - state[3]*state[8]*rc.k4(state[9])\
             - state[3]*state[8]*rc.k5(state[9]) + state[4]*state[8]*rc.k6(state[9])
    
    dnHeIIIdt = state[3]*state[8]*rc.k5(state[9]) - state[4]*state[8]*rc.k6(state[9])
    
    dnH2Idt = - state[5]*state[0]*rc.k13(state[9]) - state[5]*state[1]*rc.k11(state[9]) \
              - state[5]*state[8]*rc.k12(state[9]) + state[6]*state[8]*rc.k10(state[9]) \
              + state[6]*state[7]*rc.k19(state[9]) + state[0]*state[7]*rc.k8(state[9])
    
    dnH2IIdt = state[5]*state[1]*rc.k11(state[9]) - state[6]*state[0]*rc.k10(state[9]) \
             - state[6]*state[7]*rc.k19(state[9]) - state[6]*state[8]*rc.k18(state[9]) \
             + state[0]*state[1]*rc.k9(state[9]) + state[1]*state[7]*rc.k17(state[9])
    
    dnHmdt = - state[6]*state[7]*rc.k19(state[9]) - state[0]*state[7]*rc.k15(state[9]) \
             - state[0]*state[7]*rc.k8(state[9]) + state[0]*state[8]*rc.k7(state[9]) \
             - state[1]*state[7]*rc.k17(state[9]) - state[7]*state[8]*rc.k14(state[9])

    dnemdt = 0        
    dTdt = 0
    
    return np.array([ dnHIdt, dnHIIdt, dnHeIdt, dnHeIIdt, dnHeIIIdt, 
            dnH2Idt, dnH2IIdt, dnHmdt, dnemdt, dTdt])
  

#import ipywidgets
#
#@ipywidgets.interact(n_total = (-3.0, 9.0), e_frac = (-8.0, 0.0), rho = (1.0, 1000), \
#                     T = (0., 6.), final_t = (0.0, 25.0),
#                    safety_factor = (0.1, 10000))

def evolve(n_total = 2, H_frac = 0.76, He_frac = 0.24, H2_frac = 0.2, rho = 1e-20, final_t = 3, 
           e_frac = -4, T = np.log10(15000), safety_factor = 0.01):
    final_t = 10**final_t * 365*24*3600
    
    n_HI_initial = 10**n_total * (1.0 - 10**e_frac) *H_frac
    n_HII_initial = 10**e_frac*H_frac*10**n_total
    n_HeI_initial = 10**n_total * 10**e_frac * He_frac
    n_HeII_initial = 10**-7*10**n_total
    n_HeIII_initial = 10**-7*10**n_total
    n_H2I_initial = 10**n_total * H2_frac *10**-7
    n_H2II_initial = 10**-7*10**n_total
    n_HM_initial = 10**-7*10**n_total
    n_em_initial = 10**e_frac*10**n_total
    
    k = 1.38e-23
    gamma = 5/3
    m_H = 1.0/(1.67e-27)
    state_vector = np.array([n_HI_initial, n_HII_initial, n_HeI_initial, 
                             n_HeII_initial, n_HeIII_initial, n_H2I_initial,
                             n_H2II_initial, n_HM_initial, n_em_initial, 10**T])
    
    integrator = sint.ode(rhs)
    integrator.set_initial_value(state_vector, t=0)
    integrator.set_integrator('lsoda', method='bdf',
                             atol=state_vector*0.0, rtol=state_vector*0+0.1,)
                             #tolsf = 2.0)
    N = integrator.y[0] + integrator.y[1] + integrator.y[2] \
        + integrator.y[3] + integrator.y[4] + integrator.y[5] \
        + integrator.y[6] + integrator.y[7] + integrator.y[8]        
        
    mu = rho/N
    e = k*T/((gamma - 1)*mu*m_H) 
    
    state_vector_values = []
    nem_values = [n_em_initial]
    T_values = [T]
    ts = []
    dt = final_t / safety_factor
    ts.append(integrator.t)
    state_vector_values.append(integrator.y)
    #num= 0
    t1 = timeit.timeit()
    while integrator.t < final_t: #and num < 100:
        #S = integrator.y
        #dS = 10**(-4) + rhs(integrator.t, integrator.y)
        #dt = safety_factor * np.min(S/dS)
        
        integrator.integrate(integrator.t + dt)
        ts.append(integrator.t)
        n = integrator.y[0] + integrator.y[1] + integrator.y[2] \
            + integrator.y[3] + integrator.y[4] + integrator.y[5] \
            + integrator.y[6] + integrator.y[7] + integrator.y[8]
        mu = rho/n
        T = e * mu * (gamma - 1) / k
        integrator.y[9] = T
        state_vector_values.append(integrator.y)
        n_em = integrator.y[1] + integrator.y[3] + 2*integrator.y[4] \
             + integrator.y[6] - integrator.y[7]
        nem_values.append(n_em)
        T_values.append(T)
    t2 = timeit.timeit()
    print("Time: %.6f"%(abs(t2-t1)))
    
    state_vector_values = np.array(state_vector_values)
    ts = np.array(ts)
    plt.figure(figsize=(20,10))
    plt.loglog(ts, state_vector_values[:,0], linewidth=5, label='HI')
    plt.loglog(ts, state_vector_values[:,1], linewidth=5, label='HII')
    plt.loglog(ts, state_vector_values[:,2], linewidth=5, label='HeI')
    plt.loglog(ts, state_vector_values[:,3], linewidth=5, label='HeII')
    plt.loglog(ts, state_vector_values[:,4], linewidth=5, label='HeIII')
    plt.loglog(ts, state_vector_values[:,5], linewidth=5, label='H2I')
    plt.loglog(ts, state_vector_values[:,6], linewidth=5, label='H2II')
    plt.loglog(ts, state_vector_values[:,7], linewidth=5, label='H-')    
    plt.loglog(ts, nem_values[:], linewidth=5, label='e-')    
    
    plt.xlabel("Time [s]", size = 20)
    plt.ylabel("n", size = 20)
    plt.xticks(size = 20)
    plt.yticks(size = 20)
    plt.legend(loc='upper center', bbox_to_anchor=(1.07, 0.8), shadow=True, ncol=1, 
               fontsize = 20)
    plt.show()

import argparse

parser = argparse.ArgumentParser(description='Evolve.')