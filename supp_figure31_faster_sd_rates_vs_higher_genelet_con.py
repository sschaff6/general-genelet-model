# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 20:37:10 2020

@author: sscha
"""
import numpy as np 
import matplotlib.pyplot as plt
import GeneralGeneletModel as GGM
from matplotlib.ticker import MaxNLocator

''' 
###############################################################################
REPRESSION
###############################################################################
'''
    
kpr = 0.02 # repressor production rates
kpc = 0.02# coactivator production rates
kpi = 0.02
kd_H = 0.0003 # RNaseH degradation rates
kd_A = 0.0003 # RNaseA degradation rates
   #    All nodes
kga = [1e4,1e5,1e6] # activation rates
kgar = [5e3,5e4,5e5] # repression rates
kar = [1e4,1e5,1e6] # activator inhibition rates
kgb = [1e4,1e5,1e6] # free blocking rates
kgbc = [5e3,5e4,5e5] # coactivation rates
kbc = [1e4,1e5,1e6] # blocker inhibition rates 
kgab = [5e3,5e4,5e5] # active blocking rates
kir = [1e4,1e5,1e6] # inducer binding rates

# OG Sam method of network definition
#           G1 G2 G3     
act_vec =  [1, 2]
blk_vec =  [0, 0]
prod_vec = [-2,0]
indc_vec = [0, 0]

''' initializing topology '''
NN1 = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
#NN1.plot_topology(show_rnas=0)

# Define initial conditions
#                  dA1  dA2  dA3
dA_tot = np.array([250, 250]) # total activator added
#                 G1  G2  G3
G_tot = np.array([15, 25]) # total genelet added

# initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
#            G1 G2
G_int_vec = [1, 1]

''' initializing initial conditions '''
NN1.initial_conditions(dA_tot,G_tot,G_int_vec) # default of 0 for all other initial conditions

t_vec1 = np.linspace(0,3,1001)*3600 # seconds

csl = [[1,0,0],[0.66,0,0],[0.33,0,0]]
ls = ['-','-','-']

for n in range(len(kga)):

    my_rc = [kpr,kpc,kpi,kd_H,kd_A,kga[n],kgar[n],kar[n],kgb[n],kgbc[n],kbc[n],kgab[n],kir[n]]
    
    ''' simulating the IFFL'''
    NN1.simulate(t_vec1,1,rate_constants=my_rc)
    
    # pulling out the desired concentrations for plotting
    G1 = NN1.output_concentration['GdA1']
    G2 = NN1.output_concentration['GdA2']
    
    sim_t = NN1.sol.t
    
    fs = 12
    plt.figure(2)
    plt.subplot(2,4,1)
    plt.plot(sim_t/60,G2/G_tot[1],color=csl[n],linewidth=3,linestyle=ls[n])
    plt.ylabel('Fraction ON',fontsize=fs+1,weight='bold')
    plt.xlabel('time (min)',fontsize=fs+1,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,104)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.xaxis.set_tick_params(which='both', size=5, width=1.5, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=1.5, direction='in', right='on')
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(2)
    
    

kpr = 0.02 # repressor production rates
kpc = 0.02# coactivator production rates
kpi = 0.02
kd_H = 0.0003 # RNaseH degradation rates
kd_A = 0.0003 # RNaseA degradation rates
   #    All nodes
kga = 1e4 # activation rates
kgar = 5e3 # repression rates
kar = 1e4 # activator inhibition rates
kgb = 1e4 # free blocking rates
kgbc = 5e3 # coactivation rates
kbc = 1e4 # blocker inhibition rates 
kgab = 5e3 # active blocking rates
kir = 1e4 # inducer binding rates

my_rc = [kpr,kpc,kpi,kd_H,kd_A,kga,kgar,kar,kgb,kgbc,kbc,kgab,kir]

# OG Sam method of network definition
#           G1 G2 G3     
act_vec =  [1, 2]
blk_vec =  [0, 0]
prod_vec = [-2,0]
indc_vec = [0, 0]

G1_tot = [15,50,150]

''' initializing topology '''
NN2 = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
#NN1.plot_topology(show_rnas=0)

# Define initial conditions
#                  dA1  dA2  dA3
dA_tot = np.array([250, 250]) # total activator added

# initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
#            G1 G2
G_int_vec = [1, 1]

t_vec1 = np.linspace(0,3,1001)*3600 # seconds

csl = [[1,0,0],[0.66,0,0],[0.33,0,0]]
ls = ['-','-','-']

for n in range(len(G1_tot)):

    #                 G1  G2  G3
    G_tot = np.array([G1_tot[n], 25]) # total genelet added
    
    ''' initializing initial conditions '''
    NN2.initial_conditions(dA_tot,G_tot,G_int_vec) # default of 0 for all other initial conditions
        
    ''' simulating the IFFL'''
    NN2.simulate(t_vec1,1,rate_constants=my_rc)
    
    # pulling out the desired concentrations for plotting
    G1 = NN2.output_concentration['GdA1']
    G2 = NN2.output_concentration['GdA2']
    
    sim_t = NN1.sol.t
    
    fs = 12
    plt.figure(2)
    plt.subplot(2,4,2)
    plt.plot(sim_t/60,G2/G_tot[1],color=csl[n],linewidth=3,linestyle=ls[n])
    plt.ylabel('Fraction ON',fontsize=fs+1,weight='bold')
    plt.xlabel('time (min)',fontsize=fs+1,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,104)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.xaxis.set_tick_params(which='both', size=5, width=1.5, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=1.5, direction='in', right='on')
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(2)
    
'''
###############################################################################
COACTIVATION
###############################################################################
'''
    
kpr = 0.02 # repressor production rates
kpc = 0.02# coactivator production rates
kpi = 0.02
kd_H = 0.0003 # RNaseH degradation rates
kd_A = 0.0003 # RNaseA degradation rates
   #    All nodes
kga = [1e4,1e5,1e6] # activation rates
kgar = [5e3,5e4,5e5] # repression rates
kar = [1e4,1e5,1e6] # activator inhibition rates
kgb = [1e4,1e5,1e6] # free blocking rates
kgbc = [5e3,5e4,5e5] # coactivation rates
kbc = [1e4,1e5,1e6] # blocker inhibition rates 
kgab = [5e3,5e4,5e5] # active blocking rates
kir = [1e4,1e5,1e6] # inducer binding rates

# OG Sam method of network definition
#           G1 G2 G3     
act_vec =  [1, 2]
blk_vec =  [0, 2]
prod_vec = [2, 0]
indc_vec = [0, 0]

''' initializing topology '''
NN1 = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
#NN1.plot_topology(show_rnas=0)

# Define initial conditions
#                  dA1  dA2  dA3
dA_tot = np.array([250, 250]) # total activator added
#                 G1  G2  G3
G_tot = np.array([15, 25]) # total genelet added

# initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
#            G1 G2
G_int_vec = [1,-1]

''' initializing initial conditions '''
NN1.initial_conditions(dA_tot,G_tot,G_int_vec,dB_added=[0,212.5]) # default of 0 for all other initial conditions

t_vec1 = np.linspace(0,3,1001)*3600 # seconds

csl = [[1,0,0],[0.66,0,0],[0.33,0,0]]
ls = ['-','-','-']

for n in range(len(kga)):

    my_rc = [kpr,kpc,kpi,kd_H,kd_A,kga[n],kgar[n],kar[n],kgb[n],kgbc[n],kbc[n],kgab[n],kir[n]]
    
    ''' simulating the IFFL'''
    NN1.simulate(t_vec1,1,rate_constants=my_rc)
    
    # pulling out the desired concentrations for plotting
    G1 = NN1.output_concentration['GdA1']
    G2 = NN1.output_concentration['GdA2']
    
    sim_t = NN1.sol.t
    
    fs = 12
    plt.figure(2)
    plt.subplot(2,4,5)
    plt.plot(sim_t/60,G2/G_tot[1],color=csl[n],linewidth=3,linestyle=ls[n])
    plt.ylabel('Fraction ON',fontsize=fs+1,weight='bold')
    plt.xlabel('time (min)',fontsize=fs+1,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,104)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.xaxis.set_tick_params(which='both', size=5, width=1.5, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=1.5, direction='in', right='on')
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(2)
    
    

kpr = 0.02 # repressor production rates
kpc = 0.02# coactivator production rates
kpi = 0.02
kd_H = 0.0003 # RNaseH degradation rates
kd_A = 0.0003 # RNaseA degradation rates
   #    All nodes
kga = 1e4 # activation rates
kgar = 5e3 # repression rates
kar = 1e4 # activator inhibition rates
kgb = 1e4 # free blocking rates
kgbc = 5e3 # coactivation rates
kbc = 1e4 # blocker inhibition rates 
kgab = 5e3 # active blocking rates
kir = 1e4 # inducer binding rates

my_rc = [kpr,kpc,kpi,kd_H,kd_A,kga,kgar,kar,kgb,kgbc,kbc,kgab,kir]

# OG Sam method of network definition
#           G1 G2 G3     
act_vec =  [1, 2]
blk_vec =  [0, 2]
prod_vec = [2, 0]
indc_vec = [0, 0]

G1_tot = [15,50,150]

''' initializing topology '''
NN2 = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
#NN1.plot_topology(show_rnas=0)

# Define initial conditions
#                  dA1  dA2  dA3
dA_tot = np.array([250, 250]) # total activator added

# initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
#            G1 G2
G_int_vec = [1,-1]

t_vec1 = np.linspace(0,3,1001)*3600 # seconds

csl = [[1,0,0],[0.66,0,0],[0.33,0,0]]
ls = ['-','-','-']

for n in range(len(G1_tot)):

    #                 G1  G2  G3
    G_tot = np.array([G1_tot[n], 25]) # total genelet added
    
    ''' initializing initial conditions '''
    NN2.initial_conditions(dA_tot,G_tot,G_int_vec,dB_added=[0,212.5]) # default of 0 for all other initial conditions
        
    ''' simulating the IFFL'''
    NN2.simulate(t_vec1,1,rate_constants=my_rc)
    
    # pulling out the desired concentrations for plotting
    G1 = NN2.output_concentration['GdA1']
    G2 = NN2.output_concentration['GdA2']
    
    sim_t = NN1.sol.t
    
    fs = 12
    plt.figure(2)
    plt.subplot(2,4,6)
    plt.plot(sim_t/60,G2/G_tot[1],color=csl[n],linewidth=3,linestyle=ls[n])
    plt.ylabel('Fraction ON',fontsize=fs+1,weight='bold')
    plt.xlabel('time (min)',fontsize=fs+1,weight='bold')
    plt.xticks(fontsize=fs,weight='bold')
    plt.yticks(fontsize=fs,weight='bold')
    ax1 = plt.gca()
    ax1.set_xlim(0,104)
    ax1.set_ylim(-0.1,1.1)
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.xaxis.set_tick_params(which='both', size=5, width=1.5, direction='in', top='on')
    ax1.yaxis.set_tick_params(which='both', size=5, width=1.5, direction='in', right='on')
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(2)