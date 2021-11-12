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
3-node Ring oscillator
###############################################################################
'''
    
kpr = 0.1 # repressor production rates
kpc = 0.02# coactivator production rates
kpi = 0.02
kd_H = 0.008 # RNaseH degradation rates
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
act_vec =  [1, 2, 3]
blk_vec =  [0, 0, 0]
prod_vec = [-2,-3,-1]
indc_vec = [0, 0, 0]

''' initializing topology '''
NN1 = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
NN1.plot_topology(show_rnas=0)

# Define initial conditions
#                  dA1  dA2  dA3
dA_tot = np.array([250, 250, 250]) # total activator added
#                 G1  G2  G3
G_tot = np.array([25,  25,  25]) # total genelet added

# initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
#            G1 G2 G3
G_int_vec = [1,0,0]

''' initializing initial conditions '''
NN1.initial_conditions(dA_tot,G_tot,G_int_vec) # default of 0 for all other initial conditions

t_vec1 = np.linspace(0,24,1001)*3600 # seconds

''' simulating the IFFL'''
NN1.simulate(t_vec1,1,rate_constants=my_rc)

# pulling out the desired concentrations for plotting
G1 = NN1.output_concentration['GdA1']
G2 = NN1.output_concentration['GdA2']
G3 = NN1.output_concentration['GdA3']

sim_t = NN1.sol.t

fs = 12
plt.figure(2)
plt.subplot(3,4,1)
plt.plot(sim_t/3600,G1/G_tot[0],color=[0,0,1],linewidth=2,linestyle='-')
plt.plot(sim_t/3600,G2/G_tot[1],color=[0.6,0,0.6],linewidth=2,linestyle='-')
plt.plot(sim_t/3600,G3/G_tot[2],color=[1,0,0],linewidth=2,linestyle='-')
plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
plt.xlabel('time (hr)',fontsize=fs,weight='bold')
plt.xticks(fontsize=fs,weight='bold')
plt.yticks(fontsize=fs,weight='bold')
ax1 = plt.gca()
ax1.set_xlim(0,8)
ax1.set_ylim(-0.1,1.1)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
ax1.xaxis.set_tick_params(which='both', size=3, width=1, direction='in', top='on')
ax1.yaxis.set_tick_params(which='both', size=3, width=1, direction='in', right='on')

''' 
###############################################################################
5-node Ring oscillator
###############################################################################
'''

kpr = 0.1 # repressor production rates
kpc = 0.02# coactivator production rates
kpi = 0.02
kd_H = 0.012 # RNaseH degradation rates
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
#           G1 G2 G3 G4 G5    
act_vec =  [1, 2, 3, 4, 5]
blk_vec =  [0, 0, 0, 0, 0]
prod_vec = [-2,-3,-4,-5,-1]
indc_vec = [0, 0, 0, 0, 0]

''' initializing topology '''
NN1 = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
NN1.plot_topology(show_rnas=0)

# Define initial conditions
#                  dA1  dA2  dA3  dA4  dA5
dA_tot = np.array([250, 250, 250, 250, 250]) # total activator added
#                 G1  G2  G3  G4  G5
G_tot = np.array([25, 25, 25, 25, 25]) # total genelet added

# initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
#            G1 G2 G3 G4 G5
G_int_vec = [1, 0, 0, 0, 0]

''' initializing initial conditions '''
NN1.initial_conditions(dA_tot,G_tot,G_int_vec) # default of 0 for all other initial conditions

t_vec1 = np.linspace(0,24,1001)*3600 # seconds

''' simulating the IFFL'''
NN1.simulate(t_vec1,1,rate_constants=my_rc)

# pulling out the desired concentrations for plotting
G1 = NN1.output_concentration['GdA1']
G2 = NN1.output_concentration['GdA2']
G3 = NN1.output_concentration['GdA3']
G4 = NN1.output_concentration['GdA4']
G5 = NN1.output_concentration['GdA5']

sim_t = NN1.sol.t

fs = 12
plt.figure(2)
plt.subplot(3,4,2)
plt.plot(sim_t/3600,G1/G_tot[0],color=[0,0,1],linewidth=2,linestyle='-')
plt.plot(sim_t/3600,G2/G_tot[1],color=[0.6,0,0.6],linewidth=2,linestyle='-')
plt.plot(sim_t/3600,G3/G_tot[2],color=[1,0,0],linewidth=2,linestyle='-')
plt.plot(sim_t/3600,G4/G_tot[3],color=[0,0.6,1],linewidth=2,linestyle='-')
plt.plot(sim_t/3600,G5/G_tot[4],color=[0.7,0.7,0],linewidth=2,linestyle='-')
plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
plt.xlabel('time (hr)',fontsize=fs,weight='bold')
plt.xticks(fontsize=fs,weight='bold')
plt.yticks(fontsize=fs,weight='bold')
ax1 = plt.gca()
ax1.set_xlim(0,8)
ax1.set_ylim(-0.1,1.1)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
ax1.xaxis.set_tick_params(which='both', size=3, width=1, direction='in', top='on')
ax1.yaxis.set_tick_params(which='both', size=3, width=1, direction='in', right='on')


''' 
###############################################################################
7-node Ring oscillator
###############################################################################
'''

kpr = 0.1 # repressor production rates
kpc = 0.02# coactivator production rates
kpi = 0.02
kd_H = 0.013 # RNaseH degradation rates
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
#           G1 G2 G3 G4 G5 G6 G7   
act_vec =  [1, 2, 3, 4, 5, 6, 7]
blk_vec =  [0, 0, 0, 0, 0, 0, 0]
prod_vec = [-2,-3,-4,-5,-6,-7,-1]
indc_vec = [0, 0, 0, 0, 0, 0, 0]

''' initializing topology '''
NN1 = GGM.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
NN1.plot_topology(show_rnas=0)

# Define initial conditions
#                  dA1  dA2  dA3  dA4  dA5  dA6  dA7
dA_tot = np.array([250, 250, 250, 250, 250, 250, 250]) # total activator added
#                 G1  G2  G3  G4  G5  G6  G7
G_tot = np.array([25, 25, 25, 25, 25, 25, 25]) # total genelet added

# initial genelet states (0 = OFF, 1 = ON, -1 = BLK)
#            G1 G2 G3 G4 G5 G6 G7
G_int_vec = [1, 0, 0, 0, 0, 0, 0]

''' initializing initial conditions '''
NN1.initial_conditions(dA_tot,G_tot,G_int_vec) # default of 0 for all other initial conditions

t_vec1 = np.linspace(0,24,1001)*3600 # seconds

''' simulating the IFFL'''
NN1.simulate(t_vec1,1,rate_constants=my_rc)

# pulling out the desired concentrations for plotting
G1 = NN1.output_concentration['GdA1']
G2 = NN1.output_concentration['GdA2']
G3 = NN1.output_concentration['GdA3']
G4 = NN1.output_concentration['GdA4']
G5 = NN1.output_concentration['GdA5']
G6 = NN1.output_concentration['GdA6']
G7 = NN1.output_concentration['GdA7']

sim_t = NN1.sol.t

fs = 12
plt.figure(2)
plt.subplot(3,4,3)
plt.plot(sim_t/3600,G1/G_tot[0],color=[0,0,1],linewidth=2,linestyle='-')
plt.plot(sim_t/3600,G2/G_tot[1],color=[0.6,0,0.6],linewidth=2,linestyle='-')
plt.plot(sim_t/3600,G3/G_tot[2],color=[1,0,0],linewidth=2,linestyle='-')
plt.plot(sim_t/3600,G4/G_tot[3],color=[0,0.6,1],linewidth=2,linestyle='-')
plt.plot(sim_t/3600,G5/G_tot[4],color=[0.7,0.7,0],linewidth=2,linestyle='-')
plt.plot(sim_t/3600,G6/G_tot[5],color=[0,0.8,0],linewidth=2,linestyle='-')
plt.plot(sim_t/3600,G7/G_tot[6],color=[0,0.5,0],linewidth=2,linestyle='-')
plt.ylabel('Fraction ON',fontsize=fs,weight='bold')
plt.xlabel('time (hr)',fontsize=fs,weight='bold')
plt.xticks(fontsize=fs,weight='bold')
plt.yticks(fontsize=fs,weight='bold')
ax1 = plt.gca()
ax1.set_xlim(0,8)
ax1.set_ylim(-0.1,1.1)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
ax1.xaxis.set_tick_params(which='both', size=3, width=1, direction='in', top='on')
ax1.yaxis.set_tick_params(which='both', size=3, width=1, direction='in', right='on')