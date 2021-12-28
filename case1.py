#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 13:09:17 2021

@author: jakravit
"""
from __future__ import division
import os
import pandas as pd
import numpy as np
from build.lognorm import lognorm_params,lognorm_random
import matplotlib.pyplot as plt
import pickle
import shortuuid

# Case
sname_title = 'Case1'

# chl
sigma,scale = lognorm_params(.07,.8)
chlaData = lognorm_random(sigma, scale, 10000)

# lambda
l = np.arange(400, 902.5, 2.5)  

# run names
snames = []
runlist = {}

# indices
idx440 = np.where(l==440)
idx700 = int(np.where(l==700)[0])

# how many runs to build
runs = 10

for k in range(runs):
    
    iops = {'Phyto': {},
            'CDOM': {},
            'Det': {},
            'Min': {}}
    
    
#################### PHYTOPLANKTON ##############################################    
    
    # phyto data
    path = '/Users/jakravit/pyProjects/EAPbuild/build/classes.p'
    with open(path,'rb') as fp:
        phytodata = pickle.load(fp)

    # chl concentration
    chl = round(np.random.choice(chlaData), 3); print ('\nTot chl: {}'.format(chl))
    
    # phyto class fractions
    phyto_class_frxn = {   'Haptophytes': .3, 
                           'Diatoms': .2,
                           'Dinoflagellates': .1,
                           'Cryptophytes': .2,
                           'Green_algae': .15,
                           'Cyano_blue': .05
                           }
    
    # phyto component
    classIOPs = {}
    for c, f in phyto_class_frxn.items():
        print ('\n' + c)
        classIOPs[c] = {}
        class_chl = f * chl
        sps = []
        for i, sp in enumerate(phytodata[c].keys()):
            print (sp)
            sps.append(sp)
            info = phytodata[c][sp]
            angles = info['VSF_angles']
            rem_list = ['astar','bstar','bbstar','cstar','lambda','Qc',
                        'Sigma_c','Qb','Sigma_b','Qa','Sigma_a','Qbb','Sigma_bb',]
            classIOPs[c][sp] = info
            if  i == 0:
                frxn = .7
            else:
                frxn = .3
            classIOPs[c][sp]['sp_frxn'] = frxn
            idx = np.random.choice(len(info['astar']), 1) 
            deff = info['astar'].iloc[idx,:].index.values
            sp_chl = class_chl * frxn
            classIOPs[c][sp]['a'] = sp_chl * info['astar'].iloc[idx,:].values[0]
            classIOPs[c][sp]['b'] = sp_chl * info['bstar'].iloc[idx,:].values[0]
            classIOPs[c][sp]['c'] = sp_chl * info['cstar'].iloc[idx,:].values[0]
            classIOPs[c][sp]['bb'] = sp_chl * info['bbstar'].iloc[idx,:].values[0]
            classIOPs[c][sp]['VSF'] = sp_chl * info['VSF'][idx,:,:]
            classIOPs[c][sp]['Deff'] = deff
            classIOPs[c][sp]['sp_chl_conc'] = sp_chl
            # remove unnecessary info
            [classIOPs[c][sp].pop(x) for x in rem_list]
        
        # phyto class tot IOPs
        for p in ['a','b','c','bb','VSF']:
            class_tot_iop = 0
            for sp in sps:
                class_tot_iop = class_tot_iop + classIOPs[c][sp][p]
            classIOPs[c]['{}_tot'.format(p)] = class_tot_iop
            classIOPs[c]['class_frxn'] = f
                  
        # phyto class chl contribution
        class_tot_chl = 0
        for sp in sps:
            class_tot_chl = class_tot_chl + classIOPs[c][sp]['sp_chl_conc']            
        classIOPs[c]['class_chl'] = class_tot_chl 

    # phyto component total iops
    for p in ['a','b','c','bb','VSF']:
        comp_tot_iop = 0
        for c in phyto_class_frxn.keys():
            comp_tot_iop = comp_tot_iop + classIOPs[c]['{}_tot'.format(p)]
        classIOPs['{}_tot'.format(p)] = comp_tot_iop
    
    # phyto component total chl
    comp_tot_chl = 0
    for c in phyto_class_frxn.keys():
        comp_tot_chl = comp_tot_chl + classIOPs[c]['class_chl']
    classIOPs['TotChl'] = comp_tot_chl
            
    classIOPs['lambda'] = l
    classIOPs['VSF_angles'] = angles
    print ('\nModel tot aphy440: {}'.format(classIOPs['a_tot'][idx440][0]))
    print ('Bricaud 2004 tot aphy440: {}'.format(.06*chl**.728))
    print ('Matthews 2013 inland tot aphy440: {}'.format(.031*chl**.89))
    
    # phyto component
    iops['Phyto'] = classIOPs
            
            
########################## MINERALS ###############################################
          
    minpath = '/Users/jakravit/pyProjects/EAPbuild/build/minerals.p'
    with open(minpath, 'rb') as fp:
        datamin = pickle.load(fp)  
    rem_list = ['astar','bstar','bbstar','cstar','Qc', 'VSF_angles','VSF_theta',
                'Sigma_c','Qb','Sigma_b','Qa','Sigma_a','Qbb','Sigma_bb',]
                
    
    aphy440 = iops['Phyto']['a_tot'][idx440]
    
    # from lee 2002 for case1 waters
    r1 = np.arange(0,1.05,.05)
    r2 = np.arange(0.05,.1,.005)
    p1 = .1 + (0.5 * np.random.choice(r1) * aphy440) / (np.random.choice(r2) + aphy440) 
    sf = np.random.choice(np.arange(0.3,.8,.1))
    amin440 = sf * p1 * aphy440
    print ('\nTot amin440: {}'.format(amin440[0]))
    
    import random
    run_mins = random.sample(datamin.keys(), 2)
    min_frxn = {run_mins[0]: .7,
                run_mins[1]: .3}
    
    # mineral component
    minIOPs = {}
    sps = []
    for c, f in min_frxn.items():
        print ('\n' + c)
        minIOPs[c] = {}
        info = datamin[c]
        idx1 = np.random.choice(list(info.keys()))
        print (idx1)
        sps.append(c)
        minIOPs[c] = info[idx1] 
        cminl = (amin440 * f) / info[idx1]['astar'][0][idx440]
        astar = info[idx1]['astar'][0]
        minIOPs[c]['a'] = cminl * astar
        minIOPs[c]['b'] = cminl * info[idx1]['bstar'][0]
        minIOPs[c]['bb'] = cminl * info[idx1]['bbstar'][0]
        minIOPs[c]['VSF'] = cminl * info[idx1]['VSF'][0]
        # calculate sp slopes
        minIOPs[c]['class_slope'] = np.polyfit(l[:idx700], np.log(astar[:idx700]),1)[0]
        # remove unnecessary info
        [minIOPs[c].pop(x) for x in rem_list]
        minIOPs[c]['class_frxn'] = f
        minIOPs[c]['class_conc'] = cminl[0]
        
    # mineral component total iops
    for p in ['a','b','bb','VSF']:
        comp_tot_iop = 0 # iops
        for c in min_frxn.keys():
            comp_tot_iop = comp_tot_iop + minIOPs[c]['{}'.format(p)]
            minIOPs['{}_tot'.format(p)] = comp_tot_iop
            #comp_tot_C = comp_tot_C + minIOPs[c]['class_conc']
    
    # mineral component total conc
    comp_tot_minl = 0
    for c in min_frxn.keys():
        comp_tot_minl = comp_tot_minl + minIOPs[c]['class_conc']
    minIOPs['tot_conc'] = comp_tot_minl
    minIOPs['tot_slope'] = np.polyfit(l[:idx700], np.log(astar[:idx700]),1)[0]    
    minIOPs['lambda'] = l
    minIOPs['VSF_angles'] = angles  
    print ('\nTot min conc: {}'.format(minIOPs['tot_conc'])) 
    
    # mineral component
    iops['Min'] = minIOPs 


#################### DETRITUS ##################################################

    detpath = '/Users/jakravit/pyProjects/EAPbuild/build/det.p'
    adet440 = (1-sf) * p1 * aphy440
    print ('\nTot adet440: {}'.format(adet440[0]))
    detIOPs = {}
    with open(detpath, 'rb') as fp:
        datadet = pickle.load(fp) 
    idx1 = np.random.choice(list(datadet.keys()))
    #info = idx1.split('_')
    print ('\n'+idx1)
    astar = datadet[idx1]['astar'][0]
    cdet = adet440 / astar[idx440]
    detIOPs = datadet[idx1]
    detIOPs['a_tot'] = cdet * astar
    detIOPs['b_tot'] = cdet * datadet[idx1]['bstar'][0]
    detIOPs['bb_tot'] = cdet * datadet[idx1]['bbstar'][0]
    detIOPs['VSF_tot'] = cdet * datadet[idx1]['VSF'][0]
    detIOPs['tot_conc'] = cdet[0]
    detIOPs['tot_slope'] = np.polyfit(l[:idx700], np.log(astar[:idx700]),1)[0]       
    detIOPs['lambda'] = l
    detIOPs['VSF_angles'] = angles 
    detIOPs['tot_conc'] = cdet[0]
    # remove unnecessary info
    [detIOPs.pop(x) for x in rem_list]
    
    # det component
    print ('\nTot det conc: {}'.format(detIOPs['tot_conc']))
    iops['Det'] = detIOPs

##################### CDOM ######################################################

    l = np.arange(240,900, 2.5)
    for k in range(10):
            
        domIOPs = {}
        # random slope (240-900) from normal dist.
        slopes = np.random.normal(.02,.01,5000)
        slopes = slopes[slopes > 0]
        # plt.hist(slopes,bins=100)
        slope = np.random.choice(slopes)
        # slope = .03
        print ('\nSlope: {}'.format(slope))
        r1 = np.arange(0, 1.05, .05)
        p2 = 0.3 + (5.7 * np.random.choice(r1,1) * aphy440) / (0.02 + aphy440)
        ag440 = p2 * aphy440
        ag1 = ag440 * np.exp(-slope * (l-440))
        ag350 = ag1[np.where(l==350)]
        ag265 = ag1[np.where(l==265)]
        print ('Tot ag440: {}'.format(ag440[0]))
        print ('Tot ag350: {}'.format(ag350[0]))
        print ('Tot ag265: {}'.format(ag265[0]))
        
        #%
        cdomIOPs = {}
        
        # mu =np.random.normal(320,5,1000)
        # plt.hist(mu,bins=100)
        
        comp_data = {1: {'mu': np.random.normal(269,3,1000),
                         'std': np.random.normal(20,3,1000)}, 
                     2: {'mu': np.random.normal(299,5,1000),
                         'std': np.random.normal(20,5,1000)},
                     3: {'mu': np.random.normal(320,5,1000),
                         'std': np.random.normal(30,3,1000)},
                     4: {'mu': np.random.normal(345,5,1000),
                         'std': np.random.normal(55,5,1000)},
                     5: {'mu': np.random.normal(375,5,1000),
                         'std': np.random.normal(70,5,1000)},
                     6: {'mu': np.random.normal(407,5,1000),
                         'std': np.random.normal(90,10,1000)},
                     7: {'mu': np.random.normal(440,10,1000),
                         'std': np.random.normal(115,10,1000)},
                     8: {'mu': np.random.normal(500,10,1000),
                         'std': np.random.normal(130,10,1000)},
                     }
        gc = np.random.choice(range(0,9),1)[0]
        print ('gaus comps: {}'.format(gc))
        
        
        qx = np.random.choice(np.linspace(.1,.25,20))
        cx = ag265 * qx
        comps = {'phis':[],
                 'mus':[],
                 'stds':[],
                 'gs':[]}
        
        for c in range(gc):
            c += 1
            if c == 1:
                comps['phis'].append(cx)
                comps['mus'].append(np.random.choice(comp_data[c]['mu']))
                comps['stds'].append(np.random.choice(comp_data[c]['std']))
            else:
                cx = cx*.8
                comps['phis'].append(cx)
                comps['mus'].append(np.random.choice(comp_data[c]['mu']))
                comps['stds'].append(np.random.choice(comp_data[c]['std']))
        
        
        def ag_gauss(ag0,s,l,l0,comps):
            ag = ag0 * np.exp(-s * (l-l0))
            gcomps = 0
            for i in range(len(comps['phis'])):
                # print (i)
                gs = comps['phis'][i] * np.exp(-(l-comps['mus'][i])**2 / (2*comps['stds'][i]**2))
                comps['gs'].append(gs)
                gcomps += gs
            return ag, gcomps, comps
        
        ag, gcomps, comps = ag_gauss(ag350, slope, l, 350,comps)
        agtot = ag + gcomps
                
        ## plot ag and components
        fig, ax = plt.subplots()
        ax.plot(l,ag, label='ag_exp')
        if type(gcomps) == int:
            pass
        else:
            ax.plot(l,gcomps, label='comps_tot')
        ax.plot(l,agtot, label= 'ag_tot_abs')
        for i,k in enumerate(comps['gs']):
            ax.plot(l, k, label='comp {}'.format(i))
        ax.set_xlim(240,600)
        ax.set_ylabel('ag absorption (m$^{-1})$')
        ax.set_xlabel('lambda (nm)')
        ax.legend()
        
        #%
        #ag = ag440 * np.exp(-sy * (l-440))
        import scipy as sp
        import scipy.optimize
        
        def exp(l, a0, s):
            return a0 * np.exp(-s * (l-350))
        
        def fit_exp(l, y):
            opt_parms, parm_cov = sp.optimize.curve_fit(exp, l, y, p0=[.001, .001])
            a0, s = opt_parms
            return a0, s 
        
        # slopes
        
        # 275-295
        li0 = np.where(l == 275)[0][0]
        li1 = np.where(l == 295)[0][0]
        a, s275_295 = fit_exp(l[li0:li1],agtot[li0:li1])
        
        # S240-700
        li0 = np.where(l == 240)[0][0]
        li1 = np.where(l == 700)[0][0]
        a, s240_700 = fit_exp(l[li0:li1],agtot[li0:li1])
        
        # 300-700
        li0 = np.where(l == 300)[0][0]
        li1 = np.where(l == 700)[0][0]
        a, s300_700 = fit_exp(l[li0:li1],agtot[li0:li1])
        
        # 350-400
        li0 = np.where(l == 350)[0][0]
        li1 = np.where(l == 400)[0][0]
        a, s350_400 = fit_exp(l[li0:li1],agtot[li0:li1])
        #s350_400_poly = np.polyfit(l[li0:li1], np.log(agtot[li0:li1]),1)[0]
        
        # 350-550
        li0 = np.where(l == 350)[0][0]
        li1 = np.where(l == 550)[0][0]
        a, s350_550 = fit_exp(l[li0:li1],agtot[li0:li1])
        
        # 400-450
        li0 = np.where(l == 400)[0][0]
        li1 = np.where(l == 450)[0][0]
        a, s400_450 = fit_exp(l[li0:li1],agtot[li0:li1])
        s400_450_poly = np.polyfit(l[li0:li1], np.log(agtot[li0:li1]),1)[0]
        
        # 400-700
        li0 = np.where(l == 400)[0][0]
        li1 = np.where(l == 700)[0][0]
        a, s400_700 = fit_exp(l[li0:li1],agtot[li0:li1])
        
        
        # slope ratio (Helms 2008)
        sr = s275_295 / s350_400
        
        if sr > .7:
            break
    
    print ('\nSlopes\n240-700: {}\n300-700: {}\n350-550: {}\n350-400: {}\n275-295: {}\n400-450: {}\n400-700: {}\nSR: {}'.format(
            s240_700, s300_700, s350_550, s350_400, s275_295, s400_450, s400_700, sr))
    
    cdomIOPs['a_tot'] = agtot
    cdomIOPs['gaus_comps'] = gcomps
    cdomIOPs['num_gcomps'] = gc
    cdomIOPs['S275_295'] = s275_295
    cdomIOPs['S240_700'] = s240_700
    cdomIOPs['S300_700'] = s300_700
    cdomIOPs['S350_400'] = s350_400
    cdomIOPs['S350_550'] = s350_550
    cdomIOPs['S400_450'] = s400_450
    cdomIOPs['S400_700'] = s400_700
    cdomIOPs['slope_ratio'] = sr
    
    # cdom component
    iops['CDOM'] = cdomIOPs

############### SAVE ###########################################################
    
    uid = shortuuid.ShortUUID().random(length=10)
    uid = '{:.2f}_{:.2f}_{:.2f}_{:.2f}_{}'.format(chl,comp_tot_minl,cdet[0],ag440[0],uid)
    snames.append(uid)
    runlist[uid] = iops

with open('/Users/jakravit/pyProjects/EAPbuild/build/runlists/{}.p'.format(sname_title), 'wb') as f:
    pickle.dump(iops, f)      

            
            
                       
  
            
