from __future__ import print_function
import numpy as np
import os, sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import random
import json
import seaborn as sns


def giraldo_circvar(alpha,axis=None):
#### vector strength = 1 - circvar
    N = len(alpha)
    R = np.sqrt(np.sum(np.sin(alpha),axis)**2 + np.sum(np.cos(alpha),axis)**2)/N
    circvar = 1-R
    return circvar

def adjust_spines(ax_handle, spines):
    for loc, spine in ax_handle.spines.items():
        if loc in spines:
            spine.set_position(('outward', 10))  # outward by 10 points
            #spine.set_smart_bounds(True)
        else:
            spine.set_color('none')  # don't draw spine
    # turn off ticks where there is no spine
    if 'left' in spines:
        ax_handle.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax_handle.yaxis.set_ticks([])
    if 'bottom' in spines:
        ax_handle.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax_handle.xaxis.set_ticks([])

fig = plt.figure(figsize=(8,8))
ax = plt.subplot(1,1,1)
with open('./sideslip_and_simple_plume_trap_counts.json') as f:
    dictionary = json.load(f)

all_windspeeds_along_trajectories=[]
all_groundspeeds_along_trajectories=[]


for windspeed in dictionary:
    for wind_dir in dictionary[windspeed]:
        trapping_dictionary = {}
        trapping_info = dictionary[windspeed][wind_dir]['trapping_info']
        for trap_record in trapping_info:
            if trap_record == 'not trapped':
                continue
            else:
                trap_angle = trap_record["trap angle"]
                if trap_angle in trapping_dictionary:
                    trapping_dictionary[trap_angle] +=1
                else:
                    trapping_dictionary[trap_angle] =1
        trap_angle_hist_list = []
        for trap_angle in trapping_dictionary:
            trap_angle_hist_list.extend([float(trap_angle)]*(trapping_dictionary[trap_angle]))
        circvar = giraldo_circvar(trap_angle_hist_list)
        ax.scatter(windspeed, circvar, color = 'k')
ax.set_ylim([0,1])
ax.set_xlim([-0.1,3.0])
ax.set_ylabel('trap count circular variance')
ax.set_xlabel('simulated wind speed, m/s')
adjust_spines(ax, spines = ['bottom', 'left'])
plt.savefig('./circvar_100_meters.png', transparent = True)
plt.show()
