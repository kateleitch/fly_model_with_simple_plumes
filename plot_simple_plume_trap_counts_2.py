from __future__ import print_function
import numpy as np
import os, sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import random
import json
import seaborn as sns
sys.path.append('/home/kate/Documents/coyote_lake_field_data')
import trap_layout as trap_layout

degrees_to_rotate = 0

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

def draw_wind_and_traject_vectors(ax_handle, wind_vector,perp_vector, par_vector,heading_unit_vector, trajectory_vector, draw_legend, graph_number):
    l = 100
    if draw_legend:
        l = 0
    ax_handle.plot(np.linspace(0,wind_vector[0],l),np.linspace(0,wind_vector[1],l), 'deepskyblue', linewidth = 2, label = 'wind')
    ax_handle.scatter(wind_vector[0], wind_vector[1], s = 15, color = 'deepskyblue')
    ax_handle.plot(np.linspace(0,perp_vector[0],l),np.linspace(0,perp_vector[1],l), 'darkturquoise', linewidth = 2,label = 'wind perp')
    ax_handle.scatter(perp_vector[0], perp_vector[1], s = 15, color = 'darkturquoise')
    ax_handle.plot(np.linspace(0,par_vector[0],l),np.linspace(0,par_vector[1],l), 'cadetblue',linewidth = 2,label = 'wind par')
    ax_handle.scatter(par_vector[0], par_vector[1], s = 15, color = 'cadetblue')

    ax_handle.plot(np.linspace(0,trajectory_vector[0],l),np.linspace(0,trajectory_vector[1],l), 'k', linewidth =2,label = 'trajectory')
    ax_handle.scatter(trajectory_vector[0], trajectory_vector[1], s = 15, color = 'k')
    ax_handle.plot(np.linspace(0,3*heading_unit_vector[0],l),np.linspace(0,3*heading_unit_vector[1],l), '--k',label = 'heading')
    if draw_legend:
        ax_handle.legend(fontsize = 10, loc =2)
    else:
        y_position = -1* np.sign(trajectory_vector[1])
        ax_handle.text(-1,y_position,('%.2f m/s') %(np.dot(par_vector,heading_unit_vector)/np.linalg.norm(heading_unit_vector)),color = 'cadetblue')
        ax_handle.text(-1.9,1.6,str(graph_number))
    ax_handle.set_ylim(-2.5,2.5)
    ax_handle.set_xlim(-2.5,2.5)

def color_and_size(trap_catch_dictionary):
    return_dict = {}
    dotcolor = (0.75, 0.75, 0.75, 1)
    dotcolor = (0.5,0.5,0.5,1)
    for trap_name in trap_catch_dictionary:
        #dotsize = 5* trap_catch_dictionary[trap_name]
        offset = 0
        if trap_catch_dictionary[trap_name]>0:
            offset = 25
        dotsize = 10 * trap_catch_dictionary[trap_name] + offset
        return_dict[trap_name]={'dotcolor':dotcolor, 'dotsize':dotsize}
    return return_dict

with open('./sideslip_and_simple_plume_trap_counts.json') as f:
    dictionary = json.load(f)


fig = plt.figure(figsize=(15,15))
plt.title('plume space constant 100 meters')
row_number = 5
gs = gridspec.GridSpec(nrows=row_number, ncols=5)

subplot_count = 0
for windspeed in dictionary:
    for wind_dir in dictionary[windspeed]:

        trapping_dictionary = {}
        try:
            trapping_info = dictionary[windspeed][wind_dir]['trapping_info'][0]
        except:
            continue
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
        r, c = divmod(int(float(windspeed)*10), row_number)
        print (r,c)
        if int(float(windspeed)*10) < row_number*5:
            ax = plt.subplot(gs[r,c], polar=True)
            plt.xticks([0,np.pi/2, np.pi, 3*np.pi/2], ['E', 'N', 'W', 'S'],fontsize = 6)
            ax.set_yticks([])
            ax.grid(False)
            ax.axis("off")
            ax.set_ylim([0,1800])
            ax.set_xlim([0,2*np.pi])
            ax.text(0,0,str(windspeed))
            for key in trapping_dictionary:
                ax.scatter(float(key), 1000, color = 'gray',  edgecolor = 'white', s = trapping_dictionary[key]*25)
    subplot_count +=1

#plt.savefig('./modeled_groundspeed_vs_windspeed.svg', transparent = True)
plt.savefig('./1.png', transparent = True)
plt.show()
