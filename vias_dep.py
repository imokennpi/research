#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 21:58:55 2016

@author: kousuke
"""

inp_name = ["20161128-203422-E431_low_-700mv_100-Major",
            "20161128-193852-E431_low_-500mv_100-Major",
            "20161128-192506-E431_low_-300mv_100-Major",
            "20161128-191106-E431_low_-100mv_100-Major",
            "20161128-200618-E431_low_-50mv_100-Major",
            "20161128-195234-E431_low_50mv_100-Major",
            "20161128-162543-E431_low_100mv_100-Major",
            "20161128-153811-E431_low_300mv_110-Major",
            "20161128-151223-E431_low_500mv_110-Major",
            "20161128-202043-E431_low_700mv_100-Major",
            ]
vias_vol = [-700,
            -500,
            -300,
            -100,
            -50,
            50,
            100,
            300,
            500,
            700]

class_list = []
TMR_list = []
TMR_min = []
TMR_max = []
import MR
for i in range(len(inp_name)):
    #class_list.append(MR.MR(inp_name[i],vias_vol[i]))
    #class_list[i].MR_curve_plot()
    TMR_list.append(MR.MR(inp_name[i],vias_vol[i]).TMR)
    TMR_min.append(MR.MR(inp_name[i],vias_vol[i]).TMR_min)
    TMR_max.append(MR.MR(inp_name[i],vias_vol[i]).TMR_max)
    #TMR_list.append(MR.MR(inp_name[i],vias_vol[i]).TMR)


from matplotlib import pyplot as plt
plt.figure(figsize=(8, 6), dpi=80)
plt.xlim([-900,900])
plt.xlabel(r'Vias Voltage (V)')
plt.ylabel(r'MR ratio (%)')
plt.axvline(x=0,linestyle='dashed')
plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
plt.errorbar(vias_vol,TMR_list,yerr=[TMR_min,TMR_max],fmt='ro',ecolor='b')
plt.plot(vias_vol,TMR_list,'o') 
plt.savefig("vias_dep.png")
plt.show()
plt.close()
