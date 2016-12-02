#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 19:05:34 2016

@author: kousuke
"""

import numpy as np
import MR

major_list = ["20161128-162543-E431_low_100mv_100-Major",
              "20161128-164452-E431_low_100mv_100-Major",
              "20161128-171000-E431_low_100mv_100-Major",
              "20161128-210548-E431_low_100mv_100-Major",
              "20161128-223156-E431_low_100mv_100-Major",
              "20161128-225336-E431_low_100mv_100-Major",
              "20161128-231517-E431_low_100mv_100-Major",
              "20161128-233702-E431_low_100mv_100-Major",
              "20161128-235846-E431_low_100mv_100-Major"]
              
minor_list = ["20161128-234523-E431_low_100mv_100-minor",
              "20161128-232339-E431_low_100mv_100-minor",
              "20161128-230155-E431_low_100mv_100-minor",
              "20161128-224016-E431_low_100mv_100-minor",
             ]
              
num = len(major_list)
MR_y_list =[]
for i in range(num):
    MR_y_list.append(np.loadtxt(major_list[i] + ".csv",delimiter=",", skiprows=1,usecols=(1,)))
MR_y_array = np.zeros(len(MR_y_list[0]),np.float64)
for j in range(num):
    MR_y_array = MR_y_list[j] / num + MR_y_array
MR_x_array = np.loadtxt(major_list[0] + ".csv",delimiter=",", skiprows=1,usecols=(0,))

MR_xy = [[0 for n1 in range(2)] for n2 in range(len(MR_x_array))]

for l in range(len(MR_x_array)):
    MR_xy[l][0] = MR_x_array[l]
    MR_xy[l][1] = MR_y_array[l]

out_name ="ave_major"
np.savetxt(out_name + ".csv",MR_xy,delimiter=",")

num_m = len(minor_list)
MR_y_list_m =[]
for i1 in range(num_m):
    MR_y_list_m.append(np.loadtxt(minor_list[i1] + ".csv",delimiter=",", skiprows=1,usecols=(1,)))
MR_y_array_m = np.zeros(len(MR_y_list_m[0]),np.float64)
for j1 in range(num_m):
    MR_y_array_m = MR_y_list_m[j1] / num_m + MR_y_array_m
MR_x_array_m = np.loadtxt(minor_list[0] + ".csv",delimiter=",", skiprows=1,usecols=(0,))

MR_xy_m = [[0 for n1 in range(2)] for n2 in range(len(MR_x_array_m))]

for l1 in range(len(MR_x_array_m)):
    MR_xy_m[l1][0] = MR_x_array_m[l1]
    MR_xy_m[l1][1] = MR_y_array_m[l1]

out_name_m ="ave_minor"
np.savetxt(out_name_m + ".csv",MR_xy_m,delimiter=",")


MR_C = MR.MR(out_name,100,out_name_m)
MR_C.MR_curve_plot()
