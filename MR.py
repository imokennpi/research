#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 27 16:29:48 2016

@author: kousuke
"""
import numpy as np

class MR:
    ######################
    def __init__(self,major_filename,vias_vol,minor_filename = ""):
        self.major_filename = major_filename 
        self.vias_vol = vias_vol
        self.minor_filename = minor_filename 
        self.input_data_major = np.loadtxt(major_filename + ".csv",delimiter = ",",skiprows=1)
        i = 0
        row_len = self.input_data_major.size /self.input_data_major[0].size
        self.output_mag_field_go=[]
        self.output_resistance_go =[]
        self.output_mag_field_back=[]
        self.output_resistance_back =[]
        go_or_not = True
        self.vias_vol = self.vias_vol * 0.001
        #get datas
        while i < row_len:
            if(i > row_len / 2 ):
                go_or_not = False
            if(go_or_not):
                self.output_mag_field_go.append(self.fix_mag_field(self.input_data_major[i][0],go_or_not))
                self.output_resistance_go.append((self.vias_vol ) / self.input_data_major[i][1])
            else:
                self.output_mag_field_back.append(self.fix_mag_field(self.input_data_major[i][0],go_or_not))
                self.output_resistance_back.append((self.vias_vol ) / self.input_data_major[i][1])
            i += 1   

        self.Rp = (sum(self.output_resistance_go[0:10]) / 10 + \
                       sum(self.output_resistance_go[len(self.output_resistance_go)-10:len(self.output_resistance_go)]) / 10 +\
                           sum(self.output_resistance_back[0:10]) / 10 +\
                               sum(self.output_resistance_back[len(self.output_resistance_back)-10:len(self.output_resistance_back)]) / 10) / 4
        self.Rap = (max(self.output_resistance_go) + max(self.output_resistance_back)) / 2
        
        # to get errror bar
        Rp_min = min([min(self.output_resistance_go[0:10]),\
                          min(self.output_resistance_go[len(self.output_resistance_go)-10:len(self.output_resistance_go)]),\
                              min(self.output_resistance_back[0:10]),\
                                  min(self.output_resistance_back[len(self.output_resistance_back)-10:len(self.output_resistance_back)])])
        Rp_max = max([max(self.output_resistance_go[0:10]),\
                          max(self.output_resistance_go[len(self.output_resistance_go)-10:len(self.output_resistance_go)]),\
                              max(self.output_resistance_back[0:10]),\
                                  max(self.output_resistance_back[len(self.output_resistance_back)-10:len(self.output_resistance_back)])])  
        Rap_min = min([max(self.output_resistance_go),max(self.output_resistance_back)])
        Rap_max = max([max(self.output_resistance_go),max(self.output_resistance_back)])
        self.TMR_min = (Rap_min / Rp_max - 1) * 100
        self.TMR_max = (Rap_max / Rp_min - 1) * 100
        #
        
        self.output_ratio_go = (self.output_resistance_go / self.Rp - 1) * 100
        self.output_ratio_back = (self.output_resistance_back / self.Rp - 1) * 100
        # minor loop 
        if(minor_filename == ""):
            self.minor_exist_or_not = False
        else:
            self.minor_exist_or_not = True
            self.input_data_minor = np.loadtxt(minor_filename + ".csv",delimiter = ",",skiprows=1)
            i = 0
            row_len = self.input_data_minor.size /self.input_data_minor[0].size
            self.output_mag_field_minor=[]
            self.output_resistance_minor =[]
            while i < row_len:    
                self.output_mag_field_minor.append(self.fix_mag_field(self.input_data_minor[i][0],True))
                self.output_resistance_minor.append((self.vias_vol ) / self.input_data_minor[i][1])
                i += 1
            self.output_ratio_minor = (self.output_resistance_minor / self.Rp - 1) * 100        
        print(self.calc_TMR_ratio())
        
        
    def fix_mag_field(self,x,go_or_not):
        if(go_or_not):
            return (0.0011 * x * 1000 + 0.0627) * 0.1
        else:
            return (0.0011 * x * 1000 - 0.0546) * 0.1

    def calc_TMR_ratio(self):
        self.TMR = (self.Rap / self.Rp - 1) * 100
        return (self.Rap / self.Rp - 1) * 100

    def MR_curve_plot(self):

        #self.TMR = self.calc_TMR_ratio(output_resistance_go,output_resistance_back)
        from matplotlib import pyplot as plt
        plt.rcParams["font.size"] = 18
        
        plt.figure(figsize=(8, 6), dpi=80)
        plt.xlim([-0.3,0.3])
        plt.xlabel(r'$\mu_0 H$(T)')
        plt.ylabel(r'Resistance($\Omega$)')
        plt.axvline(x=0,linestyle='dashed')
        plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
        plt.plot(self.output_mag_field_go,self.output_resistance_go) 
        plt.plot(self.output_mag_field_back,self.output_resistance_back) 
        plt.savefig("MR_curve_wo_minor_" + self.major_filename + "_" + self.minor_filename + ".png")
        plt.show()
        plt.close()
        
        plt.figure(figsize=(8, 6), dpi=80)
        plt.xlim([-0.3,0.3])
        plt.xlabel(r'$\mu_0 H$(T)')
        plt.ylabel(r'MR ratio (%)')
        plt.axvline(x=0,linestyle='dashed')
        plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
        plt.plot(self.output_mag_field_go,self.output_ratio_go) 
        plt.plot(self.output_mag_field_back,self.output_ratio_back) 
        plt.savefig("MR_ratio_curve_wo_minor_" + self.major_filename + "_" + self.minor_filename + ".png")
        plt.show()
        plt.close()
        
        if(self.minor_exist_or_not):            
            plt.figure(figsize=(8, 6), dpi=80)
            plt.xlim([-0.3,0.3])
            plt.xlabel(r'$\mu_0 H$(T)')
            plt.ylabel(r'Resistance($\Omega$)')
            plt.axvline(x=0,linestyle='dashed')
            plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
            plt.plot(self.output_mag_field_go,self.output_resistance_go) 
            plt.plot(self.output_mag_field_back,self.output_resistance_back)
            plt.plot(self.output_mag_field_minor,self.output_resistance_minor)
            plt.savefig("MR_curve_w_minor_" + self.major_filename + "_" + self.minor_filename + ".png")
            plt.show()
            plt.close()
            
            plt.figure(figsize=(8, 6), dpi=80)
            plt.xlim([-0.3,0.3])
            plt.xlabel(r'$\mu_0 H$(T)')
            plt.ylabel(r'MR ratio (%)')
            plt.axvline(x=0,linestyle='dashed')
            plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
            plt.plot(self.output_mag_field_go,self.output_ratio_go) 
            plt.plot(self.output_mag_field_back,self.output_ratio_back)
            plt.plot(self.output_mag_field_minor,self.output_ratio_minor)
            plt.savefig("MR_ratio_curve_w_minor_" + self.major_filename + "_" + self.minor_filename + ".png")
            plt.show()
            plt.close()
            
            plt.figure(figsize=(8, 6), dpi=80)
            plt.xlim([-0.03,0.03])
            plt.xlabel(r'$\mu_0 H$(T)')
            plt.ylabel(r'Resistance($\Omega$)')
            plt.axvline(x=0,linestyle='dashed')
            plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
            plt.plot(self.output_mag_field_go,self.output_resistance_go) 
            plt.plot(self.output_mag_field_back,self.output_resistance_back)
            plt.plot(self.output_mag_field_minor,self.output_resistance_minor)
            plt.savefig("MR_curve_w_minor_ex_" + self.major_filename + "_" + self.minor_filename + ".png")
            plt.show()
            plt.close()
            
            plt.figure(figsize=(8, 6), dpi=80)
            plt.xlim([-0.03,0.03])
            plt.xlabel(r'$\mu_0 H$(T)')
            plt.ylabel(r'MR ratio (%)')
            plt.axvline(x=0,linestyle='dashed')
            plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
            plt.plot(self.output_mag_field_go,self.output_ratio_go) 
            plt.plot(self.output_mag_field_back,self.output_ratio_back)
            plt.plot(self.output_mag_field_minor,self.output_ratio_minor)
            plt.savefig("MR_ratio_curve_w_minor_ex_" + self.major_filename + "_" + self.minor_filename + ".png")
            plt.show()
            plt.close()
        
        
#class_MR = MR("20161128-164452-E431_low_100mv_100-Major",100,"20161128-165522-E431_low_100mv_100-minor")
#class_MR.MR_curve_plot()
    
    """
    def __init__(self,major_filename,vias_vol,minor_filename = ""):
        self.major_filename = major_filename 
        self.vias_vol = vias_vol
        self.minor_filename = minor_filename 
        if(minor_filename == ""):
            self.minor_exist_or_not = False
        else:
            self.minor_exist_or_not = True
            self.input_data_minor = np.loadtxt(minor_filename + ".csv",delimiter = ",",skiprows=1)
        self.input_data_major = np.loadtxt(major_filename + ".csv",delimiter = ",",skiprows=1)
         
        
        
        
    def fix_mag_field(self,x,go_or_not):
        if(go_or_not):
            return (0.0011 * x * 1000 + 0.0627) * 0.1
        else:
            return (0.0011 * x * 1000 - 0.0546) * 0.1

    def calc_TMR_ratio(self,out_go, out_back):
        Rp = (sum(out_go[0:10]) / 10 + sum(out_go[len(out_go)-10:len(out_go)]) / 10 + sum(out_back[0:10]) / 10 + sum(out_back[len(out_back)-10:len(out_back)]) / 10) / 4
        Rap = (max(out_go) + max(out_back)) / 2
        self.TMR = (Rap / Rp - 1) * 100
        return (Rap / Rp - 1) * 100

    def MR_curve_plot(self):
        i = 0
        row_len = self.input_data_major.size /self.input_data_major[0].size
        output_mag_field_go=[]
        output_resistance_go =[]
        output_mag_field_back=[]
        output_resistance_back =[]
        go_or_not = True
        vias_vol = self.vias_vol * 0.001
        while i < row_len:
            if(i > row_len / 2 ):
                go_or_not = False
            if(go_or_not):
                output_mag_field_go.append(self.fix_mag_field(self.input_data_major[i][0],go_or_not))
                output_resistance_go.append((vias_vol + 0.00003) / self.input_data_major[i][1])
            else:
                output_mag_field_back.append(self.fix_mag_field(self.input_data_major[i][0],go_or_not))
                output_resistance_back.append((vias_vol - 0.00003) / self.input_data_major[i][1])
            i += 1   
        print(self.calc_TMR_ratio(output_resistance_go,output_resistance_back))
        #self.TMR = self.calc_TMR_ratio(output_resistance_go,output_resistance_back)
        from matplotlib import pyplot as plt
        
        Rp = (sum(output_resistance_go[0:10]) / 10 + sum(output_resistance_go[len(output_resistance_go)-10:len(output_resistance_go)]) / 10 + sum(output_resistance_back[0:10]) / 10 + sum(output_resistance_back[len(output_resistance_back)-10:len(output_resistance_back)]) / 10) / 4
        output_ratio_go = (output_resistance_go / Rp - 1) * 100
        output_ratio_back = (output_resistance_back / Rp - 1) * 100
        
        plt.figure(figsize=(8, 6), dpi=80)
        plt.xlim([-0.3,0.3])
        plt.xlabel(r'$\mu_0 H$(T)')
        plt.ylabel(r'Resistance($\Omega$)')
        plt.axvline(x=0,linestyle='dashed')
        plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
        plt.plot(output_mag_field_go,output_resistance_go) 
        plt.plot(output_mag_field_back,output_resistance_back) 
        plt.savefig("MR_curve_wo_minor_" + self.major_filename + "_" + self.minor_filename + ".png")
        plt.show()
        plt.close()
        
        plt.figure(figsize=(8, 6), dpi=80)
        plt.xlim([-0.3,0.3])
        plt.xlabel(r'$\mu_0 H$(T)')
        plt.ylabel(r'MR ratio (%)')
        plt.axvline(x=0,linestyle='dashed')
        plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
        plt.plot(output_mag_field_go,output_ratio_go) 
        plt.plot(output_mag_field_back,output_ratio_back) 
        plt.savefig("MR_ratio_curve_wo_minor_" + self.major_filename + "_" + self.minor_filename + ".png")
        plt.show()
        plt.close()
        
        if(self.minor_exist_or_not):
            i = 0
            row_len = self.input_data_minor.size /self.input_data_minor[0].size
            output_mag_field_minor=[]
            output_resistance_minor =[]
            while i < row_len:    
                output_mag_field_minor.append(self.fix_mag_field(self.input_data_minor[i][0],True))
                output_resistance_minor.append((vias_vol + 0.00003) / self.input_data_minor[i][1])
                i += 1
            output_ratio_minor = (output_resistance_minor / Rp - 1) * 100
            
            plt.figure(figsize=(8, 6), dpi=80)
            plt.xlim([-0.3,0.3])
            plt.xlabel(r'$\mu_0 H$(T)')
            plt.ylabel(r'Resistance($\Omega$)')
            plt.axvline(x=0,linestyle='dashed')
            plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
            plt.plot(output_mag_field_go,output_resistance_go) 
            plt.plot(output_mag_field_back,output_resistance_back)
            plt.plot(output_mag_field_minor,output_resistance_minor)
            plt.savefig("MR_curve_w_minor_" + self.major_filename + "_" + self.minor_filename + ".png")
            plt.show()
            plt.close()
            
            plt.figure(figsize=(8, 6), dpi=80)
            plt.xlim([-0.3,0.3])
            plt.xlabel(r'$\mu_0 H$(T)')
            plt.ylabel(r'MR ratio (%)')
            plt.axvline(x=0,linestyle='dashed')
            plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
            plt.plot(output_mag_field_go,output_ratio_go) 
            plt.plot(output_mag_field_back,output_ratio_back)
            plt.plot(output_mag_field_minor,output_ratio_minor)
            plt.savefig("MR_ratio_curve_w_minor_" + self.major_filename + "_" + self.minor_filename + ".png")
            plt.show()
            plt.close()
            
            plt.figure(figsize=(8, 6), dpi=80)
            plt.xlim([-0.03,0.03])
            plt.xlabel(r'$\mu_0 H$(T)')
            plt.ylabel(r'Resistance($\Omega$)')
            plt.axvline(x=0,linestyle='dashed')
            plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
            plt.plot(output_mag_field_go,output_resistance_go) 
            plt.plot(output_mag_field_back,output_resistance_back)
            plt.plot(output_mag_field_minor,output_resistance_minor)
            plt.savefig("MR_curve_w_minor_ex_" + self.major_filename + "_" + self.minor_filename + ".png")
            plt.show()
            plt.close()
            
            plt.figure(figsize=(8, 6), dpi=80)
            plt.xlim([-0.03,0.03])
            plt.xlabel(r'$\mu_0 H$(T)')
            plt.ylabel(r'MR ratio (%)')
            plt.axvline(x=0,linestyle='dashed')
            plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
            plt.plot(output_mag_field_go,output_ratio_go) 
            plt.plot(output_mag_field_back,output_ratio_back)
            plt.plot(output_mag_field_minor,output_ratio_minor)
            plt.savefig("MR_ratio_curve_w_minor_ex_" + self.major_filename + "_" + self.minor_filename + ".png")
            plt.show()
            plt.close()
    """
