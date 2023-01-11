#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 10:56:33 2022

@author: neerbos

need to run a loop and save several different files.
format of files: shock-dens-vel-crir-uv)
"""
import uclchem
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from itertools import product

#set a parameter dictionary for phase 1 collapse model
pc,year=30.9E17,365.25*86400 #cm, s resp
name,vs,density="C",45,100000.0
av,crir,uv=1,1,1

alldensity=[1e3,1e4,1e5,1e6]
#allv=[10,15,20,25,30,35,40,45,50,55,60]
allv=[5,6,7,8,9,10,11,12,13,14,15]

allcrir=[0.1,1,10,100] #crir=cosmic ray ion rate, visex=visual extinction
alluv=[0.1,1,10,100]
elements=["H","C","O","N","S"]
#uclchem.run_model(param_dict)
def element_check(output_file, density, printit=False):
    df=uclchem.analysis.read_output_file(output_file)
    idx1=df["Density"].sub(df["Density"][0]+0.5*max(df["Density"]-df["Density"][0])).abs().idxmin() # finds index where T is at max, close to dissipation time

    #print(idx1)
    idx2=df["Time"].sub(df["Time"][idx1]+1E8/density).abs().idxmin()
    #print(df["Time"][idx2])
    df=df[:idx2]
    #get conservation values
    conserves=uclchem.analysis.check_element_conservation(df,element_list=elements)
    #check if any error is greater than 1%
    if printit:
        print(conserves)
    return all([float(x[:-1])<2 for x in conserves.values()])
    
    
def checkallfiles(): #checks whether all files reached at least 5E5years 
    for shock in ["C","J"]:
        for dens in alldensity:
            for vel in allv:
                file="../output/shocks_output/{}shock{}-{}-1-1\
.dat".format(shock, dens, vel)
                if not element_check(file,dens,printit=True):
                    print("for {} elements not conserved".format(file))
                """try:
                    df=uclchem.analysis.read_output_file(file)
                    try:
                        hond=np.max(np.where((df["Time"])>=5E5))
                    except ValueError:
                        print(file)
                except FileNotFoundError:
                    print(file, "not found")"""
                            
def abundens(vs,density,crir=1,uv=1,save=False,log=True,time=1E5,cm=False):
    for name in ["C","J"]:
        if cm:
            pc=1
            label="cm"
        else:
            pc=30.9E17
            label="zn/pc"
        file="../output/shock_outputs(3)/{}shock{}-{}-{}-{}.dat".format(name,density,vs,crir,uv)
        phase2_df=uclchem.analysis.read_output_file(file)
        species=["CH3OH","H2O","SO","SO2","OCS","CO2","HCO","C3H2"]
        try:
            print(np.max(np.where((phase2_df["Time"])>=5E5)))
        except:
            print("dtime {} dlength {}".format(12*pc/(1E5*density*year),12*vs/density))
        if log:
            scale="log"
            beginat=1
        else:
            scale="linear"
            beginat=-20
        fig,[ax,ax2]=plt.subplots(2,1,figsize=(16,9))
        ax=uclchem.analysis.plot_species(ax,phase2_df,species)
        settings=ax.set(yscale=scale,xscale=scale,
                    xlabel="time in yr", 
                    ylabel="X$_{Species}$",
                    xlim=(beginat,time))
        ax4=ax.twiny()
        ax.grid()
        halfmax=np.min(np.where(phase2_df["Density"]>=0.5*(np.max(phase2_df["Density"]-density))+density))
        time1=phase2_df["Time"][halfmax]
        time2=phase2_df["Time"][np.min(np.where(phase2_df["Time"]>=5*time1))]
        ax.axvline(time1,c="black",alpha=0.5)
        ax2.axvline(time1,c="black",alpha=0.5)
        ax.axvline(time2,c="black",alpha=0.5)
        ax2.axvline(time2,c="black",alpha=0.5)
        ax4.plot(phase2_df["Time"]*year*vs*1E5/pc,phase2_df["HCL"],alpha=0)
        ax4.set(xlabel=label,xscale=scale,yscale="log",xlim=(beginat*year*vs*1E5/pc,time*year*vs*1E5/pc))
        ax2.plot(phase2_df["Time"],phase2_df["Density"],color="black")
        ax3=ax2.twinx()
        ax3.plot(phase2_df["Time"],phase2_df["gasTemp"],color="red")
        ax2.set(xlabel="Time / year",ylabel="Density")
        ax3.set(ylabel="Temperature",facecolor="red",xscale=scale,yscale = "linear",xlim=(beginat,time))
        ax3.tick_params(axis='y', colors='red')
        ax2.grid()
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        text="shock:{} \nspeed:{} km/s \ndensity:10^{} $n_H$/cm^3\nvisex:{} \ncrir:{} * standard\nuv:{} Ha".format(name,vs,np.log10(density),av,crir,uv)
        ax2.text(0.7,-0.5,text,transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)
        if save:
            plt.savefig("../output/figures/abundens_{}:{}-{}-{}-{}-{}.png".format(name,vs,int(density),av,crir,uv))
        else:
            plt.show()

def plotall(save = False):
    for density in alldensity:
        for vs in allv:
            abundens(vs,density,av,crir,uv,log=True,time=1E6,save=save)#6*uclchem.utils.cshock_dissipation_time(vs, density))
av=1
#plotav(name,vs,density,av,crir,uv)
dens,vel=1e6,15

file="../output/shocks_output/{}shock{}-{}-1-1.dat".format("C", dens, vel)
#element_check(file,dens, printit=True)
#checkallfiles()
abundens(10,1E5,crir=1,uv=1, time=1E6)
#    abundens(vs,density,av,crir,uv,save=False)
#plotall(save=True)
