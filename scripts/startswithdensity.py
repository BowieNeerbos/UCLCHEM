#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 13:06:58 2022

@author: neerbos

see speciestable(2).py for my original version. where i define my own regions
version used for latest tables 18-10-22
"""
import uclchem
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy as np
import pandas as pd
import seaborn
from scipy import integrate
from tabulate import tabulate as tbl
from collections import defaultdict

jon=False
#set a parameter dictionary for phase 1 collapse model
pc,year=30.9E12,3.15e7 #in km, s
dens,vel=1E5,40
alldensity=[1e3,1e4,1e5,1e6]
#allv=[5,10,15,20,25,30,35,40,45,50,55,60]
allv=[5,6,7,8,9,10,11,12,13,14,15]
allcrir,alluv=[0.1,1,10,100],[0.1,1,10,100] #crir=cosmic ray ion rate uv in Habing
notspecies=["Time","Density","gasTemp","av","zeta","point","BULK","SURFACE","radfield"]
tomaslist=["SO","CH3OH","H2O","HCN","SO2"]
jonlist=["H2O","HCL","C3H2","CS","HS","H2CN","HCN","H2CS","CH3OH","NH3","NS","H2S","CH","CN","HCO","O2","OH","OCS","HNO","HNC","CO2","SO","SO2"]
elementList = [
    "H",
    "D",
    "HE",
    "C",
    "N",
    "O",
    "F",
    "P",
    "S",
    "CL",
    "LI",
    "NA",
    "MG",
    "SI",
    "PAH",
    "15N",
    "13C",
    "18O",
    "SURFACE",
    "BULK",
    "Time",
    "E-"
]
vel,dens=7,1e5
output_file="../output/shock_outputs(3)/Jshock{}-{}-1-1.dat".format(dens,vel)
df=uclchem.analysis.read_output_file(output_file=output_file)

def enhanced(df,dens,vel,printit=False,shock="C",mylist=False): #now i follow tomas' paper to define region boundaries
    """calculates enhancement factors and destruction factors
    of molecules
    """

    enhs,enhslong,dests,maxs,starts,avgs1,avgs2,={},{},{},{},{},{},{} 
    #lists of enhancement factors, destruction and max abundances respectively
    disstime=uclchem.utils.cshock_dissipation_time(vel, dens)
    timeshock=df["Time"][df["Density"].idxmax()]
    #print (timeshock)
    if shock == "C":
        idx1=df["Time"].sub(disstime).abs().idxmin()
    else:#shock ="J"
        idxlist=df[df["gasTemp"]==10].index.values
        idx1=idxlist[1] #the first value where T=10K, excluding t=0
    if mylist: #my cutoff values, at half max density
        idx1=df["Density"].sub(dens+0.5*max(df["Density"]-dens)).abs().idxmin() 
        # finds index where T is at halfmax, close to dissipation time

    #print(idx1)
    idx2=df["Time"].sub(1E8/dens+df["Time"][idx1]).abs().idxmin()
    #returns the idx where time is 100 years after shock for high density
    #and where time is 100000 years after shock at low density
    #print(df["Time"][idx1],df["Time"][idx2],dens, shock)
        
    #print(idx1,idx2,dens,vel,df["Time"].size)
    for column in df:
        if jon:
            if column in jonlist:
                avg1=integrate.trapz(df[column][0:idx1],x=df["Time"][0:idx1])/df["Time"][idx1]
                #average abundance throughout the shock
                avg2=integrate.trapz(df[column][idx1:idx2],x=df["Time"][idx1:idx2])/(df["Time"][idx2]-df["Time"][idx1])
                Tmax=np.max(df["gasTemp"])
                #average abundance after the shock
                start=df[column][0]
                change=avg1/start
                changelong=avg2/start
                destroy=change/changelong
                #if change >=1:
                avgs1[column]=avg1
                avgs2[column]=avg2
                enhs[column]=change
                enhslong[column]=changelong
                dests[column]=destroy
                maxs[column]=max(df[column])
                enhs=dict(sorted(enhs.items(),key=lambda item:item[1],reverse=True))
                starts[column]=start
        elif column not in notspecies and not column in elementList:
            if not "@" in column and not "#" in column and not "+" in column:
                avg1=integrate.trapz(df[column][0:idx1],x=df["Time"][0:idx1])/df["Time"][idx1]
                #average abundance throughout the shock
                avg2=integrate.trapz(df[column][idx1:idx2],x=df["Time"][idx1:idx2])/(df["Time"][idx2]-df["Time"][idx1])
                #average abundance after the shock
                #plt.plot(df["Time"],df[column])
                #plt.yscale("log")
                #plt.title(column)
                #plt.show()
                start=df[column][0]
                change=avg1/start
                destroy=avg1/avg2
                #if change >=1:
                avgs1[column]=avg1
                avgs2[column]=avg2
                enhs[column]=change
                dests[column]=destroy
                maxs[column]=max(df[column])
                enhs=dict(sorted(enhs.items(),key=lambda item:item[1],reverse=True))
                starts[column]=start
    if printit:      
        print("enhanced:",enhs,"\n","destroyed:",dests,"\n","starts:",starts)       
    return enhs,enhslong,dests,maxs,starts,avgs1,avgs2





output_file="../output/shocks_output/{}shock{}-5-1-1.dat".format("J",1E3)
df=uclchem.analysis.read_output_file(output_file=output_file)
enhs,long,dests,maxs,startslow,avgs1,avgs2=enhanced(df,dens,vel,shock="J",mylist=True)
output_file="../output/shocks_output/{}shock{}-5-1-1.dat".format("J",1E6)
df=uclchem.analysis.read_output_file(output_file=output_file)
enhs,long,dests,maxs,startshigh,avgs1,avgs2=enhanced(df,dens,vel,shock="J",mylist=True)
output_file="../output/shocks_output/{}shock{}-5-1-1.dat".format("J",1E4)
df=uclchem.analysis.read_output_file(output_file=output_file)
enhs,long,dests,maxs,startsmid1,avgs1,avgs2=enhanced(df,dens,vel,shock="J",mylist=True)
output_file="../output/shocks_output/{}shock{}-5-1-1.dat".format("J",1E5)
df=uclchem.analysis.read_output_file(output_file=output_file)
enhs,long,dests,maxs,startsmid2,avgs1,avgs2=enhanced(df,dens,vel,shock="J",mylist=True)
for column in startslow:
    plt.plot([3,4,5,6],[startslow[column],startsmid1[column],startsmid2[column],startshigh[column]])
    plt.yscale("log")
    plt.title(column)
    plt.show()