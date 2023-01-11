#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 13:06:58 2022

@author: neerbos

see speciestable(2).py for my original version. where i define my own regions
version used for latest tables 18-10-22

this version corresponds to the one used for the red/blue/purple/black table in the paper!
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
vel,dens=7,1e4
output_file="../output/shocks_output/Jshock{}-{}-1-1.dat".format(dens,vel)
df=uclchem.analysis.read_output_file(output_file=output_file)

def enhanced(lol,dens,vel,printit=False,shock="C",mylist=False): #now i follow tomas' paper to define region boundaries
    """calculates enhancement factors and destruction factors
    of molecules
    """
    output_file="../output/shocks_output/Jshock{}-{}-1-1.dat".format(dens,vel)
    df=uclchem.analysis.read_output_file(output_file=output_file)
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
        # finds index where density is at halfmax, close to dissipation time

    #print(idx1)
    idx2=df["Time"].sub(1E8/dens+df["Time"][idx1]).abs().idxmin()
    #returns the idx where time is 100 years after shock for high density
    #and where time is 100000 years after shock at low density
    #print(df["Time"][idx1],df["Time"][idx2],dens, shock]  
    #print(idx1,idx2,dens,vel,df["Time"].size)
    #print("time1:", df["Time"][idx1], "Time2: ", df["Time"][idx2])
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
        elif column not in notspecies and column not in elementList:
            if not "@" in column and not "#" in column and not "+" in column:
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
                print(column, maxs[column])
    if printit:      
        print("enhanced:",enhs,"\n","destroyed:",dests,"\n","maxima:",maxs)       
    return enhs,enhslong,dests,maxs,starts,avgs1,avgs2,Tmax


def table(shock="C",minenh=100,mindest=10,minabd=1E-9,writefile=False):
    """ gives enhanced species in blue and destroyed species in red
    minenh, mindest and minabd give the minimum values of the enhancement factors, destruction factors and abundances
    to consider
    returns a nested dictionary of format eg {'1000.0-10':{'blue':list, 'red':list}}"""
    allshocks=defaultdict(dict)
    count=0
    if writefile:
        f=open("../output/speciestables/{}shock-jon{}(4).tex".format(shock,jon),"w")
        r"""\usepackage{tabularx}
        \usepackage{ltablex}
        \usepackage{seqsplit}
        \usepackage[margin=1cm]{geometry} % Changing page margin
        \usepackage{array}
        \usepackage{color}
        \usepackage{collcell}"""
        f.write("\\begin{tabularx}{500pt}{s|s|s|s|s|} \n     vel\\textbackslash dens &1E3&1E4&1E5&1E6\\\\ \n\\hline \n")
    for vel in allv:
        if writefile:
            f.write("{}".format(vel))
        for dens in alldensity:

            output_file="../output/shocks_output/{}shock{}-{}-1-1.dat".format(shock,dens,vel)
            df=uclchem.analysis.read_output_file(output_file=output_file)
            enhs,enhslong,dests,maxs,starts,avgs1,avgs2,maxT=enhanced(df,dens,vel,shock=shock)
            blue,red,black,purple=[],[],[],[]
            for species in enhs:
                if avgs1[species]>=minabd: #the species that start at minimum observable
                    if enhs[species]>=minenh and dests[species]<mindest:
                        blue.append(species)
                        count+=1
                    elif enhs[species]>=minenh and dests[species]<0.1*enhs[species]:
                        purple.append(species)
                        count+=1
                    elif dests[species]>=mindest:
                        red.append(species)
                        count+=1
                elif avgs2[species] >= minabd:
                    black.append(species) #species enhanced from invisible to visible
                    count+=1
            allshocks["{}-{}".format(dens,vel)]["blue"]=blue
            allshocks["{}-{}".format(dens,vel)]["red"]=red
            allshocks["{}-{}".format(dens,vel)]["black"]=black
            allshocks["{}-{}".format(dens,vel)]["purple"]=purple
            #print(type(maxT),maxT)
            try:
                if maxT >= 10000:
                    maxT = "\\textgreater 10000"
            except TypeError:
                maxT = "\\textgreater 10000"
            if writefile:
                f.write("&\\textcolor{{gray}}{{{3}K}}  \
                        \\textcolor{{blue}}{{{0}}}  \
                        \\textcolor{{violet}}{{{4}}}  \
                        \\textcolor{{red}}{{{1}}} {2}".format(blue,red,black,maxT,purple))#prints it in latex readable format
        if writefile:
            f.write("\\\\ \n\\hline \n")
    if writefile:
        f.write("\\end{tabularx}")
        f.close()
    print(count)
    return allshocks

def JandCratio(species,mylist=True,save=False,longrun=False):
    ratio={}
    values=[]
    region="mid"
    for dens in alldensity:
        ratio[dens]={}
        for vel in allv:
            output_file="../output/shocks_output/{}shock{}-{}-1-1.dat".format("C",dens,vel)
            df=uclchem.analysis.read_output_file(output_file=output_file)
            Cenhs,Clong,Cdests,maxs,starts,avgs1,avgs2,maxT=enhanced(df,dens,vel,shock="C",mylist=mylist)

            output_file="../output/shocks_output/{}shock{}-{}-1-1.dat".format("J",dens,vel)
            df=uclchem.analysis.read_output_file(output_file=output_file)
            Jenhs,Jlong,Jdests,maxs,starts,avgs1,avgs2,maxT=enhanced(df,dens,vel,shock="J",mylist=mylist)
            if longrun:
                Jenhs=Jlong
                Cenhs=Clong
                region="post"
            ratio[dens][vel]=Jenhs[species]/Cenhs[species]
            values.append(Jenhs[species]/Cenhs[species])

    #print(ratio)
    #print(max(values))
    f=plt.figure()
    f.set_figwidth(11)
    f.set_figheight(2.5)
    for dens in alldensity:
        for vel in allv:
            plt.scatter(vel,np.log10(dens),c=ratio[dens][vel],
                        norm=clr.LogNorm(vmin=min(values),vmax=max(values)),
                        cmap="rainbow", marker='s',s=1e3)
    plt.colorbar().set_label("X(J)/X(C)")
    plt.title("J/C enhancement factors for {} - {}-shock".format(species, region))
    plt.ylabel("log density")
    plt.xlabel("shock velocity")
    if mylist:
        name="bowie"
    else:
        name="tomas"
    if save:
        plt.savefig("../output/figures/{}{}Enh_{}.png".format(name,region,species))
    plt.show()
    plt.close()
            
def onelist(dens,vel,name="C"):
    output_file="../output/shocks_output/{}shock{}-{}-1-1.dat".format(name,dens,vel)
    df=uclchem.analysis.read_output_file(output_file=output_file)
    enhanced(df,dens,vel,printit=True,shock=name)
def fulloutput():
    species=table()
    for i in species:
        print(i,species[i]["blue"])
#onelist(1E5,50,name="J")
enhanced(df,1e4,10,shock="J")
#for molecule in jonlist:
#    JandCratio(molecule,mylist=True,longrun=True,save=True)
#    JandCratio(molecule,mylist=True,longrun=False,save=True)
#table=table(shock="J",writefile=True)
#print(table)
#table=table(shock="J",writefile=True)
#print(table)
