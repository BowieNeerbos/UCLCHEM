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
import os

jon=True
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
vel,dens=7,1e6
output_file="../output/shocks_output/Jshock{}-{}-1-1.dat".format(dens,vel)
df=uclchem.analysis.read_output_file(output_file=output_file)

def enhanced(df,dens,vel,printit=False,shock="C",mylist=False, species="all"): #now i follow tomas' paper to define region boundaries
    """calculates enhancement factors and destruction factors
    of molecules
    """
    global jon, jonlist
    enhs,enhslong,dests,maxs,starts,avgs1,avgs2,avgsfull={},{},{},{},{},{},{},{}
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
    if species != "all":
        jon = True
        jonlist = [species]
    for column in df:
        if jon:
            if column in jonlist:
                avg1=integrate.trapz(df[column][0:idx1],x=df["Time"][0:idx1])/df["Time"][idx1]
                #average abundance throughout the shock
                avg2=integrate.trapz(df[column][idx1:idx2],x=df["Time"][idx1:idx2])/(df["Time"][idx2]-df["Time"][idx1])
                avg3=integrate.trapz(df[column][0:idx2],x=df["Time"][0:idx2])/df["Time"][idx2]
                Tmax=np.max(df["gasTemp"])
                #average abundance after the shock
                start=df[column][0]
                change=avg1/start
                changelong=avg3/start
                destroy=change/changelong
                #if change >=1:
                avgs1[column]=avg1
                avgs2[column]=avg2
                avgsfull[column]=avg3
                #print (shock,vel,dens,avg3,(avg1*df["Time"][idx1]+avg2*(df["Time"][idx2]-df["Time"][idx1]))/(df["Time"][idx2]))
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
                avg3=integrate.trapz(df[column][0:idx2],x=df["Time"][0:idx2])/df["Time"][idx2]
                Tmax=np.max(df["gasTemp"])
                #average abundance after the shock
                start=df[column][0]
                change=avg1/start
                changelong=avg3/start
                destroy=change/changelong
                #if change >=1:
                avgs1[column]=avg1
                avgs2[column]=avg2
                avgsfull[column]=avg3
               # print (avg3 - (avg1*df["Time"][idx1]+avg2*(df["Time"][idx2]-df["Time"][idx1]))/(df["Time"][idx2]))
                enhs[column]=change
                enhslong[column]=changelong
                dests[column]=destroy
                maxs[column]=max(df[column])
                enhs=dict(sorted(enhs.items(),key=lambda item:item[1],reverse=True))
                starts[column]=start
    if printit:      
        print("enhanced:",enhs,"\n","destroyed:",dests,"\n","maxima:",maxs)       
    return enhs,enhslong,dests,maxs,starts,avgs1,avgs2,Tmax,avgsfull


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
            enhs,enhslong,dests,maxs,starts,avgs1,avgs2,maxT,avgsfull=enhanced(df,dens,vel,shock=shock)
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

def Jtable(species,shock="J",mylist=True,save=False,longrun=False):
    ratio={}
    fullavg={}
    values=[]
    region="mid"
    for dens in alldensity:
        ratio[dens]={}
        fullavg[dens]={}
        for vel in allv:
            output_file="../output/shocks_output/{}shock{}-{}-1-1.dat".format("C",dens,vel)
            df=uclchem.analysis.read_output_file(output_file=output_file)
            Cenhs,Clong,Cdests,Cmaxs,starts,avgs1,avgs2,maxT,Cavgsfull=enhanced(df,dens,vel,shock="C",mylist=mylist)

            output_file="../output/shocks_output/{}shock{}-{}-1-1.dat".format("J",dens,vel)
            df=uclchem.analysis.read_output_file(output_file=output_file)
            Jenhs,Jlong,Jdests,Jmaxs,starts,avgs1,avgs2,maxT,Javgsfull=enhanced(df,dens,vel,shock="J",mylist=mylist,species=species)
            if longrun:
                Jenhs=Jlong
                Cenhs=Clong
                region="post"
            if shock == "J":
                ratio[dens][vel]=Jenhs[species]#so no longer a ratio
                fullavg[dens][vel]=Javgsfull[species]#gives the avg abundance achieved in j shock
                values.append(Jenhs[species])
            else:
                ratio[dens][vel]=Cenhs[species]#so no longer a ratio
                fullavg[dens][vel]=Cavgsfull[species]#gives the avg abundance achieved in j shock
                values.append(Cenhs[species])              

    #print(ratio)
    #print(max(values))
    f=plt.figure()
    f.set_figwidth(11)
    f.set_figheight(2.5)
    clrRange=10**max(abs(np.log10(min(values))),abs(np.log10(max(values))))
    for dens in alldensity:
        for vel in allv:
            plt.scatter(vel,np.log10(dens),c=ratio[dens][vel],
                        norm=clr.LogNorm(vmin=clrRange**-1,vmax=clrRange),
                        cmap="seismic", marker='s',s=1e3)
            plt.annotate((np.log10(fullavg[dens][vel])).round(2),(vel-0.25,np.log10(dens)),c="gold")
    plt.colorbar().set_label("Enhancement factor")
    plt.title("{} enhancement factors for {} - {}-shock".format(shock, species, region))
    plt.ylabel(r"$log(\rho_{pre})\; (cm^{-3})$")
    plt.xlabel("shock velocity (km/s)")
    if mylist:
        name="bowie"
    else:
        name="tomas"
    if save:
        plt.tight_layout()
        plt.savefig("../output/{0}_enhs/{0}{1}Enh_{2}.png".format(shock,region,species))
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

#for molecule in jonlist:
#    JandCratio(molecule,mylist=True,longrun=True,save=True)
#    JandCratio(molecule,mylist=True,longrun=False,save=True)
jon=False  
for molecule in df:
    if molecule not in elementList and molecule not in notspecies:
        if not "@" in molecule and not "#" in molecule and not "+" in molecule:
            if not os.path.exists("../output/J_enhs/JpostEnh_{}.png".format(molecule)):
                Jtable(molecule,mylist=True,longrun=True,save=False)
                #Jtable(molecule,mylist=True,longrun=True,save=True,shock="C")
            else:
                print("skipped", molecule)
Jtable("NH",mylist=True,longrun=True,save=False,shock="J")
#Jtable("H2",shock="J",mylist=True, longrun=True)