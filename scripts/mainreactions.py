#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 12:19:13 2022

@author: neerbos
"""
import sys
import os



dirlist=["/net/harbinger/data2/bowie/UCLCHEM2/scripts",
"/usr/lib64/python310.zip",
"/software/local/lib64/python3.10/site-packages",
"/usr/lib64/python3.10",
"/usr/lib64/python3.10/lib-dynload",
"/home/neerbos/.local/lib/python3.10/site-packages",
"/usr/lib64/python3.10/site-packages",
"/usr/lib/python3.10/site-packages",
"/usr/lib/python3.10/site-packages/IPython/extensions",
"/home/neerbos/.local/lib/python3.10/site-packages/uclchem/__init__.py",
]
for i in dirlist:
    sys.path.append(i)
#import uclchem
import uclchem


import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy as np
import pandas as pd
import seaborn
from scipy import integrate
#from tabulate import tabulate as tbl
from collections import defaultdict


def analyse(species, shock = "J", vel=8, dens=1E3):
    uclchem.analysis.analysis(species, 
                              "../output/shocks_output/{}shock{}-{}-1-1.dat".format(shock, dens, vel),
                              "../output/temp/{}-{}-{}-{}".format(shock,dens,vel,species))
    
analyse("HNCO", shock = "J")