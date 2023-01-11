#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 12:18:33 2022

@author: neerbos
"""

import sys
import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Pool
import itertools
import os
import uclchem

alldens, allv=[1E3,1E4,1E5,1E6],[5,6,7,8,9,10,11,12,13,14,15]


for dens in alldens:
    for vel in allv:
        check=False
        file1=uclchem.analysis.read_output_file(
            "../output/shocks_output/Cshock{}-{}-1-1.dat".format(dens,vel))
        file2=uclchem.analysis.read_output_file(
            "../output/shocks_output/Jshock{}-{}-1-1.dat".format(dens,vel))
        for column in file1:
            if file1[column][0] != file2[column][0]:
                check=True
        if check:
            print(dens,vel)
