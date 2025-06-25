
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 11:40:53 2024
@author: Joe Strummer 
"""
import numpy as np
import pandas as pd
import itertools
import random
import csv
import os
from pathlib import Path
import Boid_Class_PDF as boid

N = 10000 ## pop size
Time_steps = 151
start_time = 50
n_sims = 1000

cwd = Path.cwd()
mod_path = Path(__file__).parent

export_path = '../Sim/Results/'

for sim in range(542, n_sims):
    Newborns = []
    population = []
    max_ID = 0    
    ## initaialise population for each simualtion: let sex be defined through 50 50 split and ID's unique additative
    for i in range(1,N):
        population.append(boid.pop_member(ID = i))
    # record the max ID
    max_ID = max([i.ID for i in population])
    # iterate over timesteps on the ecological stage...
    for kk in range(0,Time_steps):
        population_newborns = []
        population_dead = []
################################### Sample a cohort of Focals ##################################################
        if kk == start_time:
            Newborns = [i for i in population if i.age == 0]
            filename =  str(export_path) + "Time" + str(kk) + "replication" + str(sim) + ".csv"
            with open(filename, 'w', newline = '') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(['Year', 'Kin', 'Kin ID', 'Alive', 'Age Kin',  'Age Focal', 'Alive Focal', 'Focal ID', 'Focals mother ID'])
                for ss in Newborns:
                    ss.extract_kin(kk, population)
                    writer.writerows(np.array(ss.kin_network))
################################### At later times we follow the newborn cohort ##################################################
        elif kk > start_time:
            filename = str(export_path) + "Time" + str(kk) + "replication" + str(sim) + ".csv"
            with open(filename, 'w', newline = '') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(['Year', 'Kin', 'Kin ID', 'Alive', 'Age Kin', 'Age Focal','Alive Focal', 'Focal ID', 'Focals mother ID'])
                for ss in Newborns:
                    if ss.alive == 1:
                        ss.extract_kin(kk, population)
                        writer.writerows(np.array(ss.kin_network))
     
        # loop through all population members and allow 1) aging, 2) deaths, 3) marriage events -- if married (newly or already) then can reproduce
        for ss in population:  
            ss.die() ## select those who die, but don't remove just yet -- instead make a list of those who have died
            if ss in Newborns: ## we require the cohort of Focal's to be immortal (sampling ease)
                ss.alive = 1
            if ss.alive == 0:
                population_dead.append(ss)         
        
        for ss in [i for i in population if i.alive == 1]:
            baby_prob = np.random.poisson(ss.fert, 1)[0]
            if baby_prob != 0 and ss.alive == 1:
                for i in range(0, baby_prob): ## Note that baby_prob from Poisson draw and can be >1
                    baby = boid.sex(ss, max_ID + 1)
                    max_ID = max_ID + 1
                    population_newborns.append(baby)
        
        for ss in population:  
            ss.advance_age() 
            if ss.age > 100:
                population.remove(ss)            
        
        for i in population_newborns:
            population.append(i)

        max_ID = max([i.ID for i in population])    
        if len(population_dead)>0:  
            #population.remove(population_dead) 
            for j in [i for i in population if i in population_dead]:
                    population.remove(j) 
 
        




