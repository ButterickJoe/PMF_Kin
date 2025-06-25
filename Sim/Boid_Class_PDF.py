# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 15:34:11 2024

@author: Joe Strummer
"""

#Libraries
import random
import numpy as np
import random
import csv
import pandas as pd
from pathlib import Path

cwd = Path.cwd()
mod_path = Path(__file__).parent


relative_path_fert = '../data/1974_female_fert.csv'
fert_females = pd.read_csv(relative_path_fert)
fert_females = pd.DataFrame(fert_females)
frates_f = fert_females.filter(items = ['Age' , 'Fert'])

relative_path_mort = '../data/1974_female_mort.csv'
mort_females = pd.read_csv(relative_path_mort)
mort_females = pd.DataFrame(mort_females)
mrates_f = mort_females.filter(items = ['Age' , 'Mort'])

relative_path_popstruct = '../data/1974_female_ps.csv'
stable_female_pop = pd.read_csv(relative_path_popstruct)
fv1 = pd.DataFrame(stable_female_pop)
fv1 = fv1.filter(items = ['x'])
fv1 = np.array(fv1)
fv1 = fv1.flatten()
fv1 = fv1/sum(fv1)


class pop_member(object):

    def __init__(self, **kwargs):
        
        ID = "NA"
        mother_ID = "NA"
        age = random.choices(range(0,101), weights= fv1)[0]
        kin_network = "NA"
        alive = 1
        
        if "alive" in kwargs: alive = kwargs["alive"]
        if "age" in kwargs: age = kwargs["age"]
        if "fert" in kwargs: fert = kwargs["fert"]
        if "mort" in kwargs: mort = kwargs["mort"]
        if "ID" in kwargs: ID = kwargs["ID"]
        if "mother_ID" in kwargs: mother_ID = kwargs["mother_ID"]
        if "kin_network" in kwargs: kin_network = kwargs["kin_network"]
        
        self.alive = alive
        self.age = age
        self.ID = ID
        self.mother_ID = mother_ID
        self.mort = mrates_f[mrates_f["Age"] == self.age].Mort[self.age]
        self.fert = frates_f[frates_f["Age"] == self.age].Fert[self.age]
        self.kin_network = kin_network

    def advance_age(self):
        
        self.age = self.age + 1
        ## These values change after aging...
        if self.age > 100:
            self.fert = 0
            self.mort = 1
            self.alive = 0
        else:
            self.mort = mrates_f[mrates_f["Age"] == self.age].Mort[self.age]
            self.fert = frates_f[frates_f["Age"] == self.age].Fert[self.age]
                

    
    def die(self):
        prob_die = np.random.binomial(1, 1-self.mort)
        if prob_die == 1:
            self.alive = 0   
            
        else:
            self.alive = 1

                    
    def extract_kin(self, time, population):
        
        self.kin_network = pd.DataFrame(np.zeros([4, 8]), 
                                        columns=['Year', 'Kin', 'Kin ID', 'Alive', 'Age Kin',  'Age Focal',  'Alive Focal', 'Focal ID'])
        
        time = time 

            
        ## count where to add rows in data frame using n_kin as counter
        n_kin = 0
        
        ################ Ancestors ###########################################
        
        ### Mother ............................................................
        if self.mother_ID in [i.ID for i in population]:
            mum = [i for i in population if i.ID == self.mother_ID][0]
            age = mum.age
            alive = mum.alive
            ID = mum.ID
            self.kin_network.loc[n_kin] = [time, "Mother" , ID , alive , age, self.age , self.alive, self.ID]
        else:
            self.kin_network.loc[n_kin] = [time, "Mother" , "dead ID" , 0 , "NA" , self.age, self.alive, self.ID]
            
        n_kin = n_kin + 1 # move to next row
        
        
        ################### Descendantds ###################################
        
        if self.ID in [i.mother_ID for i in population]:
            n_kin1 = 0
            alive = 1
            offspring = [i for i in population if i.mother_ID == self.ID]
            for j,i in enumerate(offspring):
                age = i.age
                alive = i.alive
                ID = i.ID
                self.kin_network.loc[j + n_kin] = [time, "Offspring" , ID, alive , age ,  self.age, self.alive, self.ID]
                n_kin1 = n_kin1 + 1 
        else:
            self.kin_network.loc[n_kin] = [time, "Offspring", "NA" , 0 ,"NA" , self.age, self.alive, self.ID]
            n_kin1 = 1

        n_kin = n_kin + n_kin1
        
        ############### Older Siblings ###############################################
        
        if len([i for i in population if i.mother_ID == self.mother_ID and i.mother_ID != "NA" and i.ID != self.ID and i.age > self.age]) != 0:
            n_kin1 = 0
            siblings_older = [i for i in population if i.mother_ID == self.mother_ID and i.ID != self.ID and i.age > self.age]
            for j,i in enumerate(siblings_older):
                age = i.age
                alive = i.alive
                ID = i.ID
                self.kin_network.loc[j + n_kin] = [time, "Older siblings" , ID , alive , age , self.age, self.alive, self.ID]
                n_kin1 = n_kin1 + 1 
        else:
            self.kin_network.loc[n_kin] = [time, "Older siblings" , "NA", 0 , "NA" , self.age, self.alive, self.ID]
            n_kin1 = 1
       
        n_kin = n_kin + n_kin1
        
        ################# Younger siblings #######################################
            
        if len([i for i in population if i.mother_ID == self.mother_ID and i.mother_ID != "NA" and i.ID != self.ID and i.age < self.age]) != 0:
            n_kin1 = 0
            siblings_younger = [i for i in population if i.mother_ID == self.mother_ID and i.ID != self.ID and i.age < self.age]
            for j,i in enumerate(siblings_younger):
                age = i.age
                alive = i.alive
                ID = i.ID
                self.kin_network.loc[j + n_kin] = [time, "Younger siblings" , ID , alive , age , self.age,  self.alive, self.ID]
                n_kin1 = n_kin1 + 1 
        else:
            self.kin_network.loc[n_kin] = [time, "Younger siblings" , "NA", 0 , "NA" , self.age, self.alive, self.ID,]
            n_kin = n_kin + 1
            
            
        return self.kin_network



def sex(parent, maxID): 
    #parents is a list of 2
    baby = pop_member(age = 0, ID = maxID)
    baby.mother_ID = parent.ID

        
    return baby
        
        
        
        
        