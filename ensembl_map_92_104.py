#!/bin/python

########################################################################################################################
# Elizabeth Hutchins, June 2022
# input: ensembl ID converter output file, lncipedia to ensembl 92 conversion table
# output: conversion table from lncipedia/ensembl 92 to ensembl 104 (matches gencode 40)
########################################################################################################################

import pandas as pd

f1 = "ensembl92.idmapper.txt"
f2 = "lncipedia_5_2_ensembl_92_genes.txt"

#read in lncipedia conversion table to ensembl 92
idmapper = pd.read_csv(f1, sep = ",", header = 0)
conv = pd.read_csv(f2, sep = "\t", header = 0)

#remove rows with repeats of the header
idmapper = idmapper[idmapper[" Release"].str.contains("Release") == False]
#remove rows > version 104 to match gencode
idmapper = idmapper[idmapper[" Release"].str.contains("105") == False]
idmapper = idmapper[idmapper[" Release"].str.contains("106") == False]
#make Release column numeric before we sort
idmapper[[' Release']] = idmapper[[' Release']].apply(pd.to_numeric)
#add column with base ensembl ID
idmapper['Old_base_ID'] = idmapper['Old stable ID'].str.replace('\.[0-9]', '').replace('\.[0-9]', '')
#sort by release column in descending order and drop duplicates
#this should result in the most recent ensembl ID for version 104 or older
idmapper = idmapper.sort_values(' Release', ascending = False).drop_duplicates(['Old_base_ID', ' New stable ID'])

#merge with lncipedia conversion table
mappings = pd.merge(conv, idmapper, left_on = 'ensemblGeneID', right_on='Old stable ID')

#output conversion between lncipedia 5.2 and ensembl 104
mappings.to_csv("ensembl104GeneUpdate.csv", sep = ",", header = True, index = False)