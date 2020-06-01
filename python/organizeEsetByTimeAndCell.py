import csv
import pandas as pd
from math import floor
from statistics import median
from collections import OrderedDict
import json

def switch(argument):
    switcher = {
        'ADA_ADE_grandparent': 'ADE',
        'ADL_parent': 'ADL',
        'AIN_parent': 'AIN',
        'AMso_parent': 'AMso',
        'ASE_parent': 'ASE',
        'ASI_parent': 'ASI',
        'ASK_parent': 'ASK',
        'AVH_parent': 'AVH',
        'DVC_parent': 'DVC',
        'IL1_parent': 'IL1',
        'M1_parent': 'M1_parent',
        'M4_parent': 'M4_parent',
        'NSM_parent': 'NSM_parent',
        'Neuroblast_ADE_ADA': 'ADE',
        'Neuroblast_ADF_AWB': ('ADF', 'AWB'),
        'Neuroblast_AFD_RMD': ('AFD', 'RMD'),
        'Neuroblast_AIZ_FLP': ('AIZ', 'FLP'),
        'Neuroblast_AIZ_FLP_RMG': ('AIZ', 'FLP'),
        'Neuroblast_ALA_RMED': ('ALA', 'RMED'),
        'Neuroblast_ALM_BDU': 'ALM_BDU',
        'Neuroblast_ASE_ASJ_AUA': ('ASE', 'ASJ', 'AUA'),
        'Neuroblast_ASG_AWA': ('ASG', 'AWA'),
        'Neuroblast_ASH_RIB': ('ASH', 'RIB'),
        'Neuroblast_ASJ_AUA': ('ASJ', 'AUA'),
        'Neuroblast_AVG_RIR': 'AVG',
        'Neuroblast_AWC_SAAVx': 'AWC',
        'Neuroblast_BAG_SMDVx': 'BAG',
        'Neuroblast_HSN_PHB': 'PHB_and_possibly_PHA',
        'Neuroblast_I6_M5': 'Neuroblast_I6_M5',
        'Neuroblast_IL1_IL2': ('IL1', 'IL2'),
        'Neuroblast_M2_M3': 'Neuroblast_M2_M3',
        'Neuroblast_PVC_LUA': 'Neuroblast_PVC_LUA',
        'Neuroblast_URX_CEPDx': 'URX',
        'OLL_parent': 'OLL',
        'OLQ_grandparent': 'OLQ',
        'OLQ_parent': 'OLQ',
        'PLM_ALN_grandparent': ('PLM', 'ALN'),
        'PLM_ALN_great_grandparent': ('PLM', 'ALN'),
        'PVQ_parent': 'PVQ_and_possibly_PVC',
        'Parent_of_AMsh_URB': ('AMsh', 'URB_and_possibly_URA'),
        'Parent_of_MI_pm1DR': 'Parent_of_MI_pm1DR',
        'Parent_of_PVP_and_rect_V': 'PVP',
        'Parents_of_PHsh_hyp8_hyp9': 'Parents_of_PHsh_hyp8_hyp9',
        'Parents_of_Y_DA6_DA7_DA9': 'Parents_of_Y_DA6_DA7_DA9',
        'RIA_parent': 'RIA',
        'RIC_parent': 'RIC',
        'RID_parent': 'RID',
        'RIM_parent': 'RIM',
        'RME_LR_parent': 'RME',
        'URA_parent': 'URB_and_possibly_URA'
    }

    return (switcher.get(argument, argument))


f = open('./Embryo_scRNA/data/namesOfNeuronalCells.txt', 'r')
namesList = f.read().splitlines()
f.close()

f = open('./Embryo_scRNA/data/usefulTimes.txt', 'r')
timesList = list(map(int, f.read().splitlines()))
f.close()

for i in range(len(timesList)):
    if timesList[i] % 10 != 0:
        timesList[i] = floor(timesList[i] / 10) * 10

f = open('./Embryo_scRNA/data/listOfGenes.txt', 'r')
genesList = f.read().splitlines()
f.close()

newNamesList = []
for i in namesList:
    newNamesList.append(switch(i.strip()))

eset = pd.read_csv('./Embryo_scRNA/data/esetWithoutUnexpressedGenes.csv')

finalDict = {}

for index, row in eset.iterrows():
    row = row.tolist()
    gene = row.pop(0).strip()

    for i, j in enumerate(row):
        geneDict = finalDict.setdefault(timesList[i], {})
        cellDict = geneDict.setdefault(newNamesList[i], {})
        geneList = cellDict.setdefault(genesList[index], [])
        geneList.append(j)

for time, v in list(finalDict.items()):
    for cellName, v1 in list(v.items()):
        if type(cellName) == tuple:
            for cellElement in cellName:
                if cellElement.strip() in v:
                    singleCell = v.get(cellElement.strip())
                    for kx, vx in v1.items():
                        singleCell[kx] = singleCell[kx] + vx
                else:
                    v[cellElement.strip()] = dict(v1)

for time, v in list(finalDict.items()):
    for cellName, v1 in list(v.items()):
        for geneName, v2 in list(v1.items()):
            v2 = [i for i in v2 if i != 0]
            medianValue = median(v2) if v2 else 0.0
            finalDict[time][cellName][geneName] = medianValue

for time, v in list(finalDict.items()):
    for cellName, v1 in list(v.items()):
        if type(cellName) == tuple:
            del(finalDict[time][cellName])

finalOrderedDict = OrderedDict(sorted(finalDict.items()))

with open('./Embryo_scRNA/data/finalDict.txt', 'w') as outfile:
    json.dump(finalOrderedDict, outfile, sort_keys=False, indent = 2)

outfile.close()