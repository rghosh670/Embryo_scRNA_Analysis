import csv
import pandas as pd

simplemine = open("./Embryo_scRNA/data/simplemine_results.txt", "r")
simplemine_reader = csv.reader(simplemine, delimiter = '\t')

idList = []

for s in simplemine_reader:
    try:
        if s[1][0] == 'W':
            idList.append(s[1])
    except IndexError:
        pass

esetIDs = open("./Embryo_scRNA/data/IDs.txt", "r")

indexList = [0]

for i in enumerate(esetIDs):
    if i[1].strip() in idList:
        indexList.append(i[0])

parsedEsetByGene = pd.read_csv("./Embryo_scRNA/data/parsedEset.csv", skiprows=lambda x: x not in indexList, usecols = lambda x: x not in indexList)
parsedEsetByGene.drop(['Unnamed: 0', 'Unnamed: 0.1'], inplace = True, axis = 1)
parsedEsetByGene.to_csv('./Embryo_scRNA/data/parsedEsetByGene.csv', index = False)