import csv
import pandas as pd

simplemine = open("./Embryo_scRNA/data/simplemine_results.txt", "r")
simplemine_reader = csv.reader(simplemine, delimiter = '\t')

idList = []

for s in simplemine_reader:
    try: # probably not the best way to avoid error caused by 'multiple genes' line
        if s[1][0] == 'W':
            idList.append(s[1])
    except IndexError:
        pass

esetIDs = open("./Embryo_scRNA/data/IDs.txt", "r")

indexList = [0] # to include header

for i in enumerate(esetIDs):
    if i[1].strip() in idList:
        indexList.append(i[0] + 1) # For some reason, have to add one to fix issue of header being treated as row

parsedEsetByGene = pd.read_csv("./Embryo_scRNA/data/parsedEset.csv", header = 0, skiprows=lambda x: x not in indexList)
parsedEsetByGene.drop(['Unnamed: 0', 'Unnamed: 0.1'], inplace = True, axis = 1)

parsedEsetByGene.to_csv('/Users/rohitghosh/Downloads/parsedEsetByGene.csv', index = False)