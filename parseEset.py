import csv
import pandas as pd

unique_names = pd.read_csv("./Embryo_scRNA/data/UniqueNames.csv")
u_names = unique_names['X']

uniqueNamesList = []
for i in u_names:
    uniqueNamesList.append(i.strip())

names = open("./Embryo_scRNA/data/names.txt", "r")

namesList = []

for i in names:
    namesList.append(i.strip())

names.close()

indexList = []

for i in range(len(namesList)):
    if namesList[i] in uniqueNamesList:
        indexList.append(i)

eset = pd.read_csv("./Embryo_scRNA/data/eset.csv", chunksize = 5500)

dfList = []

for chunk in eset:
    df = pd.DataFrame()

    for i in indexList:
        df[chunk.columns[i]] = chunk.iloc[:, i]

    dfList.append(df)

finalOutput = pd.DataFrame()

for i in dfList:
    finalOutput = pd.DataFrame.append(finalOutput, i, ignore_index = True)

finalOutput.to_csv('./Embryo_scRNA/parsedEset.csv', index = False)
