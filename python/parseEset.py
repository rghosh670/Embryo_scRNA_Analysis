import csv
import pandas as pd


f = open('./Embryo_scRNA/data/indicesOfUsefulCells.txt', 'r')
indexList = f.read().splitlines()
f.close()

indexList = [int(x.strip())+1 for x in indexList]

eset = pd.read_csv("./Embryo_scRNA/data/eset.csv", usecols=indexList)
eset.to_csv('./Embryo_scRNA/data/parsedEset.csv')