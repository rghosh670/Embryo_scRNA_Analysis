import csv
import pandas as pd

numCellsGeneExpression = pd.read_csv("./Embryo_scRNA/data/numCellsGeneExpression.csv", header = None)
zeroList = []

for i in range(len(numCellsGeneExpression)):
    if numCellsGeneExpression.iloc[i, 1] == 0:
        zeroList.append(i + 1)


parsedEsetByGene = pd.read_csv("./Embryo_scRNA/data/parsedEsetByGene.csv", skiprows = zeroList)
parsedEsetByGene.to_csv('./Embryo_scRNA/data/esetWithoutUnexpressedGenes.csv', index = False)

