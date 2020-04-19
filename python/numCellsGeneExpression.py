import csv
import pandas as pd

numCellsExpressedPerGene = {}

with open('./Embryo_scRNA/data/parsedEsetByGene.csv', newline='') as f:
    reader = csv.reader(f)
    data = list(reader)

f.close()


for i in range(1, len(data)):
    numCells = 0
    for j in range(1, len(data[i])):
        if float(data[i][j]) > 0.0:
            numCells = numCells + 1
    numCellsExpressedPerGene.update({data[i][0]: numCells})

w = csv.writer(open("./Embryo_scRNA/data/numCellsGeneExpression.csv", "w"))
for key, val in numCellsExpressedPerGene.items():
    w.writerow([key, val])
