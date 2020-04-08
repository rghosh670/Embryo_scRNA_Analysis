import csv
import pandas as pd
import matplotlib.pyplot as plt

f = open('./Embryo_scRNA/data/usefulTimes.txt', 'r')
timeList = f.read().splitlines()
f.close()

userows = list(range(0, 20))
esetWithoutUnexpressedGenes = pd.read_csv("./Embryo_scRNA/data/esetWithoutUnexpressedGenes.csv", skiprows=lambda x: x not in userows)

geneExpression = list(esetWithoutUnexpressedGenes.iloc[0, :])
print(geneExpression.pop(0))
nonZeroExpression = []
nonZeroTime = []

for i in range(len(geneExpression)):
    if geneExpression[i] != 0:
        nonZeroExpression.append(geneExpression[i])
        nonZeroTime.append(timeList[i])

nonZeroTime, nonZeroExpression = zip(*sorted(zip(nonZeroTime, nonZeroExpression)))

plt.scatter(nonZeroTime, nonZeroExpression, alpha = 0.25)
plt.show()