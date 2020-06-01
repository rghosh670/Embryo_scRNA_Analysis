from sklearn.neighbors import NearestNeighbors
import matplotlib
import matplotlib.pyplot as plt
import json
import numpy as np
import pandas as pd
import csv
import statistics
import random
from math import floor

def knnSearchFunction(timeList, expressionList, cellList, r=30):
    length = range(len(timeList))

    def my_metric(arr1, arr2):
        return abs(arr1[0] - arr2[0])

    matDict = {}
    neighborsDict = {}

    for i in length:
        print(i)
        if cellList[i] not in matDict:
            matDict[cellList[i]] = np.zeros(shape=(len(timeList), 2))

        mat = matDict[cellList[i]]
        mat[i][0] = timeList[i]
        mat[i][1] = expressionList[i]

    newExpression = expressionList[:]

    maxK = 0
    for i in length:
        print(i)
        mat = matDict[cellList[i]]
        k = len(list(x for x in mat[:, 0] if mat
                    [i][0] - r <= x <= mat[i][0] + r))
        maxK = max(k, maxK)


    for key in matDict.keys():
        print(key)
        nbrs = NearestNeighbors(n_neighbors=maxK, algorithm='ball_tree', metric=my_metric).fit(matDict[key])
        neighborsDict[key] = nbrs.kneighbors(matDict[key])

    for i in length:
        print(i)
        mat = matDict[cellList[i]]
        k = len(list(x for x in mat[:, 0] if mat
                     [i][0] - r <= x <= mat[i][0] + r))

        distances, indices = neighborsDict[cellList[i]]
        newIndex = (indices[i])[0:k]

        expressionValues = [expressionList[i] for i in newIndex]
        expressionValues = [num for num in expressionValues if num]
        newExpression[i] = statistics.median(
            expressionValues) if expressionValues else 0.0

    return newExpression


# main routine

# userGene = input('Enter Gene Name: ').lower().strip()
# userCell = input('All cells (all)? If not, which cells? ').split()
# userLabel = input('Label? Y/N: ')
# userLabel = (userLabel.lower() == 'y')
# userExcludeZeroValues = input('Exclude Zero Values? Y/N: ')
# userExcludeZeroValues = (userExcludeZeroValues.lower() == 'y')

userGene = 'rab-3'
userCell = ['IL1', 'AVH']
userLabel = True
userExcludeZeroValues = False

matplotlib.use('wxAgg')
with open('./Embryo_scRNA/data/finalDict.txt') as json_file:
    finalDict = json.load(json_file)

json_file.close()

timeList = []
expressionList = []
cellList = []

for time, v in list(finalDict.items()):  # v is dict with cellName as key
    for cellName, v1 in list(v.items()):  # v1 is dict with geneName as key
        for geneName, expression in list(v1.items()):  # v2 is actual expression value
            if geneName == userGene:
                if cellName in userCell or 'all' in userCell:
                    timeList.append(time)
                    expressionList.append(expression)
                    cellList.append(cellName)


colorDict = {}
for i in cellList:
    if i not in colorDict:
        color = [random.random(), random.random(), random.random()]
        colorDict[i] = color

plt.title(userGene, fontsize=20)
plt.xlabel('Embryo Time', fontsize=15)
plt.ylabel('Expression Level', fontsize=15)

expressionList = knnSearchFunction(timeList, expressionList, cellList)

for i in range(len(timeList)):
    if not(expressionList[i] == 0 and userExcludeZeroValues):
        plt.scatter(timeList[i], expressionList[i], marker='o',
                    color=colorDict.get(cellList[i]), s=15)
        if userLabel:
            plt.annotate(xy=(timeList[i], expressionList[i]), s=cellList[i])


mng = plt.get_current_fig_manager()
mng.frame.Maximize(True)
plt.ylim(bottom=0)
plt.show()
