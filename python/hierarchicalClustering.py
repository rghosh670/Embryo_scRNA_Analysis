from scipy.cluster.hierarchy import fcluster
from scipy import stats
import scipy.cluster.hierarchy as hac
import statistics
from math import floor
import random
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
import json
import matplotlib
matplotlib.use('wxAgg')


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
        nbrs = NearestNeighbors(
            n_neighbors=maxK, algorithm='ball_tree', metric=my_metric).fit(matDict[key])
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


with open('./Embryo_scRNA/data/finalDict.txt') as json_file:
    finalDict = json.load(json_file)

json_file.close()
myDict = {}

timeList = []
expressionList = []
cellList = []

for time, v in list(finalDict.items()):  # v is dict with cellName as key
    for cellName, v1 in list(v.items()):  # v1 is dict with geneName as key
        for geneName, expression in list(v1.items()):
            if geneName == 'rab-3':
                timeList.append(time)
                expressionList.append(expression)
                cellList.append(cellName)


cellListWithoutDuplicates = list(dict.fromkeys(cellList))
timeListWithoutDuplicates = list(dict.fromkeys(timeList))

df = pd.DataFrame(index=cellListWithoutDuplicates,
                  columns=timeListWithoutDuplicates)
df = df.fillna(1e-10)


# expressionList = knnSearchFunction(timeList, expressionList, cellList)
f = open('/Users/rohitghosh/Downloads/rab-3CorrectedExpression.txt', 'r') # generated using the function commented above
expressionList = f.read().splitlines()
expressionList = [float(i) for i in expressionList]
f.close()


def plot(indicesOfCellsToPlot, label = False):
    colorDict = {}
    for i in cellList:
        if i not in colorDict:
            color = [random.random(), random.random(), random.random()]
            colorDict[i] = color

    plt.title('rab-3', fontsize=20)
    plt.xlabel('Embryo Time', fontsize=15)
    plt.ylabel('Expression Level', fontsize=15)

    cells = [cellListWithoutDuplicates[i] for i in indicesOfCellsToPlot]

    for i in range(len(timeList)):
        if cellList[i] in cells:
            plt.scatter(timeList[i], expressionList[i], marker='o',
                        color=colorDict.get(cellList[i]), s=15)
            if label or expressionList[i] > 3.5:
                plt.annotate(xy=(timeList[i], expressionList[i]), s=cellList[i])


    mng = plt.get_current_fig_manager()
    mng.frame.Maximize(True)
    plt.ylim(bottom=0)
    plt.show()
    return
# indicesOfCellsToPlot = range(len(cellListWithoutDuplicates))
# plot(indicesOfCellsToPlot)

def dendogram(plot):
    for i in range(len(cellList)):
        df.at[cellList[i], timeList[i]] = expressionList[i]

    def my_metric(x, y):
        r = stats.pearsonr(x, y)[0]
        return 1 - r  # correlation to distance: range 0 to 2


    Z = hac.linkage(df,  method='single', metric=my_metric)

    plt.figure(figsize=(12.5, 5))
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('sample index')
    plt.ylabel('distance')
    hac.dendrogram(
        Z,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=8.,  # font size for the x axis labels
    )
    if plot:
        plt.show()
    return Z
Z = dendogram(True)

def print_clusters(timeSeries, Z, k, plot=False):
    # k Number of clusters I'd like to extract
    results = fcluster(Z, k, criterion='maxclust')

    # check the results
    s = pd.Series(results)
    clusters = s.unique()

    for c in clusters:
        cluster_indeces = s[s == c].index
        print(cluster_indeces)
        if plot:
            timeSeries.T.iloc[:, cluster_indeces].plot()
            plt.show()
    return
print_clusters(df, Z, 60, plot=False)
