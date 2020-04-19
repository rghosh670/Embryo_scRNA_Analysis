
import csv
import pandas as pd

def processIndices(csvRow):
    myList = csvRow.tolist()
    myList.pop(0)
    myList = myList[0]
    myList = myList[1:-1]
    myList = list(myList.split(", "))

    for i in range(len(myList)):
        myList[i] = (int(myList[i].strip()))

    return myList

def readCSV(filePath, indices):
    indices = [x+1 for x in indices]
    indices.insert(0, 0)
    eset = pd.read_csv(filePath, usecols = indices)

    indices.pop(0)
    indices = [x-1 for x in indices]
    return eset

def addRow(df, myList):
    myList.insert(0,-1)
    df.loc[-1] = myList
    df.index = df.index + 1
    df.sort_index(inplace=True)
    myList.pop(0)
    return df

def sortEset(df):
    new_columns = df.columns[df.iloc[df.first_valid_index()].argsort()]
    df = df[new_columns]
    return df


def generateCSVs(filePath, indices, timesList):
    eset = readCSV('./Embryo_scRNA/data/esetWithoutUnexpressedGenes.csv', indices)
    eset = addRow(eset, timesList)

    eset = sortEset(eset)
    eset.set_index('gene', inplace = True)

    chunksize = int(((len(timesList) + 1) / 10))

    timeBinChunks = []
    timeBinTuples = []

    newDataFrame = pd.DataFrame(index=eset.index)
    newDataFrame.drop([-1], inplace = True)

    for i in range(chunksize,len(timesList), chunksize):
        previous = i - chunksize
        if i + chunksize > len(timesList):
            i = len(timesList)
            break
        esetChunk = eset.iloc[:, previous:i]
        timeBinTuple = (esetChunk.iloc[0,0],esetChunk.iloc[0,-1])
        timeBinTuples.append(timeBinTuple)
        timeBinChunks.append(esetChunk)
        esetChunk.drop([-1], inplace = True)


    for i in range(len(timeBinChunks)):
        zeroProportionList = []

        for index, row in timeBinChunks[i].iterrows():
            rowList = row.tolist()
            numZeroValues = 0
            for j in range(len(rowList)):
                if rowList[j] == 0:
                    numZeroValues = numZeroValues + 1


            zeroProportionList.append(float(numZeroValues) / len(rowList))

        newDataFrame[str(timeBinTuples[i])] = zeroProportionList

    newDataFrame.to_csv(filePath)
    return None

for i in range(0, 138): # I'll be honest, I shouldn't have added this - but it worked I guess
    indexCSV = pd.read_csv("/Embryo_scRNA/data/indicesPerCellType.csv", skiprows=lambda x: x not in [i], header = None)

    f = open('./Embryo_scRNA/data/usefulTimes.txt', 'r')
    timesList = f.read().splitlines()
    f.close()

    for index, row in indexCSV.iterrows():
        filePath = './Embryo_scRNA/zeroExpressionRatiosPerCellType/' + row[0] +'.csv'
        indices = processIndices(row)
        timesList = [int(i) for j, i in enumerate(timesList) if j in indices]
        generateCSVs(filePath, indices, timesList)



