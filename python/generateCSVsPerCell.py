import csv
import pandas as pd

# I never actually ended up reading any of the python documentation so I don't really know how methods work
# I'll be honest, I'm not sure if parameters are passed by value, or how objects work, but this seemed to be fine

def processIndices(csvRow): # Reading in the indices per cell type csv file required some processing
    myList = csvRow.tolist()
    myList.pop(0)
    myList = myList[0]
    myList = myList[1:-1]
    myList = list(myList.split(", "))

    for i in range(len(myList)):
        myList[i] = (int(myList[i].strip()))

    return myList

def readCSV(filePath, indices): # I had to manipulate the indices list in order to read the eset properly
    indices = [x+1 for x in indices]
    indices.insert(0, 0)
    eset = pd.read_csv(filePath, usecols = indices)

    indices.pop(0)
    indices = [x-1 for x in indices]
    return eset

def addRow(df, myList): # copied this from stack overflow to add embryo times as the top row of eset, is removed later
    myList.insert(0,-1)
    df.loc[-1] = myList
    df.index = df.index + 1
    df.sort_index(inplace=True)
    myList.pop(0)
    return df

def sortEset(df): # sorts the eset based on the first row - in this case should always be embryo times
    new_columns = df.columns[df.iloc[df.first_valid_index()].argsort()]
    df = df[new_columns]
    return df


def generateCSVs(filePath, indices, timesList):
    eset = readCSV('./Embryo_scRNA/data/esetWithoutUnexpressedGenes.csv', indices)
    eset = addRow(eset, timesList)

    eset = sortEset(eset)
    eset.set_index('gene', inplace = True)

    eset.to_csv(filePath)
    return None

for i in range(0, 138): # I'll be honest, I shouldn't have added this - but it worked I guess. Really shouldn't be hardcoding.
    indexCSV = pd.read_csv('./Embryo_scRNA/data/indicesPerCellType.csv', skiprows=lambda x: x not in [i], header = None)

    f = open('./Embryo_scRNA/data/usefulTimes.txt', 'r')
    timesList = f.read().splitlines()
    f.close()

    for index, row in indexCSV.iterrows():
        filePath = './Embryo_scRNA/CSVsPerCellType/' + row[0] +'.csv'
        indices = processIndices(row)
        timesList = [int(i) for j, i in enumerate(timesList) if j in indices]
        generateCSVs(filePath, indices, timesList)



