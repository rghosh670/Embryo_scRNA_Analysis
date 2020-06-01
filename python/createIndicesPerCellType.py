import csv
import pandas as pd

f = open('./Embryo_scRNA/data/indicesOfUsefulCells.txt', 'r') # read in list of "relevant" indices corresponding to the ~37000 cells
indexList = f.read().splitlines()
f.close()

f = open('./Embryo_scRNA/data/names.txt', 'r') # read in full list of ~96000 names
namesList = f.read().splitlines()
f.close()

indexList = [int(x.strip()) for x in indexList] # turn all indices to ints, probably not worth the extra for loop

myDict = {}

k = 0 # should take care of mapping the 96000 names to the 37000 cells
for i in indexList:
    myDict.setdefault(namesList[i], []).append(k) # make dictionary with values corresponding to lists of indices
    k = k + 1

w = csv.writer(open("./Embryo_scRNA/data/indicesPerCellType.csv", "w")) # print dictionary out to csv
for key, val in myDict.items():
    w.writerow([key, val])
