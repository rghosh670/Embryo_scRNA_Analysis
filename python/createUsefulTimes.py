import csv
import pandas as pd

times = open("./Embryo_scRNA/data/times.txt", "r") # read in 96000 times
times_reader = csv.reader(times, delimiter = ' ')

f = open('./Embryo_scRNA/data/indicesOfUsefulCells.txt', 'r') # read in indices of 37000 cells
indexList = f.read().splitlines()
f.close()

timeList = []

k = 0 # should map the 96000 times to the 37000 cells
for i in times_reader:
    for j in i:
        if str(k) in indexList:
            timeList.append(j)
        k = k + 1


with open('./Embryo_scRNA/data/usefulTimes.txt', 'w') as f: # write out results to text file
    for item in timeList:
        f.write("%s\n" % item)

f.close()






