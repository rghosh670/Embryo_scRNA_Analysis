import csv
import pandas as pd

times = open("./Embryo_scRNA/data/times.txt", "r")
times_reader = csv.reader(times, delimiter = ' ')

f = open('./Embryo_scRNA/data/indicesOfUsefulCells.txt', 'r')
indexList = f.read().splitlines()
f.close()

timeList = []

k = 0
for i in times_reader:
    for j in i:
        if str(k) in indexList:
            timeList.append(j)
        k = k + 1


with open('./Embryo_scRNA/data/usefulTimes.txt', 'w') as f:
    for item in timeList:
        f.write("%s\n" % item)

f.close()






