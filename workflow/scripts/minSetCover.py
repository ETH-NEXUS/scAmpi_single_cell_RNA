# minimum set cover algorithm to identify therapy combinations to treat all subclones
# roche tumor profiler single cell, Franziska Singer March  2018
# greedy algorithm based on following post: https://stackoverflow.com/questions/7942312/how-do-i-make-my-implementation-of-greedy-set-cover-faster

import sys
import argparse

'''
function definitions
'''
# find the drug with minimum cost to append to cluster set cover
def findMin(sets,tempUV,weights):
    minCost = sys.float_info.max
    minElement = -1

    for i,s in enumerate(sets):
        try:
            cost = (weights[i])/(len(s.intersection(tempUV)))
            if cost < minCost:
                minCost = cost
                minElement = i
        except:
            # intersection with tempUV empty, ignore
            pass

    return (minElement, minCost)

# find the drug with minimum weight to append to cluster set cover, weight dependent on the size of clusters covered by the drug
def findMinPerc(sets,tempUV,weights,dictPerc):
    if len(dictPerc.keys()) == 0:
        print("Error. No tumor clusters identified. Rerun min set cover without option considerClusterSize.")
        sys.exit(1)

    minCost = sys.float_info.max
    minElement = -1

    for i,s in enumerate(sets):
        covered = s.intersection(tempUV)
        sumBenefit = 0.0
        for element in covered:
           # benefit = float(element) * float(dictPerc[element])
            sumBenefit = sumBenefit + float(dictPerc[element])
        try:
            cost = (weights[i])/sumBenefit
            if cost < minCost:
                minCost = cost
                minElement = i
                #print(sumBenefit)
        except:
            # intersection with tempUV empty, ignore
            pass

    return (minElement, minCost)



parser = argparse.ArgumentParser(description='Minimum set cover over a set of drugs and clusters.')
parser.add_argument('--input', dest='inTable', required=True, help='Tab separated input table with entries of form drug,clusters,weight.')
parser.add_argument('--outFile', dest='outName', required=True, help='Name of the output file.')
parser.add_argument('--percentageTable', dest='percTable', required=True, help='Tab separated input table with entries of form clusterID,number_cells_in_cluster,total_number_malignant_cells,percentage_of_cells_in_cluster.')
parser.add_argument('--considerClusterSize', dest='considerClusterSize', required=True, help='If set "yes" the cluster size is considered during calculation of minSetCover, if set "no" each covered cluster has the same value')

args = parser.parse_args()

# read in file with percentage of cells in each cluster
percTable = open(args.percTable, 'r')
percHeader = percTable.readline()

dictPerc = {}

for line in percTable:
    lineSplit = line.strip().split('\t')
    clusterID = lineSplit[0]
    percentageCells = lineSplit[3]
    dictPerc[clusterID] = percentageCells

percTable.close()
#sumPercentage = 0


infile = open(args.inTable,'r')
header = infile.readline()

# check what happens if multiple equal sets
universeArr = []  # fill first, then covert to set; contains cluster ids
setArr = [] # contains drug names - cluster IDs?
weights = []  

setToDrug = {}  # set index to drug

for line in infile:
    #print line
    lineArr = line.strip().split('\t')
    drug = lineArr[0]
    clusters = lineArr[1].split(',')
    weight = float(lineArr[2])

    for cluster in clusters:
        if cluster not in universeArr:
            universeArr.append(cluster)

    #myClustList = ','.join(str(x) for x in clusters)
    setArr.append(set(clusters))
    weights.append(weight)
    setToDrug[len(setArr)-1] = drug

infile.close()

uv = set(universeArr)
tempUV = uv
#print(setArr,weights,uv,setToDrug)

setCover = []
setCover_drugs = []
costs = []

outfile = open(args.outName,'w')
outfile.write("Drugs\tTargeted_clusters\tTargeted_percentage_of_malignant_cells\n")

# while loop ok, because: (i) it is ensured that eventually tempUV will be empty because the union of setArr will always yield the whole universe
# (ii) we do not have to remove a selected set from setarr, because it will not be chosen again (intersection with tempUV is empty)

if args.considerClusterSize == 'yes':
    while len(tempUV) != 0:
        print("Consider cluster size is set to yes.")
        (setPos, minCost) = findMinPerc(setArr,tempUV,weights,dictPerc)
        setCover.append(setArr[setPos])
        setCover_drugs.append(setToDrug[setPos])
        tempUV = tempUV.difference(setArr[setPos])
        costs.append(minCost)
        # additionally give out percentage of cells that is targeted with given set
        sumPercentage = 0
        for element in list(uv.difference(tempUV)):
            if element in dictPerc.keys():
                sumPercentage = float(sumPercentage) + float(dictPerc[element])
        # write outfile
        outfile.write("%s\t%s\t%.2f\n" %(','.join(list(setCover_drugs)),','.join(list(uv.difference(tempUV))),sumPercentage))

elif args.considerClusterSize == 'no':
    while len(tempUV) != 0:
        print("Consider cluster size is set to no.")
        (setPos, minCost) = findMin(setArr,tempUV,weights)
        setCover.append(setArr[setPos])
        setCover_drugs.append(setToDrug[setPos])
        tempUV = tempUV.difference(setArr[setPos])
        costs.append(minCost)
        # additionally give out percentage of cells that is targeted with given set
        sumPercentage = 0
        for element in list(uv.difference(tempUV)):
            if element in dictPerc.keys():
                sumPercentage = float(sumPercentage) + float(dictPerc[element])
        # write outfile
        if (sumPercentage == 0) and (len(dictPerc.keys()) == 0):
            sumPercentage = "n.a"
        outfile.write("%s\t%s\t%.2f\n" %(','.join(list(setCover_drugs)),','.join(list(uv.difference(tempUV))),sumPercentage))

else:
    print("Please indicate with 'yes' or 'no' whether to respect size of clusters or not.")

outfile.close()

print("Set cover: ", setCover)
print("Set cover drugs: ", setCover_drugs)
print("Total costs: ", sum(costs))
