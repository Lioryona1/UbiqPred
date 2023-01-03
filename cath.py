import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components


def listCreation(filename):
    """
    :param filename: PSSM file
    :return: tuple(namesList,sizesList)
    namesList = list of all the chains's name in the file
    sizesList = list of all the chains's number of amino acids in the file
    """
    namesList = []
    sizesList = []
    file1 = open(filename, 'r')
    line = file1.readline().split()
    cnt = 0
    while len(line) != 0:  # end of file
        cnt += 1
        if len(line) == 1:  # in chain header line
            sizesList.append(cnt)
            namesList.append(line[0][1:5] + line[0][-1])
            cnt = -1
        line = file1.readline().split()
    sizesList.append(cnt)
    sizesList = sizesList[1:]  # first one is redundent
    file1.close()
    return namesList, sizesList


def make_cath_df(filename, columns_number):
    """
    :param filename: cath-domain-list file
    :param columns_number: the number of columns to consider with the cath classification not include the cath domain name
    :return: dataframe of all the chains in the file and their cath classification divide to 4 different columns
    """

    df = pd.read_csv(filename, skiprows=16, header=None, delimiter=r"\s+")
    df = df.iloc[:, 0:columns_number + 1]
    cath_columns = ["n" + str(i) for i in range(1, columns_number + 1)]
    df.columns = ['chain'] + cath_columns
    df['chain'] = df['chain'].apply(lambda x: x[0:5])
    return df


def neighbor_mat(df, lst, columns_number):
    """
    :param df: cath data frame as it return from the func make_cath_df
    :param lst: list of chains
    :param columns_number: the number of columns to consider with the cath classification not include the cath domain name
    :return: matrix. mat[i][j] == 1 if there is connection between chain i and chain j
    """
    cath_columns = ["n" + str(i) for i in range(1, columns_number + 1)]
    n = len(lst)
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            similarity = df[df['chain'].isin([lst[i], lst[j]])].groupby(by=cath_columns)
            for name, group in similarity:
                if len(group['chain'].unique()) == 2:
                    mat[i][j] = mat[j][i] = 1
                    break
    return mat


def createRelatedChainslist(numberOfComponents, labels):
    """
    :param numberOfComponents: number of component = x => 0<=label values<x
    :param labels: labels
    :return: RelatedChainslist: RelatedChainslist[i] = list of all the chain index's which has the label i
    """
    relatedChainsLists = [[] for _ in range(numberOfComponents)]
    for i in range(len(labels)):
        relatedChainsLists[labels[i]].append(i)
    return relatedChainsLists


def createClusterSizesList(relatedChainslists, sizeList):
    """
    :param relatedChainslists:  relatedChainslist[i] = list of all the chain index's which has the label i
    :param relatedChainslists:  sizeList- list of all the chains's size
    :return: list of tuples (clusterIndex,size)
    """
    clusterSizes = []
    for i in range(len(relatedChainslists)):
        my_sum = 0
        for index in relatedChainslists[i]:
            my_sum += sizeList[index]
        clusterSizes.append((i, my_sum))
    return clusterSizes


def divideClusters(clusterSizes):
    """
    :param clusterSizes: list of tuples (clusterIndex,size)
    :return:  sublists,sublistsSum
    divide the list into 5 sublists such that the sum of each cluster sizes in the sublist is as close as possible
    """
    sublists = [[] for i in range(5)]
    sublistsSum = [0 for i in range(5)]
    clusterSizes.sort(reverse=True, key=lambda x: x[1])  # Sort the clusters by size descending order.
    for tup in clusterSizes:
        min_cluster_index = sublistsSum.index(min(sublistsSum))  # find the cluster with the minimal sum
        sublistsSum[min_cluster_index] += tup[1]
        sublists[min_cluster_index].append(tup[0])
    return sublists, sublistsSum


def clusterToChainList(clusterId, relatedChainsLists, nameList):
    """
    :param clusterid: list of chain indexs
    :param relatedChainslists: relatedChainslist[i] = list of all the chain index's in cluster i
    :param  namesList = list of all the chains's name in the file
    :return: chainList = list of all chain names in the cluster
    """
    cluster = relatedChainsLists[clusterId]  # get the chains in the cluster
    chainList = [nameList[i] for i in cluster]
    return chainList


def sublistsToChainLists(sublists, relatedChainsLists, nameList):
    """
    :param sublists: sublists[i] = all the clusters in sublist i
    :param relatedChainsLists: relatedChainslist[i] = list of all the chain index's in cluster i
    :return: chainLists: ChainLists[i] = list of all the chains in cluster i
    """
    chainLists = [[] for i in range(len(sublists))]
    for i in range(len(sublists)):
        for clusterId in sublists[i]:
            chainLists[i] += clusterToChainList(clusterId, relatedChainsLists, nameList)
    return chainLists


def chainListsToChainIndexDict(chainLists):
    """
    :param chainLists: chainLists[i] = list of all the chains in cluster i
    :return: chainDict: chainDict[chainName] = index of chain cluster(i if chain in ChainLists[i])
    """
    chainDict = {}
    for i in range(len(chainLists)):
        for chainName in chainLists[i]:
            chainDict[chainName] = i
    return chainDict


def dividePSSM(chainDict):
    """
    :param chainDict: chainDict[chainName] = index of chain cluster(i if chain in ChainLists[i])
    create len(chainLists) txt files. the i txt file contains the chains in chainLists[i]
    """
    filesList = [open("PSSM{}.txt".format(i), 'w') for i in range(5)]
    pssmFile = open("PSSM.txt", 'r')
    lines = pssmFile.readlines()
    fillIndex = -1  # fillIndex = i -> we now write to PSSMi.txt
    for line in lines:
        if line[0] == '>':  # header line
            fillIndex = chainDict[line[1:5] + line[-2]]
        filesList[fillIndex].write(line)
    for i in range(5):
        filesList[i].close()
    pssmFile.close()


cath_df = make_cath_df("cath-domain-list.txt", 4)
nameList, sizeList = listCreation("PSSM.txt")
print("Done1")
matHomologous = neighbor_mat(cath_df, nameList, 4)
graphHomologous = csr_matrix(matHomologous)
topology_components, topologyLabels = connected_components(csgraph=graphHomologous, directed=False,
                                                           return_labels=True)

print(nameList)
print(sizeList)
print(sum(sizeList))
print(topology_components)
print(topologyLabels)
print("Done2")
relatedChainsLists = createRelatedChainslist(topology_components, topologyLabels)
print("Done3")
clusterSizes = createClusterSizesList(relatedChainsLists, sizeList)
print("Done4")
sublists, sublistsSum = divideClusters(clusterSizes)
print("Done5")

print(relatedChainsLists)
print(clusterSizes)
print(sublists)
print(sublistsSum)


# nameList = ['1nbfA', '1nbfB', '1nbfE', '1p3qQ', '1p3qR', '1s1qA', '1s1qC', '1uzxA', '1wr6A', '1wr6B', '1wr6C', '1wr6D',
#             '1wrdA', '1xd3A', '1xd3C', '1yd8G', '1yd8H', '2ayoA', '2c7mA', '2d3gP', '2dx5A', '2fifB', '2fifD', '2fifF',
#             '2g45A', '2g45D', '2gmiA', '2gmiB', '2hd5A', '2hthB', '2ibiA', '2j7qA', '2j7qC', '2oobA', '2qhoB', '2qhoD',
#             '2qhoF', '2qhoH', '2wdtA', '2wdtC', '2wwzC', '2xbbA', '2xbbB', '3a33A', '3a9kC', '3by4A', '3c0rA', '3c0rC',
#             '3cmmA', '3cmmC', '3i3tA', '3i3tC', '3i3tE', '3i3tG', '3ifwA', '3ihpA', '3ihpB', '3jsvC', '3jsvD', '3jvzA',
#             '3jvzB', '3jvzC', '3jvzD', '3k9pA', '3kvfA', '3kw5A', '3ldzA', '3ldzD', '3ldzB', '3ldzC', '3mhsA', '3mhsB',
#             '3mhsC', '3mhsE', '3mtnA', '3mtnC', '3nheA', '3o65A', '3o65C', '3o65E', '3o65G', '3ofiA', '3ofiB', '3oj3I',
#             '3oj3J', '3oj3K', '3oj3L', '3oj3M', '3oj3N', '3oj3O', '3oj3P', '3olmA', '3phwA', '3phwC', '3phwE', '3phwG',
#             '3prmA', '3prmC', '3pt2A', '3ptfA', '3ptfB', '3tblA', '3tblB', '3tblC', '3tmpA', '3tmpC', '3tmpE', '3tmpG',
#             '3vhtA', '3vhtB']
#
# relatedChainsLists = [[0, 1, 2], [3], [4], [5, 6], [7], [8, 9, 10, 11, 15, 16], [12], [13, 14, 54, 64, 65], [17],
#                       [18, 21, 22, 23], [19], [20, 29], [24, 25], [26, 43, 59, 60, 63, 99, 100], [27],
#                       [28, 30, 50, 51, 52, 53, 74, 75, 76], [31, 32], [33], [34, 35, 36, 37], [38, 39], [40],
#                       [41, 42, 61, 62, 91], [44], [45, 46, 47], [48, 49], [55, 56], [57, 58], [66, 67, 68, 69], [70],
#                       [71], [72], [73], [77, 78, 79, 80], [81, 82], [83], [84], [85], [86], [87], [88], [89], [90],
#                       [92, 93, 94, 95, 96, 97, 98], [101, 102, 103], [104, 105, 106, 107], [108, 109]]
#
# sublists = [[15, 32, 45, 9, 4, 17, 34, 38, 10], [24, 0, 5, 16, 12, 26, 2, 20, 41], [21, 13, 27, 28, 3, 6, 31, 36, 39],
#             [33, 42, 44, 19, 11, 30, 29, 37, 22], [25, 7, 43, 23, 8, 18, 14, 1, 40, 35]]

chainLists = sublistsToChainLists(sublists, relatedChainsLists, nameList)
chainDict = chainListsToChainIndexDict(chainLists)
print(chainLists)
print(chainDict)
dividePSSM(chainDict)


