import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components


# class relatedChains:
#     def __init__(self,clusterIndex):
#         self.clusterIndex = clusterIndex
#         self.chainsIndexList = []
#
#     def addChainToRelatedChains(self):
# def listCreation(filename):
#     """
#     :param filename: PSSM file
#     :return: list of all the chains in the file
#     """
#     chains_list = []
#     file1 = open(filename, 'r')
#     line = file1.readline().split()
#     while len(line) != 0:
#         if len(line) == 1:
#             chains_list.append(line[0][1:5] + line[0][-1])
#         line = file1.readline().split()
#     file1.close()
#     return chains_list

def listCreation(filename):
    """
    :param filename: PSSM file
    :return: tuple(namesList,sizesList)
    listNames = list of all the chains's name in the file
    listSizes = list of all the chains's number of amino acids in the file
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
    df = pd.read_fwf(filename, delimiter='', skiprows=17, header=None)
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


cath_df = make_cath_df("cath-domain-list.txt", 3)
nameList, sizeList = listCreation("PSSM.txt")
print("Done1")
matTopology = neighbor_mat(cath_df, nameList, 3)
graphTopology = csr_matrix(matTopology)
topology_components, topologyLabels = connected_components(csgraph=graphTopology, directed=False,
                                                           return_labels=True)

print(nameList)
print(sizeList)
print(sum(sizeList))
print(topology_components)
print(topologyLabels)
print("Done2")
relatedChainslists = createRelatedChainslist(topology_components, topologyLabels)
print("Done3")
clusterSizes = createClusterSizesList(relatedChainslists,sizeList)
print("Done4")
sublists, sublistsSum = divideClusters(clusterSizes)
print("Done5")

print(relatedChainslists)
print(clusterSizes)
print(sublists)
print(sublistsSum)
