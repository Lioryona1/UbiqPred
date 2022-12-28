import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

def listCreation(filename):
    """
    :param filename: PSSM file
    :return: list of all the chains in the file
    """
    chains_list = []
    file1 = open(filename, 'r')
    line = file1.readline().split()
    while len(line) != 0:
        if len(line) == 1:
            chains_list.append(line[0][1:5] + line[0][-1])
        line = file1.readline().split()
    file1.close()
    return chains_list


def make_cath_df(filename):
    """
    param filename: cath-domain-list file
    return: dataframe of all the chains in the file and their cath classification divide to 4 different columns
    """
    df = pd.read_fwf(filename, delimiter='', skiprows=17, header=None)
    df = df.iloc[:, 0:5]
    df.columns = ['chain', 'n1', 'n2', 'n3', 'n4']
    df['chain'] = df['chain'].apply(lambda x: x[0:5])
    return df


def neighbor_mat(df, lst, homologous=True):
    """
    :param df: cath data frame as it return from the func make_cath_df
    :param lst: list of chains
    :param homologous: if true, make the matrix base on the homologous match, else, base on topology
    :return: matrix. mat[i][j] == 1 if there is homologous connection between chain i and chain j
    """
    n = len(lst)
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            if homologous:
                similarity = df[df['chain'].isin([lst[i], lst[j]])].groupby(by=['n1', 'n2', 'n3', 'n4'])
            else:  # topology
                similarity = df[df['chain'].isin([lst[i], lst[j]])].groupby(by=['n1', 'n2', 'n3'])
            for name, group in similarity:
                if group.shape[0] == 2:
                    mat[i][j] = mat[j][i] = 1
                    break
    return mat


cath_df = make_cath_df("cath-domain-list.txt")
lst = listCreation("PSSM.txt")
matTopology = neighbor_mat(cath_df, lst, homologous=False)
matHomologous = neighbor_mat(cath_df, lst, homologous=True)
graphTopology = csr_matrix(matTopology)
graphHomologous = csr_matrix(matHomologous)
topology_components, topologyLabels = connected_components(csgraph=graphTopology, directed=False, return_labels=True)
homologous_components, homologousLabels = connected_components(csgraph=graphHomologous, directed=False, return_labels=True)
print(topology_components)
print(topologyLabels)
print(homologous_components)
print(homologousLabels)
