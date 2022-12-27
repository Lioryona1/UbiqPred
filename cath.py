from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import pandas as pd
import numpy as np

def listCreation(filename):
    """
    :param filename: PSSM file
    :return: list of all the complexes in the file
    """
    complexes_list = []
    file1 = open(filename, 'r')
    line = file1.readline().split()
    while (len(line)!=0):
        if (len(line)==1):
            complexes_list.append(line[0][1:5] + line[0][-1])
        line = file1.readline().split()
    file1.close()
    return complexes_list


def make_cath_df(filename):
    """
    param filename: cath-domain-list file
    return: dataframe of all the complex in the file and their cath classification divide to 4 different columns
    """
    df = pd.read_fwf(filename,delimiter ='',skiprows=17, header = None)
    df = df.iloc[:, 0:5]
    df.columns = ['complex', 'n1', 'n2', 'n3', 'n4']
    df['complex'] = df['complex'].apply(lambda x: x[0:5])
    return df


def neighbor_mat(df, lst, homologous=True):
    """
    :param df: cath data frame as it return from the func make_cath_df
    :param lst: list of complexes
    :param homologous: if true, make the matrix base on the homologous match, else, base on topology
    :return: matrix. mat[i][j] == 1 if there is homologous connection between complex i and complex j
    """
    n = len(lst)
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1,n):
            if homologous:
                similarity = df[df['complex'].isin([lst[i],lst[j]])].groupby(by = ['n1','n2','n3','n4'])
            else: #topology
                similarity = df[df['complex'].isin([lst[i], lst[j]])].groupby(by=['n1', 'n2', 'n3'])
            for name, group in similarity:
                if group.shape[0] == 2:
                    mat[i][j] = mat[j][i] = 1
                    break
    return mat

graph = [
[0, 1, 1, 0, 0],
[0, 0, 1, 0, 0],
[0, 0, 0, 0, 0],
[0, 0, 0, 0, 1],
[0, 0, 0, 0, 0]
]
graph = csr_matrix(graph)

print(graph)
n_components, labels = connected_components(csgraph=graph, directed=False, return_labels=True)
print(n_components)
print(labels)

