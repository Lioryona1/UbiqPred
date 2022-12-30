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


def make_cath_df(filename, columns_number):
    """
    :param filename: cath-domain-list file
    :param columns_number: the number of columns to consider with the cath classification not include the cath domain name
    :return: dataframe of all the chains in the file and their cath classification divide to 4 different columns
    """
    df = pd.read_csv(filename, skiprows = 16, header=None, delimiter=r"\s+")
    df = df.iloc[:, 0:columns_number+1]
    cath_columns = ["n"+str(i) for i in range(1,columns_number+1)]
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


cath_df = make_cath_df("cath-domain-list.txt",3)
lst = listCreation("PSSM.txt")
matHomologous = neighbor_mat(cath_df, lst,3)
graphHomologous = csr_matrix(matHomologous)
homologous_components, homologousLabels = connected_components(csgraph=graphHomologous, directed=False, return_labels=True)
print(homologous_components)
print(homologousLabels)
