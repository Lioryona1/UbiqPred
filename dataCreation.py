import subprocess

import numpy as np
from Bio.PDB import PDBList
from Bio.PDB.MMCIFParser import MMCIFParser


# def get_alphafold_download_link(uniprot_id):
#     link_pattern = 'https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v2.pdb'
#     return link_pattern.format(uniprot_id)
#
#
# def download_alphafold_prediction(uniprot_id):
#     url = get_alphafold_download_link(uniprot_id)
#
#     result = subprocess.run(['wget', url, '-O', 'result.pdb'])
#     return result  # Result will be 0 if operation was successful.


# note : 2D3G and 3A9K in the table twice with different ligands

PDB_names_list = ['1NBF', '1P3Q', '1S1Q', '1UZX', '1WR6', '1WRD', '1XD3', '1YD8', '2AYO', '2C7M', '2D3G', '2DX5',
                  '2FIF', '2G45', '2GMI', '2HD5', '2HTH', '2IBI', '2J7Q', '2OOB', '2QHO', '2WDT', '2WWZ', '2XBB',
                  '3A33', '3A9K', '3BY4', '3C0R', '3CMM', '3I3T','3IFW', '3IHP', '3JSV', '3JVZ', '3K9P', '3KVF',
                  '3KW5', '3LDZ', '3MHS', '3MTN', '3NHE', '3O65',
                  '3OFI', '3OJ3', '3OLM', '3PHW', '3PRM', '3PT2', '3PTF', '3TBL', '3TMP', '3VHT']

pdb1 = PDBList()
#pdb1.download_pdb_files(pdb_codes=PDB_names_list, overwrite=True,
#                                 pdir='C:/Users/liory/YearC/workshop_proteins/UbiqPred/pdbs')


# fileNames = [
#     pdb1.retrieve_pdb_file(PDB_names_list[i], pdir='C:/Users/liory/YearC/workshop_proteins/UbiqPred/pdbs') for
#     i in range(len(PDB_names_list))]
#
# fileNames = ['C:/Users/liory/YearC/workshop_proteins/UbiqPred/pdbs/{}.cif'.format(PDB_names_list[i]) for i in
#              range(len(PDB_names_list))]
parser = MMCIFParser()  # create parser object

# print("file name is", fileNames[0])
# structures = [parser.get_structure(PDB_names_list[i], fileNames[i]) for i in range(len(PDB_names_list))]
#
# chains = structures[0].get_chains()
# for chain in chains:
#     print(chain.get_id(), chain.get_full_id())

threeLetterToSingelDict = {'GLY': 'G', 'ALA': 'A', 'VAL': 'V', 'LEU': 'L', 'ILE': 'I', 'THR': 'T', 'SER': 'S',
                           'MET': 'M', 'CYS': 'C', 'PRO': 'P', 'PHE': 'F', 'TYR': 'Y', 'TRP': 'W', 'HIS': 'H',
                           'LYS': 'K', 'ARG': 'R', 'ASP': 'D', 'GLU': 'E',
                           'ASN': 'N', 'GLN': 'Q'}

ubiq_names = ['UBIQ_', 'RS27A_MOUSE', 'UBC_HUMAN', 'RS27A_HUMAN', 'UBB_HUMAN', 'Q5U5U6_HUMAN', 'UBI4P_YEAST',
              'UBC_HUMAN', 'P62988', 'UBB_BOVIN', 'Q24K23_BOVIN']

spacial_cif = ['C:/Users/liory/YearC/workshop_proteins/UbiqPred/pdbs/3IFW.cif', "C:/Users/liory/YearC/workshop_proteins/UbiqPred/pdbs/3KVF.cif","C:/Users/liory/YearC/workshop_proteins/UbiqPred/pdbs/3KW5.cif"]
def findUbiqChains(filename):
    """
    param filename: pdb file (or other format of structure)
    return: chain ID's of the UBIQUITIN proteins
    """
    file1 = open(filename, 'r')
    chains = []
    found = False
    while True:
        bool = True  # for the case in 3ifw that there are ubiq_ before we expect
        line = file1.readline().split()
        for word in line:
            for ubiq in ubiq_names:
                if ubiq in word:
                    if filename in spacial_cif:
                        if bool:
                            bool = False
                            continue
                    ubiqLine = line
                    found = True
                    break
            if found:
                break
        if found:
            break
    while True:
        line = file1.readline().split()
        if "loop_" == line[0]:  # start loop_ phase after UNP phase
            break
    line = file1.readline().split()
    while "_struct" in line[0]:
        line = file1.readline().split()
    while line[0].isnumeric():
        id = line[8]
        if id in ubiqLine:
            chains.append(line[3])
        line = file1.readline().split()
    return chains


# def writeHeaderForChain(file1,s)

def atomDist(atom1, atom2):
    """
    param atom1: atom object
    param atom2: atom object
    return: the euclidian distance between the atoms
    """

    vector1 = atom1.get_vector()
    vector2 = atom2.get_vector()
    temp = vector1 - vector2  # subtracting vector
    sum_sq = np.dot(temp, temp)  # sum of the squares
    return np.sqrt(sum_sq)


def aaOutOfChain(chain):
    """
    :param chain: chain object
    :return: list of aa (not HOH molecule)
    """
    my_list = []
    amino_acids = chain.get_residues()
    for aa in amino_acids:
        name = str(aa.get_resname())
        if name in threeLetterToSingelDict.keys():  #amino acid and not other molecule
            my_list.append(aa)
    return my_list


def getLabelForAA(aa, ubiq_atoms, threshold):
    """
    param aa: amino acid object
    :param ubiq_atoms: the ubiquitin atoms
    :return: 1 if there exists an atom that is within 4 Angstrom to a ubiquitin atom else 0
    """
    for atom in aa.get_atoms():
        for ubiq_atom in ubiq_atoms:
            dist = atomDist(atom, ubiq_atom)
            if dist < threshold:
                return 1
    return 0


def structurePPBSFormat(file1, structure, structure_filename):
    """
    param filename: file to write to
    :param structure: pdb structure
    :param structure_filename:  the pdb filename of the structure
    The function write the structure into filename in PPBS format
    """
    # file1 = open(filename, 'w')
    ubiq_chains_id = findUbiqChains(structure_filename)
    print("found chains\n", ubiq_chains_id)
    # print(ubiq_chains_id)
    chains = structure.get_chains()
    ubiq_chains = [structure[0][id] for id in ubiq_chains_id]
    ubiq_amino_acids = []
    for ubiq_chain in ubiq_chains:
        ubiq_amino_acids += aaOutOfChain(ubiq_chain)
    ubiq_atoms = []
    for ubiq_aa in ubiq_amino_acids:
        ubiq_atoms += ubiq_aa.get_atoms()
    for chain in chains:
        if str(chain.get_id()) not in ubiq_chains_id:  # not a ubiquitin chain
            file1.write(">" + str(structure.get_id()).lower() + "_0-" + str(chain.get_id()) + "\n")
            chainPPBSFomrat(file1, chain, ubiq_atoms)


def chainPPBSFomrat(file1, chain, ubiq_atoms):
    """
    :param file1: file to write to
    :param chain: chain structure
    :param ubiq_atoms: list contains the ubiquitin atoms
    The function write the chain into filename in PPBS format (Without header)
    """
    amino_acids = aaOutOfChain(chain)
    for aa in amino_acids:
        name = str(aa.get_resname())
        chain_id = str(chain.get_id())
        aa_id = str(aa.get_id()[1])
        aa_type = threeLetterToSingelDict[name]
        line = [chain_id, aa_id, aa_type]
        threshold = 4
        label = str(getLabelForAA(aa, ubiq_atoms, threshold))
        line.append(label)
        file1.write(" ".join(line)+"\n")



PBBS_file = open('PSSM.txt', 'a')
structures = [parser.get_structure(PDB_names_list[i],
                                   r'C:/Users/liory/YearC/workshop_proteins/UbiqPred/pdbs/{}.cif'.format(
                                       PDB_names_list[i])) for i in range(len(PDB_names_list))]
# structure1 = parser.get_structure('1NBF', r'C:\Users\omriy\WorkshopProteins\final_project\UBIPred\UBDs\1nbf.cif')
for i in range(len(structures)):

    print(structures[i])
    structurePPBSFormat(PBBS_file, structures[i],
                            r'C:/Users/liory/YearC/workshop_proteins/UbiqPred/pdbs/{}.cif'.format(PDB_names_list[i]))
PBBS_file.close()