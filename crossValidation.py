

def listCreation(filename):
    complexes_list = []
    file1 = open(filename, 'r')
    line = file1.readline().split()
    while (len(line)!=0):
        if (len(line)==1):
            complexes_list.append(line[0][1:5] + line[0][-1])
        line = file1.readline().split()
    return complexes_list



a = listCreation('PSSM.txt')
print(a)


