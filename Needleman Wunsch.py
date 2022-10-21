# Needleman Whunsh Global Alignment Algorithm implementation

#_________________________________________________________________
#Importing Numpy for better matrix multiplication
import numpy as np
#Importing Bio.Align.substitution_matrices to provide the BLOSUM62 matrix needed for the implementation
from Bio.Align import substitution_matrices


def MAX(set):
    return (max(set))

#Cleansing the matrix data and extracting the alphabet and the modified matrix
def matrix_extraction():
    #Temporary matrix to store the BLOSUM in
    temp_matrix = substitution_matrices.load("BLOSUM62")
    
    #Empty dictionary on which storing the alphabet
    alphabet = {}
    alph = temp_matrix.alphabet
    #Filling the dictionary
    for i in range(len(alph) - 1):
        alphabet[alph[i]] = i 

    #Extracting the matrix
    blosum62 = np.asarray(temp_matrix)
    
    #Deleting the last row and the last column
    blosum62 = np.delete(blosum62, -1, 0)
    blosum62 = np.delete(blosum62, -1, 1)
    
    #Getting the minimum
    blosum62_min = blosum62.min()

    #Subtracting (or adding the absolute value of ) the minimum from the matrix
    blosum62 = blosum62 - blosum62_min

    
    #Return statement
    return alphabet, blosum62

def matrix_filling(str1, str2, gap_penalty, extension_penalty, matrix, alphabet):
    mat = np.zeros((len(str1), len(str2)))

    for i in range(len(str1)):
        for j in range(len(str2)):

            #Particular case for first element
            if i == 0 and j == 0:
                mat[i][j] = matrix[alphabet.get(str1[i]), alphabet.get(str2[j])]

            #Particular case for first row
            elif i == 0 and j != 0:
                mat[i][j] = mat[i][j-1] - gap_penalty

            #Particular case for first column
            elif i != 0 and j == 0:
                mat[i][j] = mat[i-1][j] - gap_penalty
            else:
            #General case
                a = mat[i-1][j-1] + matrix[alphabet.get(str1[i]), alphabet.get(str2[j])]
                b = mat[i-1][j] - gap_penalty
                c = mat[i][j-1] - gap_penalty
            
                mat[i][j] = MAX([a, b, c])

    print(mat)
    input()
    return max(max(mat[:,len(str2) -1]), max(mat[len(str1) -1,:])), path_finder(mat, i = len(str1) -1, j = len(str2) -1)


#Traceback function
def path_finder(matrix, i, j, str=""):

    #Particular case for first element
    if i == 0 and j == 0:
        str = str + 'P'
        str = str[::-1]
        return str
    #Particular case for first row
    if not i > 0 and j > 0:
        str = str + 'R'
        return path_finder(matrix = matrix, i = i, j = j-1, str=str)
    #Particular case) for first column
    if i > 0 and not j > 0:
        str = str + 'D'
        return path_finder(matrix = matrix, i = i-1, j = j, str=str)

    #General Case
    max = np.max([matrix[i-1][j-1], matrix[i-1][j], matrix[i][j-1] ])
    if max == matrix[i-1, j-1]:
        str = str + 'P'
        return path_finder(matrix = matrix, i = i-1, j = j-1, str=str)
    elif max == matrix[i-1, j]:
        str = str + 'D'
        return path_finder(matrix = matrix, i = i-1, j = j, str=str)
    elif max == matrix[i, j-1]:
        str = str + 'R'
        return path_finder(matrix = matrix, i = i, j = j-1, str=str)
    

def pretty_print(str1, str2, alignement, score, gap_penalty, extension_penalty):
    new_alignment = ""
    for i in range(len(alignement)):
        if i < len(str1) and i < len(str2):
            if alignement[i] == 'P':
                if str1[i] == str2[i]:
                    new_alignment = new_alignment + "|"
                else:
                    new_alignment = new_alignment + " "
            elif alignement[i] == 'D':
                str2 = str2[:i] + "-" + str2[i:]
            elif alignement[i] == 'R':
                str1 = str1[:i] + "-" + str1[i:]
        elif i < len(str1) and i == len(str2):
            for j in range(len(str2)-len(str1)):
                str1 = str1 + "-"
        elif i < len(str2) and i == len(str1):
            for j in range(len(str1)-len(str2)):
                str2 = str2 + "-"


    print(str1)
    print(new_alignment)
    print(str2)

    print("The alignment score is:",  score ,".")
    if gap_penalty == extension_penalty:
        print("The extension penalty is equal to the gap penalty which is", gap_penalty , ".")
    else:
        print("The extension penalty is" , extension_penalty , ", while the gap penalty is", gap_penalty ,".")
    pass

def main():
    print("Input first string")
    str1 = "TFDERILGVQTYWAECLA"
    print("Input second string")
    str2 = "QTFWECIKGDNATY"

    print("Gap Penalty")
    gap_penalty = input()
    print("Expension Penalty")
    extension_penalty = input()

    alphabet, matrix = matrix_extraction()
    score , alignment = matrix_filling(str1, str2, float(gap_penalty), float(extension_penalty), matrix = matrix, alphabet = alphabet)
    pretty_print(str1, str2, alignment, score, gap_penalty, extension_penalty)

main()