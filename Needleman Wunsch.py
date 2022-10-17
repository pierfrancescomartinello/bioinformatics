# Needleman Whunsh Global Alignment Algorithm implementation

#_________________________________________________________________
#Importing Numpy for better matrix multiplication
import numpy as np
#Importing Bio.Align.substitution_matrices to provide the BLOSUM62 matrix needed for the implementation
from Bio.Align import substitution_matrices


#Cleansing the matrix data and extracting the alphabet and the modified matrix
def matrix_extraction():
    #Temporary matrix to store the BLOSUM in
    temp_matrix = substitution_matrices.load("BLOSUM62")
    
    #Empty dictionary on which storing the alphabet
    alphabet = {}
    alph = temp_matrix.alphabet
    #Deleting the character used for gap (we don't use it)
    alph = alph[:-1]
    #Filling the dictionary
    for i in range(len(alph)):
        alphabet[alph[i] ] = i 

    #Extracting the matrix
    blosum62 = np.asmatrix(temp_matrix)
    
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
    mat = np.zeros([len(str1) + 1, len(str2) + 1])

    for i in range(len(str1)):
        for j in range(len(str2)):

            #Particular case for first element
            if i == 0 and j == 0:
                mat[i][j] = matrix[alphabet[str1[i]]][alphabet[str2[j]]]

            #Particular case for first row
            if i == 0 and j != 0:
                mat[i][j] = mat[i][j-1]

            #Particular case for first column
            if i != 0 and j == 0:
                mat[i][j] = mat[i-1][j]

            #General case
            mat[i][j] = max(mat[i-1][j-1] + matrix[alphabet[str1[i]]][alphabet[str2[j]]], mat[i-1][j] - gap_penalty, mat[i][j-1] - gap_penalty)
    
    return mat[len(str1)][len(str2)] , path_finder(matrix, i = len(str1), j = len(str2))


#Traceback function
def path_finder(matrix, i, j, str=''):

    #Particular case for first element

    if i == 0 and j == 0:
        str = str + 'P'
        return str
    #Particular case for first row
    if i == 0 and j != 0:
        str = str + 'R'
        path_finder(matrix = matrix, i = i-1, j = j-1, str=str)
    #Particular case for first column
    if i != 0 and j == 0:
        str = str + 'D'
        path_finder(matrix = matrix, i = i-1, j = j-1, str=str)

    #General Case
    if max(matrix[i-1][j-1], matrix[i-1][j], matrix[i][j-1]) == matrix[i-1][j-1]:
        str = str + 'P'
        path_finder(matrix = matrix, i = i-1, j = j-1, str=str)
    elif max(matrix[i-1][j-1], matrix[i-1][j], matrix[i][j-1]) == matrix[i-1][j]:
        str = str + 'D'
        path_finder(matrix = matrix, i = i-1, j = j-1, str=str)
    elif max(matrix[i-1][j-1], matrix[i-1][j], matrix[i][j-1]) == matrix[i][j-1]:
        str = str + 'R'
        path_finder(matrix = matrix, i = i-1, j = j-1, str=str)
    

def pretty_print(str1, str2, alignement, score, gap_penalty, extension_penalty):
    new_alignment = ""
    for i in range(len(alignement)):
        if alignement[i] == 'P':
            if str1[i] == str2[i]:
                new_alignment = new_alignment + "-"
            else:
                new_alignment = new_alignment + " "
        elif alignement[i] == 'D':
            str2 = str2[:i] + "-" + str2[i:]
        elif alignement[i] == 'R':
            str1 = str1[:i] + "-" + str1[i:]


    print(str1)
    print(new_alignment)
    print(str2)

    print("The alignment score is: " + score +".")
    if gap_penalty == extension_penalty:
        print("The extension penalty is equal to the gap penalty which is " + gap_penalty + ".")
    else:
        print("The extension penalty is " + extension_penalty + ", while the gap penalty is " + gap_penalty +".")
    pass

def main():
    print("Input first string")
    str1 = input()
    print("Input second string")
    str2 = input()

    print("Gap Penalty")
    gap_penalty = input()
    print("Expension Penalty")
    extension_penalty = input()

    alphabet, matrix = matrix_extraction()

    score , alignment = matrix_filling(str1, str2, gap_penalty, extension_penalty, matrix, alphabet)

    pretty_print(str1, str2, alignement, score, gap_penalty, extension_penalty)


main()