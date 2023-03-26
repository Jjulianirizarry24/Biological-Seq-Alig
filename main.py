import csv
import sys


def read_csvFile():
    file_name = sys.argv[1]
    with open(file_name, 'r') as file:
        reader = csv.reader(file)
        next(reader)
        seqs = []
        
        for row in reader:
            seqs.append(row)
            
        return seqs

    # Biological Sequence Alignment

def BioSeq(seq1,seq2):
    # Creating the num matrix

    row = len(seq1) # = 7
    tuple = len(seq2) # = 9

    scoreMatrix = [[0 for i in range(tuple + 1)] for j in range(row + 1)]
    movMatrix = [[0 for i in range(tuple)] for j in range(row)]

    # creating Score matrix
    for i in range(tuple + 1):
        scoreMatrix[0][i] = i * -2

    for i in range(row + 1):
        scoreMatrix[i][0] = i * -2



    for i in range(1,row+1): #8 belonging to seq2
        for j in range(1,tuple+1): # 10 belonging to seq1

            left = scoreMatrix[i][j-1]
            dia = scoreMatrix[i-1][j-1]
            top = scoreMatrix[i-1][j]

            # Checks if the current character matches

            if seq1[i-1] == seq2[j-1]:
                left = left - 2
                dia = dia + 1
                top = top - 2

                maxVal = max(left, dia, top)

                # Sequence of conditions that check which to add into both matrices
                if maxVal == left:
                    scoreMatrix[i][j] = left
                    movMatrix[i-1][j-1] = 'l'

                elif maxVal == top:
                    scoreMatrix[i][j] = top
                    movMatrix[i-1][j-1] = 't'

                else:
                    scoreMatrix[i][j] = dia
                    movMatrix[i-1][j-1] = 'd'
            # In case of mismatch
            else:
                left = left - 2
                dia = dia - 1
                top = top - 2

                maxVal = max(left, dia, top)

                if maxVal == left:
                    scoreMatrix[i][j] = left
                    movMatrix[i-1][j-1] = 'l'

                elif maxVal == top:
                    scoreMatrix[i][j] = top
                    movMatrix[i-1][j-1] = 't'

                else:
                    scoreMatrix[i][j] = dia
                    movMatrix[i-1][j-1] = 'd'


    #Now for backtracking...
    bsa1 = ""
    bsa2 = ""
    i = row-1 
    j = tuple-1

    while i >= 0 and j >= 0:
        
        if movMatrix[i][j] == "d":
            bsa1 += seq1[i]
            bsa2 += seq2[j]
            i = i-1
            j = j-1

        elif movMatrix[i][j] == "l":
            bsa1 += "-"
            bsa2 += seq2[j]
            j = j-1

        else:
            bsa2 += "-"
            bsa1 += seq1[i]
            i = i-1
       
    
    if i >= 0: 
        while i >= 0:
            bsa1 +=seq1[i]
            bsa2 += '-'
            i -= 1
    
    if j >= 0:    
        while j >= 0:
            bsa2 += seq2[j]
            bsa1 += '-'
            j -= 1

    return bsa1[::-1],bsa2[::-1],scoreMatrix[row][tuple]

if len(sys.argv) > 1:
    seqs = read_csvFile()
    
    for i in seqs:   
        ag1,ag2,mat = BioSeq(i[0],i[1])
        print(ag1,ag2,mat)