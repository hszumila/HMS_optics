# -*- coding: utf-8 -*-

from __future__ import division, print_function


columnNames = [
    'Xp', 'Y', 'Yp', 'D',
    'e_xFp', 'e_xpFp', 'e_yFp', 'e_ypFp', 'e_xTar'
]


def readOld(path):
    header = ''
    matrix = []

    with open(path, 'r') as fi:
        # Read header.
        line = fi.next()
        while line[0] == '!':
            header += line

            line = fi.next()

        separator = line

        # Read the rest of the matrix.
        line = fi.next()
        while line != separator:
            mLine = [float(line[1+i*16:1+(i+1)*16]) for i in xrange(4)]
            mLine.extend([int(char) for char in line[66:71]])
            matrix.append(mLine)

            line = fi.next()

    return matrix, header, separator


def read(path):
    header = ''
    matrix = []

    with open(path, 'r') as fi:
        # Read header.
        line = fi.next()
        while line[0] == '!':
            header += line

            line = fi.next()

        separator = line

        # Read the rest of the matrix.
        line = fi.next()
        while line[1:4] != '---':
            sLine = line.split()
            mLine = [float(val) for val in sLine[:4]]
            mLine.extend([int(char) for char in sLine[4]])
            matrix.append(mLine)

            line = fi.next()

    return matrix, header, separator


def write(path, matrix, header, separator):
    with open(path, 'w') as fo:
        fo.write(header)
        fo.write(separator)

        for mLine in matrix:
            fo.write(formatLine(mLine) + '\n')

        fo.write(separator)


def formatLine(mLine):
    return '{:17.9e}{:17.9e}{:17.9e}{:17.9e}  {:.0f}{:.0f}{:.0f}{:.0f}{:.0f}'.format(
        *mLine
    )


def compareExponents(rowA, rowB):
    orderA = sum(rowA[4:])
    orderB = sum(rowB[4:])

    if orderA < orderB:
        return -1
    elif orderA > orderB:
        return 1
    else:
        numA = rowA[4]*10000 + rowA[5]*1000 + rowA[6]*100 + rowA[7]*10 + rowA[8]
        numB = rowB[4]*10000 + rowB[5]*1000 + rowB[6]*100 + rowB[7]*10 + rowB[8]

        if numA < numB:
            return 1
        elif numA > numB:
            return -1
        else:
            return 0
