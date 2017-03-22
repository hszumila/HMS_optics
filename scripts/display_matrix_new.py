# -*- coding: utf-8 -*-

from __future__ import division, print_function

from pprint import pprint

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np


recMatrixPath = '../build/hms_recon_coeff_z1cm_noxbeam_fit4.dat'

matrix = []
columns = [
    'Xp', 'Y', 'Yp', 'D',
    'e_xFp', 'e_xpFp', 'e_yFp', 'e_ypFp', 'e_xTar'
]


with open(recMatrixPath, 'r') as fi:
    # Skip header.
    while fi.next()[0] == '!':
        pass

    # Read the rest of the matrix.
    line = fi.next()
    while line[1:4] != '---':
        sLine = line.split()
        mLine = [float(val) for val in sLine[:4]]
        mLine.extend([int(char) for char in sLine[4]])
        matrix.append(mLine)

        line = fi.next()


matrixIndep = np.array([mLine for mLine in matrix if mLine[-1] == 0])
matrixDep = np.array([mLine for mLine in matrix if mLine[-1] != 0])


for name, matrix in zip(['indep', 'dep'], [matrixIndep, matrixDep]):
    with PdfPages('terms_new_{}.pdf'.format(name)) as pp:
        for iCol in xrange(4):
            fig, ax = plt.subplots()
            ax.plot(matrix[:, iCol], 'o')
            ax.grid(True)
            ax.set_title('{} terms'.format(columns[iCol]))
            ax.set_xlabel('$\#_\mathrm{term}$')
            ax.set_ylabel('${}$'.format(columns[iCol]))
            ax.set_xlim(-1, len(matrix))
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            fig.tight_layout()
            pp.savefig(fig)
            plt.close(fig)
