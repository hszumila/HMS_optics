# -*- coding: utf-8 -*-

from __future__ import division, print_function

import numpy as np

import matrix_utils as mutil


recMatrixPath = 'hms_recon_coeff_z1cm_noxbeam_fit4.dat'


matrix, header, separator = mutil.read(recMatrixPath)


matrixIndep = [mLine for mLine in matrix if mLine[-1] == 0]
matrixDep = [mLine for mLine in matrix if mLine[-1] != 0]
matrixIndep.sort(cmp=mutil.compareExponents)
matrixDep.sort(cmp=mutil.compareExponents)


for name, matrix in zip(['indep', 'dep'], [matrixIndep, matrixDep]):
    fNameOut = '{}__{}.dat'.format(recMatrixPath[:-4], name)
    print(fNameOut)
    mutil.write(fNameOut, matrix, header, separator)
