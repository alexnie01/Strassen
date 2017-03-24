# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 00:58:11 2017

@author: anie
"""

dimension = 10
import math
import time
import numpy as np
import random
import sys

"""
Generate Matrix

Was used for testing purposes originally, kept commented for reference
"""

def gen_mat(dim, filename = 'temp'):
     """ generate two random matrices for testing """
     mat_1 = np.random.randint(0,2,size=(dim,dim))
     mat_2 = np.random.randint(0,2,size=(dim,dim))
    
     # reshape and save to test text file
     writ_1 = mat_1.flatten()
     writ_2 = mat_2.flatten()
     save_file = open(filename+'.txt', 'w')
     np.savetxt(save_file, writ_1, fmt = '%i')
     np.savetxt(save_file, writ_2, fmt = '%i')
     save_file.close()
     
     return mat_1, mat_2
#
#"""
#Taking in inputs
#"""
#"""
#./strassen 0 dimension inputfile
#python strassen.py arg1 arg2 arg3
#"""
#cutoff = 31
#dimension = int(sys.argv[2])
#inputfile = sys.argv[3]
#
#f = open(inputfile)
#file_numbers = [int(line.strip()) for line in f.readlines()]
#
#matrix_1_list = file_numbers[:dimension**2]
#matrix_1 = np.array(matrix_1_list).reshape(dimension,dimension)
#matrix_2_list = file_numbers[dimension**2:]
#matrix_2 = np.array(matrix_2_list).reshape(dimension,dimension)
#
#
#
#def multiply_matrices_standard(matrix_1, matrix_2, dimension):
#    matrix = []
#    for i in range(0, dimension):
#        matrix_row = []
#        for k in range(0, dimension):
#            sum_j = 0
#            for j in range(0, dimension):
#                sum_j = sum_j + matrix_1[i][j]*matrix_2[j][k]
#            matrix_row.append(sum_j)
#        matrix.append(matrix_row)
#    return matrix
#            
#
#
## #### Straussen Helper Functions
#
#def matrix_addition(matrix_1, matrix_2, dimension):
#    matrix = [[matrix_1[i][j] + matrix_2[i][j] for j in range(dimension)] for i in range(dimension)]
#    return np.array(matrix).reshape(dimension, dimension)
#
#
#def matrix_subtraction(matrix_1, matrix_2, dimension):
#    matrix = [[matrix_1[i][j] - matrix_2[i][j] for j in range(dimension)] for i in range(dimension)]
#    return np.array(matrix).reshape(dimension, dimension)
#
#def pad_zeros(matrix, dimensions_to_pad):
#    """ pads to next power of two """
#    zeros = np.matrix([np.zeros(len(matrix)) for i in range(0, dimensions_to_pad)])
#    matrix = np.vstack((matrix, zeros))
#    zeros = np.matrix([np.zeros(len(matrix)) for i in range(0, dimensions_to_pad)]).T
#    matrix = np.hstack((matrix, zeros))
#    return np.array(matrix).reshape(len(matrix), len(matrix))
#
#
#def remove_zeros(multiplied_padded, dimensions_to_remove):
#    """ removes padding to original dimension
#        must pass in number of dimensions to remove """
#    new_dimension = len(multiplied_padded)-dimensions_to_remove
#    return np.array(multiplied_padded).reshape(len(multiplied_padded), len(multiplied_padded))[0:new_dimension, 0:new_dimension]
#
#
## #### Strassens
#
#def multiply_matrices_straussens(matrix_1, matrix_2, dimension):
#    # first check if small enough to switch to standard
#    if dimension <= cutoff:
#        return multiply_matrices_standard(matrix_1, matrix_2, dimension)
#    
#
#    # deal with it if it's not a power of two
#    is_power_two = True
#    dimension_difference = 0
#    closest_power_of_two = math.log(dimension,2)
#    if math.log(dimension, 2).is_integer() == False:
#        closest_power_of_two = 2**(int(math.log(dimension,2)) + 1)
#        dimension_difference = closest_power_of_two - dimension
#        # pad matrix with zeros as needed
#        matrix_1 = pad_zeros(matrix_1, dimension_difference)
#        matrix_2 = pad_zeros(matrix_2, dimension_difference)
#        dimension = closest_power_of_two
#        
#    new_dimension = dimension / 2
#    # matrix 1 blocks
#    a = matrix_1[0:new_dimension, 0:new_dimension]
#    b = matrix_1[0:new_dimension, new_dimension:dimension]
#    c = matrix_1[new_dimension:dimension, 0:new_dimension]
#    d = matrix_1[new_dimension:dimension, new_dimension:dimension]
#    # matrix 2 blocks
#    e = matrix_2[0:new_dimension, 0:new_dimension]
#    f = matrix_2[0:new_dimension, new_dimension:dimension]
#    g = matrix_2[new_dimension:dimension, 0:new_dimension]
#    h = matrix_2[new_dimension:dimension, new_dimension:dimension]
#    
#    # calculate intermediaries
#    P1 = multiply_matrices_straussens(a, matrix_subtraction(f, h, new_dimension), new_dimension)
#    P2 = multiply_matrices_straussens(matrix_addition(a,b, new_dimension), h, new_dimension)
#    P3 = multiply_matrices_straussens(matrix_addition(c, d, new_dimension), e, new_dimension)
#    P4 = multiply_matrices_straussens(d, matrix_subtraction(g, e, new_dimension), new_dimension)
#    P5 = multiply_matrices_straussens(matrix_addition(a, d, new_dimension), matrix_addition(e, h, new_dimension), 
#                                      new_dimension)
#    P6 = multiply_matrices_straussens(matrix_subtraction(b, d, new_dimension), matrix_addition(g, h, new_dimension), 
#                                      new_dimension)
#    P7 = multiply_matrices_straussens(matrix_subtraction(a, c, new_dimension), matrix_addition(e, f, new_dimension), 
#                                      new_dimension)
#    
#    # indices refer to blocks of product matrix
#    ae_bg = matrix_addition(matrix_subtraction(matrix_addition(P5, P4, new_dimension), P2, new_dimension), P6, new_dimension)
#    af_bh = matrix_addition(P1, P2, new_dimension)
#    ce_dg = matrix_addition(P3, P4, new_dimension)
#    cf_dh = matrix_subtraction(matrix_subtraction(
#            matrix_addition(P5, P1, new_dimension), P3, new_dimension), P7, new_dimension)
#    final_matrix = np.vstack((np.hstack((ae_bg, af_bh)), np.hstack((ce_dg, cf_dh))))
#    final_matrix = final_matrix.reshape(dimension, dimension)
#    
#    # check to see if we need to remove padding in product matrix
#    if dimension_difference != 0:
#        final_matrix = remove_zeros(final_matrix, dimension_difference)
#    return final_matrix
#
#
## # ### Testing
## matrix_1, matrix_2 = generate_matrix(dimension)
## print matrix_1
## print matrix_2
#
#
#start = time.clock()
#straussens = multiply_matrices_straussens(matrix_1, matrix_2, dimension)
#
## print diagonals
#for i in range(0,len(straussens)):
#    print straussens[i][i]
## print straussens
## ## Timings for Strassen's 
## end = time.clock()
## print end - start
#
#
## ## Standard Matrix for timings
## start = time.clock()
## standard = np.matrix(multiply_matrices_standard(matrix_1, matrix_2, len(matrix_2)))
## print standard
## end = time.clock()
## print end - start
#
#
## ## Original Testing
## matrix_np = np.matmul(matrix_1, matrix_2)
## assert standard.any() == straussens.any()
## assert np.matrix(matrix_addition(matrix_1, matrix_2, len(matrix_2))).all() == np.add(matrix_1, matrix_2).all()
## assert np.matrix(matrix_subtraction(matrix_1, matrix_2, len(matrix_2))).all() == np.subtract(matrix_1, matrix_2).all()
#