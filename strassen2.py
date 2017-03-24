#!/usr/bin/python
"""
Created on Fri Mar 24 00:58:11 2017

@author: anie
"""

import math
import time
import numpy as np
import random
import sys
import pdb

""" 
Parameters and inputs
"""
"""
Taking in inputs

./strassen 0 dimension inputfile
python strassen.py arg1 arg2 arg3
"""
# switchover dimension for Strassen to classical matrix multiplication
cutoff = 4

"""
Generate Matrix
"""

def gen_mat(dim, filename = 'temp'):
     """ generate two random matrices for testing """
     mat_1 = np.random.randint(0,2,size=(dim,dim))
     mat_2 = np.random.randint(0,2,size=(dim,dim))
#     sys.stdout.write(str(mat_1))
#     sys.stdout.write(str(mat_2))
    
     # reshape and save to test text file in flat format desired by assignment
     writ_1 = mat_1.flatten()
     writ_2 = mat_2.flatten()
     # open in write mode: clear previous contents, but append as needed
     save_file = open(filename+'.txt', 'w')
     np.savetxt(save_file, writ_1, fmt = '%i')
     np.savetxt(save_file, writ_2, fmt = '%i')
     save_file.close()
     
     return mat_1, mat_2

def test():
    dim = 10
    mat_1, mat_2 = gen_mat(10)
    return (dim,'temp.txt')
def main():
    # dimension of matrices to multiply
    dim = int(sys.argv[2])
    # read in matrices single digit per line
    inputfile = sys.argv[3]
    return dim,inputfile
    
def std_matmult(mat_1, mat_2, dim):
    """ standard matrix multiplication """
    # container for product matrix
    prod_mat = np.zeros((dim,dim), dtype = int)
    for i in np.arange(dim):
        for k in np.arange(dim):
            # running sum of components for (i,k) index in product matrix
            sum_j = 0
            for j in np.arange(dim):
                sum_j += mat_1[i,j] * mat_2[j,k]
            prod_mat[i,k] = sum_j 
    return prod_mat
    
"""
Strassen's algorithm and helper functions
"""

def add_mat(mat_1, mat_2, dim):
    """ perform matrix addition element-wise """
    sum_mat = np.zeros((dim,dim), dtype=int)
    for i in np.arange(dim):
        for j in np.arange(dim):
            sum_mat[i,j] = mat_1[i,j]+mat_2[i,j]
    return sum_mat

def sub_mat(mat_1, mat_2, dim):
    """ perform matrix subtraction element-wise, mat_1-mat_2 """
    diff_mat = np.zeros((dim,dim), dtype=int)
    for i in np.arange(dim):
        for j in np.arange(dim):
            diff_mat[i,j] = mat_1[i,j]-mat_2[i,j]
    return diff_mat

    
def str_matmult(mat_1, mat_2, dim, cutoff):
    """ modified Strassen's algorithm. if input dim > switchover cutoff, use
    standard multiplication. if odd, pad one line of zeros on bottom and right,
    then recurse on dim+1. if even, use strassen's """
    if dim <= cutoff:
        return std_matmult(mat_1, mat_2, dim)
    else:
        if dim%2==1:
            # pad bottom and right with one line of zeros
            mat_10 = np.hstack((np.vstack((mat_1, np.zeros((1, dim), dtype=int))), 
                              np.zeros((dim+1,1), dtype=int)))
            mat_20 = np.hstack((np.vstack((mat_2, np.zeros((1, dim), dtype=int))), 
                              np.zeros((dim+1,1), dtype=int)))
            return str_matmult(mat_10, mat_20, dim+1, cutoff)[:-1, :-1]
        else:
            half = dim/2
            # blocks of first matrix
            a11 = mat_1[:half, :half]
            a12 = mat_1[:half, half:]
            a21 = mat_1[half:, :half]
            a22 = mat_1[half:, half:]
            
            # blocks of second matrix
            b11 = mat_2[:half, :half]
            b12 = mat_2[:half, half:]
            b21 = mat_2[half:, :half]
            b22 = mat_2[half:, half:]
            
            # strassen intermediates
            p1 = str_matmult(add_mat(a11,a22, half), add_mat(b11,b22, half), 
                             half, cutoff)
            p2 = str_matmult(add_mat(a21, a22, half), b11, half, cutoff)
            p3 = str_matmult(a11, sub_mat(b12,b22, half), half, cutoff)
            p4 = str_matmult(a22, sub_mat(b21, b11, half), half, cutoff)
            p5 = str_matmult(add_mat(a11, a12, half), b22, half, cutoff)
            p6 = str_matmult(sub_mat(a21, a11, half), add_mat(b11, b12, half),
                             half, cutoff)
            p7 = str_matmult(sub_mat(a12, a22, half), add_mat(b21, b22, half),
                             half, cutoff)
            
            # product subblocks
            c11 = add_mat(sub_mat(add_mat(p1, p4, half), p5, half), p7, half)
            c12 = add_mat(p3, p5, half)
            c21 = add_mat(p2, p4, half)
            c22 = add_mat(add_mat(sub_mat(p1, p2, half), p3, half), p6, half)
            
            # combine into product matrix
            return np.vstack((np.hstack((c11,c12)), np.hstack((c21,c22))))
            
if __name__ == "__main__":
    dim, inputfile = test()
    f = open(inputfile)
    data = np.loadtxt(inputfile, dtype = int)
    f.close()
    mat_1 = data[:dim**2].reshape(dim,dim)
    mat_2 = data[dim**2:].reshape(dim,dim)
#    sys.stdout.write(str(mat_1))
#    sys.stdout.write(str(mat_2))

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
