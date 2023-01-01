# Program that solves a system of linear equations using the matrix method with 
# Jacobi's method or Gauss-Seidel's method

import numpy as np
import argparse

# Function that solves a system of linear equations using the matrix method with Jacobi's method
def jacobi(A, b, tol, maxiter):
    k = 0
    D = np.diag(np.diag(A))
    L = np.tril(A) - D
    U = np.triu(A) - D
    x = np.zeros(len(b)).astype(np.float64)
    error = np.linalg.norm(np.dot(A, x) - b)
    success = False

    while error > tol and k < maxiter:
        # D^-1 * (b - (L + U) * x)
        x = np.dot(np.linalg.inv(D).astype(np.float64), b - np.dot(L + U, x).astype(np.float64)).astype(np.float64)
        # ||(A * x - b)|| / ||(b)|| - where ||.|| is the Euclidean norm 
        error = np.linalg.norm(np.dot(A, x) - b) / np.linalg.norm(b)
        k += 1

    if error <= tol:
        success = True

    return x, error, success

# Function that solves a system of linear equations using the matrix method with Gauss-Seidel's method
def gauss_seidel(A, b, tol, maxiter):
    k = 0
    D = np.diag(np.diag(A))
    L = np.tril(A) - D
    U = np.triu(A) - D
    x = np.zeros(len(b)).astype(np.float64)
    error = np.linalg.norm(np.dot(A, x) - b)
    success = 0

    while error > tol and k < maxiter:
        # (L + D)^-1 * (b - U * x)
        x = np.dot(np.linalg.inv(L + D).astype(np.float64), b - np.dot(U, x).astype(np.float64)).astype(np.float64)
        # ||(A * x - b)|| / ||(b)|| - where ||.|| is the Euclidean norm 
        error = np.linalg.norm(np.dot(A, x) - b) / np.linalg.norm(b)
        k += 1

    if error <= tol:
        success = 1

    return x, error, success   

# Main program
if __name__ == '__main__':
    arg = argparse.ArgumentParser()
    arg.add_argument('-g', '--gamma', type = float, default=10, 
                    help = 'Parameter which is used in construction of matrix A and vector b')
    arg.add_argument('-m', '--method', type = str, default='jacobi', choices=['jacobi', 'gauss_seidel'], 
                    help = 'Method to solve the system of linear equations')
    args = vars(arg.parse_args())

    gamma = args['gamma']

    A = np.zeros((20, 20)).astype(np.float64)

    # Fill diagonal
    for i in range(20):
        A[i, i] = gamma

    # Fill upper diagonal
    for i in range(19):
        A[i, i+1] = -1

    # Fill lower diagonal
    for i in range(1, 20):
        A[i, i-1] = -1

    # Fill vector b with x
    b = np.full(20, gamma - 2).astype(np.float64)
    b[0] = gamma - 1
    b[19] = gamma - 1
    
    tol = 1e-5
    maxiter = 300

    if args['method'] == 'jacobi':
        x, error, success = jacobi(A, b, tol, maxiter)
    elif args['method'] == 'gauss_seidel':
        x, error, success = gauss_seidel(A, b, tol, maxiter)
    else:
        print('The method is not valid')

    if success:
        print('The solution is:')
        print(x)
        print('The error is:')
        print(error)
    else:
        print('The solution was not found')