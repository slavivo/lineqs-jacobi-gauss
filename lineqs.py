# Program that solves a system of linear equations using the matrix method with 
# Jacobi's method or Gauss-Seidel's method

import numpy as np
import argparse

# Function that solves a system of linear equations using the matrix method with Jacobi's method
def jacobi(A, b, tol, maxiter, ref):
    k = 0
    D = np.diag(np.diag(A))
    L = np.tril(A) - D
    U = np.triu(A) - D
    x = np.zeros(len(b)).astype(np.float64)
    error = np.linalg.norm(np.dot(A, x) - b).astype(np.float64)
    flag = 0

    # Check if spectral radius of E - D^-1 * A < 1 to guarantee convergence
    if np.max(np.abs(np.linalg.eigvals(np.eye(20) - np.dot(np.linalg.inv(D), A)))) < 1:
        flag = 2

    # Check if A is sharply diagonally dominant to guarantee convergence
    if np.all(np.abs(np.diag(A)) > np.sum(np.abs(A), axis=1) - np.abs(np.diag(A))):
        flag = 2

    while error > tol and k < maxiter:
        # Iterative refinement
        if ref:
            r = np.dot(A, x) - b
            y = np.dot(np.linalg.inv(A), r)
            x = x - y

        # D^-1 * (b - (L + U) * x)
        x = np.dot(np.linalg.inv(D), b - np.dot(L + U, x))
        # ||(A * x - b)|| / ||(b)|| - where ||.|| is the Euclidean norm 
        error = (np.linalg.norm(np.dot(A, x) - b) / np.linalg.norm(b))
        k += 1

    if error <= tol:
        flag = 1

    return x, k, error, flag

# Function that solves a system of linear equations using the matrix method with Gauss-Seidel's method
def gauss_seidel(A, b, tol, maxiter, ref):
    k = 0
    D = np.diag(np.diag(A))
    L = np.tril(A) - D
    U = np.triu(A) - D
    x = np.zeros(len(b)).astype(np.float64)
    error = np.linalg.norm(np.dot(A, x) - b)
    flag = 0

    # Check if spectral radius of E - (D + L)^-1 * A < 1 to guarantee convergence
    if np.max(np.abs(np.linalg.eigvals(np.eye(20) - np.dot(np.linalg.inv(D + L), A)))) < 1:
        flag = 2

    # Check if A is symmetric and positive definite and has positive diagonal to guarantee convergence
    if np.allclose(A, A.T) and np.all(np.linalg.eigvals(A) > 0) and np.all(np.diag(A) > 0):
        flag = 2

    while error > tol and k < maxiter:
        # Iterative refinement
        if ref:
            r = np.dot(A, x) - b
            y = np.dot(np.linalg.inv(A), r)
            x = x - y

        # (L + D)^-1 * (b - U * x)
        x = np.dot(np.linalg.inv(L + D), b - np.dot(U, x))
        # ||(A * x - b)|| / ||(b)|| - where ||.|| is the Euclidean norm 
        error = np.linalg.norm(np.dot(A, x) - b) / np.linalg.norm(b)
        k += 1

    if error <= tol:
        flag = 1

    return x, k, error, flag   

# Main program
if __name__ == '__main__':
    arg = argparse.ArgumentParser()
    arg.add_argument('-g', '--gamma', type = float, default=10, 
                    help = 'Parameter which is used in construction of matrix A and vector b')
    arg.add_argument('-m', '--method', type = str, default='jacobi', choices=['jacobi', 'gauss_seidel'], 
                    help = 'Method to solve the system of linear equations')
    arg.add_argument('-r', '--refinement', type = bool, default=False, help='Use iterative refinement')
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
    maxiter = 1000
    ref = args['refinement']

    if args['method'] == 'jacobi':
        x, k, error, flag = jacobi(A, b, tol, maxiter, ref)
    elif args['method'] == 'gauss_seidel':
        x, k, error, flag = gauss_seidel(A, b, tol, maxiter, ref)
    else:
        print('The method is not valid')

    if flag == -1:
        print('The method did not converge due to non-convergence')
    elif flag == 1:
        print(f'The solution is: {x}\nThe error is: {error}\nThe number of iterations is: {k}')
    elif flag == 2:
        print('The method did not converge due to precision issues, but it is theoretically convergent')
    else:
        print('The method did not converge due to unknown reasons')