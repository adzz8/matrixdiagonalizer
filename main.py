import numpy as np
import sympy as sp

def pretty_print_mat(mat):
    return sp.pretty(mat)

def diagonalize_matrix(A):
    n = A.shape[0]
    A = sp.Matrix(A)
    lam = sp.symbols('λ')
    char_poly = (A - lam*sp.eye(n)).det()
    eigenvalues = sp.solve(char_poly, lam)
    print(f"The characteristic polynomial is: {char_poly}")
    print(f"The eigenvalues are: {eigenvalues}")
    C = []
    eigenvalues_in_C = []
    for val in eigenvalues:
        eigenvectors = (A - val*sp.eye(n)).nullspace()
        print(f"\nFor λ = {val}, we get:")
        for vec in eigenvectors:
            print(f"RREF leads to:\n{pretty_print_mat(vec)}")
            C.append(vec)
            eigenvalues_in_C.append(val)
    if len(C) < n:
        print("\nThe matrix is not diagonalizable because it does not have n linearly independent eigenvectors.")
        return
    C = sp.Matrix.hstack(*C)
    D = sp.diag(*eigenvalues_in_C)
    C_inv = C.inv()
    print("\nThen we have that A is diagonalizable with A = CDC⁻¹, with:")
    print(f"C:\n{pretty_print_mat(C)}")
    print(f"D:\n{pretty_print_mat(D)}")
    result = A == C*D*C_inv
    if result:
        print("\nSo, it holds that A = CDC⁻¹.")
    else:
        print("\nHowever, A = CDC⁻¹ does not hold.")


n = int(input("Enter the size of your matrix: "))
print(f"Please enter the elements of the {n}x{n} matrix row by row, separating numbers with a space.")
A = []
for i in range(n):
    row = list(map(float, input().split()))
    A.append(row)
A = np.array(A)
diagonalize_matrix(A)
