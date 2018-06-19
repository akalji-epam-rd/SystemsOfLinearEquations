package com.akalji.SystemsOfLinearEquations;

import com.akalji.matrix.*;
import com.akalji.matrix.exceptions.*;
import com.akalji.SystemsOfLinearEquations.exceptions.*;


public class PLUDecomposition {
    public static Solution PLUC(Matrix M) {
        int swaps = 0;
        int n = M.getVsize();
        int m = M.getHsize();
        int rank = n;
        Matrix A = new Matrix(M);
        Matrix P = new Matrix(M.getVsize(), M.getHsize());
        Matrix Q = new Matrix(M.getVsize(), M.getHsize());
        Matrix L = new Matrix(M.getVsize(), M.getHsize());
        Matrix U = new Matrix(M.getVsize(), M.getHsize());

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                if (i == j) {
                    P.set(i, j, 1);
                    Q.set(i, j, 1);
                } else {
                    P.set(i, j, 0);
                    Q.set(i, j, 0);
                }

        for (int i = 0; i < n; ++i) {
            double max = Math.abs(A.getElement(i, i));
            int max_row = i, max_col = i;
            for (int j = i; j < n; ++j)
                for (int k = i; k < n; ++k)
                    if (Math.abs(A.getElement(j, k)) > max) {
                        max = Math.abs(A.getElement(j, k));
                        max_row = j;
                        max_col = k;
                    }
            if (Math.abs(max) < Math.pow(10, -15)) {
                --rank;
                continue;
            }

            if (i != max_row) {
                A.swapRows(i, max_row);
                P.swapRows(i, max_row);
                swaps++;
            }

            if (i != max_col) {
                A.swapColumns(i, max_col);
                Q.swapColumns(i, max_col);
                swaps++;
            }

            for (int j = i + 1; j < n; ++j) {
                A.set(j, i, A.getElement(j, i) / A.getElement(i, i));
            }

            for (int j = i + 1; j < n; ++j) {
                for (int k = i + 1; k < n; ++k) {
                    A.set(j, k, A.getElement(j, k) - A.getElement(j, i) * A.getElement(i, k));
                }
            }
        }

        for (int i = 0; i < n; ++i) {
            for (int j = i; j < n; ++j) {
                if (i == j) {
                    L.set(j, i, 1);
                } else {
                    L.set(j, i, A.getElement(j, i));
                }
                U.set(i, j, A.getElement(i, j));
            }
        }
        return new Solution(P, Q, L, U, rank, swaps);
    }

    public static double det(Solution solution) {
        double det = 1;
        Matrix U = solution.getU();
        if (U.getVsize() != U.getHsize()) {
            throw new NotSquareMatrix();
        }
        int N = U.getVsize();
        for (int i = 0; i < N; i++) {
            det *= U.getElement(i, i);
        }

        if (solution.getSwaps() % 2 == 0)
            return det;
        else
            return -det;
    }

    public static Matrix SOLE(Solution PLUCSolution, Matrix b) {
        Matrix U = PLUCSolution.getU();
        Matrix L = PLUCSolution.getL();
        int rank = PLUCSolution.getRank();
        Matrix x = new Matrix(b.getVsize(), 1);

        int n = L.getVsize();
        Matrix y = new Matrix(n, 1);
        for (int i = 0; i < n; ++i) {
            double tmp = b.getElement(i, 0);
            for (int j = 0; j < i; ++j)
                tmp = tmp - L.getElement(i, j) * y.getElement(j, 0);
            y.set(i, 0, tmp / L.getElement(i, i));
        }
        if (rank == n) {
            for (int i = n - 1; i >= 0; --i) {
                double tmp = y.getElement(i, 0);
                for (int j = n - 1; j > i; --j)
                    tmp = tmp - U.getElement(i, j) * x.getElement(j, 0);
                x.set(i, 0, tmp / U.getElement(i, i));
            }
        } else {
            boolean res = false;
            for (int i = n - 1; i >= rank; --i)
                if (Math.abs(y.getElement(i, 0)) >= Math.pow(10, -15)) res = true;
            if (res) {
                throw new LinearEquationIncompatibleException();
            } else {
                for (int i = n - 1; i >= rank; --i) x.set(i, 0, 1);
                for (int i = rank - 1; i >= 0; --i) {
                    double tmp = y.getElement(i, 0);
                    for (int j = n - 1; j > i; --j)
                        tmp = tmp - U.getElement(i, j) * x.getElement(j, 0);
                    x.set(i, 0, tmp / U.getElement(i, i));
                }
            }
        }

        return x;
    }
}

