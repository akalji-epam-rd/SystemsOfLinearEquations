package com.akalji.SystemsOfLinearEquations;

import com.akalji.SystemsOfLinearEquations.exceptions.LinearEquationIncompatibleException;
import com.akalji.matrix.Matrix;

public class SystemOfLinearEquation {
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
