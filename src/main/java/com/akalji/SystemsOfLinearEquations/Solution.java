package com.akalji.SystemsOfLinearEquations;

import com.akalji.matrix.Matrix;

public class Solution {
    private Matrix P;
    private Matrix Q;
    private Matrix L;
    private Matrix U;
    private double det;
    private int rank;
    private int swaps;

    public Solution(Matrix p, Matrix q, Matrix l, Matrix u, int rank, int swaps) {
        P = p;
        Q = q;
        L = l;
        U = u;
        this.rank = rank;
        this.swaps = swaps;
    }

    public int getSwaps() {
        return swaps;
    }

    public Matrix getP() {
        return P;
    }

    public Matrix getQ() {
        return Q;
    }

    public Matrix getL() {
        return L;
    }

    public Matrix getU() {
        return U;
    }

    public int getRank() {
        return rank;
    }
}
