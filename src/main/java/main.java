/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

import static java.lang.Math.pow;

/**
 * @author иисусе
 */
public class main {


    double e = pow(10, -10);
    double K = 0.01;
    int RADIUS = 2;
    int T0 = 0, T1 = 1;
    double alpha = 0.01;
    int c = 2;
    int sumT = T0 + T1;
    double SQRT_K = pow(K, 0.5);
    double K1_5 = pow(K, 1.5);

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
    }

    public double T(double t, double eta) {
        return Math.exp(-t * pow(eta, 2) / this.c);
    }

    public double X(double r, double eta) {
        return Math.sin(r * eta / this.SQRT_K);
    }

    public double R(double r, double eta) {
        return X(r, eta);
    }

    public double chi(double eta) {
        return 0.5 * (RADIUS + (SQRT_K / (2 * eta)) * Math.sin(2 * RADIUS * eta / SQRT_K));
    }

    public double C(double eta) {
        double s = sumT * K / (pow(eta, 2) * chi(eta));
        double s1 = Math.sin(RADIUS * eta / SQRT_K) - RADIUS * eta * Math.cos(RADIUS * eta / SQRT_K) / SQRT_K;
        return s * s1;
    }

    public double v(double r, double t, double[] eigens) {
        double s = 0;
        for (int i = 0; i < eigens.length; i++)
            s += (C(eigens[i]) * R(r, eigens[i]) * T(t, eigens[i]));
        return s;
    }

    public double u(double r, double t, double[] eigens) {
        double s = 0;
        for (int i = 0; i < eigens.length; i++)
            s += (C(eigens[i]) * X(r, eigens[i]) * T(t, eigens[i]));
        return s / r - T1;

    }

}
