import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.solvers.BisectionSolver;
import org.apache.commons.math3.analysis.solvers.NewtonRaphsonSolver;
import org.apache.commons.math3.exception.DimensionMismatchException;

import java.sql.SQLClientInfoException;
import java.util.ArrayList;

import static java.lang.Math.*;

/**
 * Created by dodler on 22/12/15.
 @author diman
 */
public class Solver {

    private static final double H = pow(10, -5);
    double C, K, R, T0, T1, alpha;
    double sumT = T1 - T0;
    double SQRT_K;
    public final UnivariateFunction foo = new UnivariateFunction() {
        public double value(double eta) {
            double var1 = SQRT_K / pow(R, 2) - SQRT_K * alpha / R,
                    var2 = R / SQRT_K;

            return var1 / eta - 1 / tan(var2 * eta);
        }
    };
    double K1_5 = pow(K, 1.5);
    double[] eigens;
    private double eigMul = SQRT_K * (1 / pow(R, 2) - alpha / R);

    public Solver(double C, double K, double R, double T0, double T1, double alpha) {
        this.C = C;
        this.K = K;
        this.R = R;
        this.T0 = T0;
        this.T1 = T1;
        this.alpha = alpha;

        SQRT_K = pow(K, 0.5);
        K1_5 = pow(K, 1.5);
        sumT = T1 - T0;
    }

    public double[] getEigens(int i) {
        System.out.println("getEigens start");
//        UnivariateDifferentiableFunction foo = new UnivariateDifferentiableFunction()
        UnivariateFunction foo = new UnivariateFunction() {
            public double value(double eta) {
                double var1 = SQRT_K *(1/R - alpha),
                        var2 = R / SQRT_K;

                return var1 / eta - 1/ tan(eta*var2);
            }
        };

        int cnt = 0;

        BisectionSolver solver = new BisectionSolver(pow(10, -8));
        double[] result = new double[i];
//        double step = pow(10, -5), a = (Math.PI / 2) / (R / SQRT_K) - step, b = a + step, period = Math.PI / (R/ SQRT_K);
//        while (cnt < i) {
//            if (foo.value(a) * foo.value(b) < 0) {
////                System.out.print("solving a=" + a + ",b=" + b);
//                result[cnt++] = solver.solve(10000, foo, a, b);
////                if (abs(foo.value(result[cnt - 1])) > 1)
////                {
////                    System.out.print(" failed ");
////                } else
////                {
////                    System.out.print(" solved ");
////                }
////                System.out.println(" v=" + result[cnt - 1] + ", f(v)=" + foo.value(result[cnt - 1]));
//                b += (period + step);
//                a += period;
//            }
//            b += step;
//        }

        double a = 0, period = Math.PI * SQRT_K/R, b = a + period;
        while (cnt < i) {
                result[cnt++] = solver.solve(1000000, foo, a+0.01, b- 0.01);
                if (cnt < 10){
//                    System.out.println("f(a) = " + foo.value(a - 1/(cnt*0.5)));
//                    System.out.println("f(b) = " + foo.value(b));
                    System.out.println(result[cnt-1] + ",f=" + foo.value(result[cnt-1]));
                }
                b += period;
                a += period;
        }
        return result;
    }

    protected double T(double t, double eta) {
        return exp(-t * pow(eta, 2) / C);
    }

    protected double X(double r, double eta) {
        return (sin(r * eta / SQRT_K));
    }

    protected double R(double r, double eta) {
        if (r<pow(10,-8)){
            return eta/SQRT_K;
        }
        return sin(r * eta / SQRT_K) / r;
    }

    protected double chi(double eta) {
        return R / 2 - (SQRT_K / (4 * eta)) * sin(2 * R * eta / SQRT_K);
    }

    protected double C(final double eta) {
        return alpha*R*sin(R*eta/SQRT_K)*sumT*K/(pow(eta,2)*chi(eta));
//        double arg = (R * eta / SQRT_K), n1 = sumT * K / pow(eta, 2), n2 = sin(arg) - arg * cos(arg);
//        double value = n1 * n2 / chi(eta);
//        return n1 * n2 / chi(eta);
//        System.out.println("vaue={" + value + "}");
//        return value;
    }

    public double u(double r, double t, double[] eigens) {
        double s = 0;
        for (int i = 0; i < eigens.length; i++) {
//            double value = (C(eigens[i]) * R(r, eigens[i]) * T(t, eigens[i]));
//            s += value;
//            System.out.println(s);
            s += (C(eigens[i]) * R(r, eigens[i]) * T(t, eigens[i]));
        }

        return T1 - s;
    }

    /**
     * numeric experiment
     * @param eps accuracy to achieve
     * @param r value of radius - should be less then 10^-8
     * @param t value of time - can be any, better 0
     */
    public void testEps(double eps, double r, double t) {
        double eigens[] = getEigens(100000);
        System.out.println("test with params:"+ eps + ":r=" + r +":t="  + t +" started");
        double theorSum = 0;
        for (int i = 0; i < eigens.length; i++) {
            theorSum += (C(eigens[i]) * R(r, eigens[i]) * T(t, eigens[i]));
        }
        System.out.println("theorSum=" + (T1-theorSum));

        double realSum = 0;
        int i = 0;
//        while(Math.abs(realSum - theorSum) > eps){
//            realSum += (C(eigens[i]) * R(r, eigens[i]) * T(t, eigens[i]));
//            i++;
//        }

        double value = C(eigens[i]) * R(r, eigens[i]) * T(t, eigens[i]);
        while(Math.abs(value) > eps/2.0){
            value = C(eigens[i]) * R(r, eigens[i]) * T(t, eigens[i]);
            System.out.println("eigen estimation:" + Math.PI*(i+1)*SQRT_K/R + ",eigens[i]=" + eigens[i] + ",C=" + C(eigens[i]) + ",R=" + R(r,eigens[i]) + ", T=" + T(t,eigens[i]));
            i++;
        }

        System.out.println("test finished. counter=" + i);
    }
}
