package calculnumerictema3;

import static java.lang.Math.*;

/**
 *
 * @author Gabriel Budau
 */
public class CalculNumericTema3 {
    /*
     Sa se rezolve aceeasi formula A*x=b
     Matricea A: a(i,i) = 2; a(i+1,i) = 1; a(i,i+1) = 1; restul elementelor sunt 0
     Vectorul b: b(i) = i;
     Sa se rezolve cu metodele: Jacobi relaxata, Gauss-Seidel relaxata si Metoda gradientului conjugat
     Pentru primele 2 metode sa gasim parametrul de relaxare optim cum am facut la curs. Se va afisa numarul de iteratii pentru a aproxima solutia cu eroarea epsilon.
     Se va rula pentru matricile de marimea n=20, 50, 100, iar epsilon = 10^(-10)
     */

    final String DOUBLE_DISPLAY_FORMAT = "%+g";
    final String RESULT_DISPLAY_FORMAT = "%+15g";
    final Integer SYST_SIZE = 4;
    final double epsilon = 10e-10;

    double A[][];
    double b[];
    double x[];
    double x_JR[];//Jacobi relaxata
    double x_GS[];//Gauss-Seidel relaxata
    double x_GC[];//Metoda gradientului conjugat
    int n;

    protected double getAElement(int i, int j) {
        if (i == j) {
            return 2;
        }
        if ((i - 1) == j) {
            return 1;
        }
        if (i == (j - 1)) {
            return 1;
        } else {
            return 0;
        }
    }

    protected double getBelement(int i) {
        return (double) i;
    }

    protected void initArrays() {
        n = SYST_SIZE + 1;

        A = new double[n][n];
        b = new double[n];
        x = new double[n];
        x_GC = new double[n];
        x_GS = new double[n];
        x_JR = new double[n];
        for (int i = 1; i < n; i++) {
            for (int j = 1; j < n; j++) {
                A[i][j] = getAElement(i, j);
            }
            b[i] = getBelement(i);
        }
    }

    protected void printSystem(double X[]) {
        for (int i = 1; i < n; i++) {
            System.out.print("( ");
            for (int j = 1; j < n; j++) {
                System.out.printf(" " + DOUBLE_DISPLAY_FORMAT + " ", A[i][j]);
            }
            System.out.printf(") * ( " + DOUBLE_DISPLAY_FORMAT + " ) = ( " + DOUBLE_DISPLAY_FORMAT + " ) \n", X[i], b[i]);
        }
    }

    protected void JR(int iteratii) {
        int q = 0;
        int p = iteratii;
        double sigma;
        //Norma Matricii ||A|| = t;
        double t = 0.0;
        for (int i = 1; i < n; i++) {
            double suma = 0.0;
            for (int j = 1; j < n; j++) {
                suma += abs(A[i][j]);
            }
            if (suma > t) {
                t = suma;
            }
        }

        for (int k = 1; k < p - 1; k++) {
            sigma = (2 * k) / (p * t);
            double B[][] = new double[n][n];
            double temp_b[] = new double[n];

            for (int i = 1; i < n; i++) {
                for (int j = 1; j < n; j++) {
                    if (i == j) {
                        B[i][j] = 1 - sigma * A[i][j];
                    } else {
                        B[i][j] = -(sigma * A[i][j]);
                    }
                }
                temp_b[i] = sigma * b[i];
            }
            double er = 0.0;

            do {
                double y[] = new double[n];

                for (int i = 1; i < n; i++) {
                    double suma = 0.0;
                    for (int j = 1; j < n; j++) {
                        suma += B[i][j] * x_JR[j];
                    }
                    y[i] = suma + temp_b[i];
                }

                double norma_suma = 0.0;
                for (int i = 1; i < n; i++) {
                    double suma = 0.0;
                    for (int j = 1; j < n; j++) {
                        suma += A[i][j] * (y[i] - x_JR[i]) * (y[j] - x_JR[j]);
                    }
                    norma_suma += suma;

                }
                er = norma_suma;
                for (int i = 1; i < n; i++) {
                    x_JR[i] = y[i];
                }

                q++;

            } while (er >= epsilon);
            //System.out.println("N: " + n);

        }
        System.out.printf("\n%30s (", "METODA JACOBI RELAXATA");
        for (int j = 1; j < n; j++) {
            System.out.printf(" " + RESULT_DISPLAY_FORMAT + " ", x_JR[j]);
        }
        System.out.printf(")");
    }

    protected void GS(int iteratii) {
        double sigma = 0.0, sigmao = 0.0;
        double xo[] = new double[n];
        double y[] = new double[n];
        double xi[][] = new double[n][n];
        int p = iteratii;
        int m, mo = 0;
        double suma = 0.0, suma2 = 0.0, temp = 0.0, er = 0.0;
        for (int i = 1; i < n; i++) {
            x[i] = 0;
            xo[i] = 0;
            y[i] = 0;
            xi[0][i] = 0;
        }
        for (int k = 1; k <= p - 1; k++) {
            sigma = 2 * k;
            sigma = sigma / p;
            m = 0;
            do {
                y[1] = (1 - sigma) * x[1];
                suma = 0;
                for (int j = 2; j < n; j++) {
                    suma = suma + A[1][j] * x[j];
                }
                temp = sigma / A[1][1];
                y[1] = y[1] + temp * (b[1] - suma);
                for (int i = 2; i < n; i++) {
                    y[i] = (1 - sigma) * x[i];
                    suma = 0;
                    for (int j = 1; j <= i - 1; j++) {
                        suma = suma + A[i][j] * y[j];
                    }
                    suma2 = 0;
                    for (int j = i + 1; j < n; j++) {
                        suma2 = suma2 + A[i][j] * x[j];
                    }
                    temp = (double) (sigma / A[i][i]);
                    y[i] = y[i] + temp * (b[i] - suma - suma2);
                }
                suma = 0;
                for (int i = 1; i < n; i++) {
                    for (int j = 1; j < n; j++) {
                        suma = suma + A[i][j] * (y[i] - x[i]) * (y[j] - x[j]);
                    }
                }
                er = (sqrt(suma));
                for (int i = 1; i < n; i++) {
                    x[i] = y[i];
                }
                m++;
            } while (er >= epsilon);
            if (k == 1) {
                mo = m;
                sigmao = sigma;
                for (int i = 1; i < n; i++) {
                    xo[i] = x[i];
                }
            } else if (m < mo) {
                mo = m;
                sigmao = sigma;
                for (int i = 1; i < n; i++) {
                    xo[i] = x[i];
                }
            }
        }
        System.out.printf("\n%30s (", "METODA GAUSS-SEIDEL RELAXATA");
        for (int j = 1; j < n; j++) {
            System.out.printf(" " + RESULT_DISPLAY_FORMAT + " ", xo[j]);
        }
        System.out.printf(") n0 = " + mo + " | sigma = " + sigmao);

    }

    protected void GC() {
        double alfa[] = new double[n];
        double ci[] = new double[n];
        double xi[][] = new double[n][n];
        double ri[][] = new double[n][n];
        double vi[][] = new double[n][n];
        double suma = 0.0, suma2 = 0.0, suma3 = 0.0, er = 0.0;

        for (int i = 1; i < n; i++) {
            suma = 0;
            for (int j = 1; j < n; j++) {
                suma = suma + A[i][j] * xi[0][i];
            }
            ri[0][i] = b[i] - suma;
            vi[0][i] = ri[0][i];
            xi[0][i] = 0;
        }

        int k = 0;
        do {
            suma = 0;
            for (int j = 1; j < n; j++) {
                suma += ri[k][j] * ri[k][j];
            }
            suma2 = 0;
            for (int i = 1; i < n; i++) {
                suma3 = 0;
                for (int j = 1; j < n; j++) {
                    suma3 = suma3 + A[i][j] * vi[k][j];
                }
                suma2 += suma3 * vi[k][i];
            }
            alfa[k] = (suma / suma2);
            for (int i = 1; i < n; i++) {
                xi[k + 1][i] = xi[k][i] + alfa[k] * vi[k][i];
            }
            suma = 0;
            for (int j = 1; j < n; j++) {
                suma += (xi[k + 1][j] - xi[k][j]) * (xi[k + 1][j] - xi[k][j]);
            }
            er = (sqrt(suma));
            for (int i = 1; i < n; i++) {
                suma = 0;
                for (int j = 1; j < n; j++) {
                    suma += A[i][j] * xi[k + 1][j];
                }
                ri[k + 1][i] = b[i] - suma;
            }
            suma = 0;
            suma2 = 0;
            for (int j = 1; j < n; j++) {
                suma += ri[k + 1][j] * ri[k + 1][j];
                suma2 += ri[k][j] * ri[k][j];
            }
            ci[k] = suma / suma2;
            for (int i = 1; i < n; i++) {
                vi[k + 1][i] = ri[k + 1][i] + ci[k] * vi[k][i];
            }

        } while (er >= epsilon && k++ < n-2);

        System.out.printf("\n%30s (", "METODA GRADIENTULUI GONJUGAT");
        for (int j = 1; j < n; j++) {
            System.out.printf(" " + RESULT_DISPLAY_FORMAT + " ", xi[k - 1][j]);
        }
        System.out.printf(") ");
    }

    public static void main(String[] args) {
        CalculNumericTema3 cnt3 = new CalculNumericTema3();
        cnt3.initArrays();
        cnt3.printSystem(cnt3.x);
        cnt3.JR(cnt3.SYST_SIZE * 1000);
        cnt3.GS(cnt3.SYST_SIZE * 1000);
        cnt3.GC();
        System.out.print("\n\n");
    }
}
