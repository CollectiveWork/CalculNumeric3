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

    final String DOUBLE_DISPLAY_FORMAT = "%+10g";
    final String RESULT_DISPLAY_FORMAT = "%+10g";;
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
                        suma += A[i][j]*(y[i] - x_JR[i])*(y[j] - x_JR[j]);
                    }
                    norma_suma += suma;

                }
                er = norma_suma;
                for(int i = 1; i < n ; i++)
                    x_JR[i] = y[i];
                
                q++;
                
            } while (er >= epsilon);
            //System.out.println("N: " + n);
            
        }
        System.out.printf("\n%50s\n\n(", "METODA JACOBI RELAXATA");
            for (int j = 1; j < n; j++) {
                System.out.printf(" " + RESULT_DISPLAY_FORMAT + " ", x_JR[j]);
            }
            System.out.printf(")");
    }

    protected void GS(int iteratii){
        int p = iteratii;
        double sigma = 0.0;
        int nn;
        for (int k = 1; k < p-1; k ++) {
            sigma = (2*k)/p;
            nn = 0; //x = 0
            double er = 0.0;
            
            do{
                for(int i = 1; i < n; i ++){
                    if(i == 1){
                        double suma = 0.0;
                        for(int j = 1; j < nn+1; j++)
                            suma += A[i][j]*x_GS[j];
                        x_GS[i] = (1-sigma)*x_GS[i] - ((sigma/A[i][i])*(b[i] - suma));
                    }else{
                        double suma1 = 0.0;
                        double suma2 = 0.0;
                        for (int j = 1; j <= i-1; j++) {
                            suma1 += A[i][j] * x_GS[j];
                        }
                        for (int j = i + 1; j <n; j++) {
                            suma1 += A[i][j] * x_GS[j];
                        }
                         x_GS[i] = (1 - sigma) * x_GS[i] - ((sigma/A[i][i])*(b[i] - suma1 - suma2));
                    }
                        
                }
                
                double y[] = new double[n];
                for(int i = 1; i < n ; i++)
                    y[i] = x_GS[i];
                
                double norma_suma = 0.0;
                for (int i = 1; i < n; i++) {
                    double suma = 0.0;
                    for (int j = 1; j < n; j++) {
                        suma += A[i][j]*(y[i] - x_GS[i])*(y[j] - x_GS[j]);
                    }
                    norma_suma += suma;

                }
                er = norma_suma;
                for(int i = 1; i < n ; i++)
                    x_GS[i] = y[i];
                nn++;
               
                
            }while (er >= epsilon);
            
        }
         System.out.printf("\n%50s\n\n(", "METODA GAUSS-SEIDEL RELAXATA");
            for (int j = 1; j < n; j++) {
                System.out.printf(" " + RESULT_DISPLAY_FORMAT + " ", x_GS[j]);
            }
            System.out.printf(")");
    }
    
    protected void GC(int iteratii){
        
    }
    
    
    public static void main(String[] args) {
        CalculNumericTema3 cnt3 = new CalculNumericTema3();
        cnt3.initArrays();
        cnt3.printSystem(cnt3.x);
        cnt3.JR(cnt3.SYST_SIZE * 1000);
        cnt3.GS(cnt3.SYST_SIZE * 10);
    }
}
