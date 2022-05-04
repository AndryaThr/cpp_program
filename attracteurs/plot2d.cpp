#include <iostream>
#include <cmath>
#include <cstring>
#include <sstream>

using namespace std;

/* PROTOTYPES DES FONCTIONS */
// fonctions
float f(float t, float y);
float df_dz(float t, float y);
float newton(float h, float x, float y);

// méthodes de résolution
void eulerExpl(float xmin, float xmax, float y0, float h, float *y);
void eulerImpl(float xmin, float xmax, float y0, float h, float *y);
void rungeKutta(float xmin, float xmax, float y0, float h, float *y);

// affichage des résultats
void plotSolution(float xmin, float xmax, float y0, float h, float *exp, float *imp, float *rk);
void gplot(float *xd, float *exp, float *emp, float *rk, float xmin, float xmax, float ymin, float ymax, int dim);

/* --------------------------- */

int main(int argc, char *argv[]){
	cout <<"Resolution d'une equation differentielle lineaire: y'(t)=-y(t)-3.exp(-t).sin(3t)" << endl;

/// Données
    float xmin(0),
          xmax(5),
          y0(1.0),
          h(0.1);
	float euler_exp[100] = {0},
          euler_imp[100] = {0},
          runge_kutta[100] = {0};

/// Calculs
    // par la méthode d'Euler explicite
	eulerExpl(xmin, xmax, y0, h, euler_exp);
    // par la méthode d'Euler implicite
	eulerImpl(xmin, xmax, y0, h, euler_imp);
    // par la méthode de Runge-Kutta
	rungeKutta(xmin, xmax, y0, h, runge_kutta);

/// Résultats
    plotSolution(xmin, xmax, y0, h, euler_exp, euler_imp, runge_kutta);

    return 0;
}

/* --------------------------- */
/* DEFINITIONS DES FONCTIONS */

float f(float t, float y) {
    float result = -y - (3*exp(-t)*sin(3*t));
    return result;
}

float df_dz(float t, float y) {
    return (-1);
}

float newton(float h, float x, float y) {
    float z_0(10.0),
          z_k(y);
    while(0.0001<fabs(z_0 - z_k)) { // pour epsilon = 0.0001
        z_0 = z_k;
        z_k = z_0 -((z_0 - y - h*f(x, z_0))/(1 - h*df_dz(x, z_0)));
    }

    return z_k;
}

void eulerExpl(float xmin, float xmax, float y0, float h, float *y) {
    int dim = (xmax - xmin) / h;
    float x(xmin);

    y[0] = y0;
    for(int i(1); i<=dim; i++) {
        x += h;
        y[i] = y[i-1] + h*f(x, y[i-1]);
    }
}

void eulerImpl(float xmin, float xmax, float y0, float h, float *y) {
    int dim = (xmax - xmin) / h;
    float x(xmin);

    y[0] = y0;
    for(int i(1); i<=dim; i++) {
        x += h;
        y[i] = newton(h, x, y[i-1]);
    }
}

void rungeKutta(float xmin, float xmax, float y0, float h, float *y) {
    int dim = (xmax - xmin) / h;
    float x(xmin),
          k1(0), k2(0), k3(0), k4(0);

    y[0] = y0;
    for(int i(1); i<=dim; i++) {
        x += h;

        k1 = h*f(x, y[i-1]);
        k2 = h*f((x + (h/2)), (y[i-1] + (k1/2)));
        k3 = h*f((x + (h/2)), (y[i-1] + (k2/2)));
        k4 = h*f((x + h), (y[i-1] + k3));

        y[i] = y[i-1] + (k1 + 2*k2 + 2*k3 + k4)/6;
    }
}

void plotSolution(float xmin, float xmax, float y0, float h, float *exp, float *imp, float *rk) {
    float x[100] = {0},
          ymax = 1,
          ymin = -1;
    int dim = (xmax - xmin)/h;

    // valeurs de l'axe des abcisses
    for(int i(1); i<=dim; i++)
        x[i] = x[i-1] + 0.1;

    // affichage de la courbe
    gplot(x, exp, imp, rk, xmin, xmax, ymin, ymax, dim);
}

void gplot(float *x, float *exp, float *imp, float *rk, float xmin, float xmax, float ymin, float ymax, int dim){
    // FILE *GP = popen("C:\\gnuplot\\bin\\gnuplot -persist","w");
    FILE *GP = popen("gnuplot -persist","w");
    ostringstream buff;

    if (GP) {  		// gnuplot found
		cout << "Type q to close gnuplot" << endl;
/// Plot setup
        fprintf(GP, "set title 'Fonctions'\n");
        fprintf(GP, "set xlabel 'x'\n");
        fprintf(GP, "set ylabel 'y'\n");
        fprintf(GP, "set xzeroaxis\n");
        fprintf(GP, "set style data lines\n"); // if not smoothed
        buff << "set xrange [" << xmin << ":" << xmax << "]\n";
        fprintf(GP,"%s\n", buff.str().c_str());
		buff << "set yrange [" << ymin << ":" << ymax << "]\n";
        fprintf(GP,"%s\n", buff.str().c_str());

/// Reusable internal data block
        fprintf(GP, "$Fonctions << EOF\n");
        for(int i(0); i<=dim; i++)
            fprintf(GP, "%f %f %f %f\n", x[i], exp[i], imp[i], rk[i]);

        fprintf(GP, "EOF\n");
        fprintf(GP, "plot $Fonctions using 1:2 w lines t 'Euler explicite', $Fonctions using 1:3 w lines t 'Euler implicite', $Fonctions using 1:4 w lines t 'Runge-Kutta'\n");

/// Send commands to gnuplot and close pipe
        fflush(GP); // cmd vers gnuplot
        pclose(GP);
    }
    else { cout << "gnuplot not found..." << endl; }
}


