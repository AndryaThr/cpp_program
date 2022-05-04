#include <iostream>

using namespace std;

// dimension
int dim = 20000;

/* PROTOTYPES */
float f1(float y, float z);
float f2(float x, float y);
float f3(float x, float z);
void eulerExplicite(float x0, float y0, float z0, float h, float *x, float *y, float *z);
void display(float *x, float *y, float *z);

/* MAIN */
int main()
{
    cout << "ATTRACTEUR DE ROSSLER\n" << endl;
    /// DONNEES
    float y0(0), z0(0), h(0.015);
    float x0(0);
    float x[dim], y[dim], z[dim];
    cout << "\nPar la methode d'Euler explicite..." << endl;

    /// CALCULS
    eulerExplicite(x0, y0, z0, h, x, y ,z);

    /// RESULTAT
    display(x, y, z);
    return 0;
}

/* FONCTIONS */
float f1(float y, float z)
{
    return -(y+z);
}

float f2(float x, float y)
{
    return (x + (0.1)*y);
}

float f3(float x, float z)
{
    return (0.1+z*(x-14));
}

// methode d'euler explicite
void eulerExplicite(float x0, float y0, float z0, float h, float *x, float *y, float *z)
{
    int i(0);
    x[0] = x0;
    y[0] = y0;
    z[0] = z0;

    for(i=0; i<dim; i++)
    {
        x[i+1] = x[i] +h*f1(y[i], z[i]);
        y[i+1] = y[i] +h*f2(x[i], y[i]);
        z[i+1] = z[i] +h*f3(x[i], z[i]);
    }
}

void display(float *x, float *y, float *z){

/// Get/Compute data

    //FILE *GP = popen("gnuplot -persist", "w");
    FILE *GP = popen("C:\\gnuplot\\bin\\gnuplot -persist", "w");
    //FILE *GP = popen("gnuplot -persist", "w");

    if(GP){
        fprintf(GP, "set term wxt  size 800, 600\n");
        fprintf(GP, "set term wxt  title 'Rossler'\n");
        fprintf(GP, "set xlabel \"x(t)\"\n");
        fprintf(GP, "set ylabel \"y(t)\"\n");
        fprintf(GP, "set zlabel \"z(t)\"\n");
        fprintf(GP, "set style data linespoint\n");
        fprintf(GP, "set size square\n");
        fprintf(GP, "set xrange [-20:25]\n");
        fprintf(GP, "set yrange [-25:20]\n");
        fprintf(GP, "set zrange [0:55]\n");

        // données
        fprintf(GP, "$Rossler << EOD\n");
        for(int i = 0; i < dim ;i++){
            fprintf(GP, "%f %f %f \n", x[i], y[i], z[i]);
            fprintf(GP, "\n");
        }

        // traçage
        fprintf(GP, "EOD\n");
        fprintf(GP, "splot $Rossler using 1:2:3  with lines\n");

        fflush(GP);
        pclose(GP);

    }
}

