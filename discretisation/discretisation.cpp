#include <iostream>
#include <cmath>
#include <cstring>
#include <sstream>

// CETTE CONSTANTE REPRESENTE LE NOMBRE D'INTERVALLES
#define SIZE 100

using namespace std;

/* PROTOTYPES DE FONCTIONS */
float u(float x);
float u_sec(float x);

double x(0), xmin(-5), xmax(5), ymin(-2), ymax(2), dx(0);

/* Résolution par la méthode de Gauss */
class Lsolver {
    public:
        // Constructeur et destructeur
        Lsolver();
        ~Lsolver();

        // affichage de matrice et vecteurs
        void displayMatrix();
        void displayVect();
        void displayResult();
        // résolution par la méthode de Gauss
        void gaussElim();           // triangularisation de la matrice A
        void solveTriang();         // résolution du système

        int getN(){ return n; }

    private:
        float **A;          // matrice du problème A*x = B
        float *b,           // second membre de l'équation
              *x;           // inconnu (résultats)
        float U0,           // valeur de départ
              Un;
        int n;
};

int main() {
    cout << "RESOLUTION D'UNE EQUATION DIFFERENTIELLE PARTIELLE" << endl;
    cout << endl;

/// DONNEES
    Lsolver solve;
    cout << "Pour n = " << solve.getN() << " | h = " << (float)5/solve.getN() << endl;
    cout << "Voici la matrice A : " << endl;
    //solve.displayMatrix();
    cout << "Voici le vecteur b : " << endl;
    //solve.displayVect();
/// CALCULS
    //cout << "Apres pivotage partiel ... " << endl,
    solve.gaussElim();
    //cout << "A : " << endl;
    //solve.displayMatrix();
    //cout << "b : " << endl;
    //solve.displayVect();
    solve.solveTriang();

/// RESULTATS
    cout << "Comparaison des resultats : " << endl;
    solve.displayResult();
    return 0;
}

/* Définition des fonctions */
Lsolver::Lsolver() {
    // conditions initiales
    U0 = 1;     Un = 11*exp(5);
    n = SIZE;
    // allocation dynamique de la matrice A et des vecteurs b et x
    A = new float*[n-1];
    for(int i(0); i<n-1; i++)
        A[i] = new float[n-1];

    b = new float[n-1];

    x = new float[n-1];

    // remplissage de la matrice A et de b
    // A matrice diagonale
    for(int i(0); i<n-1; i++) {
        for(int j(0); j<n-1; j++) {
            if(i==j)
                A[i][j] = -2;
            else if(i==j+1 || j==i+1)
                A[i][j] = 1;
            else
                A[i][j] = 0;
        }
    }

    // b=[u''1 -u0, u''2, ... u''(n-2), u''(n-1) -un]
    float x_max = 5,                    // bornes sup et min
          x_min = 0,
          h = (x_max + x_min)/n;        // pas de chaque intervalle
    for(int i(0); i<n-1; i++) {
        x_min += h;
        b[i] = h*h*u_sec(x_min);
        if(i==0) {
            b[i] -= U0;
        } else if(i==n-2) {
            b[i] -= Un;
        }
    }
}

Lsolver::~Lsolver() {
    delete[] A;
    delete[] b;
    delete[] x;
}

void Lsolver::displayMatrix() {
    for(int i(0); i<n-1; i++) {
        cout << "[ ";
        for(int j(0); j<n-1; j++) {
            cout << " " << A[i][j] << "  ";
        }
        cout << "]" << endl;
    }
    cout << endl;
}

void Lsolver::displayVect() {
    cout << "[";
    for(int i(0); i<n-1; i++) {
        cout << " " << b[i] << " ";
    }
    cout << "]" << endl;
    cout << endl;
}

void Lsolver::gaussElim() {
    float *temp = new float[n-1];

    // On recherche le max dans une colonne
    for(int k(0); k<n-1; k++) {
        float maxNbr = fabs(A[k][k]);   // pour la comparaison du max
        int ind = k;                    // pour stocker l'indice

        for(int r(k+1); r<n-1; r++) {
            if(maxNbr < fabs(A[r][k])) {    // s'il existe un float plus grand, on stocke son indice
                maxNbr = fabs(A[r][k]);
                ind = r;
            }
        }

        // echange de lignes
        temp = A[k];
        A[k] = A[ind];
        A[ind] = temp;

        float temp2 = 0;
        temp2 = b[k];
        b[k] = b[ind];
        b[ind] = temp2;

        // pivotage
        for(int i(k+1); i<n-1; i++) {
            for(int j(k+1); j<n-1; j++) {
                A[i][j] = A[i][j] - (A[i][k] / A[k][k]) * A[k][j];
            }
            b[i] = b[i] - (A[i][k] / A[k][k]) * b[k];
            A[i][k] = 0;
        }
    }
}

void Lsolver::solveTriang() {
    for(int i(n-1-1); 0<=i; i--) {
        float s = 0;

        // calcul de la somme a_ij * x_j
        for(int j(i+1); j<n-1; j++) {
            s += A[i][j] * x[j];
        }
        x[i] = (b[i] - s) / A[i][i];
    }
}

void Lsolver::displayResult() {
    float x_max = 5,
          x_min = 0,
          h = (x_max + x_min)/n;
    cout << "  x_i | discretisation | u(x)" << endl;

    for(int i(0); i<n-1; i++) {
        x_min += h;
        cout << x_min << "      " << x[i] << "      " << u(x_min) << " " << endl;;
    }
}

float u(float x) {
    return ((2*x + 1)*exp(x));
}

float u_sec(float x) {
    return ((2*x + 5)*exp(x));
}

