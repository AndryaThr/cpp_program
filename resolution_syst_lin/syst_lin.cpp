#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

/* CLASSE Solve */
class Solve {
    public:
    // CONSTRUCTEUR & DESTRUCTEUR
        Solve(string fileName);     // constructeur
        ~Solve();                   // destructeur
    // AFFICHAGE
        void displayMat();          // matrice dim x dim
        void displayVect();         // matrice dim x 1

    // CALCUL
        void gaussElim();           // triangularisation de la matrice A
        void solveTriang();         // résolution du système

    private:
        int dim;
        float **A;
        float *b, *x;
};



/* MAIN */
int main()
{
    /// donnée
    cout << "Résolution d'un systeme linéaire A.x = b : \n" << endl;
    Solve prod1("matrice.txt");

    cout << "Voici la matrice A : " << endl;
    prod1.displayMat();
    cout << endl;

    cout << "Voici le vecteur b : " << endl;
    prod1.displayVect();
    cout << endl;

    /// calculs et résultat
    cout << "Après pivotage partiel : " << endl;
    prod1.gaussElim();
    prod1.displayMat();
    cout << endl;
    prod1.displayVect();
    cout << endl;
    cout << "La solution est : " << endl;
    prod1.solveTriang();
    cout << endl;

    return 0;
}

/* FONCTIONS */

// Constructeur
Solve::Solve(string fileName) {
    ifstream mfile(fileName);           // ouverture du fichier

    mfile >> dim;                       // récupere la dimension de la matrice

    // Allocation dynamique tableau 2 dimensions
    A = new float*;
    for(int i(0); i<dim; i++)
        A[i] = new float[dim];

    for(int i(0); i<dim; i++) {
        for(int j(0); j<dim; j++) {
            mfile >> A[i][j];
        }
    }

    x = new float[dim];         // allocation dynamique du vecteur résultat
    for(int i(0); i<dim; i++)
        x[i] = 0;               // initialisation à 0

    // allocation dynamique vecteur
    b = new float[dim];
    for(int i(0); i<dim; i++) {
        mfile >> b[i];
    }

    mfile.close();                      // fermeture du fichier
}

// Affichage de la matrice A
void Solve::displayMat() {
    for(int i(0); i<dim; i++) {
        cout << "[ ";
        for(int j(0); j<dim; j++) {
            cout << " " << A[i][j] << "  ";
        }
        cout << "]" << endl;
    }
}

// Affichage du vecteur b
void Solve::displayVect() {
    cout << "[";
    for(int i(0); i<dim; i++) {
        cout << " " << b[i] << " ";
    }
    cout << "]" << endl;
}

// Pivotage de gauss
void Solve::gaussElim() {
    float *temp = new float[dim];

    // On recherche le max dans une colonne
    for(int k(0); k<dim; k++) {
        float maxNbr = fabs(A[k][k]);   // pour la comparaison du max
        int ind = k;                    // pour stocker l'indice

        for(int r(k+1); r<dim; r++) {
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
        for(int i(k+1); i<dim; i++) {
            for(int j(k+1); j<dim; j++) {
                A[i][j] = A[i][j] - (A[i][k] / A[k][k]) * A[k][j];
            }
            b[i] = b[i] - (A[i][k] / A[k][k]) * b[k];
            A[i][k] = 0;
        }
    }
}

void Solve::solveTriang()
{
    for(int i(dim-1); 0<=i; i--) {
        float s = 0;

        // calcul de la somme a_ij * x_j
        for(int j(i+1); j<dim; j++) {
            s += A[i][j] * x[j];
        }
        x[i] = (b[i] - s) / A[i][i];
    }

    // affichage du résultat
    cout << "[";
    for(int i(0); i<dim; i++) {
        cout << " " << x[i] << " ";
    }
    cout << "]";
}

// Destructeur
Solve::~Solve() {
    delete[] A;
    delete[] b;
    delete[] x;
}
