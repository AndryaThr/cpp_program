#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

// PROTOTYPES DES FONCTIONS
float norme(float *v, int dim);
float scal(float *x, float *y, int dim);
void displayVect(float *v, int dim);
void displayMatrix(float **, int);
void displayResult(float **A, float *x, int dim);

// CLASSE POUR LA METHODE DE PUISSANCES ITEREES
class IterPow {
    public:
        IterPow(string filename);
        ~IterPow();
        // getters
        float** getMatrix() { return A; }
        float** getVectP() { return vect; }
        float* getVectX() { return x; }
        int getDim() { return dim; }
        // méthode de puissances itérées
        float powerIteration(int pos);
        // méthode de déflation
        float deflationMethod(float lbd, int pos);

    private:
        float **A,          // pour la matrice
              **vect;       // pour les vecteurs propres
        float *x,
              *x0;          // x0
        int dim;
};

/* Méthodes IterPow */
IterPow::IterPow(string filename) {
    dim = 5;

    // allocation dynamique
    x = new float[dim];
    x0 = new float[dim];
    A = new float*[dim];
    vect = new float*[dim];
    if(!A || !x || !x0 || !vect) {
        cout << "Erreur d'allocation, le programme va quitter !" << endl;
    }

    for(int i(0); i<dim; i++) {
        A[i] = new float[dim];
        vect[i] = new float[dim];
    }

    ifstream m_file(filename);

    // récupération des données depuis un fichier
    if(m_file) {
        // matrice
        for(int i(0); i<5; i++) {
            for(int j(0); j<5; j++) {
                m_file >> A[i][j];
            }
        }
        // vecteur x0
        for(int j(0); j<5; j++) {
            m_file >> x[j];
            x0[j] = x[j];
        }
    }
    m_file.close();
}

IterPow::~IterPow() {
    for(int i(0); i<dim; i++) {
        delete[] A[i];
        delete[] vect[i];
    }
    delete[] A;
    delete[] vect;
    delete[] x;
    delete[] x0;
}

float IterPow::powerIteration(int pos) {
    float l1(0), l2(1),     // lambda
          eps(1.e-6);       // epsilon

    // initialisation du vecteur x
    for(int i(0); i<dim; i++) {
        x[i] = x0[i];
    }

    float *y = new float[dim];
    if(!y) {
        cout << "Erreur !" << endl;
        return 0;
    }
    float N = norme(x, dim),    // norme
              tmp(0);

    // calcul de lambda
    while(1) {
        l2 = l1;
        N = norme(x, dim);

        // calcul du vect y
        for(int i(0); i<dim; i++) {
            y[i] = x[i] / N;
        }

        // calcul du vect x
        for(int i(0); i<dim; i++) {
            for(int k(0); k<dim; k++) {
                tmp += A[i][k]*y[k];
            }
            x[i] = tmp;
            tmp = 0;
        }
        // calcul de lambda
        l1 = scal(x, y, dim);
        if(fabs(l1 - l2)<eps)
            break;
    }

    // affectation du vecteur U1 dans x
    for(int i(0); i<dim; i++) {
        x[i] = y[i];
        vect[pos][i] = y[i];
    }

    delete[] y;
    return l1;
}

float IterPow::deflationMethod(float lbd, int pos) {
    // calcul de la matrice B
    for(int i(0); i<dim; i++) {
        for(int j(0); j<dim; j++) {
            A[i][j] -=  lbd * x[i] * x[j];
        }
    }
    // calcul de la valeur propre
    return powerIteration(pos);
}

int main()
{
    cout << "Methode des puissances iterees : " << endl;
/// DONNEES
    IterPow matrice("matrice.txt");
    float *eigenValues = new float[matrice.getDim()];

    cout << "Matrice A : " << endl;
    displayMatrix(matrice.getMatrix(), 5);
    cout << endl;

/// CALCULS
    eigenValues[0] = matrice.powerIteration(0);
    for(int i(1); i<matrice.getDim(); i++) {
        eigenValues[i] = matrice.deflationMethod(eigenValues[i-1], i);
    }

/// RESULTATS
    displayResult(matrice.getVectP(), eigenValues, matrice.getDim());

    delete[] eigenValues;
    return 0;
}

/* IMPLEMENTATION DES FONCTIONS */
// Calcul de la norme d'un vecteur
float norme(float *v, int dim) {
    float result(0);
    for(int i(0); i<dim; i++) {
        result += v[i]*v[i];
    }
    return (sqrt(result));
}

// Calcul du produit scalaire de 2 vecteurs
float scal(float *x, float *y, int dim) {
    float result(0);
    for(int i(0); i<dim; i++) {
        result += x[i]*y[i];
    }
    return result;
}

// affichage d'un vecteur
void displayVect(float *v, int dim) {
    for(int i(0); i<dim; i++) {
        cout << "[  " << v[i] << "  ]" << endl;
    }
}

// affichage d'une matrice
void displayMatrix(float **A, int dim) {
    for(int i(0); i<dim; i++) {
           cout << "[";
        for(int j(0); j<dim; j++) {
            cout << setw(4) << A[i][j] ;
        }
        cout << "   ]" << endl;
    }
}

void displayResult(float **V, float *x, int dim) {
    cout << "Voici les valeurs propres : " << endl;
    for(int i(0); i<dim; i++) {
        cout << "lambda_" << i+1 << " = " << x[i] << endl;
    }

    cout << "\nVoici les vecteurs propres : " << endl;
    for(int i(0); i<dim; i++) {
        cout << "u_" << i+1 << " = ( ";
        for(int j(0); j<dim; j++) {
            cout << V[i][j];
            if(j!=4) {
                cout << ", ";
            }
        }
        cout << ")" << endl;
    }
}
