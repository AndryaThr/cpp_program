#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;

/***************************************************
 *                    PROTOTYPES                   *
 ***************************************************/

void display(vector<float>);
void display(vector<size_t>);

/***************************************************
 *                    CLASSES                      *
 ***************************************************/

/** \brief résolution d'une équation matricielle Ax = b
 *         avec A définie positive symétrique
 *         (Utilisation du stockage profil)
 *
 *
 * \param AP: rangememt plein de la matrice A
 * \param nDiag: numéro de l'élément diagonal de chaque ligne
 * \param l: demi-largeur du profil
 * \param p: colonne de 1er élément de chaque ligne
 *
 * \param b: vecteur second membre du système
 * \param x: solution
 *
 * \param N: taille de A
 * \param K: dimension du profil de A
 */
class SolveLDLt {
private:
    vector<float> AP,
                  b,
                  x;
    vector<size_t> nDiag,
                   l,
                   p;
    size_t N,
           K;

public:
    /* constructeur et destructeur */
    SolveLDLt(string, string);
    ~SolveLDLt() {}

    /* calcul sur les index i,j */
    float Aij(size_t, size_t);                  // retourne A[i][j]
    size_t index(size_t, size_t);               // retourne index tq A[i][j] = AP[index]

    /* décomposition LDLt */
    void factor_LD();

    /* résolution d'équations */
    void solve();
};


/***************************************************
 *                       MAIN                      *
 ***************************************************/

int main() {
    cout << "+----------------------------------+" << endl;
    cout << "|        FACTORISATION LDLt        |" << endl;
    cout << "+----------------------------------+" << endl;

/// données
    SolveLDLt syst("A.txt", "b.txt");

/// traitement et résultats
    syst.factor_LD();
    syst.solve();

    return 0;
}

/***************************************************
 *                   SOLVELDLt                     *
 ***************************************************/

/** \brief constructeur
 * lecture de la matrice A dans un fichier
 *
 * \param filename: string, nom du fichier contenant la matrice A et sa dimension
 */
SolveLDLt::SolveLDLt(string m_file, string v_file) {
    /* ouverture et test d'ouverture des fichiers */
    ifstream matrix(m_file);
    ifstream vect(v_file);

    if(!matrix.is_open() || !vect.is_open()) {
        cout << "Can't find one of the files ! Exiting program" << endl;
        exit(EXIT_FAILURE);
    }

    /* lecture du fichier et initialisation des matrices */
    matrix >> this->N;         // taille de la matrice

    /* remplissage des tableaux du rangement profil de A et du vecteur v*/
    float tmp(0);
    this->K = 0;

    for(size_t i(0); i<this->N; i++) {
        // remplissage du vecteur b
        this->x.push_back(0);
        vect >> tmp;    this->b.push_back(tmp);

        // rangement profil de A
        for(size_t j(0); j<i+1; j++) {
            matrix >> tmp;
            if(tmp!=0) {
                this->p.push_back(j+1);
                this->AP.push_back(tmp);
                this->K++;
                for(size_t l(j+1); l<i+1; l++) {
                    matrix >> tmp; this->AP.push_back(tmp);
                    this->K++;
                }
                this->nDiag.push_back(this->K);
                break;
            }
        }
        if(i==0)
            this->l.push_back(0);
        else
            this->l.push_back(this->nDiag[i] - this->nDiag[i-1] - 1);
    }
    cout << "Voici le stockage profil de la matrice A : " << endl;
    cout << "\nAP = "; display(this->AP);
    cout << "l = "; display(this->l);
    cout << "p = "; display(this->p);
    cout << "nDiag = "; display(nDiag);

    matrix.close();
    vect.close();
}

/** \brief factorisation LDLt de A et stockage profil du résultat
 *
 */
void SolveLDLt::factor_LD() {
    float sum(0);

    for(size_t i(0); i<this->N; i++) {
        for(size_t j(this->p[i]-1); j<i; j++) {
            sum = 0;
            for(size_t l(this->p[i]-1); l<j; l++)
                sum += this->Aij(i,l) * this->Aij(l,l) * this->Aij(j,l);

            this->AP[this->index(i,j)] = (this->Aij(i,j) - sum) / this->Aij(j,j);

        }
        sum = 0;
        for(size_t l(this->p[i]-1); l<i; l++)
            sum += this->Aij(i,l) * this->Aij(l,l) * this->Aij(i,l);

        this->AP[this->nDiag[i]-1] -= sum;
    }
}

/** \brief résolution de l'équation matricielle Ax = b toujours en utilisant le stockage profil
 *
 */

void SolveLDLt::solve() {
    cout << "\nVoici le stockage profil des matrices L et D : " << endl;
    cout << "\nL_DP = "; display(this->AP);
    cout << "l = "; display(this->l);
    cout << "p = "; display(this->p);
    cout << "nDiag = "; display(this->nDiag);
    /* RESOLUTION DE L.D.Lt.x = b */
    float sum(0);

    /* 1. résolution de L.z = b, L matrice triangulaire inférieure */
    for(size_t i(0); i<this->N; i++) {
        sum = 0;
        for(size_t j(0); j<i; j++) {
            sum += this->Aij(i,j) * this->x[j];
        }
        this->x[i] = (this->b[i] - sum);
    }

    /* 2. résolution de D.y = z, D matrice diagonale */
    for(size_t i(0); i<this->N; i++) {
        this->x[i] /= this->Aij(i,i);
    }

    /* 3. résolution de Lt.x = y, Lt matrice triangulaire supérieure */
    for(size_t i=this->N-1; ; i--) {
            sum = 0;
            for(size_t j=this->N-1; j>i; j--) {
                sum += this->Aij(j,i) * this->x[j];
            }
            this->x[i] -= sum;
            if(i==0)
                break;
        }

    cout << "Voici la solution de l'equation : " << endl;
    cout << "x = "; display(this->x);
}

/** \brief retourne A[i][j] dans AP (Aij = 0 si (i,j) n'appartient pas au profil)
 *
 * \param i: ligne
 * \param j: colonne
 *
 * \return A[i][j]: float
 */
float SolveLDLt::Aij(size_t i, size_t j) {
    // Aij = 0 si (i,j) n'appartient pas au profil
    if(j<this->p[i]-1)
        return 0;
    return this->AP[(this->nDiag[i] - 1) - i + j];
}

/** \brief retourne index de A[i][j] dans AP
 *
 * \param i: ligne
 * \param j: colonne
 *
 * \return index: AP[index] = A[i][j]
 */
size_t SolveLDLt::index(size_t i, size_t j) {
    return (this->nDiag[i] - 1) - i + j;
}

/***************************************************
 *                    FONCTIONS                    *
 ***************************************************/

/** \brief affiche un vecteur
 *
 * \param V: Vecteur
 */
void display(vector<float> V) {
    size_t i(0);
    cout << "[";
    for(i=0; i<V.size(); i++) {
        cout << " " << V[i] << " ";
    }
    cout << "]" << endl;
    cout << endl;
}

/** \brief affiche un vecteur
 *
 * \param V: Vecteur
 */
void display(vector<size_t> V) {
    size_t i(0);
    cout << "[";
    for(i=0; i<V.size(); i++) {
        cout << " " << V[i] << " ";
    }
    cout << "]" << endl;
    cout << endl;
}

