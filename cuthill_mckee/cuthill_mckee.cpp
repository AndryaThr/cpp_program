#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iomanip>

/**
 * \author Andriamanantsoa Herilanto Tsitohaina Rasamoelina
 * \version 1.0
 */

using namespace std;

/***************************************************
 *                    PROTOTYPES                   *
 ***************************************************/

template<class T> void display(vector<T>);
template<class T> void display(vector<vector<T>>);
template<class T> void display(T**, size_t);
template<class T> void display(T*, size_t);
template<class T> vector<T> substract(vector<T>, vector<T>);
template<class T> bool contains(T, vector<T>);
bool done(vector<bool>);
size_t id(size_t, vector<size_t>);

/***************************************************
 *               		SolveLDLt                  *
 ***************************************************/
/**
 * \brief Optimisation de la résolution de l'équation A.x = b
 * Utilisation de l'algorithme de Cuthill McKee pour optimiser
 * le stockage de la matrice A par modification de son profil
 *
 * \param AP, nDiag, p: pour stocker le profil de la matrice A (vector, size_t*)
 * \param N: taille de A (size_t)
 * \param K: dimension de profil(A) (size_t)
 *
 * \param x: vecteur solution (float*)
 * \param b: vecteur second membre (float*)
 *
 * \param neighbors: tableau contenant le voisin de chaque sommet (vector)
 * \param permutation: tableau contenant la permutation des sommets du profil (vector)
 * \param edge: sommet de départ de la permutation (size_t)
 *
 * \param optim: true- utilise l'optimisation de Cuthill McKee
 *               false- factorisation LDLt et résolution sans optimisation
 */
class SolveLDLt {
/// méthodes
public:
    /* constructeur et destructeur */
    SolveLDLt(string, string, bool, bool=false);
    ~SolveLDLt();

    /* calcul sur les index i,j */
    float Aij(size_t, size_t);                                  // retourne A[i][j]
    size_t index(size_t, size_t);                               // retourne index tq A[i][j] = AP[index]
    size_t P(size_t i, size_t j);                               // retourne P[i][j]

    /* graphe et profil */
    void searchNeighbors();                                     // voisins de chaque sommet
    size_t excentricity(size_t, vector<vector<size_t>>&);       // calcul de l'excentricité d'un sommet
    size_t free_node(size_t, vector<bool>);                     // calcul du nombre de noeuds enfants non numérotés

    /* numérotation de Cuthill McKee */
    void findEdge();                                            // recherche du sommet excentrique
    void permuteEdge();                                         // modification du profil
	void permuteInverse();										// algorithme de Cuthill McKee inverse

    /* transformation de la matrice A */
    void transformProfil();                                     // nouveau stockage profil de A'

    /* optimise la matrice A */
    void optimize();                                            // étapes d'optimisation

    /* résolution de l'équation */
    void factor_LD();                                           // factorisation LDLt de AP
    void solveEquations();                                      // résolution des équations
    void transformX();                                          // si AP transformé, on transforme aussi le résultat x
    void solve();                                               // étapes de résolution

/// attributs
private:
    /* profil de la matrice A */
    vector<float> AP;
    size_t
		*nDiag,
        *p;
    size_t
		N,
        K;

    /* graphe et voisinage */
    vector<vector<size_t>> neighbors;
    vector<size_t> permutation;
    size_t edge;

    /* résolution de l'équation */
    float
		*b,
        *x;
    bool
        cuthill_optim,
        cuthill_inverse;
};

/***************************************************
 *                       MAIN                      *
 ***************************************************/

int main() {
    cout << "+----------------------------------+" << endl;
    cout << "|        FACTORISATION LDLt        |" << endl;
    cout << "+----------------------------------+" << endl;

/// données
    SolveLDLt
        syst1("A.txt", "b.txt", true),                  // utilisation de l'optimisation de Cuthill McKee -> factorisation LDLt -> Résolution
        syst2("A.txt", "b.txt", true, true),            // utilisation de l'optimisation de Cuthill McKee INVERSE -> factorisation LDLt -> Résolution
        syst3("A.txt", "b.txt", false);                 // factorisation LDLt -> Résolution

/// traitement et résultats
    cout << "--------------------------------------------------------------------------\n" << endl;
    cout << ">>>> OPTIMISATION DE CUTHILL MCKEE -> FACTORISATION LDLt -> RESOLUTION\n" << endl;
    syst1.solve();
    cout << "--------------------------------------------------------------------------\n" << endl;
    cout << ">>>> OPTIMISATION DE CUTHILL MCKEE INVERSE -> FACTORISATION LDLt -> RESOLUTION\n" << endl;
    syst2.solve();
    cout << "--------------------------------------------------------------------------\n" << endl;
    cout << ">>>> FACTORISATION LDLt -> RESOLUTION\n" << endl;
    syst3.solve();
    cout << "--------------------------------------------------------------------------\n" << endl;
    return 0;
}

/***************************************************
 *               		SolveLDLt                  *
 ***************************************************/

/**
 * \brief constructeur
 * lecture de la matrice A et le vecteur b depuis un fichier
 *
 * \param m_file: chemin du fichier de la matrice A (string)
 * \param v_file: chemin du fichier du second membre b (string)
 *
 * \return
 */
SolveLDLt::SolveLDLt(string m_file, string v_file, bool _optim, bool _inverse) {
    this->cuthill_optim = _optim;
    if(this->cuthill_optim) {
        this->cuthill_inverse = _inverse;
    }
    /* ouverture et test d'ouverture des fichiers */
    ifstream matrix(m_file);
    ifstream vect(v_file);

    if(!matrix.is_open()) {
        cout << "Can't find one of the files ! Exiting program" << endl;
        exit(EXIT_FAILURE);
    }

    /* lecture du fichier et initialisation de la matrice */
    matrix >> this->N;         // taille de la matrice

    /* remplissage des tableaux du rangement profil de A et du vecteur b */
    float tmp(0);
    this->K = 0;
    this->nDiag = new size_t[this->N];
    this->p = new size_t[this->N];

    this->b = new float[this->N];
    this->x = new float[this->N];

    for(size_t i(0); i<this->N; i++) {
        // vecteur second membre b
        vect >> tmp;
        b[i] = tmp;
        x[i] = 0;

        // rangement profil de A
        for(size_t j(0); j<i+1; j++) {
            matrix >> tmp;
            if(tmp!=0) {
                p[i] = j+1;
                this->AP.push_back(tmp);
                this->K++;
                for(size_t l(j+1); l<i+1; l++) {
                    matrix >> tmp; this->AP.push_back(tmp);
                    this->K++;
                }
                nDiag[i] = this->K;
                break;
            }
        }
    }
    matrix.close();
}

/**
 * \brief destructeur
 * free nDiag et p
 */
SolveLDLt::~SolveLDLt() {
    delete[] nDiag;
    delete[] p;
    delete[] b;
    delete[] x;
}

/**
 * \brief A[i][j]
 * retourne A[i][j] dans AP (Aij = 0 si (i,j) n'appartient pas au profil)
 *
 * \param i: ligne
 * \param j: colonne
 *
 * \return A[i][j]: float
 */
float SolveLDLt::Aij(size_t i, size_t j) {
    // Aij = 0 si (i,j) n'appartient pas au profil
    if(j<this->p[i]-1 || i<j)
        return 0;

    return this->AP[(this->nDiag[i] - 1) - i + j];
}

/**
 * \brief retourne index de A[i][j] dans AP
 *
 * \param i: ligne
 * \param j: colonne
 *
 * \return index: AP[index] = A[i][j]
 */
size_t SolveLDLt::index(size_t i, size_t j) {
    return (this->nDiag[i] - 1) - i + j;
}

/**
 * \brief retourne P[i][j]
 *
 * \param i: ligne (size_t)
 * \param j: colonne (size_t)
 *
 * \return P[i][j]
 */
size_t SolveLDLt::P(size_t i, size_t j) {
    return i==this->permutation[j] ? 1 : 0;
}

/**
 * \brief remplissage de neighbors
 * cherche le voisin de chaque sommet du graphe
 */
void SolveLDLt::searchNeighbors() {
    size_t i(0), j(0);

    for(i=0; i<this->N; i++) {
        // on parcourt chaque colonne pour trouver des valeurs non nulles
        this->neighbors.push_back(vector<size_t>());
        for(j=0; j<this->N; j++) {
            if(i<j) {
                if(Aij(j,i)!=0)
                    this->neighbors[i].push_back(j);
            } else if(i>j) {
                if(Aij(i,j)!=0)
                    this->neighbors[i].push_back(j);
            }
        }
    }
}

/**
 * \brief calcul de l'excentricité du sommet _edge
 * en utilisant les structures de niveaux
 *
 * \param _edge: sommet dont on veut calculer l'excentricité
 * \param struct_edge: structure de niveaux du sommet
 *
 * \return excentricité de _edge
 */
size_t SolveLDLt::excentricity(size_t _edge, vector<vector<size_t>> &struct_edge) {
    size_t i(0), j(0), k(0);

    for(i=0; ;i++) {
        // recherche de la structure de niveaux de _edge
        if(i==0) {
            struct_edge.push_back(vector<size_t>());
            struct_edge[0].push_back(_edge);
        } else {
            struct_edge.push_back(vector<size_t>());
            for(j=0; j<struct_edge[i-1].size(); j++) {
                size_t loc = struct_edge[i-1][j];

                for(k=0; k<neighbors[loc].size(); k++) {
                    if(!contains(neighbors[loc][k], struct_edge[i])) {
                        struct_edge[i].push_back(neighbors[loc][k]);
                    }
                }
                struct_edge[i] = substract(struct_edge[i], struct_edge[i-1]);
                if(i>1)
                    struct_edge[i] = substract(struct_edge[i], struct_edge[i-2]);
            }
        }

        if(struct_edge[i].empty())
            break;
    }
    return struct_edge.size() - 2;
}

/**
 * \brief noeuds non numérotés
 * calcul du nombre de noeuds non numérotés d'un noeud
 *
 * \param node: noeud dont on veut déterminer les enfants non numérotés (size_t)
 * \param visited: tableau de bool marquant les noeuds non numérotés (vector)
 *
 * \return le nombre de noeuds enfants non numérotés (size_t)
 */
size_t SolveLDLt::free_node(size_t node, vector<bool> visited) {
    size_t i(0), c(0);
    vector<size_t> source = this->neighbors[node];

    // on incrémente c à chaque noeud non numéroté
    for(i=0; i<source.size(); i++) {
        if(!visited[source[i]])
            c++;
    }
    return c;
}

/**
 * \brief trouve le sommet excentrique this->edge
 * en prenant un sommet au hasard pour calculer son excentricité
 */
void SolveLDLt::findEdge() {
    vector<vector<size_t>> struct_edge;
    vector<bool> visited(10, false);
    size_t i(0),
        e_i(0),
        _edge(-1);
    this->edge = 0;

    // arrêt de la recherche si tous les sommets sont parcourus ou
    // si le sommet avec la plus grande excentricité ne change pas
    while(!done(visited) && this->edge!=_edge) {
        e_i = this->excentricity(this->edge, struct_edge);
        _edge = this->edge;
        visited[this->edge] = true;

        for(i=0; i<struct_edge[e_i].size(); i++) {
            vector<vector<size_t>> tmp;
            size_t e_j = this->excentricity(struct_edge[e_i][i], tmp);
            if(e_i<e_j) {
                this->edge = struct_edge[e_i][i];
                visited[this->edge] = true;
                break;
            }
        }
        struct_edge.clear();
    }
    cout << "Le sommet peripherique final = " << this->edge << "\n" << endl;
}

/**
 * \brief permutation des sommets this->permutation
 * à partir de this->edge, on modifie le profil de A pour diminuer
 * sa dimension
 */
void SolveLDLt::permuteEdge() {
    vector<bool> visited(10, false);            // pour le marquage des noeuds non numérotés
    size_t i(0), j(0), k(0);

    /* 1er noeud */
    this->permutation.push_back(this->edge);
    visited[this->permutation[k]] = true;

    /* parcours tant que tout les noeuds ne sont pas numérotés */
    for(k=0; !done(visited); k++) {
        // noeuds enfants du sommet this->permutation[k]
        vector<size_t> child_nodes = this->neighbors[this->permutation[k]];

        // triage des noeuds enfants du sommet this->permutation[k] dans
        // l'ordre croissant selon leurs noeuds enfants non numérotés
        for(i=0; i<child_nodes.size(); i++) {
            for(j=i; j<child_nodes.size(); j++) {
                if(free_node(child_nodes[i], visited)>free_node(child_nodes[j], visited)) {
                    size_t tmp = child_nodes[i];
                    child_nodes[i] = child_nodes[j];
                    child_nodes[j] = tmp;
                }
            }
        }

        // ajout des noeuds dans le tableau contenant les permutations
        // si ils n'y sont pas encore
        for(i=0; i<child_nodes.size(); i++) {
            if(!contains(child_nodes[i], this->permutation)) {
                this->permutation.push_back(child_nodes[i]);
                visited[child_nodes[i]] = true;
            }
        }
    }
    cout << "Apres permutation = \n";
    cout << "[";
    for(i=0; i<N; i++) {
        cout << " " << i << " ";
    }
    cout << "]" << endl;
    display(this->permutation);
}

/**
 * \brief permutation inverse: permutation[i] = (N + 1) - permutation[i]
 */
void SolveLDLt::permuteInverse() {
	for(size_t i(0); i<this->N; i++) {
        this->permutation[i] = (this->N-1) - this->permutation[i];
	}
	display(this->permutation);
}

/**
 * \brief calcul de A' = tP.A.P et recherche de son nouveau stockage profil*
 */
void SolveLDLt::transformProfil() {
    size_t i(0), j(0), k(0);

    /* Calcul de la nouvelle matrice A' = tP.A.P */
    float** A = new float*[this->N];        // A[i][j] = AP[nDiag(i) - i + j)
    float** B = new float*[this->N];        // B = tP.A

    /* allocation en mémoire de A et B */
    for(i=0; i<this->N; i++) {
        A[i] = new float[this->N];
        B[i] = new float[this->N];
        for(j=0; j<i+1; j++) {
            A[i][j] = Aij(i,j);
            A[j][i] = Aij(i,j);

            B[i][j] = 0;
        }
    }

    /* B = tP.A */
    for(i=0; i<this->N; i++) {
        for(j=0; j<this->N; j++) {
            float sum(0);
            for(k=0; k<this->N; k++) {
                sum += P(k,i) * A[k][j];
            }
            B[i][j] = sum;
        }
    }

    /* A' = B.P */
    for(i=0; i<this->N; i++) {
        for(j=0; j<this->N; j++) {
            float sum(0);
            for(k=0; k<this->N; k++) {
                sum += B[i][k] * P(k,j);
            }
            A[i][j] = sum;
        }
    }

    /* création du nouveau stockage profil de A' */
    this->AP.clear();
    this->K = 0;

    for(i=0; i<this->N; i++) {
        for(j=0; j<=i; j++) {
            if(A[i][j]!=0) {
                this->p[i] = j+1;
                this->AP.push_back(A[i][j]);
                this->K++;
                for(k=j+1; k<i+1; k++) {
                    this->AP.push_back(A[i][k]);
                    this->K++;
                }
                this->nDiag[i] = this->K;
                break;
            }
        }
    }

    /* free A et B */
    for(i=0;i<N;i++) {
        delete[] A[i];
        delete[] B[i];
    }
    delete[] A;
    delete[] B;
}

/**
 * \brief optimisation de A: fait appel à toutes les méthodes nécessaires et calcul de b' = tP.b
 */
void SolveLDLt::optimize() {
    /* optimisation du profil de A */
    this->searchNeighbors();
    this->findEdge();
    this->permuteEdge();
    if(this->cuthill_inverse)
        this->permuteInverse();

    this->transformProfil();

    /* calcul de b' */
    /* copie de x */
    float *x_tmp = new float[this->N];
    for(size_t i=0; i<this->N; i++) {
        x_tmp[i] = b[i];
    }
    //display(b,10);
    for(size_t i=0; i<this->N; i++) {
        float sum(0);
        for(size_t j=0; j<this->N; j++) {
            sum += P(j,i) * x_tmp[j];
        }
        b[i] = sum;
    }

    delete[] x_tmp;
}

/**
 * \brief factorisation LDLt de A et stockage profil du résultat
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

/**
 * \brief résolution de l'équation matricielle Ax = b toujours en utilisant le stockage profil
 */
void SolveLDLt::solveEquations() {
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
}

/**
 * \brief après résolution de A'.x' = b' calcul x = P.x'
 */
void SolveLDLt::transformX() {
    size_t i(0), j(0);

    /* copie de x */
    float *x_tmp = new float[this->N];
    for(i=0; i<this->N; i++) {
        x_tmp[i] = x[i];
    }

    /* Calcul de x = P.x' */
    for(i=0; i<this->N; i++) {
        float sum(0);
        for(j=0; j<this->N; j++) {
            sum += this->P(i,j) * x_tmp[j];
        }
        x[i] = sum;
    }

    delete[] x_tmp;
}

/**
 * \brief étapes de la résolution de l'équation
 * appel à toutes les fonctions qui résolvent l'équation
 */
void SolveLDLt::solve() {
    if(this->cuthill_optim) {
        this->optimize();
    }

	if(this->cuthill_optim)
		cout << "\nApres optimisation du stockage profil de A : " << endl;
	else
		cout << "\nVoici le stockage profil de A : " << endl;

    cout << "\nL_DP = "; display(this->AP);
    cout << "p = "; display(this->p, this->N);
    cout << "nDiag = "; display(this->nDiag, this->N);
    cout << "dim = " << this->K << endl << endl;

    this->factor_LD();
    this->solveEquations();

    if(this->cuthill_optim) {
        this->transformX();
    }

    cout << "Voici la solution de l'equation : " << endl;
    cout << "x = "; display(this->x, this->N);
}

/***************************************************
 *                    FONCTIONS                    *
 ***************************************************/

/**
 * \brief affiche un vecteur (vector)
 *
 * \param V: Vecteur à afficher (vector)
 */
template<class T>
void display(vector<T> V) {
    size_t i(0);
    cout << "[";
    for(i=0; i<V.size(); i++) {
        cout << " " << V[i] << " ";
    }
    cout << "]" << endl;
}

/**
 * \brief affiche un vecteur (posize_teur)
 *
 * \param V: Vecteur à afficher (posize_teur sur tableau)
 */
template<class T>
void display(T *V, size_t dim) {
    size_t i(0);
    cout << "[";
    for(i=0; i<dim; i++) {
        cout << " " << V[i] << " ";
    }
    cout << "]" << endl;
}

/**
 * \brief affiche une matrice (vector)
 *
 * \param M: matrice à afficher (vector<vector>)
 */
template<class T>
void display(vector<vector<T>> M) {
    size_t i(0), j(0);
    for(i=0; i<M.size(); i++) {
        cout << "[";
        for(j=0; j<M[i].size(); j++) {
            cout << setw(4) << M[i][j];
        }
        cout << "]" << endl;
    }
    cout << endl;
}

/**
 * \brief affiche une matrice (posize_teur)
 *
 * \param M: matrice à afficher (posize_teur sur tableau)
 */
template<class T>
void display(T** M, size_t dim) {
    size_t i(0), j(0);
    for(i=0; i<dim; i++) {
        cout << "[";
        for(j=0; j<i+1; j++) {
            cout << setw(4) << M[i][j];
        }
        cout << "]" << endl;
    }
    cout << endl;
}

/**
 * \brief différence entre 2 ensembles
 * effectue l'opération A-B sur les ensembles
 *
 * \param A: 1er ensemble (vector)
 * \param B: 2eme ensemble (vector)
 *
 * \return l'ensemble A-B (vector)
 *
 */
template<class T>
vector<T> substract(vector<T> A, vector<T> B) {
    vector<size_t> tmp;

    for(size_t i(0); i<A.size(); i++) {
        if(!contains(A[i], B))
            tmp.push_back(A[i]);
    }

    return tmp;
}

/**
 * \brief nbr contenu dans A
 * vérifie si A contient nbr
 *
 * \param nbr: nombre à chercher (T)
 * \param A: tableau où chercher le nombre (vector)
 *
 * \return true si nbr trouvé, false sinon (bool)
 */
template<class T>
bool contains(T nbr, vector<T> A) {
    for(size_t i(0); i<A.size(); i++)
        if(A[i]==nbr)
            return true;

    return false;
}

/**
 * \brief tableau de bool true
 * vérifie si tous les éléments d'un tableau bool est true
 *
 * \param visited: tableau de bool (vector)
 *
 * \return true si tout est true, false sinon
 */
bool done(vector<bool> visited) {
    for(size_t i(0); i<visited.size(); i++) {
        if(!visited[i])
            return false;
    }
    return true;
}

/**
 * \brief retourne l'indice d'un nombre d'un tableau
 *
 * \param K: nombre dont on veut connaître l'indice (size_t)
 * \param i: tableau où l'on veut chercher le nombre (vector)
 *
 * \return indice du nombre s'il est dans le tableau et -1 sinon
 */
size_t id(size_t K, vector<size_t> t) {
    for(size_t i(0); i<t.size(); i++)
        if(t[i]==K)
            return i;
    return -1;
}
