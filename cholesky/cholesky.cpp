/*
 * Résolution de système linéaire par Méthode Gauss
 * V 1.0
 * 2021-06-22
 * Auteur(s): xxx
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <locale>
#include <iomanip>

using namespace std;

/// General helpers
float **newMat(int rows, int cols);
template <class T>
	T *newVect(size_t dim);
void displayMat(size_t dim, float **A);
void displayVec(size_t dim, float *v);
void displayTransp(size_t dim, float **A);

/* CLASSE RESOLUTION */
class Lsolver {
    public:
	Lsolver(string filename);
	~Lsolver();

	// méthode Gauss
	void displayResult(float **M);      /// pour adapter la fonction à la matrice B

	// méthode Choleski
	void triangInf();
	void invertB();
	void solveCholesky();

	// méthode Gauss
	void gaussElim();
	void solveTriangSup(float **A);     /// pour adapter la fonction à la matrice B
	void solveTriangInf();

	// accesseurs
	size_t  getdim(){ return dim;}
	float** getMat(){ return A;}
	float** getMatTri(){ return B;}
	float*  getRhs(){ return b;}

	// mutateurs
	void setdim(size_t n) {dim = n;}
	void setMat(float **M);
	void setMatTri(float **M);
	void setRhs(float *v);
private:
	size_t dim;
	float  **A; 	// matrice du problème A.x=b
	float  **B;     // matrice triangulaire inférieur
	float  *b, *x;	// second membre et inconnu du problème
	int    *p;		// pointeur utilisé pour les permutations de ligne
};

/* CLASSE MENU */
class Menu {
    public:
	    // constructeur & destructeur
		Menu();
		~Menu();
        // affiche le menu
		int displayMenu();

    private:
		string *option;
		size_t choices;
};

/* CLASSE GESTION DE MENU */
class ManageMenu {
    public:
        ManageMenu() {solver = new Lsolver("data.txt");}
        ~ManageMenu() {delete solver;}
        void manage(int ch);

    private:
        Lsolver *solver;
};

/* METHODES LSOLVER */

Lsolver::Lsolver(string filename){
/// ouverture du fichier de données
	size_t i(0), j(0);
    ifstream fichier(filename, ios::in);
    if(fichier){
        fichier >> dim;
/// allocation de la matrice, du second membre et de la solution
        A = newMat(dim, dim);
        B = newMat(dim, dim);
        b = newVect<float> (dim);
        x = newVect<float> (dim);
        p = newVect<int>  (dim);	// pointeur pour la permutation de lignes
/// remplissage des données
        for(i=0; i<dim; i++){
            for(j=0; j<dim; j++){
                fichier >> A[i][j];
                B[i][j] = 0;
            }
        }
        for(i=0; i<dim; i++){
            fichier >> b[i];
		}
        for(i=0; i<dim; i++){
            p[i] = i;
            x[i] = 0;
		}
        fichier.close();
    }
    else{
        cout << "Données non trouvées..." << endl;
    }
}

Lsolver::~Lsolver(){
	delete[] p;
	delete[] x;
	delete[] b;
	for(size_t i=0; i< dim; i++)
		delete[] A[i];
	delete[] A;
	for(size_t i=0; i< dim; i++)
		delete[] B[i];
	delete[] B;
}

void Lsolver::displayResult(float **M){
	float eps(1e-6);
    cout << "\nLe problème triangularisé:" << endl;
    for(size_t i=0; i<dim; i++){
        for(size_t j=0; j<dim-1; j++){
            if(fabs(M[p[i]][j])>eps || j>= i) cout << M[p[i]][j] << ".x" << j+1 << " + ";
            else cout << "       ";
        }
        cout << M[p[i]][dim-1] << ".x" << dim << " = " << b[p[i]] << endl;
    }
    cout << "\nLa solution:" << endl;
    for(size_t i=0; i<dim; i++)
        cout << "x" << i+1 << " = " << x[i] << endl;
}

void Lsolver::triangInf() {
    float s = 0;

    for(int i(0); i<(int) dim; i++) {
        for(int j(0); j<(int) dim; j++) {
            if(j<i) {
                for(int k(0); k<j; k++) {
                    s += B[i][k] * B[j][k];
                }
                B[i][j] = (A[i][j] - s) / B[j][j];
            } else if(i == j) {
                for(int k(0); k<i; k++) {
                    s += pow(B[i][k], 2);
                }
                B[i][i] = sqrt(A[i][i] - s);
            } else { B[i][j] = 0; }                     // matrice triangulaire inférieure
            s = 0;
        }
    }
}

void Lsolver::gaussElim(){
    size_t t(0), lpiv(0);    /// ligne de pivot courant
    float piv(0);
    for(size_t k=0; k<dim; k++){

/// Recherche du plus grand pivot
        lpiv = k;
        piv = fabs(A[p[k]][k]);
        for(size_t i=k+1; i<dim; i++){
            if(fabs(A[p[i]][k]) > piv){
                lpiv = i;
                piv = fabs(A[p[i]][k]);
            }
        }

/// Permutation de lignes
        t = p[k];
        p[k] = p[lpiv];
        p[lpiv] = t;

/// Elimination sous la ligne de pivot
        for(size_t i=k+1; i<dim; i++){
			A[p[i]][k] /= A[p[k]][k];
            for(size_t j=k+1; j<dim; j++)
                A[p[i]][j] -= (A[p[i]][k] * A[p[k]][j]);
            b[p[i]] -= (A[p[i]][k] * b[p[k]]);
            A[p[i]][k] = 0;
        }
    }
}

void Lsolver::solveTriangSup(float **M){
    float s(0);
    int i(0), j(0);
    for(i=dim-1; i>=0; i--){    /// Must go backward
        for(j=i+1, s=0; j<int(dim); j++)
            s += (M[p[i]][j]*x[j]);
        x[i] = (b[p[i]]-s)/M[p[i]][i];
    }
}

void Lsolver::solveTriangInf(){
    float s(0);
    size_t i(0), j(0);
    for(i=0; i<dim; i++){    /// on commence par la première ligne
        for(j=0, s=0; j<i; j++) {
            s += (B[i][j]*b[j]);
        }
        b[i] = (b[i] - s) / B[i][i];
    }
}

void Lsolver::invertB() {
    float temp(0);
    for(size_t i(0); i<dim; i++) {
        for(size_t j(i); j<dim; j++) {
            temp = B[i][j];
            B[i][j] = B[j][i];
            B[j][i] = temp;
        }
    }
}

void Lsolver::solveCholesky() {
    // Résolution de l'équation B * y = b, on obtient y
	solveTriangInf();
	// On transpose la matrice B
	invertB();
	// On résout Bt * x = y, on obtient le résultat x
    solveTriangSup(B);
}

void Lsolver::setMat(float **M) {
    for(size_t i(0); i<dim; i++)
        for(size_t j(0); j<dim; j++)
            B[i][j] = M[i][j];
}

void Lsolver::setMatTri(float **M) {
    for(size_t i(0); i<dim; i++)
        for(size_t j(0); j<dim; j++)
            A[i][j] = M[i][j];
}

void Lsolver::setRhs(float *v) {
    for(size_t i(0); i<dim; i++)
        b[i] = v[i];
}

/* METHODES MENU */

Menu::Menu(){
    choices = 3;

    option = newVect<string> (choices);
    option[0] = "Quit";
    option[1] = "Gauss Solve";
    option[2] = "Cholesky Solve";
}

Menu::~Menu() {}

int Menu::displayMenu() {
    int ch(0);

    cout << endl;
    for(size_t i(1); i<choices; i++)
        cout << i << " : " << option[i] << endl;
    cout << "0 : " << option[0] << endl;

    cout << "Your choice : ";   cin >> ch;

    return ch;
}

/* METHODES MANAGE MENU */

void ManageMenu::manage(int ch) {
/// Données
    if(ch == 1 || ch == 2) {
        cout << "La matrice A : " << endl;
        displayMat(solver->getdim(), solver->getMat());
        cout << "Le second membre b : " << endl;
        displayVec(solver->getdim(), solver->getRhs());
        cout << endl;
    }

/// Calculs et résultats
    switch(ch) {
    /* CALCUL PAR LA METHODE DE GAUSS */
        case 1:
            solver->gaussElim();
            solver->solveTriangSup(solver->getMat());
            solver->displayResult(solver->getMat());
            cout << endl;
            break;

    /* CALCUL PAR LA METHODE DE CHOLESKY */
        case 2:
            cout << "Recherche de la matrice triangulaire inférieure B : " << endl;
            solver->triangInf();
            solver->solveCholesky();
            solver->displayResult(solver->getMatTri());
            cout << endl;
            break;
    /* QUITTER LE PROGRAMME */
        case 0:
            cout << "Vous avez quitté le programme !" << endl;
            return;
    /* AUTRES CHOIX */
        default:
            cout << "\nRespecter le menu !!" << endl;
            return;
    }

	delete solver;
	solver = new Lsolver("data.txt");
}

/* MAIN() */

int main(){
	cout << "Locale: " << setlocale(LC_ALL,"") << endl;
    cout << "Résolution d'un système linéaire A.x=b" << endl;

/// Données
    ManageMenu prog;
    Menu menu;
    int ch(-1);

/// Calculs et résultats
    while(ch != 0) {
        ch = menu.displayMenu();
        prog.manage(ch);
        cout << "----------------------------------------" << endl;
    }

    return 0;
}

void displayMat(size_t dim, float **A){
	cout << "[";
	for(size_t i=0; i<dim; i++){
		if(i==0)cout << "[";
		else    cout << " [";
		for(size_t j=0; j<dim; j++){
			if(j<dim-1)cout << setw(8) << setprecision(5) << A[i][j] << "  ";
			else cout << setw(8) << setprecision(5) << A[i][j];
		}
		if(i<dim-1) cout << "]" << endl;
		else 		cout << "]]" << endl;
	}
}

void displayTransp(size_t dim, float **A){
	cout << "[";
	for(size_t i=0; i<dim; i++){
		if(i==0)cout << "[";
		else    cout << " [";
		for(size_t j=0; j<dim; j++){
			if(j<dim-1)cout << setw(8) << setprecision(5) << A[j][i] << "  ";
			else cout << setw(8) << setprecision(5) << A[j][i];
		}
		if(i<dim-1) cout << "]" << endl;
		else 		cout << "]]" << endl;
	}
}

void displayVec(size_t dim, float *v){
	for(size_t i=0; i<dim; i++){
		cout << " [" << v[i] << "]" << endl;
	}
}

template <class T>
T *newVect(size_t dim){
    T *v(NULL);
    v = new T[dim];		// Qu'est-ce que je fais en cas d'erreur
    return v;
}

float **newMat(int rows, int cols){
    float **mat(NULL);
    mat = newVect<float *>(rows);

    for(int i=0; i<rows; i++)
        mat[i] = newVect<float>(cols);
    return mat;
}

