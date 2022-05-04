#include <iostream>
#include <cmath>
#include <cstdlib>

#define PI 3.14159265358979

/************************************************
 *                  PROTOTYPES                  *
 ***********************************************/

double f(double, double);
double u(double, double);
double i_(int, double, double);
void four1(double *data,unsigned long nn,int isign);
void realft(double *data,unsigned long n, int isign);
void sinft(double *y, unsigned long n);
void solveGauss(double**, double*, double*, int);
void plot(double**, double, double, double, double, int, int);

using namespace std;

/************************************************
 *                    STRUCT                    *
 ***********************************************/
/* -------------------------------------------------------------------------- */
struct Vect {
    /* ----- */
    double* v;
    int dim;

    /* ----- */
    Vect(int);
    ~Vect();
    void display();
};

/* -------------------------------------------------------------------------- */
struct Mat {
    /* --- */
    double** m;
    int row,
           col;

    /* --- */
    Mat(int, int);
    ~Mat();
    void display();

    Mat T();
    //Vect operator*(Vect);
    void operator=(Mat);
};

/* -------------------------------------------------------------------------- */
class Poisson2d {
    public:
        // initialisation du problème
        Poisson2d(double*, double*, int*);
        ~Poisson2d();

        // étapes de la résolution du problème
        void eigenValues();         // delta et d
        void computeS();            // second membre S
        void computeNeta();         // neta et resolution de l'équation T*U = H
        void computeV();            // calcul de v(nu) et resolution du problème

    private:
        // nombres d'intervalles dur (0x) et (0y)
        int M, N;
        double x_min, x_max,        // bornes 0x
               y_min, y_max,        // bornes 0y
               dx, dy,              // longueur de chaque intervalle
               alpha, beta, gamma,
               delta_i;             // -1
        Vect *d,                    // vecteur contenant les diagonales de D
             *x_values,             // valeurs de x_i et y_i
             *y_values;
        Mat *S,                     // second membre
            *T,                     // matrice tridiag
            *H,                     // matrice second membre de TU = H
            *U,                     // inconnu 1
            *U_init,                // solution données par l'équation de départ
            *V,                     // matrice (nu)
            *Y;
};

/************************************************
 *                    MAIN                      *
 ***********************************************/
int main() {
    double x[2] = {-5,5},
           y[2] = {-5,5};
    int dim[2] = {128, 128};
    Poisson2d pb(x, y, dim);

    pb.eigenValues();
    pb.computeS();
    pb.computeNeta();
    pb.computeV();

    return 0;
}

/************************************************
 *                    VECT                      *
 ***********************************************/
Vect::Vect(int d) {
    dim = d;

    v = new double [dim];
    if(v == nullptr) {
        cout << "Vector cannot be created, exiting program!" << endl;
        exit(-1);
    }

    for (int i(0); i<dim; i++)
		v[i] = 0;
}

Vect::~Vect() {
    delete[] v;
}

void Vect::display() {
    cout << "[  ";
    for (int i(0); i<dim; i++)
		if(i<dim-1)
            cout << v[i] << " | ";
        else
            cout << v[i] << "  ";

	cout << "]\n" << endl;
}


/************************************************
 *                  MATRIX                      *
 ***********************************************/
Mat::Mat(int r, int c) {
    row = r;    col = c;

    m = new double* [row];
    if(m == nullptr) {
        cout << "Matrix cannot be created, exiting program!" << endl;
        exit(-1);
    }

    for (int i(0); i<row; i++) {
		m[i] = new double[col];
		if(m[i] == nullptr) {
            cout << "Matrix cannot be created, exiting program!" << endl;
            exit(-1);
        }
	}

	for (int i(0); i<row; i++)
		for (int j(0); j<col; j++)
            m[i][j] = 0;
}

Mat::~Mat() {
    for (int i(0); i<row; i++)
		delete[] m[i];

    delete[] m;
}

void Mat::display() {
	cout << "[";
	for (int i(0); i<row; i++) {
        if(i==0)
            cout << "[";
        else
            cout << " [";
		for (int j(0); j<col; j++) {
			cout << " " << m[i][j];
		}
		if(i==row-1)
            cout << " ]";
        else
            cout << " ]\n";
	}
	cout << "]\n" << endl;
}

Mat Mat::T() {
    Mat x_t(col, row);
    for(int i(0); i<row; i++) {
        for(int j(0); j<col; j++) {
            x_t.m[j][i] = m[i][j];
        }
    }
    return x_t;
}

void Mat::operator=(Mat A) {
    for (int i(0); i<row; i++)
		delete[] m[i];

    delete[] m;

    row = A.row;    col = A.col;

    m = new double* [row];
    if(m == nullptr) {
        cout << "Matrix cannot be duplicated, exiting program!" << endl;
        exit(-1);
    }

    for (int i(0); i<row; i++) {
		m[i] = new double[col];
		if(m[i] == nullptr) {
            cout << "Matrix cannot be duplicated, exiting program!" << endl;
            exit(-1);
        }
	}

	for (int i(0); i<row; i++)
		for (int j(0); j<col; j++)
            m[i][j] = A.m[i][j];
}

/************************************************
 *                  POISSON 2D                  *
 ***********************************************/
Poisson2d::Poisson2d(double* x, double* y, int* dim) {
    // bornes et dimensions
    x_min = x[0];   x_max = x[1];
    y_min = y[0];   y_max = y[1];
    M = dim[0];     N = dim[1];

    // longueur d'intervalles
    dx = ((x_max - x_min) / (M + 1));
    dy = ((y_max - y_min) / (N + 1));

    // coefficients
    alpha = -1/(dx*dx);
    gamma = -1/(dy*dy);
    beta = -2*(1 / (dx*dx) + 1 / (dy*dy));

    // allocation dynamique
    d = new Vect(M);
    x_values = new Vect(M+2);
    y_values = new Vect(M+2);

    S = new Mat(N, M);
    U_init = new Mat(M+2, N+2);
    H = new Mat(M, N);
    T = new Mat(N, N);
    U = new Mat(M, N);
    V = new Mat(N, M);
    Y = new Mat(N+2, M+2);

    if(d==nullptr || S==nullptr || T==nullptr || x_values==nullptr || y_values==nullptr || H==nullptr || U==nullptr || V==nullptr || U_init==nullptr)
        exit(-1);

    // valeurs de x_i et y_i
    for (int i(1); i <= M+2; i++) {
		x_values->v[i-1] = i_(i-1, dx, x_min);
		y_values->v[i-1] = i_(i-1, dy, y_min);
	}

	// valeurs de u(x, y)
	for(int i(0); i<M+2; i++) {
        for(int j(0); j<N+2; j++) {
            U_init->m[i][j] = u(i_(i, dx, x_min), i_(j, dy, y_min));
        }
	}
}

void Poisson2d::eigenValues() {
    double v = dy/dx;
    for(int i(0); i<M; i++) {
        // calcul de d_i
        d->v[i] = 2 + (2 * (v*v) * (1 - ( cos((PI * (i)) / (M+1) ) ) ));
    }
    // delta_i = -1
    delta_i = -1;
    //d->display();-
}

void Poisson2d::computeS() {

    for(int i(1); i<=N; i++) {
        for(int j(1); j<=M; j++) {
                S->m[i-1][j-1] = f(x_values->v[i], y_values->v[j]);
            if(i==1) {
                S->m[i-1][j-1] -= gamma * u(x_values->v[0], y_values->v[j]);
            } if(j==1) {
                S->m[i-1][j-1] -= alpha * u(x_values->v[i], y_values->v[0]);
            } if(i==N) {
                S->m[i-1][j-1] -= gamma * u(x_values->v[N+1], y_values->v[j]);
            } if(j==M) {
                S->m[i-1][j-1] -= gamma * u(x_values->v[i], y_values->v[M+1]);
            }
        }
    }
    //S->display();
}

void Poisson2d::computeNeta() {

    // S ----> neta
    for(int i(0); i<N; i++) {
        for(int j(0); j<M; j++) {
            S->m[i][j] *= dy * dy;
        }

        sinft(S->m[i], M);

        for(int j(0); j<M; j++) {
            S->m[i][j] *= 2. / (double) (M+1);
        }
    }
    //S->display();

    // H = S_Transpose => neta
    for (int i(0); i<M; i++) {
		for (int j(0); j<N; j++) {
            H->m[i][j] = S->m[j][i];
		}
	}
	//H->display();

	// resolution  T*U = H
    for(int i(0); i<M; i++) {

        for(int j(0); j<N; j++) {
            for(int k(0); k<N; k++) {
                if(j==k)
                    T->m[j][k] = d->v[i];
                else if(j==k+1 || k==j+1)
                    T->m[j][k] = -1;
            }
        }
        //T->display();

        solveGauss(T->m, H->m[i], U->m[i], N);
    }
    //U->display();
}

void Poisson2d::computeV() {
    // transposé de U = V
    for(int i(0); i<N; i++) {
        for(int j(0); j<M; j++) {
            V->m[i][j] = U->m[j][i];
        }
    }
    //V->display();
    for(int i(0); i<N; i++) {
        sinft(V->m[i], M);
    }
    //V->display();

    // solution finale
    for(int i(1); i<=N; i++) {
        for(int j(1); j<=M; j++) {
            Y->m[i][j] = V->m[i-1][j-1];
        }
    }
    for(int j(0); j<M; j++) {
        Y->m[0][j] = u(x_values->v[0], x_values->v[j]);
    }

    for(int j(0); j<M; j++) {
        Y->m[N+1][j] = u(x_values->v[N+1], x_values->v[j]);
    }

    plot(Y->m, x_min, y_min, dx, dy, M, N);
}

Poisson2d::~Poisson2d() {
    delete d;
    delete x_values;
    delete y_values;
    delete U_init;
    delete S;
    delete T;
    delete H;
    delete U;
    delete V;
}
/* -------------------------------------------------------------------------- */
double f(double x, double y) {
	return -(x*x + y*y - 2) * u(x, y);
}

double u(double x, double y) {
	return exp( (-1) * (x*x + y*y) / 2);
}

double i_(int i, double step, double start) {
    return (start + i*step);
}

void sinft(double *y, unsigned long n){
    unsigned long  j,n2=n+2;
    double sum,y1,y2;
    double theta,wi=0.0,wr=1.0,wpi,wpr,wtemp;

    theta=PI/(double) n;
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    y--;            // Décalage pour utiliser la numérotaion math
    y[1]=0.0;
    for (j=2;j<=(n>>1)+1;j++) {
        wr=(wtemp=wr)*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
        y1=wi*(y[j]+y[n2-j]);
        y2=0.5*(y[j]-y[n2-j]);
        y[j]=y1+y2;
        y[n2-j]=y1-y2;
    }
    realft(y,n,1);
    y[1]*=0.5;
    sum=y[2]=0.0;
    for (j=1;j<=n-1;j+=2) {
        sum += y[j];
        y[j]=y[j+1];
        y[j+1]=sum;
    }
}

void realft(double *data,unsigned long n, int isign){
    unsigned long i,i1,i2,i3,i4,np3;
    double c1=0.5,c2,h1r,h1i,h2r,h2i;
    double wr,wi,wpr,wpi,wtemp,theta;

    theta=PI/(double) (n>>1);
    if (isign == 1) {
        c2 = -0.5;
        four1(data,n>>1,1);
    }
    else {
        c2=0.5;
        theta = -theta;
    }
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0+wpr;
    wi=wpi;
    np3=n+3;
    for (i=2;i<=(n>>2);i++) {
        i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
        h1r=c1*(data[i1]+data[i3]);
        h1i=c1*(data[i2]-data[i4]);
        h2r = -c2*(data[i2]+data[i4]);
        h2i=c2*(data[i1]-data[i3]);
        data[i1]=h1r+wr*h2r-wi*h2i;
        data[i2]=h1i+wr*h2i+wi*h2r;
        data[i3]=h1r-wr*h2r+wi*h2i;
        data[i4] = -h1i+wr*h2i+wi*h2r;
        wr=(wtemp=wr)*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
    }
    if (isign == 1) {
        data[1] = (h1r=data[1])+data[2];
        data[2] = h1r-data[2];
    }
    else {
        data[1]=c1*((h1r=data[1])+data[2]);
        data[2]=c1*(h1r-data[2]);
        four1(data,n>>1,-1);
    }
}

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
void four1(double *data,unsigned long nn,int isign){
    unsigned long n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta;
    double tempr,tempi;

    n=nn << 1;
    j=1;
    for (i=1;i<n;i+=2) {
        if (j > i) {
            SWAP(data[j],data[i]);
            SWAP(data[j+1],data[i+1]);
        }
        m=n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax=2;
    while (n > mmax) {
        istep=mmax << 1;
        theta=isign*(2*PI/mmax);
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        for (m=1;m<mmax;m+=2) {
            for (i=m;i<=n;i+=istep) {
                j=i+mmax;
                tempr=wr*data[j]-wi*data[j+1];
                tempi=wr*data[j+1]+wi*data[j];
                data[j]=data[i]-tempr;
                data[j+1]=data[i+1]-tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
}
#undef SWAP

void solveGauss(double** A, double* b, double* x, int n) {

    double *temp = new double[n];

    for(int k(0); k<n; k++) {

        double maxNbr = fabs(A[k][k]);
        int ind = k;

        for(int r(k+1); r<n; r++) {
            if(maxNbr < fabs(A[r][k])) {
                maxNbr = fabs(A[r][k]);
                ind = r;
            }
        }

        temp = A[k];
        A[k] = A[ind];
        A[ind] = temp;

        double temp2(0);
        temp2 = b[k];
        b[k] = b[ind];
        b[ind] = temp2;

        // pivotage
        for(int i(k+1); i<n; i++) {
            for(int j(k+1); j<n; j++) {
                A[i][j] = A[i][j] - (A[i][k] / A[k][k]) * A[k][j];
            }
            b[i] = b[i] - (A[i][k] / A[k][k]) * b[k];
            A[i][k] = 0;
        }
    }

    // resolution de la matrice triangulaire supérieure
    for(int i(n-1); i>=0; i--) {
        double s(0);
        for(int j(n-1); j>=i+1; j--) {
            s += A[i][j] * x[j];
        }
        x[i] = (b[i] - s) / A[i][i];
    }
}

void plot(double **Y, double x0, double y0, double dx, double dy, int M, int N) {

    FILE *GP = popen("C:\\maxima-5.45.1\\gnuplot\\bin\\gnuplot.exe -persist", "w");

    if (GP) { // If gnuplot is found

        // fprintf(GP, "set term wxrc size 800,600\n");
        fprintf(GP, "set yrange [-2.5:2.5]\n");
        fprintf(GP, "set xrange [-2.5:2.5]\n");
        fprintf(GP, "set zrange [0:1]\n");
        fprintf(GP, "set style data points\n");
        fprintf(GP, "set decimalsign \".\"\n");

/// Reusable internal data block
        fprintf(GP, "$data << EOF\n");
        for(int i = 0; i < N+2; i++) {
            for(int j = 0; j < M+2; j++) {
                fprintf(GP, "%f %f %f\n", x0 + (j)*dx, y0 + (i)*dy, Y[i][j]);
            }
        }
        fprintf(GP, "EOF\n");
        //exp( (-1) * (x*x + y*y) / 2),
        fprintf(GP, "splot exp( (-1) * (x*x + y*y) / 2), $data using 1:2:3 with lines\n");
/// Send commands to gnuplot and close pipe
        fflush(GP);
        pclose(GP);
    } else {
        cout << "gnuplot not found..." << endl;
    }
}
