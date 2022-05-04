#include <iostream>
#include <cmath>
#include <cstdlib>

#define PI 3.14159265358979

using namespace std;

/**
 *
 * \u = x^2 * y^3
 * \f = 2*y^3 + 6*yx^2
 *
 */

double f(double, double);
double u(double, double);

void four1(double *data,unsigned long nn,int isign);
void realft(double *data,unsigned long n, int isign);
void sinft(double *y, unsigned long n);

/* --- STRUCT --- */
struct Vect {
    /* ----- */
    double* v;
    size_t dim;

    /* ----- */
    Vect(size_t);
    ~Vect();
    void display();
};

struct Mat {
    /* --- */
    double** m;
    size_t row,
           col;

    /* --- */
    Mat(size_t, size_t);
    ~Mat();
    void display();

    Mat T();
    //Mat operator*(Mat);
    Vect operator*(Vect);
    void operator=(Mat);
};

/* --- CLASS --- */
class Poisson2d {
    public:
        Poisson2d(double, double, size_t*);
        ~Poisson2d();
        void eigenValues();
        void S_i();
        void netaCompute();
        void TCompute();

        void gaussElim(Mat M, Vect v);
        void solveTriang(Mat M, Vect v);

    private:
        size_t M, N;
        double x_min, x_max,
               y_min, y_max,
               dx, dy,
               alpha, beta, gamma;
        Vect *d;
        Mat *U,
            *F,
            *S,
            *T;
};

int main() {
    size_t MN[2] = {8, 8};
    Poisson2d p(5, 5, MN);
    p.eigenValues();
    p.S_i();
    p.TCompute();
    return 0;
}

/* --- VECT --- */
Vect::Vect(size_t d) {
    dim = d;

    v = new double [dim];
    if(v == nullptr) {
        cout << "Vector cannot be created, exiting program!" << endl;
        exit(-1);
    }

    for (size_t i(0); i<dim; i++)
		v[i] = 0;
}

Vect::~Vect() {
    delete[] v;
}

void Vect::display() {
    cout << "[  ";
    for (size_t i(0); i<dim; i++)
		cout << v[i] << "  ";

	cout << "]\n" << endl;
}

/* --- MAT --- */
Mat::Mat(size_t r, size_t c) {
    row = r;    col = c;

    m = new double* [row];
    if(m == nullptr) {
        cout << "Matrix cannot be created, exiting program!" << endl;
        exit(-1);
    }

    for (size_t i(0); i<row; i++) {
		m[i] = new double[col];
		if(m[i] == nullptr) {
            cout << "Matrix cannot be created, exiting program!" << endl;
            exit(-1);
        }
	}

	for (size_t i(0); i<row; i++)
		for (size_t j(0); j<col; j++)
            m[i][j] = 0;
}

Mat::~Mat() {
    /*cout << "deleted : ( " << row << ", " << col << " )" << endl;
    display();*/
    for (size_t i(0); i<row; i++)
		delete[] m[i];

    delete[] m;
}

void Mat::display() {
	cout << "[";
	for (size_t i(0); i<row; i++) {
        if(i==0)
            cout << "[";
        else
            cout << " [";
		for (size_t j(0); j<col; j++) {
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
    for(size_t i(0); i<row; i++) {
        for(size_t j(0); j<col; j++) {
            x_t.m[j][i] = m[i][j];
        }
    }
    return x_t;
}

/* --- POISSON2D --- */
Poisson2d::Poisson2d(double x, double y, size_t* MN) {
    x_min = 0;   x_max = x;
    y_min = 0;   y_max = y;
    M = MN[0];    N = MN[1];

    dx = (x_max/(M+1));
    dy = (y_max/(N+1));

    alpha = -1/(dx*dx);
    gamma = -1/(dy*dy);
    beta = -2*(alpha + gamma);

    d = new Vect(M);
    U = new Mat(M+2,N+2);
    F = new Mat(M+2,N+2);
    T = new Mat(N, N);
    S = new Mat(N, M);

    for (size_t i = 0; i <= M+1; i++) {
		for (size_t j = 0; j <= N+1; j++) {
			/*if (i == 0 || i == M+1 || j == 0 || j == N+1)*/
				U->m[i][j] = u(x_min + i*dx, y_min + j*dy);
				F->m[i][j] = f(x_min + i*dx, y_min + j*dy);
			/*else
				U->m[i][j] = 0;
				F->m[i][j] = 0;
            */
		}
	}

	//cout << M << " ~ " << N << endl;
	//cout << dx << " ~ " << dy << endl;
	//cout << "alpha = " << alpha << " beta = " << beta << " gamma = " << gamma << endl;

	//F->display();
}

Poisson2d::~Poisson2d() {
    d->~Vect();
    U->~Mat();
    S->~Mat();
}

void Poisson2d::eigenValues() {
    // d_i and delta_i

    if(d==nullptr)
        exit(-1);

    double v = dx/dy;

    for(size_t i(0); i<M; i++) {
        d->v[i] = 2 + 2 * (v*v) * (1 - ( cos((PI * i) / (M+1) ) ) );
    }
    /*
    for(size_t i(0); i<M; i++) {
        for(size_t j(0); j<M; j++) {
            cout << " " << d->m[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    for(size_t i(0); i<M; i++) {
        for(size_t j(0); j<M; j++) {
            cout << " " << delta->m[i][j] << " ";
        }
        cout << endl;
    }
    */
}

void Poisson2d::S_i() {
    double y_tmp(0),
           x_tmp(0);

    for(size_t i(1); i<=N; i++) {
        y_tmp += dy;
        for(size_t j(1); j<=M; j++) {
            x_tmp += dx;
            if(i==1) {
                if(j==1) {
                    S->m[i-1][j-1] = f(x_tmp, y_tmp) - gamma*u(x_tmp-dx, y_tmp) - alpha*u(x_tmp, y_tmp-dy);
                } else if(j==M)
                    S->m[i-1][j-1] = f(x_tmp, y_tmp) - gamma*u(x_tmp-dx, y_tmp) - alpha*u(x_tmp, y_tmp+dy);
                else
                    S->m[i-1][j-1] = f(x_tmp, y_tmp) - gamma*u(x_tmp-dx, y_tmp);
            } else if(i==N) {
                if(j==1)
                    S->m[i-1][j-1] = f(x_tmp, y_tmp) - gamma*u(x_tmp+dx, y_tmp) - alpha*u(x_tmp, y_tmp-dy);
                else if(j==M)
                    S->m[i-1][j-1] = f(x_tmp, y_tmp) - gamma*u(x_tmp+dx, y_tmp) - alpha*u(x_tmp, y_tmp+dy);
                else
                    S->m[i-1][j-1] = f(x_tmp, y_tmp) - gamma*u(x_tmp+dx, y_tmp);
            } else {
                if(j==1)
                    S->m[i-1][j-1] = f(x_tmp, y_tmp) - alpha*u(x_tmp, y_tmp-dy);
                else if(j==M)
                    S->m[i-1][j-1] = f(x_tmp, y_tmp) - alpha*u(x_tmp, y_tmp+dy);
                else
                    S->m[i-1][j-1] = f(x_tmp, y_tmp);
            }
        }
        x_tmp = 0;
    }
}


/*
Mat Mat::operator*(Mat A) {
    if(col != A.row) {
        cout << "Cannot operate * on those two matrix!" << endl;
        exit(-1);
    }

    Mat r(row, A.col);
    for(size_t  i(0); i<r.row; i++) {
        for(size_t j(0); j<r.col; j++) {
            for(size_t k(0); k<A.row; k++) {
                r.m[i][j] += m[i][k] * A.m[k][j];
            }
        }
    }

    return r;
} */

void Poisson2d::TCompute() {

    for(size_t i(0); i<N; i++) {
        for(size_t j(0); j<M; j++) {
            S->m[i][j] *= sqrt( (M+1)/2 );
        }
    }

    for(size_t i(0); i<N; i++) {
        sinft(S->m[i], M);
    }

    S->T().display();

    for(size_t i(0); i<N; i++) {
        for(size_t j(0); j<N; j++) {
            for(size_t k(0); k<N; k++) {
                if(j==k)
                    T->m[j][k] = d->v[i];
                else if(j==k+1 || k==j+1)
                    T->m[j][k] = -1;
            }
        }
        T->display();
    }
}

Vect Mat::operator*(Vect A) {
    if(this->col != A.dim) {
        cout << "Cannot operate * !" << endl;
        exit(-1);
    }
    Vect r(row);
    for(size_t  i(0); i<row; i++) {
        for(size_t j(0); j<col; j++) {
            r.v[i] += m[i][j] * A.v[j];
        }
    }
    return r;
}

void Mat::operator=(Mat A) {
    for (size_t i(0); i<row; i++)
		delete[] m[i];

    delete[] m;

    row = A.row;    col = A.col;

    m = new double* [row];
    if(m == nullptr) {
        cout << "Matrix cannot be duplicated, exiting program!" << endl;
        exit(-1);
    }

    for (size_t i(0); i<row; i++) {
		m[i] = new double[col];
		if(m[i] == nullptr) {
            cout << "Matrix cannot be duplicated, exiting program!" << endl;
            exit(-1);
        }
	}

	for (size_t i(0); i<row; i++)
		for (size_t j(0); j<col; j++)
            m[i][j] = A.m[i][j];
}

/* --- FONCTIONS --- */
/*
double f(double x, double y) {
	return (2*y*y*y + 6*x*x*y);
}

double u(double x, double y) {
	return (x*x*y*y*y);
}
*/

double f(double x, double y) {
	return -(x*x + y*y - 2) * exp((-1) * (x*x + y*y) / 2);
}

double u(double x, double y) {
	return exp( (-1) * (x*x + y*y) / 2);
}

void Poisson2d::gaussElim(Mat M, Vect v) {
    /*float *temp = new float[n-1];

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
    }*/
}

void Poisson2d::solveTriang(Mat M, Vect v) {
    /*for(int i(n-1-1); 0<=i; i--) {
        float s = 0;

        // calcul de la somme a_ij * x_j
        for(int j(i+1); j<n-1; j++) {
            s += A[i][j] * x[j];
        }
        x[i] = (b[i] - s) / A[i][i];
    }*/
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
