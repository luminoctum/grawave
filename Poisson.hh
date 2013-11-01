#include "main.h"
#include "fftw3.h"

class PoissonSolver : FiniteMethod{
protected:
    int nrows, ncols;
public:
    PoissonSolver(){}
    PoissonSolver(int _nrows, int _ncols) :
        nrows(_nrows), ncols(_ncols){}
    Grid operator() (const Grid &var, const Boundary &varb){
        int i, j, k;
        double b1[var.rows], bN
    }
protected:
}

int _tridiagSolver(long n, double *a, double *b, double *c, double *r, double v[]){
	/* Purpose:	Solve tri-diagnol matrix
	 * Requirement: a,b,c are the sub-diag, main-diag, sub-diag, r is right hand side, v is result
	 * Reference:	wiki, tri-diagnol matrix */
	long i; double m;
	for (i=1; i<n; i++){
		m=a[i-1]/b[i-1];
		b[i]=b[i]-m*c[i-1];
		r[i]=r[i]-m*r[i-1];
	};
	v[n-1]=r[n-1]/b[n-1];
	for (i=n-2; i>=0; i--)
		v[i]=(r[i]-c[i]*v[i+1])/b[i];
};

int _PoissonSolver2(struct state domain, double **obj, double **result){
	/* Purpose: 	Perform fast poisson solver: (Laplace)result=obj
	 * Requirement: Periodic condition in x and solid wall in y
	 * Reference: 	Cheng Li */
	long i,j,k;
	double b1[domain.ny],bN[domain.ny];
	double fobj[domain.nx-2][domain.ny],fobjT[domain.ny][domain.nx-2];
	double coef[domain.ny][domain.nx-2],coefT[domain.nx-2][domain.ny];
	double diag1[domain.nx-3], diag2[domain.nx-2], diag3[domain.nx-3];

	for (i=0; i<domain.nx-2; i++){
		fftw_plan p1=fftw_plan_r2r_1d(domain.ny,obj[i+1],fobj[i],FFTW_FORWARD,FFTW_ESTIMATE);
		fftw_execute(p1);
		fftw_destroy_plan(p1);
	};
	fftw_plan p1=fftw_plan_r2r_1d(domain.ny,obj[0]     ,b1,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_plan p2=fftw_plan_r2r_1d(domain.ny,obj[domain.nx-1],bN,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p1);
	fftw_execute(p2);
	fftw_destroy_plan(p1);
	fftw_destroy_plan(p2);

	for (k=0; k<domain.ny; k++) for (i=0; i<domain.nx-2; i++) fobjT[k][i]=fobj[i][k]*domain.dx*domain.dx;
	for (k=0; k<domain.ny; k++){
		for (i=0; i<domain.nx-3; i++){
			diag1[i]=1; diag3[i]=1;
		};
		for (i=0; i<domain.nx-2; i++){
			diag2[i]=2*cos(2*PI*k/domain.ny)-4;
			if (i==0     ) fobjT[k][i]-=b1[k];
			if (i==domain.nx-3) fobjT[k][i]-=bN[k];
		};
		_tridiagSolver(domain.nx-2,diag1,diag2,diag3,fobjT[k],coef[k]);
	};
	for (k=0; k<domain.ny; k++) for (i=0; i<domain.nx-2; i++) coefT[i][k]=coef[k][i];
	for (i=0; i<domain.nx-2; i++){
		fftw_plan p1=fftw_plan_r2r_1d(domain.ny,coefT[i],coefT[i],FFTW_BACKWARD,FFTW_ESTIMATE);
		fftw_execute(p1);
		for (j=0; j<domain.ny; j++) result[i+1][j]=coefT[i][j]/domain.ny;
		fftw_destroy_plan(p1);
	};
};
