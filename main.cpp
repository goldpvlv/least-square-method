#include<iostream>
#include<math.h>
#include <conio.h>


using namespace std;


double  deviation(int k, double **A) {
	double sum = 0, div=0;

	for (int i = 0; i < k + 1; ++i) {            
		for (int j = 0; j < k + 1; ++j)
			sum += A[i][j] * A[i][j];
		div += sum;
		sum = 0;
	}

	return sqrt(div);
}

double ** zeroing(int k, double **A) {
	for (int i = 0; i < k + 1; i++) { 
		for (int j = 0; j < k + 1; j++)
			A[i][j] = 0;
	}
	return A;
}

double** product(int k, double **A, double **U) {
	
	double **c = new double*[k + 1];
	for (int i = 0; i < k + 1; i++)
		c[i] = new double[k + 1];

	for (int i = 0; i < k + 1; i++) {
		for (int j = 0; j < k + 1; j++) {
		c[i][j] = 0;
			for (int t = 0; t < k + 1; t++) 
			c[i][j] += A[i][t] * U[t][j];
		}
	}
	return c;

}


int main()
{
	int n, k, step;
	cout << "enter steps: ";
	cin >> step;
	cout << "enter number of points N: ";
	cin >> n;
	cout << "enter degree K: ";
	cin >> k;
	double *x = new double[n];
	double *y = new double[n];
	cout << "enter X: ";
	for (int i = 0; i < n; ++i) 
		cin >> x[i];

	cout << "enter Y: ";
	for (int i = 0; i < n; ++i) 
		cin >> y[i];


	double **A = new double*[k + 1]; // matrix A
	for (int i = 0; i < k + 1; i++)
		A[i] = new double[k + 1];

	double *b = new double[k + 1]; // vector b
	for (int i = 0; i < k + 1; i++)
		b[i] = 0;

	double *a = new double[k + 1]; // vector a
	for (int i = 0; i < k + 1; i++)
		a[i] = 0;
		
	double **F = new double*[k + 1]; //matrix Psi
	for (int i = 0; i < k + 1; i++)
		F[i] = new double[k + 1];

	double **E = new double*[k + 1]; //identity matrix
	for (int i = 0; i < k + 1; i++)
		E[i] = new double[k + 1];
	
	double **c = new double*[k + 1];
	for (int i = 0; i < k + 1; i++)
		c[i] = new double[k + 1];

	A = zeroing(k, A);

	for (int i = 0; i < k + 1; i++) { //enter matrix A
		for (int j = 0; j < k + 1; j++)
			for (int t = 0; t < n; t++) 
				A[i][j] += pow(x[t], i + j);		
	}


	for (int i = 0; i < k + 1; ++i) { //enter vector b
		for (int j = 0; j < n; ++j) {
			b[i] += pow(x[j], i)*y[j];
		}
	}

	double **U = new double*[k + 1]; //approximation. 
	for (int i = 0; i < k + 1; i++)
		U[i] = new double[k + 1];

	c = zeroing(k, c);

	c = product(k, A, A);

	for (int i = 0; i < k + 1; ++i) { // first approximation
		for (int j = 0; j < k + 1; ++j)
			U[i][j] = A[i][j] / deviation(k,c);
	}

	for (int i = 0; i < k + 1; i++) { // enter identity matrix
		for (int j = 0; j < k + 1; j++) {
			if (i == j)
				E[i][j] = 1;
			else
				E[i][j] = 0;
		}
	}

	double eps = 0.00001;
	int oks = 0;


	do {                                    //Shultz's method

		c = zeroing(k, c);

		c = product(k, A, U);
		
		for (int i = 0; i < k + 1; ++i) {
			for (int j = 0; j < k + 1; ++j)
				F[i][j] = E[i][j] - c[i][j];             //find Psi
		}
		
		for (int i = 0; i < k + 1; ++i) {
			for (int j = 0; j < k + 1; ++j) {
				c[i][j] = 0;
				c[i][j] = E[i][j] + F[i][j];
			}
		}

		U = product(k, U, c);

		oks++;

	} while ((oks < step ) and (deviation(k, F) > eps));

	for (int i = 0; i < k + 1; i++) {
		a[i] = 0;
		for (int j = 0; j < k + 1; j++)
			a[i] += U[i][j] * b[j];                        //find vector a by applying inverse matrix U and vector b
	}

	cout << "vector a"<<endl;
	for (int i = 0; i < k + 1; i++) {                    // output a
		cout << a[i] << " ";
	}

	_getch();
	return 0;
}
