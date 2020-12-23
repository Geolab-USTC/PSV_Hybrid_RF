#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
	Matrix A is specified as a[row][column]; that is
	a[row] is a pointer to an array that contains all
	the values for that row.
*/

#define SMALL 1.0e-100

double ** init_pivot(dimension,rows)
/*
	Initializes the system solver by providing an array to place
	the augmented matrix in.  Also holds the original row order
	for use by reaugment().
*/
int dimension;
int **rows;
{
	double **temp_ptr;
	int i;

	temp_ptr = (double **) malloc(sizeof(double *)*dimension);
	rows[0]  = (int    *) malloc(4*dimension);
	for (i=0; i<dimension; i++) {
		temp_ptr[i] = (double *) malloc(sizeof(double)*(dimension+1));
		rows[0][i] = i;
	}
	return(temp_ptr);
}

int solve_pivot(a, rows, n)
/*
	Solves a system of n linear equations with n unknowns by
	Gaussian elimination with partial pivoting.  Call init_pivot()
	first to provide appropriate storage of the augmented
	matrix.  Solve_pivot provides an upper triangular augmented
	matrix (overwriting a) to be used by the routine backsub()
	for back substitution.  The lower triangular matrix is also
	stored (in the zeroed portion of a) for use by reaugment().
*/
double **a;
int *rows, n;
{
	double len, max, *hold;
	int i, j, k, ihold, ret_val;
	void pivot();

	/* find "length" of matrix a */
	len = 0.;
	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
			len += a[i][j] * a[i][j];
	len *= SMALL / (double) (n*n);
	ret_val = 0;

	for (i=0; i<n-1; i++) {

		/* check if all values in column are much smaller than length */
		max = fabs(a[i][i]);
		k = i;
		for (j=i+1; j<n; j++) {
			if (max < fabs(a[j][i])) {
				max = fabs(a[j][i]);
				k = j;
			}
		}
		if (max < len) {
			fprintf(stderr,"Solve_pivot: Matrix looks singular, column %d\n",i);
			ret_val = -1;
		}

		/* interchange rows */
		if (i != k) {
			hold = a[k];
			a[k] = a[i];
			a[i] = hold;
			ihold = rows[k];
			rows[k] = rows[i];
			rows[i] = ihold;
		}

		pivot(a, n, i);
	}
	return(ret_val);
}

void pivot(a, n, k)
/*
	Performs pivoting operation on lower right-hand submatrix.
	(Gaussian elimination... the row interchange is done before call.)
	Element a[k][k] is the pivot.  The lower triangular matrix
	(multipliers) is stored in the zeroed portion of the upper
	triangular matrix for future reference (the unit diagonal
	elements of the lower triangular matrix are implicit).
*/
double **a;
int n, k;
{
	double ratio;
	int i, j;

	for (i=k+1; i<n; i++) {
		ratio = -a[i][k] / a[k][k];
		a[i][k] = ratio;		/* store lower tri. value */
		for (j=k+1; j<=n; j++)
			a[i][j] += ratio * a[k][j];
	}
}

void backsub(a, n)
/*
	Performs backsubstitiution to get unknowns.
*/
double **a;
int n;
{
	double sum;
	int i, j;

	a[n-1][n] /= a[n-1][n-1];
	for (i=n-2; i>=0; i--) {
		sum = 0.0;
		for (j=i+1; j<n; j++)
			sum += a[i][j] * a[j][n];
		a[i][n] -= sum;
		a[i][n] /= a[i][i];
	}
}

void reaugment(a, b, rows, n)
/*
	Uses the lower triangular matrix saved from solve_pivot to
	reconstruct a new augmented matrix for any given vector b
	in A x = b.
*/
double **a, *b;
int *rows, n;
{
	int i, j;
	double sum;

	for (i=0; i<n; i++) {
		sum = 0.0;
		for (j=0; j<i; j++)
			sum += a[i][j] * a[j][n];
		sum += b[rows[i]];
		a[i][n] = sum;
	}
}
