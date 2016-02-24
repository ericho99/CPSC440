#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

void inv_double_gs(double *a, int n, double *u, double *b);
void orthogonalize(double *u, double *b, int n, int col);
void normalize(double *u, double *b, int n, int col);
void updateSubtractVector(double *subtract, double *u, int n, double dp, int col);
double dotProductInMatrix(double *m, int n, int col1, int col2);
double dotProductTwoMatrices(double *m1, int row, double *m2, int col, int n);
double *transposeMatrix(double *m, int n);
void matrixMultiply(double *m1, double *m2, int n);
void assignIdentityMatrix(double *m, int n);
double *generateZeroVector(int n);
double *generateRandomMatrix(int n);
void copy_into(double *from, double *to, int n);
void setEntry(double *a, int n, int i, int j, double newVal);
double getEntry(double *a, int n, int i, int j);
void printMatrix(double *a, int n);
void printVector(double *a, int n);
int checkOrthogonal(double *m, int n, double precision);
int checkEqualMatrix(double *m1, double *m2, int n, double precision);

// just testing code for a sanity check
int main() {
  int n = 1200;
  // adjust this precision value based on how close numbers should be
  double precision = 0.0001;

  double *u = malloc(n * n * sizeof(double));
  double *b = malloc(n * n * sizeof(double));
  double *a = generateRandomMatrix(n);

  // call main function, wrapped by timer code
  clock_t start = clock(), diff;
  inv_double_gs(a, n, u, b);
  diff = clock() - start;

  int msec = diff * 1000 / CLOCKS_PER_SEC;
  printf("Time taken %d seconds %d milliseconds", msec / 1000, msec % 1000);


  if (checkOrthogonal(u, n, precision) == 1) {
    printf("success\n");
  }

  // multiplying a with its inverse and setting it to a
  matrixMultiply(a, b, n);
  double *checkIdentity = malloc(n * n * sizeof(double));
  assignIdentityMatrix(checkIdentity, n);
  // check that A A^-1 is equal to identity matrix
  if (checkEqualMatrix(checkIdentity, a, n, precision) == 1) {
    printf("another success yay\n");
  }

  // free variables
  free(checkIdentity);
  free(a);
  free(u);
  free(b);
}

void inv_double_gs(double *a, int n, double *u, double *b) {
  assignIdentityMatrix(b, n);
  copy_into(a, u, n);

  // orthogonalize each column of u according to gram schmidt
  int i;
  for (i = 0; i < n; ++i) {
    // orthogonalize twice for double gram schmidt
    orthogonalize(u, b, n, i);
    orthogonalize(u, b, n, i);
  }

  double *uTranspose = transposeMatrix(u, n); 
  matrixMultiply(b, uTranspose, n);
  free(uTranspose);
}

// orthogonalize column wrt previous columns, where u is the orthogonal matrix,
// b is the transformation matrix, n is the size of the matrix, and col
// is the column to orthogonalize
void orthogonalize(double *u, double *b, int n, int col) {
  // creates the vector which we subtract from the column of A
  // which is now contained in U
  int i;
  double *subtract = generateZeroVector(n);
  double *subtractT = generateZeroVector(n);
  for (i = 0; i < col; ++i) {
    double dp = dotProductInMatrix(u, n, i, col);
    updateSubtractVector(subtract, u, n, dp, i);
    updateSubtractVector(subtractT, b, n, dp, i);
  }

  // performs the subtraction
  for (i = 0; i < n; ++i) {
    setEntry(u, n, i, col, getEntry(u, n, i, col) - subtract[i]);
    setEntry(b, n, i, col, getEntry(b, n, i, col) - subtractT[i]);
  }

  free(subtract);
  free(subtractT);

  // normalize the resulting column
  normalize(u, b, n, col);
}

// normalizes specified column in u, and also divides by the same amount in b
void normalize(double *u, double *b, int n, int col) {
  int i;
  double size = sqrt(dotProductInMatrix(u, n, col, col));
  for (i = 0; i < n; ++i) {
    setEntry(u, n, i, col, (getEntry(u, n, i, col) / size));
    setEntry(b, n, i, col, (getEntry(b, n, i, col) / size));
  }
}

// updates the subtract vector with the dot product multiplied by the
// entries of the specified column of u
void updateSubtractVector(double *subtract, double *u, int n, double dp, int col) {
  int i;
  for (i = 0; i < n; ++i) {
    subtract[i] += (dp * getEntry(u, n, i, col));
  }
}

// performs the dot product of two columns in a matrix u
double dotProductInMatrix(double *m, int n, int col1, int col2) {
  int i;
  double sum = 0.0;
  for (i = 0; i < n; ++i) {
    sum += getEntry(m, n, i, col1) * getEntry(m, n, i, col2);
  }
  return sum;
}

// dot product of row and column in two matrices
double dotProductTwoMatrices(double *m1, int row, double *m2, int col, int n) {
  int i;
  double sum = 0.0;
  for (i = 0; i < n; ++i) {
    sum += getEntry(m1, n, row, i) * getEntry(m2, n, i, col);
  }
  return sum;
}

// transposes a square matrix, spits out a new matrix
double *transposeMatrix(double *m, int n) {
  int i, j;
  double *newMat = malloc(n * n * sizeof(double));
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j){
      setEntry(newMat, n, i, j, getEntry(m, n, j, i));
    }
  }
  return newMat;
}

// multiplies two matrixes, assigns value to first matrix
void matrixMultiply(double *m1, double *m2, int n) {
  int i, j;
  double *finalMat = malloc(n * n * sizeof(double));
  // multiply ith column of m2 with jth row of m1
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      setEntry(finalMat, n, j, i, dotProductTwoMatrices(m1, j, m2, i, n));
    }
  }

  copy_into(finalMat, m1, n);
  free(finalMat);
}

// assigns identity matrix to matrix m
void assignIdentityMatrix(double *m, int n) {
  int i, j;
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      if (i == j) {
        setEntry(m, n, i, j, 1.0);
      } else {
        setEntry(m, n, i, j, 0.0);
      }
    }
  }
}

// generates a vector of 0s of size n
double *generateZeroVector(int n) {
  int i;
  double *v = malloc(n * sizeof(double));
  for (i = 0; i < n; ++i) {
    v[i] = 0.0;
  }

  return v;
}

// generates matrix of random doubles of size n
double *generateRandomMatrix(int n) {
  int i;
  double *newMat = malloc(n * n * sizeof(double));
  for (i = 0; i < n * n; ++i) {
    // fill with random numbers between 0 and 10
    newMat[i] = rand() / ((double) RAND_MAX / 10);
  }
  return newMat;
}

// copies one matrix into another
void copy_into(double *from, double *to, int n){
  int i;
  for (i = 0; i < n * n; ++i){
    to[i] = from[i];
  }
}

// sets entry in matrix a
void setEntry(double *a, int n, int i, int j, double newVal){
  a[i * n + j] = newVal;
}

// gets entry in matrix a
double getEntry(double *a, int n, int i, int j){
  return a[i * n + j];
}

// prints a square matrix of size n, for debugging
void printMatrix(double *a, int n) {
  int i, j;
  printf("printing a matrix...\n");
  for (i = 0 ; i < n; ++i){
    for (j = 0; j < n; ++j){
      printf("%lf ", getEntry(a, n, i, j));
    }
    printf("\n");
  }
  printf("\n");
}

// prints a vector of size n, for debugging
void printVector(double *a, int n) {
  int i;
  printf("printing a vector...\n");
  for (i = 0; i < n; ++i){
    printf("%lf\n", a[i]);
  }
  printf("\n");
}

// checks if all a matrices columns are orthonormal, for debugging
int checkOrthogonal(double *m, int n, double precision) {
  int i, j;
  for (i = 0; i < n; ++i) {
    double dp = dotProductInMatrix(m, n, i, i);
    // checks if columns add to 1 with a certain precision
    if (dp > (1.0 + precision) || dp < (1.0 - precision)) {
      return 0;
    }

    // checks if rows are orthogonal to a certain precision
    for (j = 0; j < i; ++j) {
      double dp = dotProductInMatrix(m, n, i, j);
      if (dp > (0.0 + precision) || dp < (0.0 - precision)) {
        return 0;
      }
    }
  }

  return 1;
}

// checks if two matrices are equal to a certain precision, for debugging
int checkEqualMatrix(double *m1, double *m2, int n, double precision) {
  int i;
  for (i = 0; i < n * n; ++i) {
    double v1 = m1[i];
    double v2 = m2[i];
    if (v1 > (v2 + precision) || v1 < (v2 - precision)) {
      return 0;
    }
  }
  return 1;
}
