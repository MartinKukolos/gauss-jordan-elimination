#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX 10
#define EPSILON 1e-10

void printMatrix(double mat[][MAX], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            printf("%8.3f ", mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

int gaussJordan(double mat[][MAX], int n, double solution[]) {
    // Forward elimination to create upper triangular form
    for (int i = 0; i < n; i++) {
        // Find pivot
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(mat[k][i]) > fabs(mat[maxRow][i])) {
                maxRow = k;
            }
        }
        
        // Swap maximum row with current row
        for (int k = i; k <= n; k++) {
            double tmp = mat[maxRow][k];
            mat[maxRow][k] = mat[i][k];
            mat[i][k] = tmp;
        }
        
        // Check for singular matrix
        if (fabs(mat[i][i]) < EPSILON) {
            printf("Matrix is singular or nearly singular!\n");
            return 0;
        }
        
        // Make all rows below this one 0 in current column
        for (int k = i + 1; k < n; k++) {
            double factor = mat[k][i] / mat[i][i];
            for (int j = i; j <= n; j++) {
                mat[k][j] -= factor * mat[i][j];
            }
        }
    }
    
    // Backward elimination to create reduced row echelon form
    for (int i = n - 1; i >= 0; i--) {
        // Divide the pivot row by the pivot element
        double pivot = mat[i][i];
        for (int j = i; j <= n; j++) {
            mat[i][j] /= pivot;
        }
        
        // Make all rows above this one 0 in current column
        for (int k = i - 1; k >= 0; k--) {
            double factor = mat[k][i];
            for (int j = i; j <= n; j++) {
                mat[k][j] -= factor * mat[i][j];
            }
        }
    }
    
    // Extract solution
    for (int i = 0; i < n; i++) {
        solution[i] = mat[i][n];
    }
    
    return 1;
}

int main() {
    int n;
    double mat[MAX][MAX];
    double solution[MAX];
    
    printf("Enter the number of equations (max %d): ", MAX);
    scanf("%d", &n);
    
    if (n > MAX || n <= 0) {
        printf("Invalid number of equations!\n");
        return 1;
    }
    
    printf("\nEnter the augmented matrix [A|b]:\n");
    printf("(Enter coefficients row by row, including the constant term)\n");
    
    for (int i = 0; i < n; i++) {
        printf("Row %d: ", i + 1);
        for (int j = 0; j <= n; j++) {
            scanf("%lf", &mat[i][j]);
        }
    }
    
    printf("\nOriginal augmented matrix:\n");
    printMatrix(mat, n);
    
    if (gaussJordan(mat, n, solution)) {
        printf("Final reduced row echelon form:\n");
        printMatrix(mat, n);
        
        printf("Solution:\n");
        for (int i = 0; i < n; i++) {
            printf("x%d = %.6f\n", i + 1, solution[i]);
        }
    }
    
    return 0;
}

/* 
Example input for the system:
2x + y - z = 8
-3x - y + 2z = -11
-2x + y + 2z = -3

Enter: 3
Then enter:
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3

Expected output: x1 = 2, x2 = 3, x3 = -1
*/

/*
Enter the number of equations (max 10): 5

Enter the augmented matrix [A|b]:
(Enter coefficients row by row, including the constant term)
Row 1: 2 1 -1 0 1 3
Row 2: 1 -1 2 1 0 4
Row 3: 0 1 1 -1 2 5
Row 4: 1 0 -1 2 1 2
Row 5: -1 2 0 1 -1 1

Original augmented matrix:
    2.000    1.000   -1.000    0.000    1.000    3.000
    1.000   -1.000    2.000    1.000    0.000    4.000
    0.000    1.000    1.000   -1.000    2.000    5.000
    1.000    0.000   -1.000    2.000    1.000    2.000
    -1.000    2.000    0.000    1.000   -1.000    1.000

Final reduced row echelon form:
    1.000    0.000    0.000    0.000    0.000    1.074
    0.000    1.000    0.000    0.000    0.000    1.333
    0.000    0.000    1.000    0.000    0.000    1.778
    0.000    0.000    0.000    1.000    0.000    0.704
    0.000    0.000    0.000    0.000    1.000    1.296

Solution:
x1 = 1.074074
x2 = 1.333333
x3 = 1.777778
x4 = 0.703704
x5 = 1.296296
*/
