#include <stdlib.h>
#include <matrix.h>

void free_matrix(Matrix* m) {
    /*
        Free memory space.
    */
    free(m->data);
    free(m);
}

Matrix* init_matrix(int row, int col, double* data) {
    /*
        Initialize a matrix structure.
    */
    Matrix* m;
    int i,j;
    m = (Matrix *)malloc(sizeof(Matrix));
    m->row = row;
    m->col = col;
    m->data = (double *)malloc(row * col * sizeof(double));
    if (data == NULL)
        return m;
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            SET(m, i, j, data[i + j * row]);
        }
    }
    return m;
}