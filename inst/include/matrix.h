typedef struct {
    int row;
    int col;
    double* data;
}
Matrix;

#define SET(m, r, c, v) ((m)->data[(r) * (m)->col + (c)] = (v))
#define GET(m, r, c) ((m)->data[(r) * (m)->col + (c)])

void free_matrix(Matrix* m);
Matrix* init_matrix(int row, int col, double* data);
