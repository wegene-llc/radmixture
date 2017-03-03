typedef struct {
    int row;
    int col;
    double* data;
}
Matrix;

// Set the value of an entry.
#define SET(m, r, c, v) ((m)->data[(r) * (m)->col + (c)] = (v))
// Get the value of an entry
#define GET(m, r, c) ((m)->data[(r) * (m)->col + (c)])

void free_matrix(Matrix* m);
Matrix* init_matrix(int row, int col, double* data);
