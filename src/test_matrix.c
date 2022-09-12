#include "./s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
    if (rows > 0 && columns > 0) {
        result->columns = columns;
        result->rows = rows;
    } else {
        return ErrorMatrixpParam;
    }
    result->matrix = (double**) malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; ++i)
            result->matrix[i] = malloc(rows * columns * sizeof(double));
    fi(rows)
        fj(columns)
            result->matrix[i][j] = 0;
    return OK;
}

void s21_remove_matrix(matrix_t *A) {
    for (int i = 0; i < A->rows; ++i)
            free(A->matrix[i]);
    free(A->matrix);
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
    if (EQ(A, B)) {
        for (int i = 0; i < A->rows; ++i) {
            for (int j = 0; j < A->columns; ++j) {
                if (fabsl(A->matrix[i][j] - B->matrix[i][j]) > S21_EPS)
                    return FAILURE;
            }
        }
    } else {
        return FAILURE;
    }
    return SUCCESS;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
    int ans;
    if (EQ(A, B)) {
        ans = s21_create_matrix(A->rows, A->columns, result);
        fi(A->rows)
            fj(A->columns)
                result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
    } else {
        ans = ErrorMatrixOpetation;
    }
    return ans;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
    int ans;
    if (EQ(A, B)) {
        ans = s21_create_matrix(A->rows, A->columns, result);
        fi(A->rows)
            fj(A->columns)
                result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
    } else {
        ans = ErrorMatrixOpetation;
    }
    return ans;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
    int ans;
    ans = s21_create_matrix(A->rows, A->columns, result);
    fi(A->rows)
        fj(A->columns)
            result->matrix[i][j] = A->matrix[i][j] * number;
    return ans;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
    int ans;
    if (A->columns == B->rows) {
        ans = s21_create_matrix(A->rows, B->columns, result);
        fi(A->rows)
            fj(B->columns)
                fk(A->columns)
                    result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
    } else {
        ans =  ErrorMatrixOpetation;
    }
    return ans;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
    int ans;
    ans = s21_create_matrix(A->columns, A->rows, result);
    fi(A->rows)
        fj(A->columns)
            result->matrix[j][i] = A->matrix[i][j];
    return ans;
}
int just_determinant(matrix_t *A) {
    double summ = 0;
    if (A->rows == 1)
        summ = A->matrix[0][0];
    if (A->rows == 2) {
        summ = A->matrix[0][0] * A->matrix[1][1] -
                A->matrix[0][1] * A->matrix[1][0];
    }
    if (A->rows >= 3) {
        fi(A->rows) {
            fj(A->columns) {
                matrix_t minor;
                select_minor_matrix(A, &minor, i, j);
                summ += A->matrix[i][j] * pow(-1, i + j) *
                                        just_determinant(&minor);
                s21_remove_matrix(&minor);
            }
            return summ;
        }
    }
    return summ;
}

void print_matr(matrix_t *A) {
    fi(A->rows) {
        fj(A->columns)
            printf("%f ", A->matrix[i][j]);
        printf("\n");
    }
    printf("\n");
}

void select_minor_matrix(matrix_t *A, matrix_t *minor, int i, int j) {
    s21_create_matrix(A->rows- 1, A->columns - 1, minor);
    for (int m = 0, f = 0; m < A->rows; ++m, ++f) {
        if (m == i) ++m;
        if (m == A->rows) break;
        for (int n= 0, p = 0; n < A->columns; ++n, ++p) {
            if (n == j) ++n;
            if (n == A->columns) break;
            minor->matrix[f][p] = A->matrix[m][n];
        }
    }
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
    int ans;
    ans = s21_create_matrix(A->rows, A->columns, result);
    if (A->columns == A->rows) {
            if (A->rows == 1)
                result->matrix = A->matrix;
            if (A->rows == 2) {
                fi(A->rows) {
                    fj(A->columns) {
                        matrix_t minor;
                        select_minor_matrix(A, &minor, i, j);
                        result->matrix[i][j] = pow(-1, i + j) *
                                                minor.matrix[0][0];
                        s21_remove_matrix(&minor);
                    }
                }
            }
            if (A->rows >= 3) {
                fi(A->rows)
                    fj(A->columns) {
                        matrix_t minor;
                        select_minor_matrix(A, &minor, i, j);
                        result->matrix[i][j] = pow(-1, i + j) *
                                                just_determinant(&minor);
                        s21_remove_matrix(&minor);
                    }
            }
    } else {
        ans = ErrorMatrixOpetation;
    }
    return ans;
}

int s21_determinant(matrix_t *A, double *result) {
    int ans = OK;
    if (A->rows == A->columns)
        *result = just_determinant(A);
    else
        ans = ErrorMatrixOpetation;
    return ans;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
    int ans;
    double deter_A;
    ans = s21_determinant(A, &deter_A);
    if (A->rows == A->columns && ans == 0 && deter_A != 0) {
        matrix_t calc_matr, mult_matr;
        ans = s21_calc_complements(A, &calc_matr);
        ans = s21_mult_number(&calc_matr, 1 / deter_A, &mult_matr);
        ans = s21_transpose(&mult_matr, result);
        s21_remove_matrix(&calc_matr);
        s21_remove_matrix(&mult_matr);
    } else if (deter_A == 0) {
        ans = ErrorMatrixpParam;
    } else {
        ans = ErrorMatrixOpetation;
    }
    return ans;
}
