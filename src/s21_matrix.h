#ifndef SRC_S21_MATRIX_H_
#define SRC_S21_MATRIX_H_
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define OK 0
#define ErrorMatrixpParam 1
#define ErrorMatrixOpetation 2
#define SUCCESS 1
#define FAILURE 0
#define EQ(M, N) M->rows == N->rows && M->columns == N->columns ? 1 : 0
#define fi(line) for (int i = 0; i < line; ++i)
#define fj(line) for (int j = 0; j < line; ++j)
#define fk(line) for (int k = 0; k < line; ++k)
#define S21_EPS 1E-7


typedef struct matrix_struct {
    double** matrix;
    int rows;
    int columns;
} matrix_t;


int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);
int s21_eq_matrix(matrix_t *A, matrix_t *B);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);


int just_determinant(matrix_t *A);
void select_minor_matrix(matrix_t *A, matrix_t *minor, int i, int j);
void print_matr(matrix_t *A);

#endif  //  SRC_S21_MATRIX_H_
