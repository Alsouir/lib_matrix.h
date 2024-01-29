#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  if (result == NULL || rows <= 0 || columns <= 0) {
    return 1;
  }
  int result_create = 0;
  result->rows = rows;
  result->columns = columns;
  result->matrix = calloc(rows, sizeof(double *));
  for (int i = 0; i < rows; i++) {
    result->matrix[i] = calloc(columns, sizeof(double));
  }
  return result_create;
}

void s21_remove_matrix(matrix_t *A) {
  if (A != NULL && A->rows > 0 && A->columns > 0) {
    for (int i = 0; i < A->rows; i++) {
      free(A->matrix[i]);
    }
    free(A->matrix);
    A->matrix = NULL;
    A->rows = 0;
    A->columns = 0;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  if (A == NULL || B == NULL || A->rows != B->rows ||
      A->columns != B->columns || A->rows == 0 || A->columns == 0) {
    return FAILURE;
  }
  int result = SUCCESS;
  for (int i = 0; i < A->rows && result == SUCCESS; i++) {
    for (int j = 0; j < A->columns && result == SUCCESS; j++) {
      if ((fabs(A->matrix[i][j] - B->matrix[i][j])) >= 1e-7) {
        result = FAILURE;
      }
    }
  }
  return result;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A == NULL || B == NULL || result == NULL || A->rows == 0 ||
      A->columns == 0 || B->columns == 0 || B->rows == 0) {
    return 1;
  } else if (A->rows != B->rows || A->columns != B->columns) {
    return 2;
  }
  s21_create_matrix(A->rows, A->columns, result);
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
    }
  }
  return 0;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A == NULL || B == NULL || result == NULL || A->rows == 0 ||
      A->columns == 0 || B->columns == 0 || B->rows == 0) {
    return 1;
  } else if (A->rows != B->rows || A->columns != B->columns) {
    return 2;
  }
  s21_create_matrix(A->rows, A->columns, result);
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
    }
  }
  return 0;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  if (A == NULL || result == NULL || A->rows == 0 || A->columns == 0) {
    return 1;
  }
  s21_create_matrix(A->rows, A->columns, result);
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] * number;
    }
  }
  return 0;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A == NULL || B == NULL || result == NULL || A->rows == 0 ||
      A->columns == 0 || B->columns == 0 || B->rows == 0) {
    return 1;
  } else if (A->columns != B->rows) {
    return 2;
  }
  s21_create_matrix(A->rows, B->columns, result);
  double sum_result_ij = 0;
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < B->columns; j++) {
      sum_result_ij = 0;
      for (int k = 0; k < A->columns; k++) {
        sum_result_ij = sum_result_ij + A->matrix[i][k] * B->matrix[k][j];
      }
      result->matrix[i][j] = sum_result_ij;
    }
  }
  return 0;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  if (A == NULL || result == NULL || A->columns < 1 || A->rows < 1) {
    return 1;
  }
  s21_create_matrix(A->columns, A->rows, result);
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[j][i] = A->matrix[i][j];
    }
  }
  return 0;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  if (A == NULL || result == NULL) {
    return 1;
  } else if (A->columns != A->rows) {
    return 2;
  }
  s21_create_matrix(A->columns, A->rows, result);

  int rows_minor = 0;
  int columns_minor = 0;
  double minor_determinant = 0;
  matrix_t minor;
  int sing = 1;
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      if ((i + j) % 2 == 0) {
        sing = 1;
      } else {
        sing = -1;
      }
      s21_create_matrix(A->columns - 1, A->rows - 1, &minor);
      rows_minor = 0;
      for (int k = 0; k < A->columns; k++) {
        if (k == i) continue;
        columns_minor = 0;
        for (int t = 0; t < A->columns; t++) {
          if (t == j) continue;
          minor.matrix[rows_minor][columns_minor] = A->matrix[k][t];
          columns_minor++;
        }
        rows_minor++;
      }
      minor_determinant = 0;
      s21_determinant(&minor, &minor_determinant);
      result->matrix[i][j] =
          sing * minor_determinant;  // pow(-1, (i + j)) * minor_determinant;
      s21_remove_matrix(&minor);
    }
  }

  return 0;
}

int s21_determinant(matrix_t *A, double *result) {
  if (A == NULL || result == NULL || A->columns < 1 || A->rows < 1) {
    return 1;
  } else if (A->columns != A->rows) {
    return 2;
  }
  if (A->columns == 1) {
    *result = A->matrix[0][0];
  } else if (A->columns == 2) {
    *result =
        A->matrix[1][1] * A->matrix[0][0] - A->matrix[1][0] * A->matrix[0][1];
  } else {
    int sing = 1;
    matrix_t auxiliary_matrix;
    for (int i = 0; i < A->rows; i++) {
      s21_create_matrix(A->columns - 1, A->rows - 1, &auxiliary_matrix);
      int t = 0;
      for (int j = 1; j < A->rows; j++) {
        int k = 0;
        for (int l = 0; l < A->rows; l++) {
          if (l != i) {
            auxiliary_matrix.matrix[t][k] = A->matrix[j][l];
            k++;
          }
        }
        t++;
      }
      double result_auxiliary_matrix = 0;
      s21_determinant(&auxiliary_matrix, &result_auxiliary_matrix);
      *result += A->matrix[0][i] * result_auxiliary_matrix * sing;
      sing = -sing;
      s21_remove_matrix(&auxiliary_matrix);
    }
  }
  return 0;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  if (A == NULL || result == NULL || A->rows < 1 || A->columns < 1) {
    return 1;
  } else if (A->columns != A->rows) {
    return 2;
  }
  double determinant = 0;
  int res = 0;
  s21_determinant(A, &determinant);

  if (determinant != 0) {
    matrix_t calc_complements;
    s21_calc_complements(A, &calc_complements);
    matrix_t transpose;
    s21_transpose(&calc_complements, &transpose);
    determinant = 1 / determinant;
    s21_mult_number(&transpose, determinant, result);
    s21_remove_matrix(&calc_complements);
    s21_remove_matrix(&transpose);
  } else {
    res = 2;
  }
  return res;
}
