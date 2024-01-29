#include <check.h>

#include "s21_matrix.h"

START_TEST(add_test) {
  matrix_t a, b, c, d;
  s21_create_matrix(3, 3, &a);
  s21_create_matrix(3, 3, &b);
  s21_create_matrix(3, 3, &d);
  a.matrix[0][0] = 1;
  a.matrix[0][1] = 2;
  a.matrix[0][2] = 3;
  a.matrix[1][0] = 1;
  a.matrix[1][1] = 2;
  a.matrix[1][2] = 3;
  a.matrix[2][0] = 1;
  a.matrix[2][1] = 2;
  a.matrix[2][2] = 3;
  b.matrix[0][0] = 3;
  b.matrix[0][1] = 2;
  b.matrix[0][2] = 1;
  b.matrix[1][0] = 3;
  b.matrix[1][1] = 2;
  b.matrix[1][2] = 1;
  b.matrix[2][0] = 3;
  b.matrix[2][1] = 2;
  b.matrix[2][2] = 1;
  d.matrix[0][0] = 4;
  d.matrix[0][1] = 4;
  d.matrix[0][2] = 4;
  d.matrix[1][0] = 4;
  d.matrix[1][1] = 4;
  d.matrix[1][2] = 4;
  d.matrix[2][0] = 4;
  d.matrix[2][1] = 4;
  d.matrix[2][2] = 4;
  int res = s21_sum_matrix(&a, &b, &c);
  if (!res)
    ck_assert_int_eq(s21_eq_matrix(&c, &d), 1);
  else
    ck_assert_int_eq(res, 0);
  s21_remove_matrix(&c);
  s21_remove_matrix(&a);
  s21_remove_matrix(&b);
  s21_remove_matrix(&d);
}
END_TEST

START_TEST(sub_test) {
  matrix_t a, b, c, d;
  s21_create_matrix(3, 3, &a);
  s21_create_matrix(3, 3, &b);
  s21_create_matrix(3, 3, &d);
  a.matrix[0][0] = 1;
  a.matrix[0][1] = 2;
  a.matrix[0][2] = 3;
  a.matrix[1][0] = 1;
  a.matrix[1][1] = 2;
  a.matrix[1][2] = 3;
  a.matrix[2][0] = 1;
  a.matrix[2][1] = 2;
  a.matrix[2][2] = 3;
  b.matrix[0][0] = 3;
  b.matrix[0][1] = 2;
  b.matrix[0][2] = 1;
  b.matrix[1][0] = 3;
  b.matrix[1][1] = 2;
  b.matrix[1][2] = 1;
  b.matrix[2][0] = 3;
  b.matrix[2][1] = 2;
  b.matrix[2][2] = 1;
  d.matrix[0][0] = -2;
  d.matrix[0][1] = 0;
  d.matrix[0][2] = 2;
  d.matrix[1][0] = -2;
  d.matrix[1][1] = 0;
  d.matrix[1][2] = 2;
  d.matrix[2][0] = -2;
  d.matrix[2][1] = 0;
  d.matrix[2][2] = 2;
  int res = s21_sub_matrix(&a, &b, &c);
  if (!res)
    ck_assert_int_eq(s21_eq_matrix(&c, &d), 1);
  else
    ck_assert_int_eq(res, 0);
  s21_remove_matrix(&c);
  s21_remove_matrix(&a);
  s21_remove_matrix(&b);
  s21_remove_matrix(&d);
}
END_TEST

START_TEST(mul_num) {
  matrix_t a, c, d;
  double b = -2;
  s21_create_matrix(3, 3, &a);
  s21_create_matrix(3, 3, &d);
  a.matrix[0][0] = 1;
  a.matrix[0][1] = 2;
  a.matrix[0][2] = 3;
  a.matrix[1][0] = 1;
  a.matrix[1][1] = 2;
  a.matrix[1][2] = 3;
  a.matrix[2][0] = 1;
  a.matrix[2][1] = 2;
  a.matrix[2][2] = 3;
  d.matrix[0][0] = -2;
  d.matrix[0][1] = -4;
  d.matrix[0][2] = -6;
  d.matrix[1][0] = -2;
  d.matrix[1][1] = -4;
  d.matrix[1][2] = -6;
  d.matrix[2][0] = -2;
  d.matrix[2][1] = -4;
  d.matrix[2][2] = -6;
  int res = s21_mult_number(&a, b, &c);
  if (!res)
    ck_assert_int_eq(s21_eq_matrix(&c, &d), 1);
  else
    ck_assert_int_eq(res, 0);
  s21_remove_matrix(&c);
  s21_remove_matrix(&a);
  s21_remove_matrix(&d);
}
END_TEST

START_TEST(mul_matrix) {
  matrix_t a, b, c, d;
  s21_create_matrix(2, 3, &a);
  s21_create_matrix(3, 4, &b);
  s21_create_matrix(2, 4, &d);
  a.matrix[0][0] = 1;
  a.matrix[0][1] = 2;
  a.matrix[0][2] = 3;
  a.matrix[1][0] = 1;
  a.matrix[1][1] = 2;
  a.matrix[1][2] = 3;
  b.matrix[0][0] = 3;
  b.matrix[0][1] = 2;
  b.matrix[0][2] = 1;
  b.matrix[0][3] = 1;
  b.matrix[1][0] = 3;
  b.matrix[1][1] = 2;
  b.matrix[1][2] = 1;
  b.matrix[1][3] = 1;
  b.matrix[2][0] = 3;
  b.matrix[2][1] = 2;
  b.matrix[2][2] = 1;
  b.matrix[2][3] = 1;
  d.matrix[0][0] = 18;
  d.matrix[0][1] = 12;
  d.matrix[0][2] = 6;
  d.matrix[0][3] = 6;
  d.matrix[1][0] = 18;
  d.matrix[1][1] = 12;
  d.matrix[1][2] = 6;
  d.matrix[1][3] = 6;
  int res = s21_mult_matrix(&a, &b, &c);
  if (!res)
    ck_assert_int_eq(s21_eq_matrix(&c, &d), 1);
  else
    ck_assert_int_eq(res, 0);
  s21_remove_matrix(&c);
  s21_remove_matrix(&a);
  s21_remove_matrix(&b);
  s21_remove_matrix(&d);
}
END_TEST

START_TEST(trans_test) {
  matrix_t a, c, d;
  s21_create_matrix(2, 3, &a);
  s21_create_matrix(3, 2, &d);
  a.matrix[0][0] = 1;
  a.matrix[0][1] = 2;
  a.matrix[0][2] = 3;
  a.matrix[1][0] = 1;
  a.matrix[1][1] = 2;
  a.matrix[1][2] = 3;
  d.matrix[0][0] = 1;
  d.matrix[0][1] = 1;
  d.matrix[1][0] = 2;
  d.matrix[1][1] = 2;
  d.matrix[2][0] = 3;
  d.matrix[2][1] = 3;
  int res = s21_transpose(&a, &c);
  if (!res)
    ck_assert_int_eq(s21_eq_matrix(&c, &d), 1);
  else
    ck_assert_int_eq(res, 0);
  s21_remove_matrix(&c);
  s21_remove_matrix(&a);
  s21_remove_matrix(&d);
}
END_TEST

START_TEST(inv_test) {
  matrix_t a, c, d;
  s21_create_matrix(3, 3, &a);
  s21_create_matrix(3, 3, &d);
  a.matrix[0][0] = 1;
  a.matrix[0][1] = 2;
  a.matrix[0][2] = 3;
  a.matrix[1][0] = 1;
  a.matrix[1][1] = 5;
  a.matrix[1][2] = 3;
  a.matrix[2][0] = 7;
  a.matrix[2][1] = 2;
  a.matrix[2][2] = 3;
  d.matrix[0][0] = -1.0 / 6;
  d.matrix[0][1] = 0;
  d.matrix[0][2] = 1.0 / 6;
  d.matrix[1][0] = -1.0 / 3;
  d.matrix[1][1] = 1.0 / 3;
  d.matrix[1][2] = 0;
  d.matrix[2][0] = 11.0 / 18;
  d.matrix[2][1] = -2.0 / 9;
  d.matrix[2][2] = -1.0 / 18;
  int res = s21_inverse_matrix(&a, &c);
  if (!res)
    ck_assert_int_eq(s21_eq_matrix(&c, &d), 1);
  else
    ck_assert_int_eq(res, 0);
  s21_remove_matrix(&c);
  s21_remove_matrix(&a);
  s21_remove_matrix(&d);
}
END_TEST

START_TEST(eq_test) {
  matrix_t a, c;
  s21_create_matrix(2, 2, &a);
  s21_create_matrix(2, 2, &c);
  a.matrix[0][0] = NAN;
  a.matrix[0][1] = 0;
  a.matrix[1][0] = 0;
  a.matrix[1][1] = 0;
  c.matrix[0][0] = NAN;
  c.matrix[0][1] = 1;
  c.matrix[1][0] = 0;
  c.matrix[1][1] = 0;
  int res = s21_eq_matrix(&a, &c);
  ck_assert_int_eq(res, 0);
  s21_remove_matrix(&c);
  s21_remove_matrix(&a);
}
END_TEST

START_TEST(det_inv_test_0) {
  matrix_t a, c;
  double b = 0.0;
  s21_create_matrix(5, 5, &a);
  s21_create_matrix(5, 5, &c);
  a.matrix[0][0] = 1;
  a.matrix[0][1] = 2;
  a.matrix[0][2] = 3;
  a.matrix[0][3] = 4;
  a.matrix[0][4] = 5;
  a.matrix[1][0] = 6;
  a.matrix[1][1] = 7;
  a.matrix[1][2] = 8;
  a.matrix[1][3] = 9;
  a.matrix[1][4] = -10;
  a.matrix[2][0] = 1;
  a.matrix[2][1] = 2;
  a.matrix[2][2] = 3;
  a.matrix[2][3] = 4;
  a.matrix[2][4] = 5;
  a.matrix[3][0] = 6;
  a.matrix[3][1] = 7;
  a.matrix[3][2] = 8;
  a.matrix[3][3] = 9;
  a.matrix[3][4] = -10;
  a.matrix[4][0] = 1;
  a.matrix[4][1] = 2;
  a.matrix[4][2] = 3;
  a.matrix[4][3] = 4;
  a.matrix[4][4] = 5;
  s21_determinant(&a, &b);
  ck_assert_double_eq_tol(b, 0, 1e-6);
  int res = s21_inverse_matrix(&a, &c);
  ck_assert_int_eq(res, 2);
  s21_remove_matrix(&a);
  s21_remove_matrix(&c);
}
END_TEST

START_TEST(det_test_1) {
  matrix_t a;
  double b = 0.0;
  s21_create_matrix(5, 5, &a);
  a.matrix[0][0] = 55;
  a.matrix[0][1] = 54;
  a.matrix[0][2] = 53;
  a.matrix[0][3] = -52;
  a.matrix[0][4] = 51;
  a.matrix[1][0] = 50;
  a.matrix[1][1] = 49;
  a.matrix[1][2] = 48;
  a.matrix[1][3] = 47;
  a.matrix[1][4] = -10;
  a.matrix[2][0] = 1;
  a.matrix[2][1] = 2;
  a.matrix[2][2] = 3;
  a.matrix[2][3] = 4;
  a.matrix[2][4] = -7;
  a.matrix[3][0] = 6;
  a.matrix[3][1] = 7;
  a.matrix[3][2] = 8;
  a.matrix[3][3] = 9;
  a.matrix[3][4] = 11;
  a.matrix[4][0] = 13;
  a.matrix[4][1] = 17;
  a.matrix[4][2] = 19;
  a.matrix[4][3] = 11;
  a.matrix[4][4] = 7;
  s21_determinant(&a, &b);
  ck_assert_double_eq_tol(b, 208624.0, 1e-6);
  s21_remove_matrix(&a);
}
END_TEST

START_TEST(subnormal) {
  double **arr;
  arr = calloc(2, sizeof(double *));
  for (int i = 0; i < 2; i++) {
    arr[i] = calloc(2, sizeof(double));
  }
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      arr[i][j] = 1.1 + i;
    }
  }
  double **farr;
  farr = calloc(1, sizeof(double *));
  farr[0] = calloc(1, sizeof(double));
  farr[0][0] = 3.0;
  matrix_t f = {farr, 1, 1}, a = {0, 2, 3}, b = {0, 0, 0}, c, d = {0, 2, 2},
           e = {arr, 0, 0};
  double i = 0;
  int res = s21_create_matrix(-1, 0, &a);
  ck_assert_int_eq(res, 1);
  res = s21_sum_matrix(&a, &b, &c);
  ck_assert_int_eq(res, 1);
  res = s21_sum_matrix(&a, &b, NULL);
  ck_assert_int_eq(res, 1);
  res = s21_sum_matrix(&a, &d, &c);
  ck_assert_int_eq(res, 2);
  res = s21_sub_matrix(&a, &b, &c);
  ck_assert_int_eq(res, 1);
  res = s21_sub_matrix(&a, &b, NULL);
  ck_assert_int_eq(res, 1);
  res = s21_sub_matrix(&a, &d, &c);
  ck_assert_int_eq(res, 2);
  res = s21_mult_number(&a, i, NULL);
  ck_assert_int_eq(res, 1);
  res = s21_mult_number(&b, i, &c);
  ck_assert_int_eq(res, 1);
  res = s21_mult_matrix(&a, &b, &c);
  ck_assert_int_eq(res, 1);
  res = s21_mult_matrix(NULL, &d, &c);
  ck_assert_int_eq(res, 1);
  res = s21_mult_matrix(&a, &d, &c);
  ck_assert_int_eq(res, 2);
  res = s21_transpose(&a, NULL);
  ck_assert_int_eq(res, 1);
  res = s21_calc_complements(NULL, &a);
  ck_assert_int_eq(res, 1);
  res = s21_calc_complements(&a, &c);
  ck_assert_int_eq(res, 2);
  res = s21_transpose(NULL, &c);
  ck_assert_int_eq(res, 1);
  res = s21_inverse_matrix(&a, NULL);
  ck_assert_int_eq(res, 1);
  res = s21_inverse_matrix(&e, &c);
  ck_assert_int_eq(res, 1);
  res = s21_inverse_matrix(&a, &c);
  ck_assert_int_eq(res, 2);
  res = s21_determinant(&a, NULL);
  ck_assert_int_eq(res, 1);
  res = s21_determinant(&a, &i);
  ck_assert_int_eq(res, 2);
  res = s21_determinant(&f, &i);
  ck_assert_int_eq(res, 0);
  res = s21_eq_matrix(&a, &b);
  ck_assert_int_eq(res, 0);
  for (int i = 0; i < 2; i++) {
    free(arr[i]);
  }
  free(farr[0]);
  free(arr);
  free(farr);
  arr = NULL;
}
END_TEST

START_TEST(inverse_matrix_1) {
  matrix_t A, B, C;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &C);
  A.matrix[0][0] = 2;
  A.matrix[0][1] = 5;
  A.matrix[0][2] = 7;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 4;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = -2;
  A.matrix[2][2] = -3;
  s21_inverse_matrix(&A, &B);
  C.matrix[0][0] = 1;
  C.matrix[0][1] = -1;
  C.matrix[0][2] = 1;
  C.matrix[1][0] = -38;
  C.matrix[1][1] = 41;
  C.matrix[1][2] = -34;
  C.matrix[2][0] = 27;
  C.matrix[2][1] = -29;
  C.matrix[2][2] = 24;
  ck_assert_int_eq(s21_eq_matrix(&B, &C), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
}
END_TEST

START_TEST(inverse_matrix_2) {
  matrix_t A, B, C;
  s21_create_matrix(4, 4, &A);
  s21_create_matrix(4, 4, &C);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[0][3] = 4;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = 2;
  A.matrix[1][3] = 1;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = 2;
  A.matrix[2][2] = 1;
  A.matrix[2][3] = 3;
  A.matrix[3][0] = 3;
  A.matrix[3][1] = 7;
  A.matrix[3][2] = 2;
  A.matrix[3][3] = 4;
  s21_inverse_matrix(&A, &B);
  C.matrix[0][0] = -0.125;
  C.matrix[0][1] = 0.19166666;
  C.matrix[0][2] = 0.325;
  C.matrix[0][3] = -0.16666666;
  C.matrix[1][0] = -0.125;
  C.matrix[1][1] = 0.05833333;
  C.matrix[1][2] = -0.075;
  C.matrix[1][3] = 0.16666666;
  C.matrix[2][0] = 0.125;
  C.matrix[2][1] = 0.675;
  C.matrix[2][2] = 0.275;
  C.matrix[2][3] = -0.5;
  C.matrix[3][0] = 0.25;
  C.matrix[3][1] = -0.58333333;
  C.matrix[3][2] = -0.25;
  C.matrix[3][3] = 0.33333333;
  ck_assert_int_eq(s21_eq_matrix(&B, &C), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
}
END_TEST

START_TEST(inverse_matrix_3) {
  matrix_t A, B;
  s21_create_matrix(3, 4, &A);
  ck_assert_int_eq(s21_inverse_matrix(&A, &B), 2);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(inverse_matrix_4) {
  // success with task reference values
  matrix_t A = {};
  matrix_t result = {};
  matrix_t eq_matrix = {};
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &eq_matrix);
  A.matrix[0][0] = 2, A.matrix[0][1] = 5, A.matrix[0][2] = 7;
  A.matrix[1][0] = 6, A.matrix[1][1] = 3, A.matrix[1][2] = 4;
  A.matrix[2][0] = 5, A.matrix[2][1] = -2, A.matrix[2][2] = -3;
  eq_matrix.matrix[0][0] = 1, eq_matrix.matrix[0][1] = -1,
  eq_matrix.matrix[0][2] = 1;
  eq_matrix.matrix[1][0] = -38, eq_matrix.matrix[1][1] = 41,
  eq_matrix.matrix[1][2] = -34;
  eq_matrix.matrix[2][0] = 27, eq_matrix.matrix[2][1] = -29,
  eq_matrix.matrix[2][2] = 24;
  ck_assert_int_eq(s21_inverse_matrix(&A, &result), 0);
  ck_assert_int_eq(s21_eq_matrix(&result, &eq_matrix), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&eq_matrix);
}
END_TEST

START_TEST(calc_complements) {
  matrix_t A = {};
  matrix_t result = {};
  matrix_t eq_matrix = {};
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &eq_matrix);
  A.matrix[0][0] = 1, A.matrix[0][1] = 2, A.matrix[0][2] = 3;
  A.matrix[1][0] = 0, A.matrix[1][1] = 4, A.matrix[1][2] = 2;
  A.matrix[2][0] = 5, A.matrix[2][1] = 2, A.matrix[2][2] = 1;
  ck_assert_int_eq(s21_calc_complements(&A, &result), 0);
  eq_matrix.matrix[0][0] = 0, eq_matrix.matrix[0][1] = 10,
  eq_matrix.matrix[0][2] = -20;
  eq_matrix.matrix[1][0] = 4, eq_matrix.matrix[1][1] = -14,
  eq_matrix.matrix[1][2] = 8;
  eq_matrix.matrix[2][0] = -8, eq_matrix.matrix[2][1] = -2,
  eq_matrix.matrix[2][2] = 4;
  ck_assert_int_eq(s21_eq_matrix(&result, &eq_matrix), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&eq_matrix);
}
END_TEST

int main(void) {
  Suite *s1 = suite_create("Core");
  TCase *tc1_1 = tcase_create("Core");
  SRunner *sr = srunner_create(s1);
  int nf;

  suite_add_tcase(s1, tc1_1);

  tcase_add_test(tc1_1, add_test);
  tcase_add_test(tc1_1, sub_test);
  tcase_add_test(tc1_1, mul_num);
  tcase_add_test(tc1_1, mul_matrix);
  tcase_add_test(tc1_1, trans_test);
  tcase_add_test(tc1_1, inv_test);
  tcase_add_test(tc1_1, eq_test);
  tcase_add_test(tc1_1, det_inv_test_0);
  tcase_add_test(tc1_1, det_test_1);
  tcase_add_test(tc1_1, subnormal);
  tcase_add_test(tc1_1, inverse_matrix_1);
  tcase_add_test(tc1_1, inverse_matrix_2);
  tcase_add_test(tc1_1, inverse_matrix_3);
  tcase_add_test(tc1_1, inverse_matrix_4);
  tcase_add_test(tc1_1, calc_complements);

  srunner_set_fork_status(sr, CK_NOFORK);
  srunner_run_all(sr, CK_ENV);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);

  return nf == 0 ? 0 : 1;
}