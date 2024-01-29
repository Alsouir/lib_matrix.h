# lib_matrix.h

Реализация функций библиотеки matrix.h

В данном проекте реализована своя библиотека для обработки числовых матриц на языке программирования С. Матрицы являются одной из базовых структур данных в программировании, например, они применяются для представления табличных значений, для вычислительных задач и нейронных сетей.  


## Реализованные операции над матрицами 

Все операции (кроме сравнения матриц) возвращают результирующий код:  
- 0 - OK
- 1 - Ошибка, некорректная матрица   
- 2 - Ошибка вычисления (несовпадающие размеры матриц; матрица, для которой нельзя провести вычисления и т.д.)

1. Создание матриц (create_matrix)
2. Очистка матриц (remove_matrix)
3. Сравнение матриц (eq_matrix) (Сравнение должно происходит до 7 знака после запятой включительно.)
4. Сложение (sum_matrix) и вычитание матриц (sub_matrix)
5. Умножение матрицы на число (mult_number). Умножение двух матриц (mult_matrix)
6. Транспонирование матрицы (transpose)
7. Минор матрицы и матрица алгебраических дополнений (calc_complements)
8. Определитель матрицы (determinant)
9. Обратная матрица (inverse_matrix)


