CC=gcc
CFLAGS=-c -Wall  -Werror -Wextra -g -std=c11
CCTEST=gcc -fprofile-arcs -ftest-coverage

all: s21_matrix.a

s21_matrix.a: s21_matrix.o
	ar rc s21_matrix.a s21_matrix.o

s21_matrix.o: s21_matrix.c 
	$(CC) $(CFLAGS) s21_matrix.c

test.o: test.c 
	$(CC) $(CFLAGS) test.c

test: clean test.o s21_matrix.a
ifeq ($(shell uname),Linux)   
	$(CC) test.o s21_matrix.a -o unit_tests -lcheck -lm -lsubunit  
	./unit_tests 
else   
	$(CC) test.o s21_matrix.a -o unit_tests -lcheck  
	./unit_tests 
endif 


html_s21_matrix.a: html_s21_matrix.o 
	ar rc html_s21_matrix.a s21_matrix.o 

html_s21_matrix.o: s21_matrix.c 
	$(CCTEST) $(CFLAGS) s21_matrix.c

html_test: test.o html_s21_matrix.a
ifeq ($(shell uname),Linux)   
	$(CCTEST)  test.o html_s21_matrix.a -o unit_tests -lcheck -lm -lsubunit   
else   
	$(CCTEST)  test.o html_s21_matrix.a -o unit_tests -lcheck    
endif

gcov_report:clean html_test
	./unit_tests
	lcov -t "unit_tests" -o fizzbuzz.info -c -d .
	genhtml -o report fizzbuzz.info

clean:
	rm -rf *.o *.a *.gcno *.gcda *.info *.gcov unit_tests report