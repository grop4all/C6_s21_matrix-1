CFLAGS=-c -Wall -Werror -Wextra -std=c11

ifeq ($(shell uname), Linux)
	LINUX_LIBS=-lsubunit -lrt -lpthread -lm
endif

all: clean s21_matrix.a

s21_matrix.a:
	rm -rf *.o
	gcc $(CFLAGS) s21_matrix.c
	ar r s21_matrix.a s21_matrix.o
	ranlib s21_matrix.a

s21_matrix.o: s21_matrix.c s21_matrix.h
	gcc $(CFLAGS) s21_matrix.c

test: s21_matrix.o
	gcc test_matrix.c s21_matrix.o -o check.out -lcheck $(LINUX_LIBS)
	./check.out

gcov_report: test_matrix.c s21_matrix.c
	clear
	gcc test_matrix.c s21_matrix.c -o coverage.o -lcheck $(LINUX_LIBS) -fprofile-arcs -ftest-coverage  
	./coverage.o
	lcov -c -d . -o coverage.info
	genhtml coverage.info -o coverage
	#open coverage/index.html

clean:
	rm -rf *.o *.out *.a
	rm -rf *.gcda *.gcno
	rm -rf coverage coverage.info 
	
leaks: test
	gcc -g test_matrix.c s21_matrix.o -o check.out -lcheck $(LINUX_LIBS)
	valgrind -q --leak-check=full --show-leak-kinds=all --track-origins=yes --show-reachable=no ./check.out

leaks2: test
	CK_FORK=no leaks --atExit -- ./check.out

cppcheck:
	cppcheck --enable=all --suppress=missingIncludeSystem ./*.c

lint:
	cp ./../materials/linters/.clang-format .
	clang-format -n *.c *.h
	clang-format -i *.c *.h
	clang-format -i *.c *.h
	rm .clang-format
