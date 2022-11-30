SHELL = /bin/sh


CC = mpicxx

#DEG = -DDfinding_eigenvalues
#DEG = -DDminus_shift -DDfinding_eigenvalues
#DEG = -DDEBUG_solve -DDEBUG_reflection -DDEBUG_get_x
#FLAGS = -Wall -Wextra -pedantic -Werror -g
FLAGS = -Wall -Wextra -pedantic -Werror -O3 


all: clean prog
prog: main.o prole.o lu_matrix.o
	$(CC) $(FLAGS) $(DEG) main.o prole.o lu_matrix.o
main.o: main.cpp errors.h
	$(CC) -c $(FLAGS) $(DEG) main.cpp
prole.o: prole.cpp prole.h errors.h
	$(CC) -c $(FLAGS) $(DEG) prole.cpp
lu_matrix.o: lu_matrix.cpp lu_matrix.h errors.h
	$(CC) -c $(FLAGS) $(DEG) lu_matrix.cpp
clean: 
	rm  -f *.o *.out 
clean_txt:
	rm  -f *.txt
