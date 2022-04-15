CC = g++
CFLAGS = -g -Wall -Wextra -Wpedantic -pedantic-errors -O3

NAME = smoothing
OBJS = io.o
DEPS = *.h
INCLUDE = -I../libraries/ -I/usr/include/eigen3/ -I/afs/umbc.edu/users/a/d/adamb/pub/ 
#INCLUDE = -I eigen-3.4.0

LIB = -L/usr/lib64 -lblas


# build target executable
$(NAME): $(OBJS)
	$(CC) -o $@ $^ $(LIB)

# build object files for each cpp file
# $@ is the target, $^ are the dependencies
%.o: %.cpp $(DEPS)
	$(CC) $(INCLUDE) -c -o $@ $< $(CFLAGS)

.PHONY: clean
clean:
	rm -f *.o $(NAME) out.obj

.PHONY: run
run:
	./$(NAME) bunny.obj out.obj 1 1 50

.PHONY: implicit
implicit:
	./smoothing -i 10e-5 bunny.obj out_implicit.obj 1 1 50
