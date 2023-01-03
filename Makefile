CC = gcc
CFLAGS = -Wall -Wextra -g -O3 -lm

SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)
EXE = ApproxSchrod

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(OBJ) $(CFLAGS) -o $@ 

$(OBJ): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean

clean: 
	rm -f *.o
