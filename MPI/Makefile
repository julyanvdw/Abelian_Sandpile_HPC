# Compiler and flags
CC       = mpicc
CFLAGS   = -Wall -O2

# Number of processors for MPI run
NP       = 8

# Target executable
TARGET   = MPI_abelian_sandpile
SRC      = MPI_abelian_sandpile.c
OBJ      = $(SRC:.c=.o)

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

run: $(TARGET)
	mpirun -np $(NP) ./$(TARGET) 1024 1024 0

clean:
	rm -f $(TARGET) $(OBJ)