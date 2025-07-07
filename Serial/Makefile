CC = gcc
CFLAGS = -Wall -O2

TARGET = serial_abelian_sandpile

all: $(TARGET)

$(TARGET): $(TARGET).c
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).c

run: $(TARGET)
	./$(TARGET) 1024 1024 0

clean:
	rm -f $(TARGET)