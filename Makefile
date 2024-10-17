all: test_c


WARNINGS = -Wall
DEBUG = -ggdb -fno-omit-frame-pointer
OPTIMIZE = -O2


#test_c: Makefile main.c
#	$(CC) -o $@ $(WARNINGS) $(DEBUG) $(OPTIMIZE) main.c

test_c: Makefile main.c
	$(CC) -o $@ $(WARNINGS) $(DEBUG) main.c

clean:
	rm -f test_c

# Builder will call this to install the application before running.
install:
	echo "Installing is not supported"
	./test_c

