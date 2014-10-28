CC=cc
CFLAGS=-g -Wall
LDFLAGS=-llapack -lblas -lm
BINARY=shallow_water

SOURCE_FILES=shallow_water.c\
             f_calc.c\
             init.c\
             io.c\
             jacobian.c

OBJECTS=$(SOURCE_FILES:.c=.o)

all: $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $(BINARY)

$(OBJECTS):
	$(CC) $(CFLAGS) -o $@ -c $(@:.o=.c)

delete:
	rm -r 16x16_Jac/
	mkdir 16x16_Jac/

clean:
	rm $(OBJECTS) $(BINARY)
