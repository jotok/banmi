P = flas_example
OBJECTS = banmi_util.o banmi.o
CFLAGS = -g -Wall -pg -fPIC
LDLIBS = -lgsl -lgslcblas -lm
CC=c99

$(P): $(OBJECTS)

thit:
	$(CC) `pkg-config --cflags guile-1.8` -shared -o libthit.so \
		$(OBJECTS) -fPIC thit.c $(LDLIBS)

clean: $(OBJECTS)
	rm -f $(P) $(OBJECTS) gmon.out libbanmi-guile.so
