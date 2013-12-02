P = flas_example
OBJECTS = banmi_util.o banmi.o
CFLAGS = -g -Wall -pg -fPIC
LDLIBS = -lgsl -lgslcblas -lm
CC=c99

$(P): $(OBJECTS)

banmi.o: banmi_util.o

mico.o: banmi_util.o

mico: $(OBJECTS)

thit: $(OBJECTS)
	$(CC) `pkg-config --cflags guile-2.0` -shared -o libthit.so \
		$(OBJECTS) -fPIC thit.c $(LDLIBS)

clean:
	rm -f $(P) $(OBJECTS) gmon.out libbanmi-guile.so
