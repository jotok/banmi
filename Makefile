P = flas_example
OBJECTS = banmi_util.o banmi.o mico.o
CFLAGS = -g -Wall -pg -fPIC
LDLIBS = -lgsl -lgslcblas -lm
CC=c99

$(P): $(OBJECTS)

banmi.o: banmi_util.o

mico.o: banmi_util.o

thit:
	$(CC) `pkg-config --cflags guile-1.8` -shared -o libthit.so \
		$(OBJECTS) -fPIC thit.c $(LDLIBS)

clean: $(OBJECTS)
	rm -f $(P) $(OBJECTS) gmon.out libbanmi-guile.so
