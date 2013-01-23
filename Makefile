P=flas
OBJECTS = banmi_util.o banmi.o
CFLAGS = -g -Wall -pg -fPIC
LDLIBS = -lgsl -lgslcblas -lm
CC=c99

$(P): $(OBJECTS)

guile:
	$(CC) `pkg-config --cflags guile-1.8` -shared -o libbanmi-guile.so \
		$(OBJECTS) -fPIC banmi_guile.c $(LDLIBS)

clean: $(OBJECTS)
	rm -f $(P) $(OBJECTS) gmon.out libbanmi-guile.so
