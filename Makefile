P=flas
OBJECTS = banmi_util.o banmi.o
CFLAGS = -g -Wall
LDLIBS = -lgsl -lgslcblas -lm
CC=c99

$(P): $(OBJECTS)

clean:
	rm -f $(P) $(OBJECTS)
