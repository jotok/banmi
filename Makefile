P=banmi
OBJECTS=banmi_util.o
CFLAGS = -g -Wall
LDLIBS = -lgsl -lgslcblas -lm
CC=c99

$(P): $(OBJECTS)

clean:
	rm -f $(P) $(OBJECTS)
	rm -f paper/*.{aux,log,pdf}
