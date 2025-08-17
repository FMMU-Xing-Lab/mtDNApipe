# scripts/overlap_corrector/Makefile
# #
# # A simple Makefile to compile the overlap_corrector tool.
# # To compile, navigate to this directory and run 'make'.
# #

CC = g++
CFLAGS = -O2 -std=c++11
LDFLAGS =
SOURCES = overlap_corrector.cpp
EXECUTABLE = overlap_corrector

all: $(EXECUTABLE)

$(EXECUTABLE): $(SOURCES)
	$(CC) $(CFLAGS) $(SOURCES) -o $(EXECUTABLE) $(LDFLAGS)

clean:
	rm -f $(EXECUTABLE)
