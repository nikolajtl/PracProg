CXX      = g++
CXXFLAGS = -O2 -Wall -std=c++23
LDFLAGS  = -llapacke -llapack -lblas

SRC = main.cc matrix.cc
OBJ = $(SRC:.cc=.o)
TARGET = main.exe
OUTFILE = out.txt
DATAFILE = data.txt

all: $(TARGET) run plot

# Rule to link the object files into the final executable
$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $(TARGET)

# Rule to compile the .cc files into .o files
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET) 2>&1 | tee $(OUTFILE)

plot: $(DATAFILE)
	gnuplot -e "set terminal pdf; set output 'plot.pdf'; \
	set logscale x; set logscale y; set xlabel 'N'; set ylabel 'Relative error'; \
	set title 'Singular integral'; \
	plot 'data.txt' using 1:2 with lines title 'Plain MC, actual error', \
		'data.txt' using 1:3 with lines title 'Plain MC, estimated error', \
		'data.txt' using 1:4 with lines title 'Quasi-random MC, estimated error'"

clean:
	rm -f $(OBJ) $(TARGET) $(OUTFILE) $(DATAFILE) plot.pdf

rebuild: clean all