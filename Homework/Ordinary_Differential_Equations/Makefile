CXX      = g++
CXXFLAGS = -O2 -Wall -std=c++23
LDFLAGS  = -llapacke -llapack -lblas

TARGET   = main.exe
SRC      = main.cc matrix.cc
OBJ      = $(SRC:.cc=.o)  # Automatically replaces .cc with .o

all: $(TARGET) run plot

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

plot: funcvals.txt
	gnuplot -e "\
	set terminal pdf; set output 'plot.pdf'; \
	set xlabel 'x'; set ylabel 'y'; \
	set xrange [-8 : 8]; set yrange [-5 : 15]; set size ratio -1; \
	plot 'funcvals.txt' using 1:2 with lines linestyle 2 title 'Relativistic orbit' "

clean:
	rm -f $(TARGET) $(OBJ) funcvals.txt plot.pdf