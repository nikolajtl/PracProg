# Meste kopieret fra Linear_Equations Makefile'en

CXX      = g++
CXXFLAGS = -O2 -Wall -std=c++23
LDFLAGS  = -llapacke -llapack -lblas

TARGET   = main.exe
SRC      = main.cc matrix.cc matrix.h

all: $(TARGET) run plot

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

run: $(TARGET)
	./$(TARGET)

plot: data.txt
	gnuplot -e "\
	set terminal pdf; set output 'plot.pdf'; \
	set xlabel 'x'; set ylabel 'y'; \
	set xrange [-1 : 11]; \
	set yrange [-1.2 : 2.2]; \
	plot 'data.txt' using 1:2 with points title 'Data', \
		'splinedata.txt' with lines title 'Quadratic Spline', \
		'integraldata.txt' with lines linestyle 1 dashtype 2 title 'Integral of spline' "

clean:
	rm -f $(TARGET) data.txt splinedata.txt integraldata.txt plot.pdf