CXX      = g++
CXXFLAGS = -O2 -Wall -std=c++23
LDFLAGS  = -llapacke -llapack -lblas

TARGET   = main.exe
SRC      = main.cc
OUTPUT	 = out.txt erfdata.txt convdata.txt
ERFPLOT	 = erfplot.pdf
CONVPLOT = convplot.pdf


all: $(TARGET) run ploterf plotconv

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

run: $(TARGET)
	./$(TARGET)
	python3 compare.py

ploterf: erfdata.txt
	gnuplot -e "\
	set terminal pdf; set output 'erfplot.pdf'; \
	set xlabel 'x'; set ylabel 'erf(x)'; \
	plot 'erfdata.txt' using 1:2 with lines linestyle 2 notitle "

plotconv: convdata.txt
	gnuplot -e "\
	set terminal pdf; set output 'convplot.pdf'; \
	set xlabel 'acc'; set ylabel '|myerf(1) - erf(1)|'; \
	set xrange [1e-10 : 1]; set yrange [1e-15 : 1e-3]; \
	set logscale x; set logscale y; \
	plot 'convdata.txt' using 1:2 with points notitle "

clean:
	rm -f $(TARGET) $(OUTPUT) $(ERFPLOT) $(CONVPLOT)
