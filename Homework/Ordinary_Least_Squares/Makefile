CXX      = g++
CXXFLAGS = -O2 -Wall -std=c++23
LDFLAGS  = -llapacke -llapack -lblas

TARGET   = main.exe
SRC      = main.cc matrix.cc matrix.h
FITFILE = fit.txt
DATAFILE = data.txt
PLOTFILE = fitplot.pdf

all: $(TARGET) run plot

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

run: $(TARGET)
	./$(TARGET)

# Data:
TS = 1 2 3 4 6 9 10 13 15
NS = 117 100 88 72 53 29.5 25.2 15.2 11.1
DNS = 6 5 4 4 4 3 3 2 2

plot: $(FITFILE)
	gnuplot -e "\
	set terminal pdf; set output 'fitplot.pdf'; \
	set xlabel 'Time [days]'; set ylabel 'N [arb. units]'; \
	set xrange [0 : 16]; \
	set yrange [-10 : 160]; \
	plot 'data.txt' using 1:2:3 with errorbars pointtype 7 pointsize 0.3 title 'Data', \
		'fit.txt' with lines title 'Fit', \
		'fit_ap_lp.txt' with lines title 'Fit, ++', \
		'fit_am_lm.txt' with lines title 'Fit, --', \
		'fit_ap_lm.txt' with lines title 'Fit, +-', \
		'fit_am_lp.txt' with lines title 'Fit, -+' "

clean:
	rm -f $(TARGET) $(FITFILE) $(PLOTFILE) fit_ap_lp.txt fit_am_lm.txt fit_ap_lm.txt fit_am_lp.txt