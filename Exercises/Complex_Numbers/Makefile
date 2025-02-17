CXX = clang++

CXXFLAGS = -std=c++20

all : main.exe run

main.exe : main.cc
	$(CXX) $(CXXFLAGS) -o main.exe main.cc

run:
	./main.exe

clean:
	$(RM) main.exe
