CXX = clang++

CXXFLAGS = -std=c++20

all : io1.exe io2.exe io3.exe

io1.exe : io1.cc
	$(CXX) $(CXXFLAGS) -o io1.exe io1.cc

io2.exe : io2.cc
	$(CXX) $(CXXFLAGS) -o io2.exe io2.cc

io3.exe : io3.cc
	$(CXX) $(CXXFLAGS) -o io3.exe io3.cc

clean :
	$(RM) io1.exe io2.exe io3.exe out.txt outputfile

test:
	echo $(CXX)
