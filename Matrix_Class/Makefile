CXXFLAGS += -std=c++23 -Wall -Wextra -DLONG_DOUBLE
LDLIBS += -lstdc++ -lm
MAXERR = 1
ifeq ($(CXX),g++)
  CXXFLAGS += -fmax-errors=$(MAXERR)
else ifeq ($(findstring clang++,$(CXX)),clang++)
  CXXFLAGS += -ferror-limit=$(MAXERR)
else
  CXXFLAGS += -fmax-errors=$(MAXERR)
endif

Out.txt:main ; ./$< > $@
main:main.o matrix.o
main.o:matrix.h
matrix.o:matrix.h
clean:
	$(RM) [Oo]ut* main *.o
test:
	echo $(findstring clang++,$(CXX)) 
	echo $(CXX)
	echo $(CXXFLAGS)