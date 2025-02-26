# CXXFLAGS += -Wall -Wextra -std=c++23 -ferror-limit=1
# LDLIBS += -lstdc++ -lm
# CODE = $(filter %.cs,$^)
# DLLS = $(addprefix -reference:,$(filter %.dll,$^))
# MKLIB = mcs -out: $@ -target:library $(DLLS) $(CODE)
# MKEXE = mcs -target:exe $(DLLS) $(CODE)

# Out.txt: main.exe main
# 	mono $< > $@
# 	./main > $@
# main.exe: main.cs ; $(MKEXE)

CXXFLAGS += -Wall -Wextra -std=c++23 -fmax-errors=1
LDLIBS += -lm

# Filter source files and shared libraries
CODE = $(filter %.cc,$^)  # Extracts .cc source files
LIBS = $(addprefix -l,$(basename $(filter %.so,$^)))  # Extracts .so libraries

# Compilation commands
MKLIB = g++ -shared -o $@ $(CODE) $(LIBS) $(CXXFLAGS) $(LDLIBS)  # Build shared library
MKEXE = g++ -o $@ $(CODE) $(LIBS) $(CXXFLAGS) $(LDLIBS)  # Build executable

# Run executable and save output
Out.txt: main.exe main
	./main.exe > $@
	./main > $@

# Build main executable
main.exe: main.cc
	$(MKEXE)

clean:
	rm -f *.exe *.so Out.txt