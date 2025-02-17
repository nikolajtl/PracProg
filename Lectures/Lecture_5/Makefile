# Compiler
CXX = clang++

# Compiler Flags
CXXFLAGS = -Wall -Wextra -std=c++17

# Executable Name
TARGET = main.exe

# Source Files
SRC = main.cc

# Build Target
all: $(TARGET)

# Compile the program
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Clean build files
clean:
	rm -f $(TARGET)
