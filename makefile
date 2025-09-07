# Makefile for z_calculator_thread project

CXX = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra
LIBS = -lmpc -lmpfr -lgmp
TARGET = bin/z_calculator_thread

SOURCES =                     \
  custom_y.cpp                \
  z_walk.cpp                  \
  z_calculator_thread.cpp     \
  z_calculator_small.cpp      \
  z_calculator.cpp            \
  main.cpp                    

HEADERS = config.h 

all: $(TARGET)

$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES) $(LIBS)

clean:
	rm -f $(TARGET)
