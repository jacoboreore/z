# Makefile for z_calculator_thread project

CXX = g++
CXXFLAGS = -O2 -std=c++17 
# including warnings, to use this line comment the line above 
#CXXFLAGS = -O2 -std=c++17 #-Wall -Wextra

LIBS = -lmpc -lmpfr -lgmp
TARGET = bin/z_calculator_thread

SOURCES =                     \
  custom_y.cpp                \
  z_theta.cpp                 \
  z_walk.cpp                  \
  z_calculator_thread.cpp     \
  z_calculator_small.cpp      \
  z_calculator.cpp            \
  refine_newton.cpp           \
  main.cpp                    

HEADERS = config.h            \
          logs.h

all: $(TARGET)

$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES) $(LIBS)

clean:
	rm -f $(TARGET)
