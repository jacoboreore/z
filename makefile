# Makefile for Z project

CXX      = g++
CXXFLAGS = -O2 -std=c++20 -MMD -MP
LDFLAGS  = -pthread

# add pure headers here 
LIBS     = -lmpc -lmpfr -lgmp -fopenmp

SRC_DIR  = .
OBJ_DIR  = build
BIN_DIR  = bin
TARGET   = $(BIN_DIR)/z_calculator_thread

SOURCES = \
  custom_y.cpp                \
  z_theta.cpp                 \
  z_walk.cpp                  \
  z_calculator_rs_gabcke.cpp  \
  z_calculator_thread.cpp     \
  z_calculator_small.cpp      \
  z_calculator.cpp            \
  z_ball.cpp                  \
  refine_newton.cpp           \
  settings.cpp                \
  cli.cpp                     \
  utils.cpp                   \
  main.cpp

OBJECTS = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(SOURCES))
DEPS    = $(OBJECTS:.o=.d)

all: $(TARGET)

$(TARGET): $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)

$(OBJ_DIR)/%.o: %.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@

clean:
	rm -rf $(OBJ_DIR) $(TARGET)

.PHONY: all clean
-include $(DEPS)
