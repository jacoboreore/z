# Z project
CXX      = g++
MPREAL_INC ?= mpreal

# --- compile/link flags ---
# high preciosion OS, bad performance
#CPPFLAGS += -DOS_USE_MPFR -DMPFR_USE_NO_MACRO -DMPFR_USE_INTMAX_T \
#            -I$(MPREAL_INC) -include stdint.h
# low  preciosion OS, good performance
CPPFLAGS += -DMPFR_USE_NO_MACRO -DMPFR_USE_INTMAX_T \
            -I$(MPREAL_INC) -include stdint.h
CXXFLAGS += -O3 -march=native -mtune=native -std=c++20 -MMD -MP -fopenmp $(CPPFLAGS)
LDFLAGS  += -pthread
LIBS     += -lmpc -lmpfr -lgmp -fopenmp -lstdc++fs

# Debug / tracing
DBGFLAGS = -O1 -g -fno-omit-frame-pointer -D_GLIBCXX_ASSERTIONS -D_FORTIFY_SOURCE=2 \
           -Wall -Wextra -Wconversion -Wshadow -Wpedantic -Wfloat-equal
SANFLAGS = -fsanitize=address,undefined

SRC_DIR  = .
OBJ_DIR  = build
BIN_DIR  = bin
TARGET   = $(BIN_DIR)/z

SOURCES =                     \
  custom_y.cpp                \
  z_theta.cpp                 \
  z_walk.cpp                  \
  z_calculator_rs_gabcke.cpp  \
  z_calculator_small.cpp      \
  z_calculator.cpp            \
  z_ball.cpp                  \
  refine_newton.cpp           \
  refine_newton_thread.cpp    \
  refine_stage1_request.cpp   \
  refine_stage2_decide.cpp    \
  z_distributed_prechunk.cpp  \
  settings.cpp                \
  cli.cpp                     \
  utils.cpp                   \
  main.cpp                      


OBJECTS = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(SOURCES))
DEPS    = $(OBJECTS:.o=.d)

# force rebuild when flags change 
HASH_CMD ?= md5sum          # on MSYS2 this exists; if not, use: sha1sum, e.g docker
FLAGS_ALL := $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(LIBS)
FLAGS_HASH := $(shell printf '%s' "$(FLAGS_ALL)" | $(HASH_CMD) | cut -d' ' -f1)
FLAGS_STAMP := $(OBJ_DIR)/.flags-$(FLAGS_HASH)

all: $(FLAGS_STAMP) $(TARGET)

$(FLAGS_STAMP):
	@mkdir -p $(OBJ_DIR)
	@rm -f $(OBJ_DIR)/.flags-*
	@touch $@

$(TARGET): $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBS)

$(OBJ_DIR)/%.o: %.cpp $(FLAGS_STAMP) | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@

debug: CXXFLAGS += $(DBGFLAGS) $(SANFLAGS)
debug: LDFLAGS  += $(SANFLAGS)
debug: all

clean:
	rm -rf $(OBJ_DIR) $(TARGET)

.PHONY: all clean debug
-include $(DEPS)
