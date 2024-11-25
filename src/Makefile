# Include paths and flags
toInclude = $(R_LIBRARY_DIR)/INLA/include/
PKG_CFLAGS = -I$(R_HOME)/include

# Object files
obj = cgeneric_debugger.o

# Compiler settings
CC = gcc
CXX = g++

# Compiler flags based on system
ifeq ($(shell uname -s), Linux)
    flags = -O2 -Wall -Wextra -fpic -g -fopenmp
    LINK_FLAGS = -L$(R_HOME)/lib -lR -ldl -lm
else ifeq ($(shell uname -s), Darwin)
    ARCH := $(shell uname -m)
    flags = -O2 -Wall -Wextra -fpic -g -arch $(ARCH) -fvisibility=default
    LINK_FLAGS = -L$(R_HOME)/lib -lR -ldl -lm -arch $(ARCH) -Wl,-export_dynamic
endif

# Default target
all: inla_cgeneric_debugger.so

# Check or download cgeneric.h
cgeneric.h:
	if [ ! -f cgeneric.h ]; then \
		echo "Downloading cgeneric.h..."; \
		curl -o cgeneric.h https://raw.githubusercontent.com/hrue/r-inla/devel/inlaprog/src/cgeneric.h; \
	fi

# Compile object files
%.o: %.c cgeneric.h
	$(CC) $(flags) $(PKG_CFLAGS) -I$(toInclude) -c $< -o $@

# Link shared library
inla_cgeneric_debugger.so: $(obj)
	$(CC) $(flags) -shared $(obj) -o inla_cgeneric_debugger.so $(LINK_FLAGS)

# Clean up
clean:
	rm -f *.o
	rm -f inla_cgeneric_debugger.so
	rm -f cgeneric.h

.PHONY: all clean