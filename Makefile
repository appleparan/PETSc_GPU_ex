.DEFAULT_GOAL := all
#### PROJECT SETTINGS ####
# The name of the executable to be created
TARGET := poisson
MPIDIR = /usr/

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

# Compiler used
CC := $(MPIDIR)/bin/mpicc
# source file directory
SRC := src
# object file directory
OBJ := obj

SOURCES := $(wildcard $(SRC)/*.c)
OBJECTS := $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SOURCES))
CUSPDIR = /opt/cusp/
# Add additional include paths
CFLAGS = -O2 -I$(SRC_PATH) ${PETSC_CCPPFLAGS} -I${CUSPDIR}
LFLAGS = -L${PETSC_CSH_LIB_PATH} 
LDLIBS = ${PETSC_KSP_LIB}

#### END PROJECT SETTINGS ####
all: $(OBJECTS) 
	$(CC) $^ -o $(TARGET) $(LFLAGS) $(LDLIBS)

$(OBJ)/%.o: $(SRC)/%.c
	mkdir -p $(OBJ)
	$(CC) -I$(SRC) $(CFLAGS) -c $< -o $@

# Removes all build files
.PHONY: cleanall
cleanall:
	echo "Deleting all obj files"
	rm -r $(OBJ)
	rm $(TARGET)

