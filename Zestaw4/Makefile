# Makefile
# Compiler
LIBS =

DIR = $(notdir $(CURDIR))
# Compilation target directory
TARGET_DIR = ./target
# Object directory in target
OBJ_DIR = $(TARGET_DIR)/obj
# Executables directory in target
EXE_DIR = $(TARGET_DIR)/exe
# Source code directory
SRC_DIR = ./src

# Names for object and exe files for programs
Z4_PROGRAM = main
HFILES = 

# Misc
VEC_GEN = vector_gen
BENCH = bench

# Solvers
SOLVERS_DIR = solvers/
SP_SOLVER = sparse_solver
D_SOLVER = dense_solver
TH_SOLVER = thomas_solver

SOLVERS = $(SP_SOLVER) $(D_SOLVER) $(TH_SOLVER)



# OBJECTS
Z4_PROGRAM_OBJ = $(OBJ_DIR)/$(Z4_PROGRAM).o
# VEC_GEN_OBJ = $(OBJ_DIR)/$(VEC_GEN).o
# BENCH_OBJ = $(OBJ_DIR)/$(BENCH).o
# SOLVERS_OBJ = $(addprefix $(OBJ_DIR)/$(SOLVERS_DIR), $(addsuffix .o, $(SOLVERS)))

OBJ_FILES = $(Z4_PROGRAM_OBJ)
DEP_FILES = $(OBJ_FILES:.o=.d)

# Executables path
Z4_PROGRAM_EXE = $(EXE_DIR)/$(Z4_PROGRAM).x
EXES = $(Z4_PROGRAM_EXE)

# Compilator flags
CFLAGS = -Wall -Wextra -std=c++23 -pedantic -I ./../eigen-5.0.1/ -I$(SRC_DIR) -O3 -march=native
# Linker flags
LFLAGS =
# Compiler and linker
CO = clang++
LD = $(CO)

# Compilation rules
# $(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HFILES) target_dir
# 	$(CO) $(CFLAGS) -c $< -o $@
# $(OBJ_DIR)/%.o: $(SRC_DIR)/%.c target_dir
# 	$(CO) $(CFLAGS) -c $< -o $@


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	mkdir -p $(dir $@)
	$(CO) $(CFLAGS) -c $< -o $@


$(Z4_PROGRAM_EXE): $(OBJ_FILES)
	mkdir -p $(EXE_DIR)
	$(LD) -o $@ $(LFLAGS) $^

.PHONY: run
run: $(Z4_PROGRAM_EXE) Makefile target_dir
	$(Z4_PROGRAM_EXE)

.PHONY: compile
compile: $(Z4_PROGRAM_EXE) Makefile target_dir
# 	Do nothing

.PHONY : clean
clean :
	$(RM) $(TARGET) *.o core; rm -rf $(TARGET_DIR)

.PHONY: target_dir
target_dir:
	mkdir -p $(EXE_DIR)
	mkdir -p $(OBJ_DIR)
	mkdir -p $(OBJ_DIR)/$(SOLVERS_DIR)

$(EXE_DIR):
	mkdir -p $(EXE_DIR)

.PHONY:  tar
tar: clean
	(cd ../; tar -cvzf $(DIR).tar.gz  $(DIR) )