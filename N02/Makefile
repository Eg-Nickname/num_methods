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
N01_PROGRAM = main
HFILES = 

# Misc
VEC_GEN = vector_gen
BENCH = bench

# Solvers
SOLVERS_DIR = solvers/
SP_SOLVER = sparse_solver
D_SOLVER = dense_solver
SM_SOLVER = shermanmorrison_solver

SOLVERS = $(SP_SOLVER) $(D_SOLVER) $(SM_SOLVER)



# OBJECTS
N01_PROGRAM_OBJ = $(OBJ_DIR)/$(N01_PROGRAM).o
VEC_GEN_OBJ = $(OBJ_DIR)/$(VEC_GEN).o
BENCH_OBJ = $(OBJ_DIR)/$(BENCH).o
SOLVERS_OBJ = $(addprefix $(OBJ_DIR)/$(SOLVERS_DIR), $(addsuffix .o, $(SOLVERS)))

OBJ_FILES = $(N01_PROGRAM_OBJ) $(SOLVERS_OBJ) $(VEC_GEN_OBJ) $(BENCH_OBJ)
DEP_FILES = $(OBJ_FILES:.o=.d)

# Executables path
N01_PROGRAM_EXE = $(EXE_DIR)/$(N01_PROGRAM).x
EXES = $(N01_PROGRAM_EXE)

# Compilator flags
CFLAGS = -Wall -Wextra -std=c++23 -pedantic -I ./../eigen-5.0.1/ -I$(SRC_DIR)
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


$(N01_PROGRAM_EXE): $(OBJ_FILES) 
	$(LD) -o $@ $(LFLAGS) $^

.PHONY: run
run: $(N01_PROGRAM_EXE) Makefile target_dir
	$(N01_PROGRAM_EXE)

# .PHONY: compile
# compile: $(SORTING_PROGRAM_EXE) Makefile
# # 	Do nothing

.PHONY : clean
clean :
	$(RM) $(TARGET) *.o core; rm -rf $(TARGET_DIR)

.PHONY: target_dir
target_dir:
	mkdir -p $(OBJ_DIR)
	mkdir -p $(OBJ_DIR)/$(SOLVERS_DIR)
	mkdir -p $(EXE_DIR)