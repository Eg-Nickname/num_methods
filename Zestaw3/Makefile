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
JAC_PROGRAM = jacobi
HFILES = 

# Executables objects path
JAC_PROGRAM_OBJ= $(OBJ_DIR)/$(JAC_PROGRAM).o
# Executables path
JAC_PROGRAM_EXE = $(EXE_DIR)/$(JAC_PROGRAM).x

# Compilator flags
CFLAGS = -Wall -Wextra -std=c++23 -pedantic -I ./../eigen-5.0.1/
# Linker flags
LFLAGS = 
# Compiler and linker
CO = clang++
LD = $(CO)

# Compilation rules
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HFILES) target_dir
	$(CO) $(CFLAGS) -c $< -o $@
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c target_dir
	$(CO) $(CFLAGS) -c $< -o $@

$(JAC_PROGRAM_EXE): $(JAC_PROGRAM_OBJ) 
	$(LD) -o $@ $(LFLAGS) $^

.PHONY: run
run: $(JAC_PROGRAM_EXE) Makefile target_dir
	$(JAC_PROGRAM_EXE)

# .PHONY: compile
# compile: $(SORTING_PROGRAM_EXE) Makefile
# # 	Do nothing

.PHONY : clean
clean :
	$(RM) $(TARGET) *.o core; rm -rf $(TARGET_DIR)

.PHONY: target_dir
target_dir:
	mkdir -p $(OBJ_DIR)
	mkdir -p $(EXE_DIR)