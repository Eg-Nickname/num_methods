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
CFLAGS = -Wall -Wextra -std=c++23 -pedantic -I ./../eigen-5.0.1/ -I$(SRC_DIR) -O3 -march=native -DNDEBUG
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

.PHONY: compile
compile: $(N01_PROGRAM_EXE) Makefile target_dir
# 	Do nothing

.PHONY : clean
clean :
	$(RM) $(TARGET) *.o core; rm -rf $(TARGET_DIR)

.PHONY: target_dir
target_dir:
	mkdir -p $(OBJ_DIR)
	mkdir -p $(OBJ_DIR)/$(SOLVERS_DIR)
	mkdir -p $(EXE_DIR)

.PHONY: bench_data_dir
bench_data_dir:
	mkdir -p bench_data

.PHONY: run_sparse_lu_10k
run_sparse_lu_10k: $(N01_PROGRAM_EXE) Makefile target_dir bench_data_dir
	$(N01_PROGRAM_EXE) -b ./bench_data/optimised_sparse_lu_10k.csv SPARSE_LU 30 100 10000 1 40
# 	$(N01_PROGRAM_EXE) -b ./bench_data/optimised_sparse_lu_50k.csv SPARSE_LU 100 200 50000 1 20

.PHONY: run_sparse_qr_10k
run_sparse_qr_10k: $(N01_PROGRAM_EXE) Makefile target_dir bench_data_dir
	$(N01_PROGRAM_EXE) -b ./bench_data/optimised_sparse_qr_10k.csv SPARSE_QR 30 100 10000 1 20
# 	$(N01_PROGRAM_EXE) -b ./bench_data/optimised_sparse_qr_10k.csv SPARSE_QR 20 50 50000 1 20

.PHONY: run_dense_full_lu_10k
run_dense_full_lu_10k: $(N01_PROGRAM_EXE) Makefile target_dir bench_data_dir
	$(N01_PROGRAM_EXE) -b ./bench_data/optimised_fullpiv_dense_lu_10k.csv DENSE_FULL_LU 30 100 10000 1 20

.PHONY: run_dense_full_qr_10k
run_dense_full_qr_10k: $(N01_PROGRAM_EXE) Makefile target_dir bench_data_dir
	$(N01_PROGRAM_EXE) -b ./bench_data/optimised_fullpiv_dense_qr_10k.csv DENSE_FULL_QR 30 100 10000 1 20

.PHONY: run_dense_par_lu_10k
run_dense_par_lu_10k: $(N01_PROGRAM_EXE) Makefile target_dir bench_data_dir
	$(N01_PROGRAM_EXE) -b ./bench_data/optimised_par_dense_lu_10k.csv DENSE_PAR_LU 30 100 10000 1 20

.PHONY: run_dense_pqr_qr_10k
run_dense_par_qr_10k: $(N01_PROGRAM_EXE) Makefile target_dir bench_data_dir
	$(N01_PROGRAM_EXE) -b ./bench_data/optimised_par_dense_qr_10k.csv DENSE_PAR_QR 30 100 10000 1 20

.PHONY: run_sherman_10k 
run_sherman_10k: $(N01_PROGRAM_EXE) Makefile target_dir bench_data_dir
	$(N01_PROGRAM_EXE) -b ./bench_data/optimised_sherman_10k.csv SHERMAN_MORRISON 30 100 10000 1 20

.PHONY: run_big_mat
run_big_mat: $(N01_PROGRAM_EXE) Makefile target_dir bench_data_dir
	$(N01_PROGRAM_EXE) -b ./bench_data/optimised_sparse_lu_100k.csv SPARSE_LU 100 100 100000 1 20
	$(N01_PROGRAM_EXE) -b ./bench_data/optimised_sherman_100k.csv SHERMAN_MORRISON 100 10 100000 1 20

