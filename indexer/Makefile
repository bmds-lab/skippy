# Compiler and flags
CC = g++
#CFLAGS = -O3 -std=c++11 -g 
CFLAGS = -O3 -std=c++17 -g -pthread
#CFLAGS = -O3 -std=c++11 -g -fprofile-instr-generate -fcoverage-mapping
#CFLAGS = -g -o -std=c++11

# Source directory and files
SOURCE_DIR = .
#SOURCES = main.cpp encode_kmer_2bit.cpp hash.cpp indexGenome.cpp serialiseKmersMap.cpp packGenomeBlob.cpp
SOURCES = \
    main.cpp \
    encode_kmer_2bit.cpp \
    hash.cpp \
    indexGenome.cpp \
    serialiseKmersMap.cpp \
    packGenomeBlob.cpp \
    generate_seed_matrix.cpp \
    SeedMatrix.cpp \
    MonteCarlo.cpp \
    greedyOrder.cpp \
    buildGreedyOrder.cpp \
    report_seed_matrix.cpp
# Header directory
HEADER_DIR = headers

# Object directory and files
OBJECT_DIR = objects
OBJECTS = $(SOURCES:%.cpp=$(OBJECT_DIR)/%.o)

# Executable name
EXECUTABLE = indexReference

.PHONY: all
all:
	$(MAKE) -C ..

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@

$(OBJECTS): $(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp | $(OBJECT_DIR)
	$(CC) $(CFLAGS) -I$(HEADER_DIR) -c $< -o $@

$(OBJECT_DIR):
	mkdir -p $(OBJECT_DIR)

clean:
	rm -rf $(OBJECT_DIR) $(EXECUTABLE) test_contract

# ── Test harness ──────────────────────────────────────────────────────────────
# Standalone driver for generateSeedMatrix, so a contract can be checked without
# reindexing a whole reference.
TEST_SOURCES = \
    test_contract.cpp \
    encode_kmer_2bit.cpp \
    hash.cpp \
    generate_seed_matrix.cpp \
    SeedMatrix.cpp \
    MonteCarlo.cpp \
    greedyOrder.cpp \
    buildGreedyOrder.cpp \
    report_seed_matrix.cpp

TEST_OBJECTS = $(TEST_SOURCES:%.cpp=$(OBJECT_DIR)/%.o)

test_contract: $(TEST_OBJECTS)
	$(CC) $(CFLAGS) $(TEST_OBJECTS) -o $@

$(OBJECT_DIR)/test_contract.o: $(SOURCE_DIR)/test_contract.cpp | $(OBJECT_DIR)
	$(CC) $(CFLAGS) -I$(HEADER_DIR) -c $< -o $@

.PHONY: test
test: test_contract
	./test_contract 2X3X2X4X2X3 skippy subs 500 300 20

