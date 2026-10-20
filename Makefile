CXX ?= g++
CXXFLAGS ?= -O3 -std=c++20 -pthread
CLANG_TIDY ?= clang-tidy

INDEX_SOURCES = main.cpp encode_kmer_2bit.cpp hash.cpp indexGenome.cpp \
    serialiseKmersMap.cpp packGenomeBlob.cpp generate_seed_matrix.cpp \
    SeedMatrix.cpp MonteCarlo.cpp greedyOrder.cpp buildGreedyOrder.cpp \
    report_seed_matrix.cpp
SEARCH_SOURCES = main.cpp kmerUtilities.cpp serialiseKmersMap.cpp \
    smithWaterman.cpp hash.cpp ssw_cpp.cpp findSpans.cpp SeedMatrix.cpp
# ssw_cpp.cpp is bundled third-party code.
LINT_SEARCH_SOURCES = $(filter-out ssw_cpp.cpp,$(SEARCH_SOURCES))

INDEX_OBJECTS = $(patsubst %.cpp,build/indexer/%.o,$(INDEX_SOURCES))
SEARCH_OBJECTS = $(patsubst %.cpp,build/search/%.o,$(SEARCH_SOURCES)) build/search/ssw.o
INDEX_HEADERS = $(wildcard indexer/*.h indexer/*.hpp indexer/headers/*)
SEARCH_HEADERS = $(wildcard search/*.h search/*.hpp search/headers/*)

# These names are defined by both tools. Keep each implementation distinct in
# the linked executable, including their different SeedMatrix layouts.
INDEX_RENAMES = main GENOME_DIR getBaseName murmurHash3 buildGreedyOrder \
    loadGreedyOrders read_entire_file saveGreedyOrders writeTextBlobToFile \
    readTextBlobFromFile initialiseGreedyOrders xorHash REFERENCE fileNames \
    SeedMatrix GreedyState StripePattern gGreedyOrders gPattern
INDEX_DEFINES = $(foreach name,$(INDEX_RENAMES),-D$(name)=index_$(name))

ifdef USE_EDLIB
SEARCH_OBJECTS += build/search/alignSequences_edlib.o build/search/edlib.o
SEARCH_FLAGS = -DUSE_EDLIB -Isearch/edlib-1.2.7/edlib/include
LINT_SEARCH_SOURCES += alignSequences_edlib.cpp
endif

.PHONY: all clean test lint
all: skippy

skippy: main.cpp $(INDEX_OBJECTS) $(SEARCH_OBJECTS)
	$(CXX) $(CXXFLAGS) main.cpp $(INDEX_OBJECTS) $(SEARCH_OBJECTS) -o $@

build/indexer/%.o: indexer/%.cpp $(INDEX_HEADERS) | build/indexer
	$(CXX) $(CXXFLAGS) $(INDEX_DEFINES) -Iindexer/headers -c $< -o $@

build/search/%.o: search/%.cpp $(SEARCH_HEADERS) | build/search
	$(CXX) $(CXXFLAGS) $(SEARCH_FLAGS) -Dmain=search_main -Isearch/headers -c $< -o $@

build/search/ssw.o: search/ssw.c $(SEARCH_HEADERS) | build/search
	$(CXX) $(CXXFLAGS) -Isearch/headers -c $< -o $@

build/search/edlib.o: search/edlib-1.2.7/edlib/src/edlib.cpp | build/search
	$(CXX) $(CXXFLAGS) $(SEARCH_FLAGS) -c $< -o $@

build/indexer build/search:
	mkdir -p $@

clean:
	rm -rf build skippy

test: skippy
	python3 -m tests.test_cli
	python3 -m tests.test_mito

lint:
	$(CLANG_TIDY) main.cpp -- $(CXXFLAGS)
	@for source in $(addprefix indexer/,$(INDEX_SOURCES)); do \
		$(CLANG_TIDY) "$$source" -- $(CXXFLAGS) $(INDEX_DEFINES) -Iindexer/headers || exit $$?; \
	done
	@for source in $(addprefix search/,$(LINT_SEARCH_SOURCES)); do \
		$(CLANG_TIDY) "$$source" -- $(CXXFLAGS) $(SEARCH_FLAGS) -Dmain=search_main \
			-Isearch/headers || exit $$?; \
	done
