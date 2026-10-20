# Compiler and flags
CC = g++
#CFLAGS = -pthread -O2 -std=c++20 -Iedlib-1.2.7/edlib/include
CFLAGS = -pthread -O3 -std=c++20

# edlib is opt-in: smithWaterman.cpp only calls it under #ifdef USE_EDLIB.
# Build with `make USE_EDLIB=1`, which needs edlib-1.2.7 unpacked here.
ifdef USE_EDLIB
CFLAGS += -Iedlib-1.2.7/edlib/include -DUSE_EDLIB
endif

# Source directory and files
SOURCE_DIR = .
CPP_SOURCES = \
main.cpp \
kmerUtilities.cpp \
serialiseKmersMap.cpp \
smithWaterman.cpp \
hash.cpp \
ssw_cpp.cpp \
findSpans.cpp \
SeedMatrix.cpp

C_SOURCES = ssw.c

EDLIB_SRC = edlib-1.2.7/edlib/src/edlib.cpp

ifdef USE_EDLIB
CPP_SOURCES += alignSequences_edlib.cpp
endif

SOURCES = $(CPP_SOURCES) $(C_SOURCES)

# Header directory
HEADER_DIR = headers
HEADERS = $(wildcard $(HEADER_DIR)/*.h)

# Object directory and files
OBJECT_DIR = objects
CPP_OBJECTS = $(CPP_SOURCES:%.cpp=$(OBJECT_DIR)/%.o)
C_OBJECTS   = $(C_SOURCES:%.c=$(OBJECT_DIR)/%.o)
EDLIB_OBJ   = $(OBJECT_DIR)/edlib.o
OBJECTS     = $(CPP_OBJECTS) $(C_OBJECTS)

ifdef USE_EDLIB
OBJECTS += $(EDLIB_OBJ)
endif

# Executable name
EXECUTABLE = skippy

.PHONY: all
all:
	$(MAKE) -C .. USE_EDLIB=$(USE_EDLIB)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@

$(CPP_OBJECTS): $(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp $(HEADERS) | $(OBJECT_DIR)
	$(CC) $(CFLAGS) -I$(HEADER_DIR) -c $< -o $@

$(C_OBJECTS): $(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.c $(HEADERS) | $(OBJECT_DIR)
	$(CC) $(CFLAGS) -I$(HEADER_DIR) -c $< -o $@

$(EDLIB_OBJ): $(EDLIB_SRC) | $(OBJECT_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJECT_DIR):
	mkdir -p $(OBJECT_DIR)

clean:
	rm -rf $(OBJECT_DIR) $(EXECUTABLE)
