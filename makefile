CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -pedantic
TARGET = simplify

SRCS = simplify.cpp polygon.cpp geometry.cpp simplifier.cpp
OBJS = $(SRCS:.cpp=.o)

DATA_DIR = data
OUT_DIR = outputs
TARGET_VERTICES ?= 50

INPUTS = \
	$(DATA_DIR)/input_blob_with_two_holes.csv \
	$(DATA_DIR)/input_cushion_with_hexagonal_hole.csv \
	$(DATA_DIR)/input_lake_with_two_islands.csv \
	$(DATA_DIR)/input_original_01.csv \
	$(DATA_DIR)/input_original_02.csv \
	$(DATA_DIR)/input_original_03.csv \
	$(DATA_DIR)/input_original_04.csv \
	$(DATA_DIR)/input_original_05.csv \
	$(DATA_DIR)/input_original_06.csv \
	$(DATA_DIR)/input_original_07.csv \
	$(DATA_DIR)/input_original_08.csv \
	$(DATA_DIR)/input_original_09.csv \
	$(DATA_DIR)/input_original_10.csv \
	$(DATA_DIR)/input_rectangle_with_two_holes.csv \
	$(DATA_DIR)/input_wavy_with_three_holes.csv

.PHONY: all clean test dirs

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

simplify.o: simplify.cpp polygon.h simplifier.h geometry.h
	$(CXX) $(CXXFLAGS) -c simplify.cpp

polygon.o: polygon.cpp polygon.h
	$(CXX) $(CXXFLAGS) -c polygon.cpp

geometry.o: geometry.cpp geometry.h polygon.h
	$(CXX) $(CXXFLAGS) -c geometry.cpp

simplifier.o: simplifier.cpp simplifier.h geometry.h polygon.h
	$(CXX) $(CXXFLAGS) -c simplifier.cpp

dirs:
	mkdir -p $(OUT_DIR)

test: all dirs
	for f in $(INPUTS); do \
		echo "Running $$f"; \
		./$(TARGET) "$$f" $(TARGET_VERTICES) > "$(OUT_DIR)/$$(basename $$f .csv)_out.txt"; \
	done

clean:
	rm -f $(OBJS) $(TARGET)
	rm -rf $(OUT_DIR)