### Data Structures Assignment 2: Polygon Simplification

## Setup & Compilation
1. Open a WSL terminal in this directory.
2. Run `make` to compile the source code. This will generate the `simplify` executable.

## Usage
Run the program using the following syntax:
```bash
./simplify <input_file.csv> <target_vertices>
```

**Basic Example:**
```bash
./simplify data/sample.csv 10
```
*By default, the simplified CSV data and area summaries are printed to the console (standard output), and debug messages are printed to standard error.*

## Useful Commands (Saving Output)
To keep your console clean and save the outputs to files, use shell redirection:
```bash
./simplify data/sample.csv 10 > output.csv 2> debug.txt 
```
*   `> output.csv`: Saves the simplified polygon data and final summaries.
*   `2> debug.txt`: Saves the debug logs.

## Dependencies
This project is built using standard C++17. It does not require any third-party libraries to compile or run.

## Test Results
Running the simplifier on the provided `sample.csv` with a target of 10 vertices yields the following results:

**Input:** `data/sample.csv` (12 vertices)  
**Target:** 10 vertices

**Summary Output:**
```
Total signed area in input: 3.210000e+00
Total signed area in output: 3.210000e+00
Total areal displacement: 2.071795e-01
```
