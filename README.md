# FK-Alpha: Gukov-Manolescu Knot Invariant Computation

A computational framework for computing Gukov-Manolescu two-variable series F_K(x,q) for knot complements using braid representations.

## 🚀 Quick Start

### Method 1: Complete Setup and Run
```bash
# 1. Install Python dependencies
pip install -r requirements.txt

# 2. Compile the FK computation engine for your platform
g++ -O3 -o fk_linux fk.cpp -std=c++17

# 3. Edit main.py to set your braid (lines 132-137)
#    Example: braids = [[1, 1, 1]] for trefoil knot

# 4. Run the computation
python3 main.py

# 5. Check results in Data/Output/[knot_name].json
```

### Method 2: Direct C++ Binary (Advanced)
```bash
# Use pre-computed input files
./fk_linux "Data/Input/input_file" "Data/Output/output_file"
```

### Method 3: Quick Test with Pre-computed Example
```bash
# Run trefoil knot [1,1,1] with degree 5 (no Python needed)
g++ -O3 -o fk_linux fk.cpp -std=c++17
./fk_linux "Data/Input/braid_len_3_deg_5_knot_0" "Data/Output/trefoil_example"
cat Data/Output/trefoil_example.json
```

## 🏗️ Build Requirements

### Dependencies
```bash
# Install all required packages
pip install -r requirements.txt

# Or install individually:
pip install gurobipy numpy sympy

# System requirements
- C++17 compatible compiler (g++, clang++)
- Python 3.7+
- Linux/macOS (tested on Ubuntu/Debian)
```

### Important Notes
- **Platform-specific compilation required**: The included `_` binary is for macOS and won't work on Linux
- **Gurobi license**: Academic license available for free, commercial license required for commercial use
- **Memory**: Computations scale as O(degree^components), keep degree ≤ 25 for practical use

### Compilation
```bash
# Basic compilation (recommended)
g++ -O3 -o fk_linux fk.cpp -std=c++17

# Debug version
g++ -g -O0 -o fk_debug fk.cpp -std=c++17

# With profiling
g++ -O3 -pg -o fk_profile fk.cpp -std=c++17
```

## 📁 Project Structure

### Core Files
- **`fk.cpp`** - Main C++ implementation of Gukov-Manolescu computation
- **`main.py`** - Python orchestrator for complete workflow
- **`sign_diagrams.py`** - Braid preprocessing and sign assignment
- **`braid_ilp.py`** - Integer Linear Programming setup
- **`alex.py`** - Alexander polynomial verification

### Header Files
- **`qalg.hpp`** - Quantum algebra operations (q-binomials, q-Pochhammer)
- **`solution_pool.hpp`** - ILP solver interface
- **`linalg.hpp`** - Linear algebra utilities
- **`bilvector.hpp`** - Bi-indexed vector for polynomial storage

### Data Directories
- **`Data/Input/`** - CSV files with braid and ILP data
- **`Data/Output/`** - JSON results with computed invariants
- **`Data/Caches/`** - Cached quantum computations

## 💻 Usage Examples

### Example 1: Trefoil Knot [1,1,1]
```bash
# Method 1: Use pre-computed input
./fk_linux "Data/Input/braid_len_3_deg_5_knot_0" "Data/Output/trefoil"

# Method 2: Generate input via Python
# Edit main.py:
#   braids = [[1, 1, 1]]
#   titres = ["trefoil"]
python main.py
```

### Example 2: Custom Braid
```python
# Edit main.py
def main():
    DEGREE = 10
    titres = ["my_knot"]
    braids = [[1, -2, 1, -2]]  # Figure-eight knot
    
    for title, braid in zip(titres, braids):
        main(title, braid, DEGREE)
```

### Example 3: Direct C++ Usage
```bash
# Assuming you have prepared input_file.csv
./fk_linux "path/to/input_file" "path/to/output_file"

# Check results
cat path/to/output_file.json
```

## 📊 Input/Output Format

### Input CSV Format
```csv
degree,                    # Computation degree bound
components,                # Number of link components  
writhe,                    # Writhe of the link
gen1,r1,gen2,r2,...       # Crossing generators and R-matrix types
component1,component2,...  # Closed strand components
top1,bottom1,top2,bot2,... # Crossing component assignments
criteria_line1,            # ILP criteria constraints
/                          # Separator
inequality_line1,          # ILP inequality constraints  
/                          # Separator
assignment_data            # Symbolic variable assignments
```

### Output JSON Format
```json
{
  "braid": [1,1,1],
  "inversions": [0,0,0,0,0,0],
  "coefficient_q_powers": [
    [[1,-1]],      // x^0: -q^1 
    [],            // x^1: no terms
    [[2,1]],       // x^2: q^2
    [[3,1]]        // x^3: q^3
  ]
}
```

## 🧪 Testing

```bash
# Run trefoil computation test
./fk_linux "Data/Input/braid_len_3_deg_5_knot_0" "Data/Output/test_trefoil"

# Verify against known Alexander polynomial
python main.py  # Includes verification step

# Performance testing
time ./fk_linux "Data/Input/braid_len_8_deg_10_knot_1" "Data/Output/perf_test"
```

### Platform-Specific Notes
- **Linux**: Use `g++` compiler, ensure `build-essential` package installed
- **macOS**: Use `clang++` or `g++` from Homebrew 
- **Windows**: Use WSL or Visual Studio with C++17 support