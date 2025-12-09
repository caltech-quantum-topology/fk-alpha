# FK-Alpha: Gukov-Manolescu Knot Invariant Computation

A computational framework for computing Gukov-Manolescu two-variable series F_K(x,q) for knot complements using braid representations. This implementation processes crossing data through quantum algebraic operations and outputs polynomial coefficients representing knot invariants.

## 🚀 Quick Start

### Method 1: Using Python Interface (Recommended)
```bash
# Edit main.py to set your braid and run computation
python main.py
```

### Method 2: Direct C++ Binary
```bash
# Compile the FK computation engine
g++ -O3 -o fk_linux fk.cpp -std=c++17

# Run computation with input/output files
./fk_linux "Data/Input/input_file" "Data/Output/output_file"
```

### Method 3: Using Pre-computed Examples
```bash
# Run trefoil knot [1,1,1] with degree 5
./fk_linux "Data/Input/braid_len_3_deg_5_knot_0" "Data/Output/trefoil_example"
```

## 🏗️ Build Requirements

### Dependencies
```bash
# Required Python packages
pip install sympy gurobipy

# System requirements
- C++17 compatible compiler (g++, clang++)
- Python 3.7+
- Gurobi optimizer (for ILP solving)
```

### Compilation
```bash
# Basic compilation
g++ -O3 -o fk_linux fk.cpp -std=c++17

# With debug information
g++ -g -O2 -o fk_debug fk.cpp -std=c++17

# For profiling
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

This represents: **F_K(x,q) = -q·x⁰ + q²·x² + q³·x³**

## 🔬 Mathematical Background

The Gukov-Manolescu invariant F_K(x,q) is computed through:

1. **Braid Preprocessing**: Convert knot braid to crossing data
2. **ILP Setup**: Generate angle state constraints  
3. **Quantum Algebra**: Apply q-binomial coefficients and q-Pochhammer symbols
4. **Accumulation**: Sum contributions from all feasible angle states

The computation follows the framework from:
*S. Gukov, C. Manolescu, "A two-variable series for knot complements"*

## ⚙️ Performance Notes

- **Memory**: Scales as O(degree^components)
- **Time**: Depends on ILP solver efficiency
- **Caching**: Pre-computed quantum operations stored in `Data/Caches/`
- **Degree Limits**: Practical computation typically degree ≤ 25

## 🧪 Testing

```bash
# Run trefoil computation test
./fk_linux "Data/Input/braid_len_3_deg_5_knot_0" "Data/Output/test_trefoil"

# Verify against known Alexander polynomial
python main.py  # Includes verification step

# Performance testing
time ./fk_linux "Data/Input/braid_len_8_deg_10_knot_1" "Data/Output/perf_test"
```

## 🛠️ Troubleshooting

### Common Issues

1. **Gurobi License**: Ensure Gurobi is properly licensed
2. **Binary Format**: Recompile if moving between architectures
3. **Memory Limits**: Reduce degree for large braids
4. **Missing Input**: Check CSV file format matches specification

### Debug Mode
```bash
# Compile with debug info
g++ -g -DDEBUG -o fk_debug fk.cpp -std=c++17

# Run with detailed output
./fk_debug "input" "output" 2>&1 | tee debug.log
```