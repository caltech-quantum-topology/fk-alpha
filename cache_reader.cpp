#include "cache_reader.hpp"
#include <iostream>

bilvector<bilvector<std::vector<std::vector<int>>>> read_pochhammer_cache(const std::string& filename, int max_rows) {
    bilvector<bilvector<std::vector<std::vector<int>>>> cache;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cout << "ERROR: Unable to open cache file: " << filename << std::endl;
        exit(0);
    }
    
    std::string line;
    std::getline(file, line, '\n');
    int num_entries = string_to_int(line);
    
    for (int entry = 0; entry < num_entries; entry++) {
        // Read the two indices and entry size
        std::getline(file, line, '\n');
        int comma_pos = line.find(",");
        int index1 = string_to_int(line.substr(0, comma_pos));
        int second_comma = line.find(",", comma_pos + 1);
        int index2 = string_to_int(line.substr(comma_pos + 1, second_comma - comma_pos - 1));
        int entry_length = string_to_int(line.substr(second_comma + 1, line.find(",", second_comma + 1) - second_comma - 1));
        
        // Limit rows to read based on max_rows parameter
        int rows_to_read = (max_rows > 0 && max_rows < entry_length) ? max_rows : entry_length;
        
        // Initialize cache entry - bilvector will auto-expand for negative indices
        cache[index1][index2].resize(rows_to_read);
        
        // Read the specified number of rows
        for (int row = 0; row < rows_to_read; row++) {
            std::getline(file, line, '\n');
            
            // First pass: read min_degree and max_degree to determine size
            size_t first_comma = line.find(",");
            size_t second_comma = line.find(",", first_comma + 1);
            
            int min_degree = string_to_int(line.substr(0, first_comma));
            int max_degree = string_to_int(line.substr(first_comma + 1, second_comma - first_comma - 1));
            int expected_coeffs = max_degree - min_degree + 1;
            
            // Preallocate row vector with exact size needed
            cache[index1][index2][row].resize(expected_coeffs + 2);
            cache[index1][index2][row][0] = min_degree;
            cache[index1][index2][row][1] = max_degree;
            
            // Parse remaining coefficients directly into preallocated space
            size_t pos = second_comma + 1;
            for (int coeff_idx = 0; coeff_idx < expected_coeffs; coeff_idx++) {
                comma_pos = line.find(",", pos);
                if (comma_pos == std::string::npos) comma_pos = line.length();
                
                std::string num_str = line.substr(pos, comma_pos - pos);
                cache[index1][index2][row][2 + coeff_idx] = string_to_int(num_str);
                pos = comma_pos + 1;
            }
        }
        
        // Skip remaining rows if we're not reading all of them
        for (int row = rows_to_read; row < entry_length; row++) {
            std::getline(file, line, '\n');
        }
    }
    
    file.close();
    return cache;
}