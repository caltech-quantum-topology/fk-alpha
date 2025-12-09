#include "cachealg.hpp"
#include <algorithm>
#include <stdexcept>

std::vector<bilvector<int64_t>> cache_polynomial_multiply(
    const std::vector<std::vector<int>>& cache_poly,
    const std::vector<bilvector<int64_t>>& input_poly,
    int X_MAX
) {
    // Note: Cache size can now be variable, so we handle truncation below
    
    // Check if input polynomial has correct size
    if (static_cast<int>(input_poly.size()) != X_MAX + 1) {
        throw std::runtime_error("Input polynomial size mismatch: expected " + std::to_string(X_MAX + 1) + 
                                " degrees but got " + std::to_string(input_poly.size()) + " degrees.");
    }
    
    // Initialize result with appropriate size
    std::vector<bilvector<int64_t>> result(X_MAX + 1, bilvector<int64_t>(0, 1, 20, 0));
    
    // Truncate cache to at most X_MAX + 1 entries
    int cache_size_to_use = std::min(X_MAX + 1, static_cast<int>(cache_poly.size()));
    
    // Outer loop: input polynomial degrees (0 to X_MAX)
    for (int input_deg = 0; input_deg <= X_MAX; input_deg++) {
        const bilvector<int64_t>& input_laurent = input_poly[input_deg];
        
        // Inner loop: cache degrees (0 to min(cache_size, X_MAX - input_deg))
        int max_cache_deg = std::min(cache_size_to_use - 1, X_MAX - input_deg);
        for (int cache_deg = 0; cache_deg <= max_cache_deg; cache_deg++) {
            const std::vector<int>& cache_row = cache_poly[cache_deg];
            
            // Check for malformed cache rows
            if (cache_row.size() < 3) {
                throw std::runtime_error("Malformed cache row at degree " + std::to_string(cache_deg) + 
                                       ": row has " + std::to_string(cache_row.size()) + 
                                       " elements but needs at least 3. Cache generation error.");
            }
            
            int min_laurent_deg = cache_row[0];
            int max_laurent_deg = cache_row[1];
            
            // Verify we have the right number of coefficients
            int expected_coeffs = max_laurent_deg - min_laurent_deg + 1;
            if (cache_row.size() != expected_coeffs + 2) {
                throw std::runtime_error("Malformed cache row at degree " + std::to_string(cache_deg) + 
                                       ": expected " + std::to_string(expected_coeffs + 2) + 
                                       " elements but got " + std::to_string(cache_row.size()) + 
                                       ". Cache generation error.");
            }
            
            // Calculate resulting power series degree (guaranteed <= X_MAX)
            int result_deg = input_deg + cache_deg;
            
            // Multiply cache Laurent polynomial with input Laurent polynomial
            for (int cache_laurent_idx = 0; cache_laurent_idx < expected_coeffs; cache_laurent_idx++) {
                int cache_laurent_deg = min_laurent_deg + cache_laurent_idx;
                int64_t cache_coeff = cache_row[2 + cache_laurent_idx];
                
                // Skip zero coefficients
                if (cache_coeff == 0) continue;
                
                // Multiply with all terms in input Laurent polynomial
                for (int input_laurent_deg = input_laurent.get_max_nindex(); 
                     input_laurent_deg <= input_laurent.max_pindex(); 
                     input_laurent_deg++) {
                    
                    int64_t input_coeff = input_laurent[input_laurent_deg];
                    
                    // Skip zero coefficients
                    if (input_coeff == 0) continue;
                    
                    // Calculate resulting Laurent degree
                    int result_laurent_deg = cache_laurent_deg + input_laurent_deg;
                    
                    // Add to result
                    result[result_deg][result_laurent_deg] += cache_coeff * input_coeff;
                }
            }
        }
    }
    
    return result;
}

std::vector<bilvector<int64_t>> cache_polynomial_multiply_reversed(
    const std::vector<std::vector<int>>& cache_poly,
    const std::vector<bilvector<int64_t>>& input_poly,
    int X_MAX
) {
    // Note: Cache size can now be variable, so we handle truncation below
    
    // Check if input polynomial has correct size
    if (static_cast<int>(input_poly.size()) != X_MAX + 1) {
        throw std::runtime_error("Input polynomial size mismatch: expected " + std::to_string(X_MAX + 1) + 
                                " degrees but got " + std::to_string(input_poly.size()) + " degrees.");
    }
    
    // Initialize result with appropriate size
    std::vector<bilvector<int64_t>> result(X_MAX + 1, bilvector<int64_t>(0, 1, 20, 0));
    
    // Truncate cache to at most X_MAX + 1 entries
    int cache_size_to_use = std::min(X_MAX + 1, static_cast<int>(cache_poly.size()));
    
    // Outer loop: input polynomial degrees (0 to X_MAX)
    for (int input_deg = 0; input_deg <= X_MAX; input_deg++) {
        const bilvector<int64_t>& input_laurent = input_poly[input_deg];
        
        // Inner loop: cache degrees (0 to min(cache_size, X_MAX - input_deg))
        int max_cache_deg = std::min(cache_size_to_use - 1, X_MAX - input_deg);
        for (int cache_deg = 0; cache_deg <= max_cache_deg; cache_deg++) {
            const std::vector<int>& cache_row = cache_poly[cache_deg];
            
            // Check for malformed cache rows
            if (cache_row.size() < 3) {
                throw std::runtime_error("Malformed cache row at degree " + std::to_string(cache_deg) + 
                                       ": row has " + std::to_string(cache_row.size()) + 
                                       " elements but needs at least 3. Cache generation error.");
            }
            
            int min_laurent_deg = cache_row[0];
            int max_laurent_deg = cache_row[1];
            
            // Verify we have the right number of coefficients
            int expected_coeffs = max_laurent_deg - min_laurent_deg + 1;
            if (cache_row.size() != expected_coeffs + 2) {
                throw std::runtime_error("Malformed cache row at degree " + std::to_string(cache_deg) + 
                                       ": expected " + std::to_string(expected_coeffs + 2) + 
                                       " elements but got " + std::to_string(cache_row.size()) + 
                                       ". Cache generation error.");
            }
            
            // Calculate resulting power series degree (guaranteed <= X_MAX)
            int result_deg = input_deg + cache_deg;
            
            // Multiply cache Laurent polynomial with input Laurent polynomial
            for (int cache_laurent_idx = 0; cache_laurent_idx < expected_coeffs; cache_laurent_idx++) {
                int cache_laurent_deg = min_laurent_deg + cache_laurent_idx;
                int64_t cache_coeff = cache_row[2 + cache_laurent_idx];
                
                // Skip zero coefficients
                if (cache_coeff == 0) continue;
                
                // Multiply with all terms in input Laurent polynomial
                for (int input_laurent_deg = input_laurent.get_max_nindex(); 
                     input_laurent_deg <= input_laurent.max_pindex(); 
                     input_laurent_deg++) {
                    
                    int64_t input_coeff = input_laurent[input_laurent_deg];
                    
                    // Skip zero coefficients
                    if (input_coeff == 0) continue;
                    
                    // Calculate resulting Laurent degree
                    int result_laurent_deg = -cache_laurent_deg + input_laurent_deg;
                    
                    // Add to result
                    result[result_deg][result_laurent_deg] += cache_coeff * input_coeff;
                }
            }
        }
    }
    
    return result;
}
