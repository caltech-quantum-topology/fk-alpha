#include <vector>
#include <cstddef>
#include <cstdint>
#include <string>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <cmath>
#include <array>
#include <cstdio> 
#include <stdexcept>
#include <limits>

#include "solution_pool.hpp"
#include "string_to_numeric.hpp"
#include "linalg64.hpp"
#include "bilvector.hpp"
#include "qalg64.hpp"

int main() {
    int x_deg = 21; // this value can be adjusted as needed 
    std::vector<int> blocks = {1};  // keep this value as {1}, even for links 
    int components = 1;  // keep this value as 1, even for links 
                    
    std::string data_file = "Data/Caches/caching_knots_pochhammers.txt";
    std::ofstream data_out;
    data_out.open(data_file, std::ofstream::out | std::ofstream::trunc);
    if (!data_out.is_open()) {
        std::cerr << "Error opening file: " << data_file << std::endl;
        return 1;
    }
    // Use a stringstream to collect entries, then write with correct count
    int count = 0;
    for (int k = 0; k < x_deg; k++) {
         for (int j = 0; j <= k + 1; j++) {
                count++;
            }
    }
    for (int k = -1; k >= -x_deg; k--) {
         for (int j = -x_deg; j <= k + 1; j++) {
                count++;
          }
     }
    data_out << count << "\n";
    
    for (int k = 0; k < x_deg; k++) {
        for (int j = 0; j <= k + 1; j++) {
             std::cout << "Computing Pochhammer^(+1) for k=" << k << ", j=" << j << "\r";
             
             // Initialize polynomial term storage with initial coefficient
             std::vector<int> X_MAX = {k - j + 1};
             std::vector<bilvector<int64_t>> term(X_MAX[0] + 1, bilvector<int64_t>(0, 1, 20, 0));
             term[0][0] = 1;  // Set initial coefficient to 1
             x_q_pochhammer(term, k, j, 0, components, X_MAX,  blocks);
             
             // Write entry to buffer
             data_out << k << "," << j << "," << X_MAX[0] + 1 << ",\n";
             for (int i = 0; i < term.size(); i++) {
                 data_out << term[i].get_max_nindex() << "," << term[i].get_max_pindex() << ","; // Write min and max q powers
                 // Write all coefficient-q_power pairs for this component combination
                 for (int j_inner = term[i].get_max_nindex(); j_inner <= term[i].get_max_pindex(); j_inner++) {
                     data_out << term[i][j_inner] << ",";  // [q_power, coefficient]
                 }
                 data_out << "\n";
             }
         }
    }
    for (int k = -1; k >= -x_deg; k--) {
         for (int j = -x_deg; j <= k + 1; j++) {
              std::cout << "Computing Pochhammer^(+1) for k=" << k << ", j=" << j << "\r";
              // Initialize polynomial term storage with initial coefficient
              std::vector<int> X_MAX = {k - j + 1};
              std::vector<bilvector<int64_t>> term(X_MAX[0] + 1,  bilvector<int64_t>(0, 1, 20, 0));
              term[0][0] = 1;  // Set initial coefficient to 1
              x_q_pochhammer(term, k, j, 0, components, X_MAX, blocks);
              
              // Write entry to buffer
              data_out << k << "," << j << "," << X_MAX[0] + 1 << ",\n";
              for (int i = 0; i < term.size(); i++) {
                  data_out << term[i].get_max_nindex() << "," << term[i].get_max_pindex() << ","; // Write min and max q powers
                  // Write all coefficient-q_power pairs for this component combination
                  for (int j_inner = term[i].get_max_nindex(); j_inner <= term[i].get_max_pindex(); j_inner++) {
                      data_out << term[i][j_inner] << ",";  // [q_power, coefficient]
                  }
                  data_out << "\n";
              }
          }
     }
     
    data_out.close();

    int inv_count = 0;
    for (int k = 0; k < x_deg; k++) {
         for (int j = 0; j > -x_deg; j--) {
                inv_count++;
          }
     }

    // Cache inverse Pochhammer symbols
    std::string inv_data_file = "Data/Caches/caching_knots_inv_pochhammers.txt";
    std::ofstream inv_data_out;
    inv_data_out.open(inv_data_file, std::ofstream::out | std::ofstream::trunc);
    if (!inv_data_out.is_open()) {
        std::cerr << "Error opening file: " << inv_data_file << std::endl;
        return 1;
    }
    inv_data_out << inv_count << "\n";
    for (int k = 0; k < x_deg; k++) {
        for (int j = 0; j > -x_deg; j--) {
             std::cout << "Computing Pochhammer^(-1) for k=" << k << ", j=" << j << "\r";
             
             // Initialize polynomial term storage with initial coefficient
             std::vector<int> X_MAX = {x_deg - k / 2};
             std::vector<bilvector<int64_t>> term(X_MAX[0] + 1,  bilvector<int64_t>(0, 1, 20, 0));
             term[0][0] = 1;  // Set initial coefficient to 1
             x_q_inv_pochhammer(term, k, j, 0, components, X_MAX, blocks);
             
             // Write entry
             inv_data_out << k << "," << j << "," << X_MAX[0] + 1 << ",\n";
             for (int i = 0; i < term.size(); i++) {
                 inv_data_out << term[i].get_max_nindex() << "," << term[i].get_max_pindex() << ","; // Write min and max q powers
                 // Write all coefficient-q_power pairs for this component combination
                 for (int j_inner = term[i].get_max_nindex(); j_inner <= term[i].get_max_pindex(); j_inner++) {
                     inv_data_out << term[i][j_inner] << ",";  // [q_power, coefficient]
                 }
                 inv_data_out << "\n";
             }
         }
    }
    
    inv_data_out.close();

    int binom_count = 0;
     for (int k = 0; k <= x_deg; k++) {
          for (int j = 0; j <= k; j++) {
                 binom_count++;
           }
      }
     for (int k = 0; k <= x_deg; k++) {
          for (int j = -1; j >= -x_deg; j--) {
                 binom_count++;
           }
      }
     // Cache q-binomial coefficients 
     std::string binom_data_file = "Data/Caches/caching_knots_binomials.txt";
     std::ofstream binom_data_out;
     binom_data_out.open(binom_data_file, std::ofstream::out | std::ofstream::trunc);
     if (!binom_data_out.is_open()) {
         std::cerr << "Error opening file: " << binom_data_file << std::endl;
         return 1;
     }
     binom_data_out << binom_count << "\n";
     for (int k = 0; k <= x_deg; k++) {
           for (int j = 0; j <= k; j++) {
              std::cout << "Computing Binom for k=" << k << ", j=" << j << "\n"; 
              std::vector<bilvector<int64_t>> term = std::vector<bilvector<int64_t>>(1,  bilvector<int64_t>(0, 1, 20, 0));
              term[0][0] = 1;
              pp_q_binom(term, k, j, false);

              // Write entry
              binom_data_out << k << "," << j << "," << 1 <<",\n";
              for (int i = 0; i < term.size(); i++) {
                  binom_data_out << term[i].get_max_nindex() << "," << term[i].get_max_pindex() << ","; // Write min and max q powers
                  // Write all coefficient-q_power pairs for this component combination
                  for (int j_inner = term[i].get_max_nindex(); j_inner <= term[i].get_max_pindex(); j_inner++) {
                      binom_data_out << term[i][j_inner] << ",";  // [q_power, coefficient]
                  }
                  binom_data_out << "\n";
              }
          }
     }
     for (int k = 0; k <= x_deg; k++) {
           for (int j = -1; j >= -x_deg; j--) {
               std::cout << -x_deg << "\n";
               std::cout << "Computing Binom for k=" << k << ", j=" << j << "\r";
               std::vector<bilvector<int64_t>> term = std::vector<bilvector<int64_t>>(1,  bilvector<int64_t>(0, 1, 20, 0));
               term[0][0] = 1;
               np_q_binom(term, j, k, false);

               // Write entry
               binom_data_out << j << "," << k << "," << 1 << ",\n";
               for (int i = 0; i < term.size(); i++) {
                   binom_data_out << term[i].get_max_nindex() << "," << term[i].get_max_pindex() << ","; // Write min and max q powers
                   // Write all coefficient-q_power pairs for this component combination
                   for (int j_inner = term[i].get_max_nindex(); j_inner <= term[i].get_max_pindex(); j_inner++) {
                       binom_data_out << term[i][j_inner] << ",";  // [q_power, coefficient]
                   }
                   binom_data_out << "\n";
               }
           }
      }
     binom_data_out.close();
}
