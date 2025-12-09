#include "qalg.hpp"

#include <iostream>

void print(std::vector<bilvector<int>>& term) {
    std::cout << "\n";
    for (int i = 0; i < 10 + 1; i++) {
        for (int j = term[i].get_max_nindex(); j <= term[i].get_max_pindex(); j++) {
            if (j == 0) {
                std::cout << " zero ";
            }
            std::cout << term[i][j] << " ";
        }
        std::cout << "\n";
    }
}


int main() { // tests
    std::vector<bilvector<int>> term(10 + 1,  bilvector<int>(0, 1, 20, 0));
    term[0][0] = 1;
    x_q_inv_pochhammer(term, 10, 5, 2);
    print(term);
    x_q_inv_pochhammer(term, 10, 5, 2);
    print(term);
    term = std::vector<bilvector<int>>(10 + 1,  bilvector<int>(0, 1, 20, 0));
    term[0][0] = 1;
    x_q_pochhammer(term, 10, 5, 2);
    print(term);
    term = std::vector<bilvector<int>>(10 + 1,  bilvector<int>(0, 1, 20, 0));
    term[0][0] = 1;
    np_q_binom(term, 0, -5, 7, false);
    print(term);

    std::cout << "\n";
}