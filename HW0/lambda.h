#ifndef LAMBDA_H
#define LAMBDA_H
#include <vector>

using by_ref_t = std::vector<int> (*)(std::vector<int>&);
using by_value_t = std::vector<int> (*)(std::vector<int>);
using cmp_lambda_t = bool (*)(int &, int&);

extern by_ref_t by_ref;
extern by_value_t by_value;
extern cmp_lambda_t cmp_lambda;
extern int g_counter;

void init_globals();

std::vector<int> ascending_sort(std::vector<int>& vec);
#endif
