#include "log_add.h"

#include <Rcpp.h>

// [[Rcpp::export]]
double log_sum(const std::vector<double> & x) {
    auto iter = x.begin();
    auto end = x.end();
    double  res = *(iter++);
    for(iter; iter != end; ++iter)
        res = R::logspace_add(res, *iter);
    return res;
}

// [[Rcpp::export]]
std::vector<double> log_cumsum(const std::vector<double> & x) {
    std::vector<double> res(x.size());
    auto x_it = x.begin();
    auto res_it = res.begin();
    auto res_end = res.end();
    auto res_it_lag = res_it;
    *(res_it++) = *(x_it++);
    for(; res_it != res_end; ++res_it) 
        *res_it = R::logspace_add(*(res_it_lag++), *(x_it++));
    return res;
}
