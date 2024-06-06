#ifndef FUNCTOR_H
#define FUNCTOR_H

#include <vector>

class IAccumulator {
public:
	virtual float operator() (const float &a) = 0;
};

template<typename F>
std::vector<float> cdf(const std::vector<float>& pdfv, F& acc);
#endif
