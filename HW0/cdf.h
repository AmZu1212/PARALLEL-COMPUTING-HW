#include "functor.h"
class Accumulator : public IAccumulator
{
private:
    float totalSum = 0;

public:
    // easy default constructor
    Accumulator() : totalSum(0) {}

    // overriding the opertator turning this into a functor
    float operator()(const float &a) override
    {
        totalSum += a;
        return totalSum;
    }

    // Default destructor
    //~Accumulator();
};

// PDF is a probability density function vector
template <typename F>
std::vector<float> cdf(const std::vector<float> &pdfv, F &acc)
{
    // result is a cdf, which means it is a cummulative sum
    // so with each step the value increases
    std::vector<float> result;

    // acc is the accumulation function i think
    float sum = 0;
    for (float current : pdfv)
    {
        sum = acc(current);
        result.push_back(sum);
    }
    return result;
}
