#include <functor.h>
class Accumulator : public IAccumulator
{
private:
    float totalSum = 0;
public:
   
    // easy default constructor
    Accumulator() : totalSum(0) {}
    
    // overriding the opertator turning this into a functor
    float operator() (const float &a) override
    {
        totalSum += a;
        return totalSum;
    }

    // Default destructor
    ~Accumulator();
};