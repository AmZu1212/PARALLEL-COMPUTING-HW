#include "lambda.h"
#include <iostream>
#include <algorithm>



// this function initializes the by_value, by_ref, cmp_lambda
// lambda pointers.
void init_globals()
{
    /*
     *  Excercise 1, a+b:
     *  Define a lambda function which gets a vector of integers by value
     *  and returns a new vector that contains only the even numbers from
     *  the input vector. Assign that lambda function to the extern global
     *  variable "by_value".
     */
    by_value = [](std::vector<int> numVector)
    {
        std::vector<int> evenNumsOnly;
        for (int num : numVector)
        {
            if (num % 2 == 0)
            {
                evenNumsOnly.push_back(num);
            }
        }
        return evenNumsOnly;
    };

    /*
     *  Excercise 1, c+d:
     *  Define a lambda function which gets a vector of integers by reference
     *  and returns a new vector that contains only the even numbers from
     *  the input vector. Assign that lambda function to the extern global
     *  variable "by_ref".
     */
    by_ref = [](std::vector<int> &numVector)
    {
        std::vector<int> evenNumsOnly;
        for (int num : numVector)
        {
            if (num % 2 == 0)
            {
                evenNumsOnly.push_back(num);
            }
        }
        return evenNumsOnly;
    };

    /*
     *  Excercise 1, e+f:
     *  Write a lambda function that compares two integers by their absolute values.
     *  This lambda function gets a global integer counter and increment it per each
     *  compare operation. Assign that lambda function to the extern global variable
     *  "cmp_lambda".
     */
    cmp_lambda = [](int &a, int &b) -> bool
    {

        int nums[2] = {a, b}; // backup
        // make sure they're both positive
        for(int i = 0; i < 2; i++)
        {
            if(nums[i] < 0)
            {
                nums[i] *= -1;
            }
        }

        //forgot to increment the comaprison counter whoops
        ++g_counter;
        
        // convention is "a < b ?"
        return nums[0] < nums[1];
    };
}


/*
*  Excercise 2, g+h:
*  Write the ascending_sort function to sort a vector of integers in 
*  ascending order of their absolute values, using the C++ STL sort. 
*  The ascending_sort() should return the sorted vector should print 
*  the total amount of comparison operations applied during the sort.
*/
std::vector<int> ascending_sort(std::vector<int> &vec)
{
    // initialize lambda pointers
    init_globals();
    // sort from beggining to end using out cmp_lambda
    std::sort(vec.begin(), vec.end(), cmp_lambda);
    // print the amount of comparisons
    std::cout << g_counter << std::endl;
    return vec;
}


