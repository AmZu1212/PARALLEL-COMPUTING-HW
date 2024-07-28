/*
 * A parallel implementation of a stencil computation to solve the
 * 2-dimensional heat diffusion problem by Jacobi iteration using TBB.
 *
 * Copyright (c) 2013 Australian National University. All rights reserved.
 */

/* Updated by Dr. Yariv Aridor, 2022 */

#include <iostream>
#include <algorithm>
#include <mutex>
#include <limits>
#include <chrono>
#include <tbb/tbb.h>
#include "heat.h"
using namespace std;

size_t n_x = 700;
size_t n_y = 700;
int max_iter = 1000;
double t_edge = 100.0;
double converge = 0.01;

// allocate data arrays
double *t_old = new double[n_x * n_y]();
double *t_new = new double[n_x * n_y]();

void delete_heap()
{
    delete[] t_old;
    delete[] t_new;
}

void init_heap()
{
    size_t i, j;
    // fix boundary values
    j = 0;
    for (i = 0; i < n_x; i++)
        t_new[j * n_x + i] = t_old[j * n_x + i] = t_edge;
    j = n_y - 1;
    for (i = 0; i < n_x; i++)
        t_new[j * n_x + i] = t_old[j * n_x + i] = t_edge;
    i = 0;
    for (j = 0; j < n_y; j++)
        t_new[j * n_x + i] = t_old[j * n_x + i] = t_edge;
    i = n_x - 1;
    for (j = 0; j < n_y; j++)
        t_new[j * n_x + i] = t_old[j * n_x + i] = t_edge;

    for (j = 1; j < n_x - 1; j++)
        for (i = 1; i < n_y - 1; i++)
            t_new[j * n_x + i] = t_old[j * n_x + i] = 0;
}

bool validate_heat(double *matrix)
{
    /* Printout the result */
    size_t i, j;
    const size_t kMaxPrint = 5;
    double reference[(kMaxPrint - 1) * (kMaxPrint - 1)] = {99.8727, 99.7457, 99.6192, 99.4935, 99.7457, 99.492, 99.2392, 98.988, 99.6192, 99.2392, 98.8607, 98.4845, 99.4935, 98.988, 98.4845, 97.984};
    for (j = 1; j < kMaxPrint; j++)
    {
        for (i = 1; i < kMaxPrint; i++)
        {
            double diff = fabs(matrix[j * n_x + i] - reference[(j - 1) * (kMaxPrint - 1) + (i - 1)]);
            if (diff > 1e-4)
                return false;
        }
    }
    return true;
}

void heat(int mode)
{
    int iter = 0;
    double max_diff;
    int rc = 0;

    while (iter < max_iter)
    {
        iter++;
        max_diff = 0.f;

        if (mode == 0)
        { // sequential execution
            for (size_t j = 1; j < n_y - 1; j++)
            {
                for (size_t i = 1; i < n_x - 1; i++)
                {
                    t_new[j * n_x + i] = 0.25 * (t_old[j * n_x + i + 1] + t_old[j * n_x + i - 1] + t_old[(j + 1) * n_x + i] + t_old[(j - 1) * n_x + i]);
                    double tdiff = fabs(t_old[j * n_x + i] - t_new[j * n_x + i]);
                    max_diff = max(max_diff, tdiff);
                }
            }
        }

        // Under (mode==1), add a parallel implementation using parallel_for with a global max_diff
        // variable to accumulate the convergence values.
        // parallel_for (global_diff w/lock)

        // basically just make a parallel for with such that the ranges are based on the rows and cols of the range variable given by tbb.
        if (mode == 1)
        {
            using namespace tbb;
            std::mutex mtx;
            parallel_for(blocked_range2d<size_t>(1, n_y - 1, 1, n_x - 1),
                         [&](const blocked_range2d<size_t> &range)
                         {
                             double localMaxDiff = 0.0;
                             for (size_t j = range.rows().begin(); j < range.rows().end(); ++j)
                             {
                                 for (size_t i = range.cols().begin(); i < range.cols().end(); ++i)
                                 {
                                     // Compute the new temperature at (j, i)
                                     t_new[j * n_x + i] = 0.25 * (t_old[j * n_x + i + 1] + t_old[j * n_x + i - 1] + t_old[(j + 1) * n_x + i] + t_old[(j - 1) * n_x + i]);
                                     // Calculate the temperature difference
                                     double tdiff = fabs(t_old[j * n_x + i] - t_new[j * n_x + i]);
                                     // Update the local max diff
                                     localMaxDiff = std::max(localMaxDiff, tdiff);
                                 }
                             }
                             // Lock and update the total max_diff (when the loc kgoes outof scope it unlocks)
                             std::lock_guard<std::mutex> lock(mtx);
                             max_diff = std::max(max_diff, localMaxDiff); // this makes the code go really slow because of mutex contention.
                         }

            );
        }
        // Under (mode==2), add a parallel implementation using parallel_for with optimize mode==1
        // by using local max_diff variables to accumulate the local max convergence values for each
        // thread before the global max_diff variable is updated.
        // parallel_for (local_diff)
        if (mode == 2)
        {
            using namespace tbb;
            // Create a vector to store the local max_diff for each thread
            // tbb::this_task_arena::max_concurrency() gives us the largest amount of chunks possible os we dont run into space issues :)
            vector<double> local_max_diffs(this_task_arena::max_concurrency(), 0.0);

            parallel_for(blocked_range2d<size_t>(1, n_y - 1, 1, n_x - 1),
                [&](const blocked_range2d<size_t> &r)
                {
                    // Get the running thread's index
                    int thread_index = this_task_arena::current_thread_index();
                    // Get reference to store the local diff later
                    double &local_diff = local_max_diffs[thread_index];

                    for (size_t j = r.rows().begin(); j < r.rows().end(); ++j)
                    {
                        for (size_t i = r.cols().begin(); i < r.cols().end(); ++i)
                        {
                            // Compute the new temperature at (j, i)
                            t_new[j * n_x + i] = 0.25 * (t_old[j * n_x + i + 1] + t_old[j * n_x + i - 1] + t_old[(j + 1) * n_x + i] + t_old[(j - 1) * n_x + i]);
                            // Calculate the temperature difference
                            double tdiff = fabs(t_old[j * n_x + i] - t_new[j * n_x + i]);
                            // Update the local max diff
                            local_diff = max(local_diff, tdiff);
                        }
                    }
                }
                );
            // Get the maximal diff from the vector for this run. this is better than a contended mutex.
            max_diff = *max_element(local_max_diffs.begin(), local_max_diffs.end());
        }

        // Under (mode==3), add a parallel implementation using parallel_reduce.
        // parallel_reduce
        if (mode == 3)
        {
            using namespace tbb;
            max_diff = parallel_reduce(blocked_range2d<size_t>(1, n_y - 1, 1, n_x - 1), 0.0, 
                [&](const blocked_range2d<size_t> &r, double init) -> double
                {
                    // update the local variable to initial value of 0.
                    double local_max_diff = init;
                    for (size_t j = r.rows().begin(); j < r.rows().end(); ++j) {
                        for (size_t i = r.cols().begin(); i < r.cols().end(); ++i) {
                            // Compute the new temperature at (j, i)
                            t_new[j * n_x + i] = 0.25 * (t_old[j * n_x + i + 1] + t_old[j * n_x + i - 1] + t_old[(j + 1) * n_x + i] + t_old[(j - 1) * n_x + i]);
                            // Calculate the temperature difference
                            double tdiff = fabs(t_old[j * n_x + i] - t_new[j * n_x + i]);
                            // Update the local max diff
                            local_max_diff = std::max(local_max_diff, tdiff);
                        }
                    }
                    return local_max_diff; 
                },
                // choose the maximum between 2 chunks
                [](double x, double y) -> double 
                {
                    return std::max(x, y);
                }
            );
        }

        // swap array pointers
        double *temp = t_new;
        t_new = t_old;
        t_old = temp;
        if (max_diff < converge)
            break;
    }

    cout << "iterations: " << iter << " convergence:" << max_diff << endl;

    // return (validate_heat(t_new) == false) ? -1 : 0;
}
