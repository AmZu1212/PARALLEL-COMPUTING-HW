#ifndef NBODY_H

#include <cmath>
#include <tbb/tbb.h>        // for TBB
#include <immintrin.h>      // for AVX
#define TBBGRAIN 1024       // 4 (SIZE OF FLOAT) * 6 (SIZE OF PARTICLE) * 1024(TBB GRAIN) =~ 24KB
using namespace std;

const int nParticles = 16384 * 2;
const float dt = 0.01f;
const float softening = 1e-20f;

typedef struct
{
    float x, y, z, vx, vy, vz;
} OneParticle;

float global_X[nParticles];
float global_Y[nParticles];
float global_Z[nParticles];
float global_Vx[nParticles];
float global_Vy[nParticles];
float global_Vz[nParticles];


struct ParticleType
{
    float x, y, z;
    float vx, vy, vz;
    float trash1, trash2;
};

ParticleType particles[nParticles];



// ======================= SUPPORT VECTORS ==========================================

__m256 zeroVector = _mm256_set1_ps(0.0f);
__m256 oneVector = _mm256_set1_ps(1.0f);
__m256 dtVector = _mm256_set1_ps(dt);
__m256 softVector = _mm256_set1_ps(softening);

// ========================== SERIAL FUNCTIONS ==========================================

void get_particle_serial(int i, OneParticle *p)
{
    p->x = particles[i].x;
    p->y = particles[i].y;
    p->z = particles[i].z;
    p->vx = particles[i].vx;
    p->vy = particles[i].vy;
    p->vz = particles[i].vz;
}
void init_particles_serial()
{
    for (unsigned int i = 0; i < nParticles; i++)
    {
        particles[i].x = (float)(i % 15);
        particles[i].y = (float)((i * i) % 15);
        particles[i].z = (float)((i * i * 3) % 15);
        particles[i].vx = 0;
        particles[i].vy = 0;
        particles[i].vz = 0;
    }
}
void move_particles_serial()
{

    // Loop over particles that experience force
    for (int i = 0; i < nParticles; i++)
    {

        // Components of the gravity force on particle i
        float Fx = 0, Fy = 0, Fz = 0;

        // Loop over particles that exert force: vectorization expected here
        for (int j = 0; j < nParticles; j++)
        {

            // Avoid singularity and interaction with self
            const float softening = 1e-20f;

            // Newton's law of universal gravity
            const float dx = particles[j].x - particles[i].x;
            const float dy = particles[j].y - particles[i].y;
            const float dz = particles[j].z - particles[i].z;

            const float rr1 = 1.0f / sqrt(dx * dx + dy * dy + dz * dz + softening);
            const float drPowerN32 = rr1 * rr1 * rr1;
            // Calculate the net force
            Fx += dx * drPowerN32;
            Fy += dy * drPowerN32;
            Fz += dz * drPowerN32;
        }

        // Accelerate particles in response to the gravitational force
        particles[i].vx += dt * Fx;
        particles[i].vy += dt * Fy;
        particles[i].vz += dt * Fz;
    }

    // Move particles according to their velocities
    // O(N) work, so using a serial loop
    for (int i = 0; i < nParticles; i++)
    {
        particles[i].x += particles[i].vx * dt;
        particles[i].y += particles[i].vy * dt;
        particles[i].z += particles[i].vz * dt;
    }
}

// ========================== NEW DATABASE =============================================

// this gets 30x speed up!!!!
// boring, no optimizations here.
void get_particle_parallel(int i, OneParticle *p){
    p->x = global_X[i];
    p->y = global_Y[i];
    p->z = global_Z[i];
    p->vx = global_Vx[i];
    p->vy = global_Vy[i];
    p->vz = global_Vz[i];
}
void init_particles_parallel()
{
    tbb::parallel_for( // (range, lambda)
    tbb::blocked_range<unsigned int>(0, nParticles, TBBGRAIN),
    [](tbb::blocked_range<unsigned int> &range)
    {
        for (unsigned int i = range.begin(); i < range.end(); i++)
        {
            global_X[i]  = (float)(i % 15);
            global_Y[i]  = (float)(((i * i) % 15));
            global_Z[i]  = (float)((i * i * 3) % 15);
            global_Vx[i] = 0.0f;
            global_Vy[i] = 0.0f;
            global_Vz[i] = 0.0f;
        }
    });
}
void move_particles_parallel()
{
    tbb::parallel_for(tbb::blocked_range<unsigned int>(0, nParticles, TBBGRAIN),
    [](const tbb::blocked_range<unsigned int> &range)
    {
        for (unsigned int i = range.begin(); i < range.end(); ++i)
        {

            // Components of the gravity force on particle i
            float Fx = 0, Fy = 0, Fz = 0;
            __m256 FxVector = zeroVector;
            __m256 FyVector = zeroVector;
            __m256 FzVector = zeroVector;

            // creation of x, y, z vectors for calculations. because in the j iterations they stay the same.
            __m256 PixVector = _mm256_set1_ps(global_X[i]);
            __m256 PiyVector = _mm256_set1_ps(global_Y[i]);
            __m256 PizVector = _mm256_set1_ps(global_Z[i]);

            // i will vectorize in jumps of j+8's.
            for (unsigned int j = 0; j < nParticles; j += 8)
            {
                // creation of more x, y, z vectors for dx/y/z.
                __m256 PjxVector = _mm256_loadu_ps(&global_X[j]);
                __m256 PjyVector = _mm256_loadu_ps(&global_Y[j]);
                __m256 PjzVector = _mm256_loadu_ps(&global_Z[j]);

                // Newton's law of universal gravity
                __m256 dxVector = _mm256_sub_ps(PjxVector, PixVector);
                __m256 dyVector = _mm256_sub_ps(PjyVector, PiyVector);
                __m256 dzVector = _mm256_sub_ps(PjzVector, PizVector);

                // a TON of calcs for R1 ^ 3
                __m256 dxPow2 = _mm256_mul_ps(dxVector, dxVector);
                __m256 dyPow2 = _mm256_mul_ps(dyVector, dyVector);
                __m256 dzPow2 = _mm256_mul_ps(dzVector, dzVector);
                // denominator calcs
                __m256 temp1 = _mm256_add_ps(dxPow2, softVector);
                __m256 temp2 = _mm256_add_ps(dyPow2, dzPow2);
                __m256 temp3 = _mm256_add_ps(temp1, temp2);
                __m256 denominatorVector = _mm256_sqrt_ps(temp3);

                // finally RR1
                __m256 rr1Vector = _mm256_div_ps(oneVector, denominatorVector);

                // now ot the power of 3
                __m256 drPowerN32Vector = _mm256_mul_ps(rr1Vector, _mm256_mul_ps(rr1Vector, rr1Vector));

                // Calculate the net force
                FxVector = _mm256_add_ps(FxVector, _mm256_mul_ps(dxVector, drPowerN32Vector));
                FyVector = _mm256_add_ps(FyVector, _mm256_mul_ps(dyVector, drPowerN32Vector));
                FzVector = _mm256_add_ps(FzVector, _mm256_mul_ps(dzVector, drPowerN32Vector));
            }

            // step out of vectorization for velocity calculations. need to unvectorize Fx/y/z.
            // this works 100%.
            float *TempArray = (float*) &FxVector;
            Fx = TempArray[0] + TempArray[1] + TempArray[2] + TempArray[3]
                + TempArray[4] + TempArray[5] + TempArray[6] + TempArray[7];
            
            TempArray = (float*) &FyVector;
            Fy = TempArray[0] + TempArray[1] + TempArray[2] + TempArray[3]
                + TempArray[4] + TempArray[5] + TempArray[6] + TempArray[7];

            TempArray = (float*) &FzVector;
            Fz = TempArray[0] + TempArray[1] + TempArray[2] + TempArray[3]
                + TempArray[4] + TempArray[5] + TempArray[6] + TempArray[7];

            // Accelerate particles in response to the gravitational force
            global_Vx[i] += dt * Fx;
            global_Vy[i] += dt * Fy;
            global_Vz[i] += dt * Fz;
        }
    });

    // Move particles according to their velocities
    tbb::parallel_for(tbb::blocked_range<unsigned int>(0, nParticles, TBBGRAIN),
    [](const tbb::blocked_range<unsigned int> &range)
    {
        for (unsigned int i = range.begin(); i < range.end(); ++i)
        {
            global_X[i] += global_Vx[i] * dt;
            global_Y[i] += global_Vy[i] * dt;
            global_Z[i] += global_Vz[i] * dt;
        }
    });
}

// ========================== TBB + VECTORIZATION FUNCTIONS ===============================
/*
// for now only 12x speed up...
// barely vectorized here...
void get_particle_parallel_OLD(int i, OneParticle *p)
{
    // this only works when we add 2 trash floats to avoid overwriting.
    __m256 data1 = _mm256_loadu_ps(&particles[i].x);
    _mm256_storeu_ps(&p->x, data1);
}

// cant vectorize here...
void init_particles_parallel_OLD()
{
    tbb::parallel_for( // (range, lambda)
        tbb::blocked_range<unsigned int>(0, nParticles, TBBGRAIN),
        [](tbb::blocked_range<unsigned int> &range)
        {
            for (unsigned int i = range.begin(); i < range.end(); i++)
            {
                particles[i].x = (float)(i % 15);
                particles[i].y = (float)(((i * i) % 15));
                particles[i].z = (float)((i * i * 3) % 15);
                particles[i].vx = 0.0f;
                particles[i].vy = 0.0f;
                particles[i].vz = 0.0f;
            }
        });
}

// main optimization will come from here:
void move_particles_parallel_OLD(){
    // trying only with tbb for now :(
    // Compute forces using TBB
    tbb::parallel_for(tbb::blocked_range<unsigned int>(0, nParticles, TBBGRAIN),
        [](const tbb::blocked_range<unsigned int> &range)
        {
            for (unsigned int i = range.begin(); i < range.end(); ++i)
            {

                // Components of the gravity force on particle i
                float Fx = 0, Fy = 0, Fz = 0;
                __m256 FxVector = zeroVector;
                __m256 FyVector = zeroVector;
                __m256 FzVector = zeroVector;

                // creation of x, y, z vectors for calculations. because in the j iterations they stay the same.
                __m256 PixVector = _mm256_set1_ps(particles[i].x);
                __m256 PiyVector = _mm256_set1_ps(particles[i].y);
                __m256 PizVector = _mm256_set1_ps(particles[i].z);

                // i will vectorize in jumps of j+8's.
                for (unsigned int j = 0; j < nParticles; j += 8)
                {
                    // create the data arrays to store into the vectors.
                    float jxArray[8] = 
                    {
                        particles[j + 0].x, particles[j + 1].x, particles[j + 2].x, particles[j + 3].x,
                        particles[j + 4].x, particles[j + 5].x, particles[j + 6].x, particles[j + 7].x
                    };

                    float jyArray[8] = 
                    {
                        particles[j + 0].y, particles[j + 1].y, particles[j + 2].y, particles[j + 3].y,
                        particles[j + 4].y, particles[j + 5].y, particles[j + 6].y, particles[j + 7].y
                    };
                    
                    float jzArray[8] = 
                    {
                        particles[j + 0].z, particles[j + 1].z, particles[j + 2].z, particles[j + 3].z,
                        particles[j + 4].z, particles[j + 5].z, particles[j + 6].z, particles[j + 7].z
                    };
                    // creation of more x, y, z vectors for dx/y/z.
                    __m256 PjxVector = _mm256_loadu_ps(jxArray);
                    __m256 PjyVector = _mm256_loadu_ps(jyArray);
                    __m256 PjzVector = _mm256_loadu_ps(jzArray);

                    // Newton's law of universal gravity
                    __m256 dxVector = _mm256_sub_ps(PjxVector, PixVector);
                    __m256 dyVector = _mm256_sub_ps(PjyVector, PiyVector);
                    __m256 dzVector = _mm256_sub_ps(PjzVector, PizVector);

                    // a TON of calcs for R1 ^ 3
                    __m256 dxPow2 = _mm256_mul_ps(dxVector, dxVector);
                    __m256 dyPow2 = _mm256_mul_ps(dyVector, dyVector);
                    __m256 dzPow2 = _mm256_mul_ps(dzVector, dzVector);
                    // denominator calcs
                    __m256 temp1 = _mm256_add_ps(dxPow2, softVector);
                    __m256 temp2 = _mm256_add_ps(dyPow2, dzPow2);
                    __m256 temp3 = _mm256_add_ps(temp1, temp2);
                    __m256 denominatorVector = _mm256_sqrt_ps(temp3);

                    // finally RR1
                    __m256 rr1Vector = _mm256_div_ps(oneVector, denominatorVector);

                    // now ot the power of 3
                    __m256 drPowerN32Vector = _mm256_mul_ps(rr1Vector, _mm256_mul_ps(rr1Vector, rr1Vector));

                    // Calculate the net force
                    FxVector = _mm256_add_ps(FxVector, _mm256_mul_ps(dxVector, drPowerN32Vector));
                    FyVector = _mm256_add_ps(FyVector, _mm256_mul_ps(dyVector, drPowerN32Vector));
                    FzVector = _mm256_add_ps(FzVector, _mm256_mul_ps(dzVector, drPowerN32Vector));
                }

                // step out of vectorization for velocity calculations. need to unvectorize Fx/y/z.
                // this works 100%.
                float *TempArray = (float*) &FxVector;
                Fx = TempArray[0] + TempArray[1] + TempArray[2] + TempArray[3]
                   + TempArray[4] + TempArray[5] + TempArray[6] + TempArray[7];
                
                TempArray = (float*) &FyVector;
                Fy = TempArray[0] + TempArray[1] + TempArray[2] + TempArray[3]
                   + TempArray[4] + TempArray[5] + TempArray[6] + TempArray[7];

                TempArray = (float*) &FzVector;
                Fz = TempArray[0] + TempArray[1] + TempArray[2] + TempArray[3]
                   + TempArray[4] + TempArray[5] + TempArray[6] + TempArray[7];

                // Accelerate particles in response to the gravitational force
                particles[i].vx += dt * Fx;
                particles[i].vy += dt * Fy;
                particles[i].vz += dt * Fz;
            }
        });

    // Move particles according to their velocities
    tbb::parallel_for(tbb::blocked_range<unsigned int>(0, nParticles, TBBGRAIN),
        [](const tbb::blocked_range<unsigned int> &range)
        {
            for (unsigned int i = range.begin(); i < range.end(); ++i)
            {
                particles[i].x += particles[i].vx * dt;
                particles[i].y += particles[i].vy * dt;
                particles[i].z += particles[i].vz * dt;
            }
        });
}

// =================================== CODE GRAVEYARD =========================================
// =================== TBB ONLY FUNCTIONS ==========================
// only 3.9x speed up.
// TBB Only Optimizations. 1024 is the best grain size.
void get_particle_parallel_TBB(int i, OneParticle *p)
{
    p->x = particles[i].x;
    p->y = particles[i].y;
    p->z = particles[i].z;
    p->vx = particles[i].vx;
    p->vy = particles[i].vy;
    p->vz = particles[i].vz;
}
void init_particles_parallel_TBB()
{
    tbb::parallel_for( // (range, lambda)
        tbb::blocked_range<unsigned int>(0, nParticles, TBBGRAIN),
        [](tbb::blocked_range<unsigned int> &range)
        {
            // trying only with TBB...
            for (unsigned int i = range.begin(); i < range.end(); ++i)
            {
                particles[i].x = (float)(i % 15);
                particles[i].y = (float)(((i * i) % 15));
                particles[i].z = (float)((i * i * 3) % 15);
                particles[i].vx = 0.0f;
                particles[i].vy = 0.0f;
                particles[i].vz = 0.0f;
            }
        });
}
void move_particles_parallel_TBB()
{
    // trying only with tbb for now :(
    // Compute forces using TBB
    tbb::parallel_for(tbb::blocked_range<unsigned int>(0, nParticles, TBBGRAIN),
        [](const tbb::blocked_range<unsigned int> &range)
        {
            const float softening = 1e-20f;

            for (unsigned int i = range.begin(); i < range.end(); ++i)
            {

                // Components of the gravity force on particle i
                float Fx = 0, Fy = 0, Fz = 0;

                // Loop over particles that exert force: vectorization expected here
                for (unsigned int j = 0; j < nParticles; ++j)
                {

                    // Newton's law of universal gravity
                    float dx = particles[j].x - particles[i].x;
                    float dy = particles[j].y - particles[i].y;
                    float dz = particles[j].z - particles[i].z;

                    float rr1 = 1.0f / sqrt(dx * dx + dy * dy + dz * dz + softening);
                    float drPowerN32 = rr1 * rr1 * rr1;

                    // Calculate the net force
                    Fx += dx * drPowerN32;
                    Fy += dy * drPowerN32;
                    Fz += dz * drPowerN32;
                }

                // Accelerate particles in response to the gravitational force
                particles[i].vx += dt * Fx;
                particles[i].vy += dt * Fy;
                particles[i].vz += dt * Fz;
            }
        });

    // Move particles according to their velocities
    tbb::parallel_for(tbb::blocked_range<unsigned int>(0, nParticles, TBBGRAIN),
        [](const tbb::blocked_range<unsigned int> &range)
        {
            for (unsigned int i = range.begin(); i < range.end(); ++i)
            {
                particles[i].x += particles[i].vx * dt;
                particles[i].y += particles[i].vy * dt;
                particles[i].z += particles[i].vz * dt;
            }
        });
}

*/
#endif
