/**************************************************/
/*
/*              DON'T TOUCH
/*
/**************************************************/
#ifndef NBODY_H
// #define NBODY_H
#include <cmath>
#include <tbb/tbb.h> // for TBB
#include <immintrin.h> // for AVX
using namespace std;
#define VECLEN 8
#define TBBGRAIN 1024 // 4 (SIZE OF FLOAT) * 6 * 1024 =~ 24KB


typedef struct
{
    float x;
    float y;
    float z;
    float vx;
    float vy;
    float vz;
} OneParticle;

const int nParticles = 16384 * 2;
const float dt = 0.01f;

struct ParticleType
{
    float x, y, z;
    float vx, vy, vz;
};

ParticleType particles[nParticles];


// Helper function to horizontally sum the elements of an AVX vector
float horizontal_sum(__m256 v)
{
    __m128 vlow = _mm256_castps256_ps128(v);
    __m128 vhigh = _mm256_extractf128_ps(v, 1);
    vlow = _mm_add_ps(vlow, vhigh);
    __m128 shuf = _mm_movehdup_ps(vlow);
    __m128 sums = _mm_add_ps(vlow, shuf);
    shuf = _mm_movehl_ps(shuf, sums);
    sums = _mm_add_ss(sums, shuf);
    return _mm_cvtss_f32(sums);
}



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

// ========================== TBB ONLY FUNCTIONS ==========================================

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

// ========================== TBB + VECTORIZATION FUNCTIONS ===============================

void get_particle_parallel(int i, OneParticle *p)
{
    // Load first 4 floats (x, y, z, vx) using 128-bit load
    __m128 data1 = _mm_loadu_ps(&particles[i].x);

    // Store first 4 floats
    _mm_storeu_ps(&p->x, data1);
    
    // Store last 2 floats
    p->vy = particles[i].vy;
    p->vz = particles[i].vz;
}


void init_particles_parallel()
{
    tbb::parallel_for( // (range, lambda)
        tbb::blocked_range<unsigned int>(0, nParticles, TBBGRAIN),
        [](tbb::blocked_range<unsigned int> &range)
        {

            // 8 floats in parallel is 256 bits. we'll work on 8 at a time.
            /*for (unsigned int i = range.begin(); i < range.end(); i = i + 8)
            {

                //particles[i].x = (float)(i % 15); //Convertion:
                __m256 particlesXVector = _mm256_set_ps(particles[i].x,     particles[i + 1].x, particles[i + 2].x, particles[i + 3].x,
                                                        particles[i + 4].x, particles[i + 5].x, particles[i + 6].x, particles[i + 7].x);

                __m256 indexesVector = _mm256_set_ps(i,     i + 1, i + 2, i + 3, 
                                                     i + 4, i + 5, i + 6, i + 7);
                
                _mm256_store_ps(particles[i].x);
            }*/

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

void move_particles_parallel(){
        tbb::parallel_for(tbb::blocked_range<unsigned int>(0, nParticles, TBBGRAIN),
        [](const tbb::blocked_range<unsigned int> &range)
        {
            __m256 softeningVector = _mm256_set1_ps(1e-20f);
            __m256 zeroVector = _mm256_setzero_ps();

            for (unsigned int i = range.begin(); i < range.end(); i++)
            {

                // Components of the gravity force on particle i
                //float Fx = 0, Fy = 0, Fz = 0;
                __m256 FxVector = zeroVector;
                __m256 FyVector = zeroVector;
                __m256 FzVector = zeroVector;


                // Loop over particles that exert force: vectorization expected here
                for (unsigned int j = 0; j < nParticles; j++)
                {

                    // Newton's law of universal gravity
                    //float dx = particles[j].x - particles[i].x;
                    //float dy = particles[j].y - particles[i].y;
                    //float dz = particles[j].z - particles[i].z;

                    // J vectors.
                    __m256 particlesJXVector = _mm256_set1_ps(particles[j].x);
                    __m256 particlesJYVector = _mm256_set1_ps(particles[j].y);
                    __m256 particlesJZVector = _mm256_set1_ps(particles[j].z);
                    
                    // I vectors.
                    __m256 particlesIXVector = _mm256_set1_ps(particles[i].x);
                    __m256 particlesIYVector = _mm256_set1_ps(particles[i].y);
                    __m256 particlesIZVector = _mm256_set1_ps(particles[i].z);

                    // dx, dy, dz vectors.
                    __m256 dxVector = _mm256_sub_ps(particlesJXVector, particlesIXVector);
                    __m256 dyVector = _mm256_sub_ps(particlesJYVector, particlesIYVector);
                    __m256 dzVector = _mm256_sub_ps(particlesJZVector, particlesIZVector);
                    

                    //float rr1 = 1.0f / sqrt(dx * dx + dy * dy + dz * dz + softening);
                    __m256 dxVectorp2 = _mm256_mul_ps(dxVector, dxVector);
                    __m256 dyVectorp2 = _mm256_mul_ps(dyVector, dyVector);
                    __m256 dzVectorp2 = _mm256_mul_ps(dzVector, dzVector);

                    __m256 pythagorasVector = _mm256_add_ps(_mm256_add_ps(dxVectorp2, softeningVector), _mm256_add_ps(dyVectorp2, dzVectorp2));

                    __m256 sqrtVector = _mm256_sqrt_ps(pythagorasVector);
                    __m256 rr1Vector = _mm256_div_ps(_mm256_set1_ps(1.0f), sqrtVector);

                    // float drPowerN32 = rr1 * rr1 * rr1;
                    __m256 drPowerN32Vector = _mm256_mul_ps(rr1Vector, _mm256_mul_ps(rr1Vector, rr1Vector));
                    

                    // Calculate the net force
                    //Fx += dx * drPowerN32;
                    //Fy += dy * drPowerN32;
                    //Fz += dz * drPowerN32;
                    FxVector = _mm256_add_ps(FxVector, _mm256_mul_ps(dxVector, drPowerN32Vector));
                    FyVector = _mm256_add_ps(FyVector, _mm256_mul_ps(dyVector, drPowerN32Vector));
                    FzVector = _mm256_add_ps(FzVector, _mm256_mul_ps(dzVector, drPowerN32Vector));
                }

                // Horizontal sum of the vectors
                __m256 sumFx = _mm256_hadd_ps(FxVector, FxVector);
                sumFx = _mm256_hadd_ps(sumFx, sumFx);
                float Fx = _mm256_cvtss_f32(sumFx);

                __m256 sumFy = _mm256_hadd_ps(FyVector, FyVector);
                sumFy = _mm256_hadd_ps(sumFy, sumFy);
                float Fy = _mm256_cvtss_f32(sumFy);

                __m256 sumFz = _mm256_hadd_ps(FzVector, FzVector);
                sumFz = _mm256_hadd_ps(sumFz, sumFz);
                float Fz = _mm256_cvtss_f32(sumFz);


                // Accelerate particles in response to the gravitational force
                //particles[i].vx += dt * Fx;
                //particles[i].vy += dt * Fy;
                //particles[i].vz += dt * Fz;
                __m256 particlesVXVector = _mm256_set1_ps(particles[i].vx);
                __m256 particlesVYVector = _mm256_set1_ps(particles[i].vy);
                __m256 particlesVZVector = _mm256_set1_ps(particles[i].vz);
                __m256 dtVector = _mm256_set1_ps(dt);

                particlesVXVector = _mm256_add_ps(particlesVXVector, _mm256_mul_ps(dtVector, _mm256_set1_ps(Fx)));
                particlesVYVector = _mm256_add_ps(particlesVYVector, _mm256_mul_ps(dtVector, _mm256_set1_ps(Fy)));
                particlesVZVector =_mm256_add_ps(particlesVZVector, _mm256_mul_ps(dtVector, _mm256_set1_ps(Fz)));
            }
        });

   // Vectorized movement calculation [works]
    tbb::parallel_for(tbb::blocked_range<unsigned int>(0, nParticles, TBBGRAIN),
        [](const tbb::blocked_range<unsigned int> &range) {
            __m256 vDt = _mm256_set1_ps(dt);

            for (unsigned int i = range.begin(); i < range.end(); i++) {
                __m256 vX = _mm256_set1_ps(particles[i].x);
                __m256 vY = _mm256_set1_ps(particles[i].y);
                __m256 vZ = _mm256_set1_ps(particles[i].z);

                __m256 vVx = _mm256_set1_ps(particles[i].vx);
                __m256 vVy = _mm256_set1_ps(particles[i].vy);
                __m256 vVz = _mm256_set1_ps(particles[i].vz);

                vX = _mm256_add_ps(vX, _mm256_mul_ps(vVx, vDt));
                vY = _mm256_add_ps(vY, _mm256_mul_ps(vVy, vDt));
                vZ = _mm256_add_ps(vZ, _mm256_mul_ps(vVz, vDt));

                particles[i].x = _mm256_cvtss_f32(vX); // Store first element
                particles[i].y = _mm256_cvtss_f32(vY);
                particles[i].z = _mm256_cvtss_f32(vZ);
            }
        });
}



//=================================== CODE GRAVEYARD ===================================

//void get_particle_parallel(int i, OneParticle *p){
    // Load particle data into an AVX register
    // __m256 particleData = _mm256_loadu_ps((float*)&particles[i]);
    // Store the loaded data back into the OneParticle struct
    //_mm256_storeu_ps((float*)p, particleData);
//}

// INIT_PARTICLE BUT TBB AND VECTORIZED
/*void init_particles_parallel()
{
    tbb::parallel_for( // (range, lambda)
          tbb::blocked_range<unsigned int>(0, nParticles, TBBGRAIN),
         [](tbb::blocked_range<unsigned int> &range) {
         
           
            
            // Velocities set up
            const __m256 zeroVector = _mm256_setzero_ps();
            for (unsigned int i = range.begin(); i < range.end(); i += 8)
            {
                _mm256_storeu_ps(&particles[i].vx, zeroVector);
                _mm256_storeu_ps(&particles[i].vy, zeroVector);
                _mm256_storeu_ps(&particles[i].vz, zeroVector);
            }

            // Vectorized initialization for x, y, z
            for (unsigned int i = range.begin(); i < range.end(); i += 8)
            {
                // Load the indexes into an array
                __m256 indices = _mm256_set_ps(i + 7, i + 6, i + 5, i + 4, i + 3, i + 2, i + 1, i);

                // Calculate x values
                //__m256 xValues = _mm256_fmadd_ps(indices, _mm256_set1_ps(1.0f / 15.0f), zeroVector);
                __m256 xValues = _mm256_add_ps(_mm256_mul_ps(indices, _mm256_set1_ps(1.0f / 15.0f)), zeroVector); // x = i % 15
                xValues = _mm256_round_ps(xValues, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);

                // Calculate y values
                //__m256 yValues = _mm256_fmadd_ps(indices, indices, zeroVector);
                // __m256 yValues = _mm256_add_ps(_mm256_mul_ps(indices, indices, zeroVector)); 
                __m256 yValues = _mm256_add_ps(_mm256_mul_ps(indices, indices), zeroVector); // y = i * i

                yValues = _mm256_div_ps(yValues, _mm256_set1_ps(15.0f)); // y = (i * i) % 15
                yValues = _mm256_round_ps(yValues, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);

                // Calculate z values
                __m256 zValues = _mm256_mul_ps(yValues, _mm256_set1_ps(3.0f)); // z = i * i * 3
                zValues = _mm256_div_ps(zValues, _mm256_set1_ps(15.0f)); // z = (i * i * 3) % 15
                zValues = _mm256_round_ps(zValues, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);

                // Store the calculated x, y, z values
                _mm256_storeu_ps(&particles[i].x, xValues);
                _mm256_storeu_ps(&particles[i].y, yValues);
                _mm256_storeu_ps(&particles[i].z, zValues);
            }
        });
}*/

// parallel version
/*void move_particles_parallel(){
    // Loop over particles that experience force in parallel
    tbb::parallel_for(0, nParticles, 1, [&](int i) {
        // Components of the gravity force on particle i
        __m256 Fx = _mm256_setzero_ps();
        __m256 Fy = _mm256_setzero_ps();
        __m256 Fz = _mm256_setzero_ps();

        // Loop over particles that exert force
        for (int j = 0; j < nParticles; j += 8) {
            // Load position of particle i
            __m256 pos_i_x = _mm256_set1_ps(particles[i].x);
            __m256 pos_i_y = _mm256_set1_ps(particles[i].y);
            __m256 pos_i_z = _mm256_set1_ps(particles[i].z);

            // Load positions of up to 8 particles
            __m256 pos_j_x = _mm256_loadu_ps(&particles[j].x);
            __m256 pos_j_y = _mm256_loadu_ps(&particles[j].y);
            __m256 pos_j_z = _mm256_loadu_ps(&particles[j].z);

            // Compute deltas
            __m256 dx = _mm256_sub_ps(pos_j_x, pos_i_x);
            __m256 dy = _mm256_sub_ps(pos_j_y, pos_i_y);
            __m256 dz = _mm256_sub_ps(pos_j_z, pos_i_z);

            // Compute squared distances with softening
            const float softening = 1e-20f;
            __m256 dist_sq = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(dx, dx), _mm256_mul_ps(dy, dy)), 
                                            _mm256_add_ps(_mm256_mul_ps(dz, dz), _mm256_set1_ps(softening)));
            __m256 rr1 = _mm256_rsqrt_ps(dist_sq); // 1/sqrt(dist_sq)
            __m256 drPowerN32 = _mm256_mul_ps(_mm256_mul_ps(rr1, rr1), rr1); // rr1^3

            // Calculate the net force components 
            // Fx = _mm256_fmadd_ps(dx, drPowerN32, Fx);
            // Fy = _mm256_fmadd_ps(dy, drPowerN32, Fy);
            // Fz = _mm256_fmadd_ps(dz, drPowerN32, Fz);

            // [fmadd not supported]
            Fx = _mm256_add_ps(_mm256_mul_ps(dx, drPowerN32), Fx);
            Fy = _mm256_add_ps(_mm256_mul_ps(dy, drPowerN32), Fy);
            Fz = _mm256_add_ps(_mm256_mul_ps(dz, drPowerN32), Fz);

        }

        // Reduce the forces (if using more than 8 particles, sum the results)
        float force_x = 0.0f, force_y = 0.0f, force_z = 0.0f;

        // Store the computed forces into a temporary array
        float force[8];
        _mm256_storeu_ps(force, Fx);
        for (int k = 0; k < 8; ++k) {
            force_x += force[k];
        }
        Fx = _mm256_set1_ps(force_x);

        _mm256_storeu_ps(force, Fy);
        for (int k = 0; k < 8; ++k) {
            force_y += force[k];
        }
        Fy = _mm256_set1_ps(force_y);

        _mm256_storeu_ps(force, Fz);
        for (int k = 0; k < 8; ++k) {
            force_z += force[k];
        }
        Fz = _mm256_set1_ps(force_z);


        // Accelerate particles in response to the gravitational force
        //__m256 dt_vec = _mm256_set1_ps(dt);
        // Extract the sum of forces as a scalar
        float force_x_sum = _mm256_hadd_ps(Fx, Fx)[0];
        float force_y_sum = _mm256_hadd_ps(Fy, Fy)[0];
        float force_z_sum = _mm256_hadd_ps(Fz, Fz)[0];


        particles[i].vx += force_x_sum * dt;
        particles[i].vy += force_y_sum * dt;
        particles[i].vz += force_z_sum * dt;

    });

    // Move particles according to their velocities (serial loop)
    for (int i = 0; i < nParticles; i++) {
        particles[i].x += particles[i].vx * dt;
        particles[i].y += particles[i].vy * dt;
        particles[i].z += particles[i].vz * dt;
    }
}*/

#endif
