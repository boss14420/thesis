/*
 * =====================================================================================
 *
 *       Filename:  hydrodynamic-flow.cc
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  03/18/2015 11:23:52 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Boss14420 (), firefox at gmail dot com
 *   Organization:
 *
 * =====================================================================================
 */

#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>
#include "common.hpp"

#define SWAP(x, y) (x ^= y ^= x ^= y);

#define STRIDE 1024
#define WIDTH 1023
#define HEIGHT 1023

#define OBSTACLE_MIN_X 25
#define OBSTACLE_MAX_X 60
#define OBSTACLE_MIN_Y 10
#define OBSTACLE_MAX_Y 40

//#define dx .01f
//#define dy .01f
//#define dt .01f
#define INDEXU(x, y) ((y) * (STRIDE) + (x))
#define INDEXV(x, y) ((y) * (STRIDE) + (x))
#define INDEXP(x, y) ((y)*STRIDE+ (x))
#define SQR(x) ((x) * (x))
#define nu 0.1f
#define Dtolerance 0.01f

#define PAD (256/sizeof(T) - 2)
//#define PAD 0

int *nodivergence;
//__device__ int nodivergence = 1;

#define cellPerThreadX 32
#define cellPerThreadY 32
#define boundaryCellPerThreadX 16
#define boundaryCellPerThreadY 16
#define threadPerBlockX 32
#define threadPerBlockY 8

__const__ dim3 dimBlock(threadPerBlockX, threadPerBlockY);
__const__ dim3 dimGrid( (((WIDTH + dimBlock.x - 1)/dimBlock.x) + cellPerThreadX - 1)/cellPerThreadX,
              (((HEIGHT+ dimBlock.y - 1)/dimBlock.y) + cellPerThreadY - 1)/cellPerThreadY
            );

__const__ dim3 dimGrid2( (((WIDTH + dimBlock.x - 1)/dimBlock.x) + boundaryCellPerThreadX - 1)/boundaryCellPerThreadX,
               (((HEIGHT+ dimBlock.y - 1)/dimBlock.y) + boundaryCellPerThreadY - 1)/boundaryCellPerThreadY
             );


/*
   P[INDEXP(i, j)] <-- P[i, j]
   uc[INDEXU(i, j)] <-- u[i-1/2, j]
   vc[INDEXV(i, j)] <-- v[i, j-1/2]
 */

template <class T>
struct SimpleBoundary
{
private:
public:
    int x0, x1, y0, y1;
    T InFlowU, InFlowV;

public:
//    __global__ SimpleBoundary(SimpleBoundary const &sb);
    __host__ __device__ SimpleBoundary(int x0, int x1, int y0, int y1, T InFlowU, T InFlowV)
        : x0(x0), x1(x1), y0(y0), y1(y1), InFlowU(InFlowU), InFlowV(InFlowV)
    {}

//    __host__ __device__ SimpleBoundary(SimpleBoundary const &sb)
//        : x0(sb.x0), x1(sb.x1), y0(sb.y0), y1(sb.y1), InFlowU(sb.InFlowU), InFlowV(sb.InFlowV)
//    {}

//    __host__ __device__ ~SimpleBoundary() {}

    __device__ bool isInFlowBoundary(int i, int j) const { return i == 0; }
    __device__ bool isOutFlowBoundary(int i, int j) const { return i == WIDTH; }
    __device__ bool isFloorBoundary(int i, int j) const { return j == 0; }
    __device__ bool isCeilingBoundary(int i, int j) const { return j == HEIGHT; }

    __device__ bool isObstacle(int i, int j) const {
        return x0 <= i && i <= x1 && y0 <= j && j <= y1;
    }

    __device__ T inflowU() const { return InFlowU; }
    __device__ T inflowV() const { return InFlowV; }
};

template <typename T, typename BoundaryCond>
__global__
void update_boundary(T *uc, T *vc, T *P, BoundaryCond bound)
{
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int x = blockIdx.x * blockDim.x + threadIdx.x;

    int startX = x*boundaryCellPerThreadX + 1;
    int endX = min(startX + boundaryCellPerThreadX, WIDTH);
    int startY = y*boundaryCellPerThreadY + 1;
    int endY = min(startY + boundaryCellPerThreadY, HEIGHT);
    // printf("Thread (%d, %d), startX = %d, startY = %d\n", x, y, startX, startY);

    int i, j;
    for (j = startY; j < endY; ++j) {
        uc[INDEXU(0, j)] = bound.inflowU();
        vc[INDEXV(0, j)] = bound.inflowV();
        uc[INDEXU(WIDTH, j)] = uc[INDEXU(WIDTH-1, j)];
        vc[INDEXV(WIDTH, j)] = vc[INDEXV(WIDTH-1, j)];
        P[INDEXP(WIDTH, j)] = P[INDEXP(WIDTH-1,j)];
    }
    for (i = startX; i < endX; ++i) {
        uc[INDEXU(i, 0)] = uc[INDEXU(i, HEIGHT-1)];
        vc[INDEXV(i, 0)] = vc[INDEXU(i, HEIGHT-1)];
        P[INDEXP(i, 0)] = P[INDEXP(i, HEIGHT-1)];
//        uc[INDEXU(i, 0)] = 0;//uc[INDEXU(i, 1)];
//        vc[INDEXV(i, 0)] = 0;//vc[INDEXU(i, 1)];

        uc[INDEXU(i, HEIGHT)] = uc[INDEXU(i, 1)];
        vc[INDEXV(i, HEIGHT)] = vc[INDEXV(i, 1)];
        P[INDEXV(i, HEIGHT)] = P[INDEXV(i, 1)];
//        uc[INDEXU(i, HEIGHT)] = 0;// uc[INDEXU(i, HEIGHT-1)];
//        vc[INDEXV(i, HEIGHT)] = 0;// vc[INDEXV(i, HEIGHT-1)];
    }
}

template <typename T, typename BoundaryCond>
__global__
void update_uv(T const *uc, T const *vc, T const *P, T *un, T *vn, T dt, T dx, T dy, BoundaryCond bound)
{
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int x = blockIdx.x * blockDim.x + threadIdx.x;

    int startX = x*cellPerThreadX + 1;
    int endX = min(startX + cellPerThreadX, WIDTH);
    int startY = y*cellPerThreadY + 1;
    int endY = min(startY + cellPerThreadY, HEIGHT);

    int i, j;
    for (j = startY; j < endY; ++j) {
        for (i = startX; i < endX; ++i) {
            if (bound.isObstacle(i, j)) {
                // zero velocity
                un[INDEXU(i, j)] = 0;
                vn[INDEXU(i, j)] = 0;
            } else {
                // TODO: use shared memory
                T Du, Dv;
                T uim12j    = uc[INDEXU(i-1,j)];        // u(i-1/2, j)
                T uip12j    = uc[INDEXU(i, j)];         // u(i+1/2, j)
                T uip32j    = uc[INDEXU(i+1, j)];       // u(i+3/2, j)
                T uip12jp1  = uc[INDEXU(i, j+1)];       // u(i+1/2, j+1)
                T uip12jm1  = uc[INDEXU(i, j-1)];       // u(i+1/2, j-1)
                T uim12jp1  = uc[INDEXU(i-1, j+1)];     // u(i-1/2, j+1)

                T vijp12    = vc[INDEXV(i, j)];         // v(i, j+1/2)
                T vijp32    = vc[INDEXP(i, j+1)];       // v(i, j+3/2);
                T vijm12    = vc[INDEXV(i, j-1)];       // v(i, j-1/2)
                T vip1jp12  = vc[INDEXV(i+1, j)];       // v(i+1, j+1/2)
                T vip1jm12  = vc[INDEXV(i+1, j-1)];     // v(i+1, j-1/2)
                T vim1jp12  = vc[INDEXV(i-1, j)];       // v(i-1, j+1/2);


                // Du
                T uij       = .5 * (uip12j + uim12j);
                T uip1j     = .5 * (uip32j + uip12j);

                Du = -1/dx * (uip1j*uip1j - uij*uij);

                // uv(i+1/2, j+1/2)
                T uvip12jp12 = .5*(uip12j + uip12jp1) * .5*(vip1jp12 + vijp12);
                // u(i+1/2, j+1/2)
                T uvip12jm12 = .5*(uip12jm1 + uip12j) * .5*(vip1jm12 + vijm12);

                Du += -1/dy * (uvip12jp12 - uvip12jm12);
                Du += -1/dx*(P[INDEXP(i,j)] - P[INDEXP(i-1,j)]);
                Du += nu*( (uip32j - 2*uip12j + uim12j) / (dx*dx)
                          +(uip12jp1 - 2*uip12j + uip12jm1) / (dy*dy) );

//                DFu(uc, vc, P, dx, dy, i, j, Du);
                un[INDEXU(i, j)] = uc[INDEXU(i, j)] + dt * Du;


                // Dv
                T vij       = .5 * (vijp12 + vijm12);       // v(i, j)
                T vijp1     = .5 * (vijp32 + vijp12);       // v(i, j+1)

                Dv = -1/dy * (vijp1*vijp1 - vij*vij);

                // uv(i-1/2, j+1/2)
                T uvim12jp12 = .5*(uim12j + uim12jp1) * .5*(vim1jp12 + vijp12);
                Dv += -1/dx * (uvip12jm12 - uvim12jp12);
                Dv += -1/dy * (P[INDEXP(i, j)] - P[INDEXP(i, j-1)]);
                Dv += nu*( (vijp32 - 2*vijp12 + vijm12) / (dy*dy)
                          +(vip1jp12 - 2*vijp12 + vim1jp12) / (dx*dx) );

//                DFv(uc, vc, P, dx, dy, i, j, Dv);
                vn[INDEXV(i, j)] = vc[INDEXV(i, j)] + dt * Dv;
            }
        }
    }
}

template <typename T> struct double_vec;

template <> struct double_vec<float> { typedef float2 type; };
template <> struct double_vec<double> { typedef double2 type; };

template <typename T, typename BoundaryCond>
__global__
void adjust_puv(T const *uc, T const *vc, T *P, T *un, T *vn, 
				T dt, T dx, T dy, T beta, 
                int *nodivergence,
                BoundaryCond bound, bool cellType)
{
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int x = blockIdx.x * blockDim.x + threadIdx.x;

    int shift = (y % 2) ^ cellType;

    int startX = x*cellPerThreadX + 1;
    int endX = min(startX + cellPerThreadX, WIDTH);
    int startY = y*cellPerThreadY + 1;
    int endY = min(startY + cellPerThreadY, HEIGHT);

    T D, delta_P;
    int thread_nodivergence = 1;
    //typename double_vec<T>::type u12; // coalesced access u
    T *vij, *vijp1;
    T *uij, *uip1j;

    for (int j = startY; j < endY; ++j, shift = 1-shift) {
        for (int i = startX + shift; i < endX; i += 2) {
    //int i = startX + shift;
    //int j = startY;
            if (bound.isObstacle(i, j)) {
                //P[INDEXP(i, j)] = 0;
            } else {
                uij     = un + INDEXU(i, j);
                uip1j   = un + INDEXU(i+1, j);
                //u12     = *((float2*) uij); // u12 = { u(i, j), u(i+1, j) }
                vij     = vn + INDEXV(i, j);
                vijp1   = vn + INDEXV(i, j+1);
                //D = 1/dx * (un[INDEXU(i+1,j)] - un[INDEXU(i,j)])
                //    +1/dy * (vn[INDEXV(i,j+1)] - vn[INDEXV(i,j)]);
                //D = 1/dx * (u12.y - u12.x) + 1/dy * (*vijp1 - *vij);
                D = 1/dx * (*uip1j - *uij) + 1/dy * (*vijp1 - *vij);
                //if (fabs(D) > Dtolerance) {
                    delta_P = -beta * D;
                    P[INDEXP(i,j)] += delta_P;
                    //un[INDEXU(i,j)] -= (dt/dx)*delta_P;
                    *uij            -= (dt/dx)*delta_P;
                    //un[INDEXU(i+1,j)] += (dt/dx)*delta_P;
                    *uip1j          += (dt/dx)*delta_P;
                    *vij            -= (dt/dy)*delta_P;
                    *vijp1          += (dt/dy)*delta_P;
                    //vn[INDEXV(i,j)] -= (dt/dy)*delta_P;
                    //vn[INDEXV(i,j+1)] += (dt/dy)*delta_P;
                    //thread_nodivergence = 0;
                    thread_nodivergence &= (fabs(D) <= Dtolerance);
                //}
            }
        }
    }

    //int warp_nodivergence = __all(thread_nodivergence);
    // first thread in a warp
    //if ( (threadIdx.y * blockDim.x + threadIdx.x) % warpSize == 0) {

    int bn = __syncthreads_and(thread_nodivergence);
    // first thread in a block
    if ( (threadIdx.y == 0) && (threadIdx.x == 0) ) {
        //atomicAnd(&nodivergence, bn);
        atomicAnd(nodivergence, bn);
    }
}

template <typename T, typename BoundaryCond>
//__global__
void time_step(T *uc, T *vc, T *P, T *un, T *vn, 
				T dt, T print_step, T dx, T dy, T beta, BoundaryCond bound)
{
    /*
    dim3 dimBlock(16, 16);
    dim3 dimGrid( (((WIDTH + dimBlock.x - 1)/dimBlock.x) + cellPerThreadX - 1)/cellPerThreadX,
            (((HEIGHT+ dimBlock.y - 1)/dimBlock.y) + cellPerThreadY - 1)/cellPerThreadY
            );

    dim3 dimGrid2( (((WIDTH + dimBlock.x - 1)/dimBlock.x) + boundaryCellPerThreadX - 1)/boundaryCellPerThreadX,
               (((HEIGHT+ dimBlock.y - 1)/dimBlock.y) + boundaryCellPerThreadY - 1)/boundaryCellPerThreadY
             );
    */

    int steps = print_step / dt;
    while(steps--) {
        update_uv<<<dimGrid, dimBlock>>>(uc, vc, P, un, vn, dt, dx, dy, bound);
        cudaDeviceSynchronize();
        update_boundary<<<dimGrid2, dimBlock>>>(un, vn, P, bound);
        cudaDeviceSynchronize();

        //int iteration = 0;
        int hnodivergence = 1;

        do {
            // printf("Iteration %d\r", ++iteration);

            //hnodivergence = 1;
            //*nodivergence = 1;
            cudaMemsetAsync(nodivergence, 1, sizeof(int), 0);

            adjust_puv<<<dimGrid, dimBlock>>>(uc, vc, P, un, vn, dt, dx, dy, beta, nodivergence, bound, true);
            cudaDeviceSynchronize();
            adjust_puv<<<dimGrid, dimBlock>>>(uc, vc, P, un, vn, dt, dx, dy, beta, nodivergence, bound, false);
            cudaDeviceSynchronize();

            cudaMemcpy(&hnodivergence, nodivergence, sizeof(int), cudaMemcpyDeviceToHost);
        } while (!hnodivergence);

        // printf("\n");
        update_boundary<<<dimGrid2, dimBlock>>>(un, vn, P, bound);
        cudaDeviceSynchronize();

        // swap (uc, un), (vc, vn)
        T *tmpc = uc; uc = un; un = tmpc;
        tmpc = vc; vc = vn; vn = tmpc;
    }
}

template <typename T, typename BoundaryCond>
void initialize(T* &ucurrent, T* &vcurrent, T* &unew, T* &vnew, T* &P, T* &huc, T* &hvc,
                int* &nodivergence,
                BoundaryCond bound)
{
    // aglinment on u[1], v[1]
    cudaMalloc(&ucurrent, (PAD + STRIDE*STRIDE) * sizeof(T));
    cudaMalloc(&vcurrent, (PAD + STRIDE*STRIDE) * sizeof(T));
    cudaMalloc(&unew, (PAD + STRIDE*STRIDE) * sizeof(T));
    cudaMalloc(&vnew, (PAD + STRIDE*STRIDE) * sizeof(T));
    cudaMalloc(&P, (PAD + STRIDE*STRIDE) * sizeof(T));
    ucurrent += PAD;
    vcurrent += PAD;
    unew += PAD;
    vnew += PAD;
    P += PAD;
    // cudaMemset(ucurrent, 0, STRIDE*STRIDE * sizeof(T));
    // cudaMemset(vcurrent, 0, STRIDE*STRIDE * sizeof(T));
    cudaMemsetAsync(unew, 0, STRIDE*STRIDE * sizeof(T), 0);
    cudaMemsetAsync(vnew, 0, STRIDE*STRIDE * sizeof(T), 0);
    cudaMemsetAsync(P, 0, STRIDE*STRIDE * sizeof(T), 0);

    cudaMalloc(&nodivergence, sizeof(*nodivergence));

    // inflow boundary
//    for (int j = 0; j <= HEIGHT; ++j) {
//        ucurrent[INDEXU(0, j)] = bound.inflowU();
//        vcurrent[INDEXU(0, j)] = bound.inflowV();
//    }
    // update_boundary<<<dimGrid2, dimBlock>>>(ucurrent, vcurrent, P, bound);
    // cudaDeviceSynchronize();

    huc = (T*) std::malloc((PAD + STRIDE * STRIDE) * sizeof(T));
    hvc = (T*) std::malloc((PAD + STRIDE * STRIDE) * sizeof(T));
    huc += PAD;
    hvc += PAD;
    std::memset(huc, 0, STRIDE * STRIDE * sizeof(T));
    std::memset(hvc, 0, STRIDE * STRIDE * sizeof(T));
    for (int j = 0; j <= HEIGHT; ++j) {
        huc[j*STRIDE] = 1;
        hvc[j*STRIDE] = 0;
    }

    // cudaMemcpy(huc, ucurrent, STRIDE*STRIDE*sizeof(T), cudaMemcpyDeviceToHost);
    // cudaMemcpy(hvc, vcurrent, STRIDE*STRIDE*sizeof(T), cudaMemcpyDeviceToHost);
    cudaMemcpyAsync(ucurrent, huc, STRIDE*STRIDE*sizeof(T), cudaMemcpyHostToDevice, 0);
    cudaMemcpyAsync(vcurrent, hvc, STRIDE*STRIDE*sizeof(T), cudaMemcpyHostToDevice, 0);
    // cudaDeviceSynchronize();
}

template <typename T>
void freememory(T* ucurrent, T* vcurrent, T* unew, T* vnew, T* P, T* huc, T* hvc, int *nodivergence)
{
    cudaFree(ucurrent - PAD);
    cudaFree(vcurrent - PAD);
    cudaFree(unew - PAD);
    cudaFree(vnew - PAD);
    cudaFree(P - PAD);
    cudaFree(nodivergence);
    std::free(huc - PAD);
    std::free(hvc - PAD);
}

template <typename T>
void print_velocity(T const *uc, T const *vc, int index)
{
    char name[20];
    std::printf("Step %d\n", index);

    T const *u = uc, *v = vc;
//    T *u = new T[STRIDE*STRIDE];
//    T *v = new T[STRIDE*STRIDE];
//    for (int j = 0; j < HEIGHT; ++j)
//        for (int i = 0; i < WIDTH; ++i) {
//            u[INDEXU(i, j)] = .5 * (uc[INDEXU(i, j)] + uc[INDEXU(i+1, j)]);
//            v[INDEXU(i, j)] = .5 * (vc[INDEXU(i, j)] + vc[INDEXU(i, j+1)]);
//        }

//    std::cout << "u\n";
    std::snprintf(name, 19, "uc%04d.pgm", index);
    exportPixmap(u, WIDTH, HEIGHT, STRIDE, name);
    std::snprintf(name, 19, "uc%04d.mat", index);
    exportMatlab(u, WIDTH, HEIGHT, STRIDE, name);
    //exportValue(uc);

//    std::cout << "\nv\n";
    std::snprintf(name, 19, "vc%04d.pgm", index);
    exportPixmap(v, WIDTH, HEIGHT, STRIDE, name);
    std::snprintf(name, 19, "vc%04d.mat", index);
    exportMatlab(v, WIDTH, HEIGHT, STRIDE, name);
    //exportValue(vc);
    //std::cout << "\n\n";

//    delete[] u;
//    delete[] v;
}

template <typename T, typename BoundaryCond>
void flow(T total_time, T print_step, T dt, T dx, T dy, T beta0, BoundaryCond &bound)
{
    T *uc, *vc, *un, *vn, *P;
    //int *blk_nodivergence;
    T *huc, *hvc;
    initialize(uc, vc, un, vn, P, huc, hvc, nodivergence, bound);
    T beta = beta0 / (2*dt*(1/(dx*dx) + 1/(dy*dy)));

    T accumulate_time = 0;
    T last_printed = 0;
    int index = 0;
    print_velocity(huc, hvc, index++);
    while (accumulate_time < total_time) {
        accumulate_time += print_step;
        //time_step<<<1,1>>>(uc, vc, P, un, vn, dt, print_step, dx, dy, beta, bound);
        time_step(uc, vc, P, un, vn, dt, print_step, dx, dy, beta, bound);
        cudaDeviceSynchronize();

        std::swap(uc, un);
        std::swap(vc, vn);

        // if (accumulate_time >= last_printed + print_step) {
            cudaMemcpyAsync(huc, uc, STRIDE*STRIDE*sizeof(T), cudaMemcpyDeviceToHost, 0);
            cudaMemcpyAsync(hvc, vc, STRIDE*STRIDE*sizeof(T), cudaMemcpyDeviceToHost, 0);
            // cudaDeviceSynchronize();
            last_printed += print_step;
            print_velocity(huc, hvc, index++);
        // }
    }

    freememory(uc, vc, un, vn, P, huc, hvc, nodivergence);
}

int main()
{
    typedef float T;

/////// calculate occupancy
    int numBlocks;
    int potentialBlockSize;
    int minGridSize;//, gridSize;

    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &potentialBlockSize, adjust_puv<T, SimpleBoundary<T> >, 0, WIDTH*HEIGHT);
    std::cout << "minGridSize: " << minGridSize
                << ", potentialblockSize: " << potentialBlockSize << '\n';

    int device;
    cudaDeviceProp prop;
    int activeWarps;
    int maxWarps;

    cudaGetDevice(&device);
    cudaGetDeviceProperties(&prop, device);

    int blockSize = dimBlock.x * dimBlock.y;
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &numBlocks,
            adjust_puv<T, SimpleBoundary<T> >,
            blockSize,
            0);

    activeWarps = numBlocks * blockSize/ prop.warpSize;
    maxWarps = prop.maxThreadsPerMultiProcessor / prop.warpSize;

    std::cout << "activeWarps: " << activeWarps
                << ", maxWarps: " << maxWarps << '\n';
    std::cout << "Occupancy: " << (double)activeWarps/maxWarps<< '\n';
    std::cout << "No of blocks: " << dimGrid.x * dimGrid.y << '\n';

/////////////////////////////////////////////////////////////////////////////////////
    SimpleBoundary<T> sb ( OBSTACLE_MIN_X, OBSTACLE_MAX_X,
                               OBSTACLE_MIN_Y, OBSTACLE_MAX_Y,
                               1.f, 0.0f );
    T dx = .01, dy = .01, dt = .0001;
    T total_time = 900*dt, print_step = 100*dt;
    T beta0 = 1.7f;

    flow<T, SimpleBoundary<T> >(total_time, print_step, dt, dx, dy, beta0, sb);
}
