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

#define STRIDE 128
#define WIDTH 120
#define HEIGHT 120

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

// int *nodivergence;
int nodivergence = 1;

#define cellPerThreadX 2
#define cellPerThreadY 1
#define boundaryCellPerThreadX 16
#define boundaryCellPerThreadY 16

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
    SimpleBoundary(int x0, int x1, int y0, int y1, T InFlowU, T InFlowV)
        : x0(x0), x1(x1), y0(y0), y1(y1), InFlowU(InFlowU), InFlowV(InFlowV)
    {}

//    __host__ __device__ SimpleBoundary(SimpleBoundary const &sb)
//        : x0(sb.x0), x1(sb.x1), y0(sb.y0), y1(sb.y1), InFlowU(sb.InFlowU), InFlowV(sb.InFlowV)
//    {}

//    __host__ __device__ ~SimpleBoundary() {}

    bool isInFlowBoundary(int i, int j) const { return i == 0; }
    bool isOutFlowBoundary(int i, int j) const { return i == WIDTH; }
    bool isFloorBoundary(int i, int j) const { return j == 0; }
    bool isCeilingBoundary(int i, int j) const { return j == HEIGHT; }

    bool isObstacle(int i, int j) const {
        return x0 <= i && i <= x1 && y0 <= j && j <= y1;
    }

    T inflowU() const { return InFlowU; }
    T inflowV() const { return InFlowV; }
};

template <typename T, typename BoundaryCond>
void update_boundary(T *uc, T *vc, T *P, BoundaryCond bound)
{
    // printf("Thread (%d, %d), startX = %d, startY = %d\n", x, y, startX, startY);

    int i, j;
    for (j = 1; j < HEIGHT+1; ++j) {
        uc[INDEXU(0, j)] = bound.inflowU();
        vc[INDEXV(0, j)] = bound.inflowV();
        uc[INDEXU(WIDTH, j)] = uc[INDEXU(WIDTH-1, j)];
        vc[INDEXV(WIDTH, j)] = vc[INDEXV(WIDTH-1, j)];
        P[INDEXP(WIDTH, j)] = P[INDEXP(WIDTH-1,j)];
    }
    for (i = 1; i < WIDTH+1; ++i) {
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
void update_uv(T const *uc, T const *vc, T const *P, T *un, T *vn, T dt, T dx, T dy, BoundaryCond bound)
{
    int i, j;
    #pragma omp parallel for private(j)
    for (j = 1; j < HEIGHT; ++j) {
        for (i = 1; i < WIDTH; ++i) {
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

template <typename T, typename BoundaryCond>
void adjust_puv(T const *uc, T const *vc, T *P, T *un, T *vn,
		T dt, T dx, T dy, T beta,
                BoundaryCond bound, bool cellType)
{
    T D, delta_P;
//    int shift = (1 % 2) ^ cellType;
    int i, j;

#pragma omp parallel for private(j, D, delta_P) reduction(and:nodivergence)
    for (j = 1; j < HEIGHT; ++j) {
        int shift = (j % 2) ^ cellType;
        for (i = shift + 1; i < WIDTH; i += 2) {
            if (bound.isObstacle(i, j)) {
                //P[INDEXP(i, j)] = 0;
            } else {
                D = 1/dx * (un[INDEXU(i+1,j)] - un[INDEXU(i,j)])
                    +1/dy * (vn[INDEXV(i,j+1)] - vn[INDEXV(i,j)]);
                if (fabs(D) > Dtolerance) {
                    delta_P = -beta * D;
                    P[INDEXP(i,j)] += delta_P;
                    un[INDEXU(i,j)] -= (dt/dx)*delta_P;
                    un[INDEXU(i+1,j)] += (dt/dx)*delta_P;
                    vn[INDEXV(i,j)] -= (dt/dy)*delta_P;
                    vn[INDEXV(i,j+1)] += (dt/dy)*delta_P;
                    nodivergence = 0;
                }
            }
        }
    }

}

template <typename T, typename BoundaryCond>
void time_step(T *uc, T *vc, T *P, T *un, T *vn,
				T dt, T print_step, T dx, T dy, T beta, BoundaryCond bound)
{
    int steps = print_step / dt;
    while(steps--) {
        update_uv(uc, vc, P, un, vn, dt, dx, dy, bound);
        update_boundary(un, vn, P, bound);

//        int iteration = 0;
            //int nodivergence = 1;

        do {
//            printf("Iteration %d\r", ++iteration);

            nodivergence = 1;
                            //for (int i = 0; i < dimGrid.x * dimGrid.y; ++i) blk_nodivergence[i] = 1;

            adjust_puv(uc, vc, P, un, vn, dt, dx, dy, beta, bound, true);
            adjust_puv(uc, vc, P, un, vn, dt, dx, dy, beta, bound, false);
        } while (!nodivergence);

//        printf("\n");
        update_boundary(un, vn, P, bound);

        // swap (uc, un), (vc, vn)
        T *tmpc = uc; uc = un; un = tmpc;
        tmpc = vc; vc = vn; vn = tmpc;
    }
}

template <typename T, typename BoundaryCond>
void initialize(T* &ucurrent, T* &vcurrent, T* &unew, T* &vnew, T* &P, T* &huc, T* &hvc,
                BoundaryCond bound)
{
    ucurrent = (T*) malloc(STRIDE*STRIDE * sizeof(T));
    vcurrent = (T*) malloc(STRIDE*STRIDE * sizeof(T));
    unew = (T*) malloc(STRIDE*STRIDE * sizeof(T));
    vnew = (T*) malloc(STRIDE*STRIDE * sizeof(T));
    P = (T*) malloc(STRIDE*STRIDE * sizeof(T));
    // cudaMemset(ucurrent, 0, STRIDE*STRIDE * sizeof(T));
    // cudaMemset(vcurrent, 0, STRIDE*STRIDE * sizeof(T));
    memset(unew, 0, STRIDE*STRIDE * sizeof(T));
    memset(vnew, 0, STRIDE*STRIDE * sizeof(T));
    memset(P, 0, STRIDE*STRIDE * sizeof(T));

    //nodivergence = malloc(dimGrid.x * dimGrid.y * sizeof(*nodivergence));

    // inflow boundary
//    for (int j = 0; j <= HEIGHT; ++j) {
//        ucurrent[INDEXU(0, j)] = bound.inflowU();
//        vcurrent[INDEXU(0, j)] = bound.inflowV();
//    }
    // update_boundary<<<dimGrid2, dimBlock>>>(ucurrent, vcurrent, P, bound);
    // cudaDeviceSynchronize();

    //huc = (T*) std::malloc(STRIDE * STRIDE * sizeof(T));
    //hvc = (T*) std::malloc(STRIDE * STRIDE * sizeof(T));
    std::memset(ucurrent, 0, STRIDE * STRIDE * sizeof(T));
    std::memset(vcurrent, 0, STRIDE * STRIDE * sizeof(T));
    for (int j = 0; j <= HEIGHT; ++j) {
        ucurrent[j*STRIDE] = bound.inflowU();
        vcurrent[j*STRIDE] = bound.inflowV();
    }

    // cudaMemcpy(huc, ucurrent, STRIDE*STRIDE*sizeof(T), cudaMemcpyDeviceToHost);
    // cudaMemcpy(hvc, vcurrent, STRIDE*STRIDE*sizeof(T), cudaMemcpyDeviceToHost);
    //cudaMemcpy(ucurrent, huc, STRIDE*STRIDE*sizeof(T), cudaMemcpyHostToDevice);
    //cudaMemcpy(vcurrent, hvc, STRIDE*STRIDE*sizeof(T), cudaMemcpyHostToDevice);
    // cudaDeviceSynchronize();
}

template <typename T>
void freememory(T* ucurrent, T* vcurrent, T* unew, T* vnew, T* P, T* huc, T* hvc)
{
    free(ucurrent);
    free(vcurrent);
    free(unew);
    free(vnew);
    free(P);
    //free(nodivergence);
    //std::free(huc);
    //std::free(hvc);
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
    T *huc, *hvc;
    initialize(uc, vc, un, vn, P, huc, hvc, bound);
    T beta = beta0 / (2*dt*(1/(dx*dx) + 1/(dy*dy)));

    T accumulate_time = 0;
    T last_printed = 0;
    int index = 0;
    print_velocity(uc, vc, index++);
    while (accumulate_time < total_time) {
        accumulate_time += print_step;
        time_step(uc, vc, P, un, vn, dt, print_step, dx, dy, beta, bound);
        // cudaDeviceSynchronize();

        std::swap(uc, un);
        std::swap(vc, vn);

        // if (accumulate_time >= last_printed + print_step) {
            //cudaMemcpy(huc, uc, STRIDE*STRIDE*sizeof(T), cudaMemcpyDeviceToHost);
            //cudaMemcpy(hvc, vc, STRIDE*STRIDE*sizeof(T), cudaMemcpyDeviceToHost);
            // cudaDeviceSynchronize();
            last_printed += print_step;
            print_velocity(uc, vc, index++);
        // }
    }

    freememory(uc, vc, un, vn, P, huc, hvc);
}

int main()
{
    SimpleBoundary<float> sb ( OBSTACLE_MIN_X, OBSTACLE_MAX_X,
                               OBSTACLE_MIN_Y, OBSTACLE_MAX_Y,
                               1.f, 0.0f );
    float dx = .01, dy = .01, dt = .00001;
    float total_time = 900*dt, print_step = 100*dt;
    float beta0 = 1.7f;

    flow<float, SimpleBoundary<float> >(total_time, print_step, dt, dx, dy, beta0, sb);
}
