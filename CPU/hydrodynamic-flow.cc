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
#include "common.hpp"

#define SWAP(x, y) (x ^= y ^= x ^= y);

#define STRIDE 256
#define WIDTH 200
#define HEIGHT 200
//#define dx .01f
//#define dy .01f
//#define dt .01f
#define INDEXU(x, y) ((y) * (STRIDE) + (x))
#define INDEXV(x, y) ((y) * (STRIDE) + (x))
#define INDEXP(x, y) ((y)*STRIDE+ (x))
#define SQR(x) ((x) * (x))
#define nu 0.1f
#define Dtolerance 0.01f


/*
   P[INDEXP(i, j)] <-- P[i, j]
   uc[INDEXU(i, j)] <-- u[i-1/2, j]
   vc[INDEXV(i, j)] <-- v[i, j-1/2]
 */

template <class T>
class SimpleBoundary
{
private:
public:
    int x0, x1, y0, y1;
    T InFlowU, InFlowV;

public:
    bool isInFlowBoundary(int i, int j) const { return i == 0; }
    bool isOutFlowBoundary(int i, int j) const { return i == WIDTH; }
    bool isFloorBoundary(int i, int j) const { return j == 0; }
    bool isCeilingBoundary(int i, int j) const { return j == HEIGHT; }

    bool isObstace(int i, int j) const {
        return x0 <= i && i <= x1 && y0 <= j && j <= y1;
    }

    T inflowU() const { return InFlowU; }
    T inflowV() const { return InFlowV; }
};


// finite differencing of u
template <typename T>
void DFu(T const *uc, T const *vc, T const *P, T dx, T dy, int i, int j, T &Du)
{
    Du = -1/dx*.25f*(SQR(uc[INDEXU(i+1,j)] + uc[INDEXU(i,j)]) - SQR(uc[INDEXU(i-1,j)] + uc[INDEXU(i,j)]));
    Du += -1/dy*.25f*( (uc[INDEXU(i,j+1)]+uc[INDEXU(i,j)])*(vc[INDEXV(i+1,j)]+vc[INDEXV(i,j)])
                      -(uc[INDEXU(i,j)]+uc[INDEXU(i,j-1)])*(vc[INDEXV(i+1,j-1)]+vc[INDEXV(i,j-1)]) );
    Du += -1/dx*(P[INDEXP(i,j)] - P[INDEXP(i-1,j)]);
    Du += nu*( (uc[INDEXU(i+1,j)] - 2*uc[INDEXU(i,j)] + uc[INDEXU(i-1,j)])/(dx*dx)
                +(uc[INDEXU(i,j+1)] - 2*uc[INDEXU(i,j)] + uc[INDEXU(i,j-1)])/(dy*dy) );
}

// finite differencing of v (without dP/dy)
template <typename T>
void DFv(T const *uc, T const *vc, T const *P, T dx, T dy, int i, int j, T &Dv)
{
    Dv = -1/dy*.25f*(SQR(vc[INDEXV(i,j+1)] + vc[INDEXV(i,j)]) - SQR(vc[INDEXV(i,j-1)] + vc[INDEXV(i,j)]))
          -1/dy*.25f*( (uc[INDEXU(i,j+1)]+uc[INDEXU(i,j)])*(vc[INDEXV(i+1,j)]+vc[INDEXV(i,j)])
                      -(uc[INDEXU(i-1,j)]+uc[INDEXU(i-1,j+1)])*(vc[INDEXV(i,j)]+vc[INDEXV(i-1,j)]) )
          -1/dy*(P[INDEXP(i,j)] - P[INDEXP(i,j-1)])
          + nu*( (vc[INDEXV(i,j+1)] - 2*vc[INDEXV(i,j)] + vc[INDEXV(i,j-1)])/(dx*dx)
                +(vc[INDEXV(i+1,j)] - 2*vc[INDEXV(i,j)] + vc[INDEXV(i-1,j)])/(dy*dy) );
}

template <typename T, typename BoundaryCond>
void update_boundary(T *uc, T *vc, T *P, BoundaryCond bound)
{
    int i, j;
    for (j = 0; j <= HEIGHT; ++j) {
        uc[INDEXU(0, j)] = bound.inflowU();
        vc[INDEXV(0, j)] = bound.inflowV();
        uc[INDEXU(WIDTH, j)] = uc[INDEXU(WIDTH-1, j)];
        vc[INDEXV(WIDTH, j)] = vc[INDEXV(WIDTH-1, j)];
        P[INDEXP(WIDTH, j)] = P[INDEXP(WIDTH-1,j)];
    }
    for (i = 1; i <= WIDTH-1; ++i) {
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
    for (j = 1; j < HEIGHT; ++j) {
        for (i = 1; i < WIDTH; ++i) {
            if (bound.isObstace(i, j)) {
                // zero velocity
                un[INDEXU(i, j)] = 0;
                vn[INDEXU(i, j)] = 0;
            } else {
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
void time_step(T const *uc, T const *vc, T *P, T *un, T *vn, T dt, T dx, T dy, T beta, BoundaryCond bound)
{
    update_uv(uc, vc, P, un, vn, dt, dx, dy, bound);
    update_boundary(un, vn, P, bound);

    T D, delta_P;
    bool incompressible = false;
    int iteration = 0;
    do {
        incompressible = true;
        std::cout << "Iteration " << ++iteration << '\r';
        for (int j = 1; j < HEIGHT; ++j) {
            for (int i = 1; i < WIDTH; ++i) {
                if (bound.isObstace(i, j)) {
                    //P[INDEXP(i, j)] = 0;
                } else {
                    D = 1/dx * (un[INDEXU(i+1,j)] - un[INDEXU(i,j)])
                        +1/dy * (vn[INDEXV(i,j+1)] - vn[INDEXV(i,j)]);
                    if (std::fabs(D) > Dtolerance) {
//                        if (D > 1.5) {
//                            std::cout << i << ", " << j << '\n';
//                            std::exit(-1);
//                        }
                        delta_P = -beta * D;
                        P[INDEXP(i,j)] += delta_P;
                        un[INDEXU(i,j)] -= (dt/dx)*delta_P;
                        un[INDEXU(i+1,j)] += (dt/dx)*delta_P;
                        vn[INDEXV(i,j)] -= (dt/dy)*delta_P;
                        vn[INDEXV(i,j+1)] += (dt/dy)*delta_P;
                        incompressible = false;
                    }
                }
            }
        }
    } while (!incompressible);
    update_boundary(un, vn, P, bound);
    std::cout << '\n';
}

template <typename T, typename BoundaryCond>
void initialize(T* &ucurrent, T* &vcurrent, T* &unew, T* &vnew, T* &P, BoundaryCond bound)
{
    ucurrent = (T*) std::calloc(STRIDE*STRIDE, sizeof(T));
    vcurrent = (T*) std::calloc(STRIDE*STRIDE, sizeof(T));
    unew = (T*) std::calloc(STRIDE*STRIDE, sizeof(T));
    vnew = (T*) std::calloc(STRIDE*STRIDE, sizeof(T));
    P = (T*) std::calloc(STRIDE*STRIDE, sizeof(T));

    // inflow boundary
    for (int j = 0; j <= HEIGHT; ++j) {
        ucurrent[INDEXU(0, j)] = bound.inflowU();
        vcurrent[INDEXU(0, j)] = bound.inflowV();
    }
}

template <typename T>
void print_velocity(T const *uc, T const *vc, int index)
{
    char name[20];
    std::cout << "Step " << index << '\n';

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
void flow(T total_time, T print_step, T dt, T dx, T dy, T beta0, BoundaryCond bound)
{
    T *uc, *vc, *un, *vn, *P;
    initialize(uc, vc, un, vn, P, bound);
    T beta = beta0 / (2*dt*(1/(dx*dx) + 1/(dy*dy)));

    T accumulate_time = 0;
    T last_printed = 0;
    int index = 0;
    print_velocity(uc, vc, index++);
    while (accumulate_time < total_time) {
        accumulate_time += dt;
        time_step(uc, vc, P, un, vn, dt, dx, dy, beta, bound);

        std::swap(uc, un);
        std::swap(vc, vn);
        if (accumulate_time >= last_printed + print_step) {
            last_printed += print_step;
            print_velocity(uc, vc, index++);
        }
    }
}

int main()
{
    SimpleBoundary<float> sb { 25, 50, 0, 50, 1.f, 0.0f };
    float dx = .01, dy = .01, dt = .00001;
    float total_time = 4000*dt, print_step = 1000*dt;
    float beta0 = 1.7f;

    flow<double, decltype(sb)>(total_time, print_step, dt, dx, dy, beta0, sb);
}
