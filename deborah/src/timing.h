/*
 * timing.h
 *
 *  Created on: 31/lug/2010
 *      Author: Luca
 */

#ifndef DEBORAH_TIMING_H_
#define DEBORAH_TIMING_H_

#include "types.h"
#include <time.h>

#define READTSC(buf) __asm__ __volatile__ ( \
  "rdtsc\n\t" \
  "movl %%eax,0(%0)\n\t" \
  "movl %%edx,4(%0)" \
        : : "r" (buf) : "eax", "edx" )


// Coarse-precision timer based on clock() API
class BasicTimer {
public:
        BasicTimer() { t0 = t1 = clock(); };
        double time() { return (double)(t1-t0)/(double)CLOCKS_PER_SEC; };
        double mark(bool t) {if(t) t0=clock(); else t1=clock(); return time();};
private:
        clock_t t0, t1;
};


#ifdef __MINGW_H

class NanoTimer {
public:
  NanoTimer();
  double time();
  double mark(bool t);
private:
  int64 t0, t1, prec;
};

#else

// Accurate timer based on clock_gettime() API
class NanoTimer {
public:
        NanoTimer() {
                struct timespec tres;
                ctype = CLOCK_MONOTONIC;
                clock_gettime(ctype, &t0);
                t1 = t0;
                clock_getres(ctype, &tres);
                prec = tres.tv_sec * 1000000000LL + tres.tv_nsec;
        };
        double time() {
                int64 dt;
                dt = t1.tv_sec * 1000000000LL + t1.tv_nsec - (t0.tv_sec * 1000000000LL + t0.tv_nsec);
                return ((double)(dt/prec))/1000000000.0;
        };
        double mark(bool t) {
                if(t) clock_gettime(ctype, &t0);
                else clock_gettime(ctype, &t1);
                return time();
        };
        int64 precision() {
                return prec;
        }
        int64 timestamp(bool which) {
                return which ? t0.tv_sec * 1000000000LL + t0.tv_nsec : t1.tv_sec * 1000000000LL + t1.tv_nsec;
        }
        void set_clockid(int t) {
                ctype = (clockid_t) t;
        }
private:
        struct timespec t0, t1;
        int64 prec;
        clockid_t ctype;
};
#endif


#endif /* TIMING_H_ */
