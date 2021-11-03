/*
 * timing.cpp
 *
 *  Created on: 31/lug/2010
 *      Author: Luca
 */

#include "timing.h"

#ifdef __MINGW_H

#include <windows.h>
NanoTimer::NanoTimer() {
    LARGE_INTEGER pff;
    QueryPerformanceFrequency(&pff);
    prec = pff.QuadPart;
};

double NanoTimer::time() {
    return (double)(t1 - t0) / (double)prec;
};

double NanoTimer::mark(bool t) {
    LARGE_INTEGER pfc;
    QueryPerformanceCounter(&pfc);
    if(t) t0 = pfc.QuadPart;
    else t1 = pfc.QuadPart;
    return time();
};
#endif
