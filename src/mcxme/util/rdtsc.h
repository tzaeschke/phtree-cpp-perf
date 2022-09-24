/*
 * rdtsc.h
 * 
 * Created on: Jan 17, 2009
 *
 * Darius Sidlauskas (dariuss@cs.au.dk)
 * Copyright © 2008-2012 Aalborg University
 * Copyright © 2012-2013 Aarhus University
 *
 * Uses a dedicated ReaD Time-Stamp Counter (RDTSC) in order
 * to obtain very fine time measurements with least overhead.
 *
 */

#ifndef SRC_RDTSC_H_
#define SRC_RDTSC_H_

#include <stdlib.h>
#include <stdint.h>

/*
 * ReaD Time-Stamp Counter (RDTSC) is used for very fine
 * time measurements.
 */
extern "C" {
inline uint64_t RDTSC() {
    uint32_t lo, hi;
    __asm__ __volatile__ (  // serialize
            "xorl %%eax,%%eax \n        cpuid"
            ::: "%rax", "%rbx", "%rcx", "%rdx");
     // We cannot use "=A", since this would use %rax on x86_64
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return (uint64_t) hi << 32 | lo;
}
}

/*
 * In order to use the following methods, cycles per second (CPS)
 * have to be defined based on your CPU!
 */
#define SANDYBRIDGE_CPS 3400000000ULL  // 1x 4-core i7-3770 Intel (3.40GHz)

#ifdef CPS
  #define CYCLES_PER_SECOND CPS
#else
#define CYCLES_PER_SECOND SANDYBRIDGE_CPS
#endif

/*
 * Warning: depends on CPU frequency, see above
 */
inline double get_nanoseconds( uint64_t start, uint64_t end ) {
    return ( end - start ) / ( CYCLES_PER_SECOND / 1000000000.0 );
}

/*
 * Warning: depends on CPU frequency, see above
 */
inline double get_microseconds( uint64_t start, uint64_t end ) {
    return ( end - start ) / ( CYCLES_PER_SECOND / 1000000.0 );
}

/*
 * Warning: depends on CPU frequency, see above
 */
inline double get_milliseconds( uint64_t start, uint64_t end ) {
    return ( end - start ) / ( CYCLES_PER_SECOND / 1000.0 );
}

/*
 * Warning: depends on CPU frequency, see above
 */
inline double get_seconds( uint64_t start, uint64_t end ) {
    return ( end - start ) / (CYCLES_PER_SECOND / 1.0);
}

#endif  // SRC_RDTSC_H_
