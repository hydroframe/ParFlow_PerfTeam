/* Generic wrapper header file for managing different parallel implementations */

#ifndef _PF_PARALLEL_H
#define _PF_PARALLEL_H

#define PlusEquals(a, b) a += b
#define AMPS_RANK_0 !amps_Rank(amps_CommWorld)
#define CU_CALL(x)
#define CU_NORET(x)

#ifdef HAVE_OMP

/* Utility macros for inserting OMP pragmas in macros */
#define EMPTY()
#define DEFER(x) x EMPTY()
#define PRAGMA(args) _Pragma( #args )

#include "pf_omploops.h"

#endif // HAVE_OMP

#ifdef HAVE_CUDA

#include "pf_cuda.h"

#undef CU_CALL
#define CU_CALL( x ) if(AMPS_RANK_0) { return x; }
#undef CU_NORET
#define CU_NORET( x ) if(AMPS_RANK_0) { x; return; }
#endif // HAVE_CUDA

#endif // _PF_PARALLEL_H
