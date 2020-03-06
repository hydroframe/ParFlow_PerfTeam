/* Macro redefinitions for OMP enabled loops */

#ifndef _PF_OMPLOOPS_H
#define _PF_OMPLOOPS_H

#ifndef HAVE_OMP

#define MASTER(body) body
#define SINGLE(body, ...) body
#define SINGLE_REGION(body, ...) body
#define BARRIER
#define BEGIN_REGION
#define END_REGION

#else

#include <omp.h>
#include <stdarg.h>

#define NO_OMP_PARALLEL()																								\
	if (omp_get_active_level() != 0)																			\
	{																																			\
		fprintf(stderr,																											\
						"(%d) Error: Hit Parallel region in %s when not allowed\n", \
						amps_Rank(amps_CommWorld), __FUNCTION__);										\
		exit(-1);																														\
	}

/* Utility macros for inserting OMP pragmas in macros */
#define EMPTY()
#define DEFER(x) x EMPTY()
#define DEFER2(x) x EMPTY EMPTY() ()

#undef PRAGMA
#define PRAGMA(args) _Pragma( #args )

#undef PlusEquals
#define PlusEquals(a, b) OMP_PlusEquals(&(a), b)

#undef LOCALS
#define LOCALS(...) DEFER(_LOCALS)(__VA_ARGS__)
#define _LOCALS(...) ,__VA_ARGS__

#undef NO_LOCALS
#define NO_LOCALS DEFER(_NO_LOCALS)
#define _NO_LOCALS

#define MASTER(body) PRAGMA(omp master)         \
  {                                             \
    body;                                       \
  }

#define SINGLE(body, ...) PRAGMA(omp single __VA_ARGS__)  \
  {                                                       \
    body;                                                 \
  }

/* @MCB:
   By default, OpenMP nested parallel regions are disabled, and
   will only spawn one thread in the inner regions.  Because the
   regions are dynamically scoped, this means the "active" region
   will only have a thread count of one.  This ensures that any
   calls within the body region that encounter omp pragmas such as
   loop distributions, reductions, or barriers only apply as a
   single threaded context.

   A common reason for this is the InitVector calls.  If there are
   sections of code that need to be executed serially within a parallel
   region, such as in clustering.c, but have InitVector calls within
   said sections, then the region context needs to understand it only has
   one thread executing in order to handle the for loops properly.
*/
#define SINGLE_REGION(body, ...)                \
  PRAGMA(omp single __VA_ARGS__)                \
  {                                             \
    PRAGMA(omp parallel)                        \
    {                                           \
      body;                                     \
    }                                           \
  }

// @MCB: Note, the trailing {} is to deal with potential semicolons
#define BARRIER PRAGMA(omp barrier) {}

#define BEGIN_REGION PRAGMA(omp parallel) {
#define END_REGION }

/*
  Keeps the BoxLoop macros much tidier
  For GrGeomInLoop macros, the macro name will be used, not its expanded value
  i.e. In BoxLoopI# NoWait will expand to "for nowait", but in GrGeomInLoop it will
  be concatinated with _GrGeomInLoop resulting in NoWait_GrGeomInLoop
  This is necessary because going over the box clusters in a NewParallel section
  wraps the entire body of the macro inside the region for better performance,
  and so simply adjusting the pragma on the for loops is not sufficient.

  @MCB: TODO: Refer to the BCLoop implementations.
  The GrGeom versions can be converted similarly for better reuse.
*/
#define NewParallel parallel for
#define InParallel for
#define NoWait for nowait

/* @MCB: Include other macro files AFTER we redefine the LOCALS pragmas so it doesn't break */
#include "pf_omp_grgeom.h"
#include "pf_omp_bcloops.h"

/* Helper Functions */
extern "C++"{

#pragma omp declare simd
  template<typename T>
  inline void OMP_PlusEquals(T *arr_loc, T value)
  {
#pragma omp atomic update
    *arr_loc += value;
  }

#pragma omp declare simd
  template<typename T>
  inline T _HarmonicMean(T a, T b)
  {
    return (a + b) ? (2.0 * a * b) / (a + b) : 0;
  }

#pragma omp declare simd
  template<typename T>
  inline T _HarmonicMeanDZ(T a, T b, T c, T d)
  {
    return ((c * b) + (a * d)) ?  (((c + d) * a * b) / ((b * c) + (a * d))) : 0;
  }

#pragma omp declare simd
  template<typename T>
  inline T _UpstreamMean(T a, T b, T c, T d)
  {
    return (a - b) >= 0 ? c : d;
  }

#pragma omp delcare simd
  template<typename T>
  inline T _ArithmeticMean(T a, T b)
  {
    return (0.5 * (a + b));
  }

#undef AtomicSet
#define AtomicSet(a, b) OMP_AtomicSet(&(a), b)
  template<typename T>
  inline void OMP_AtomicSet(T *var, T val)
  {
#pragma omp atomic write
    *var = val;
  }

}

inline
void omp_Printf(const char *fmt, ...)
{
#pragma omp master
  {
    va_list argp;
    va_start(argp, fmt);
    vfprintf(stderr, fmt, argp);
    va_end(argp);
  }
}

inline
void omp_AllPrintf(const char *fmt, ...)
{
#pragma omp critical
  {
    va_list argp;
    va_start(argp, fmt);
    vfprintf(stderr, fmt, argp);
    va_end(argp);
  }

#pragma omp barrier
}

inline
void omp_AnyPrintf(const char *fmt, ...)
{
#pragma omp single nowait
  {
    va_list argp;
    va_start(argp, fmt);
    vfprintf(stderr, fmt, argp);
    va_end(argp);
  }
}


/*------------------------------------------------------------------------------
  BoxLoop Macro Redefinitions
  NOTE: The private pragmas will look odd.  This is because the NO_LOCALS macro
  expands to nothing, and the LOCALS(...) macro expands to ,__VA_ARGS__
  and so will insert the necessary comma itself.
  TODO: Come up with something that doesn't look like a crazy person wrote it
  -----------------------------------------------------------------------------*/

/* Util function to calculate the desired step in BoxLoopIX macros.
   Not tested with strides other than 1, but the math seems right.
   Pragma is to enable SIMD function calls from SIMD loops */
#pragma omp declare simd uniform(idx, nx, ny, sx, jinc, kinc)
inline int
INC_IDX(int idx, int i, int j, int k,
        int nx, int ny, int sx,
        int jinc, int kinc)
{
  return (k * kinc + (k * ny + j) * jinc +
          (k * ny * nx + j * nx + i) * sx) + idx;
}

#if 1

#undef InitMatrixLoop
#define InitMatrixLoop(s, stencil, ptr, mat_sub,                        \
                       i, j, k,                                         \
                       ix, iy, iz, nx, ny, nz,                          \
                       i1, nx1, ny1, nz1, sx1, sy1, sz1,                \
                       body)                                            \
  {                                                                     \
    for (s = 0; s < StencilSize(stencil); s++)                          \
    {                                                                   \
      ptr = SubmatrixElt(mat_sub, s, ix, iy, iz);                       \
      i1 = 0;                                                           \
      DeclareInc(PV_jinc_1, PV_kinc_1, nx, ny, nz, nx1, ny1, nz1, sx1, sy1, sz1); \
      PRAGMA(omp for collapse(2) private(i, j, k, i1))                  \
        for (k = iz; k < iz + nz; k++)                                  \
        {                                                               \
          for (j = iy; j < iy + ny; j++)                                \
          {                                                             \
            PRAGMA(omp simd)                                            \
              for (i = ix; i < ix + nx; i++)                            \
              {                                                         \
                i1 = INC_IDX((i - ix), (j - iy), (k - iz),              \
                             nx, ny, sx1, PV_jinc_1, PV_kinc_1);        \
                body;                                                   \
              }                                                         \
          }                                                             \
        }                                                               \
    }                                                                   \
  }

#endif

#undef _BoxLoopI0
#define _BoxLoopI0(pragma, locals,                          \
                   i, j, k,                                 \
                   ix, iy, iz,                              \
                   nx, ny, nz,                              \
                   body)                                    \
  {                                                         \
    PRAGMA(omp pragma collapse(3) private(i, j, k locals))  \
      for (k = iz; k < iz + nz; k++)                        \
      {                                                     \
        for (j = iy; j < iy + ny; j++)                      \
        {                                                   \
          for (i = ix; i < ix + nx; i++)                    \
          {                                                 \
            body;                                           \
          }                                                 \
        }                                                   \
      }                                                     \
  }

#undef _BoxLoopI1
#define _BoxLoopI1(pragma, locals,                                      \
                   i, j, k,                                             \
                   ix, iy, iz, nx, ny, nz,                              \
                   i1, nx1, ny1, nz1, sx1, sy1, sz1,                    \
                   body)                                                \
  {                                                                     \
    int i1_start = i1;                                                  \
    DeclareInc(PV_jinc_1, PV_kinc_1, nx, ny, nz, nx1, ny1, nz1, sx1, sy1, sz1); \
    PRAGMA(omp pragma collapse(3) private(i, j, k, i1 locals))          \
      for (k = iz; k < iz + nz; k++)                                    \
      {                                                                 \
        for (j = iy; j < iy + ny; j++)                                  \
        {                                                               \
          for (i = ix; i < ix + nx; i++)                                \
          {                                                             \
            i1 = INC_IDX(i1_start, (i - ix), (j - iy), (k - iz),        \
                         nx, ny, sx1, PV_jinc_1, PV_kinc_1);            \
            body;                                                       \
          }                                                             \
        }                                                               \
      }                                                                 \
  }

#undef _BoxLoopI2
#define _BoxLoopI2(pragma, locals,                                      \
                   i, j, k,                                             \
                   ix, iy, iz, nx, ny, nz,                              \
                   i1, nx1, ny1, nz1, sx1, sy1, sz1,                    \
                   i2, nx2, ny2, nz2, sx2, sy2, sz2,                    \
                   body)                                                \
  {                                                                     \
    int i1_start = i1;                                                  \
    int i2_start = i2;                                                  \
    DeclareInc(PV_jinc_1, PV_kinc_1, nx, ny, nz, nx1, ny1, nz1, sx1, sy1, sz1); \
    DeclareInc(PV_jinc_2, PV_kinc_2, nx, ny, nz, nx2, ny2, nz2, sx2, sy2, sz2); \
    PRAGMA(omp pragma collapse(3) private(i, j, k, i1, i2 locals))      \
      for (k = iz; k < iz + nz; k++)                                    \
      {                                                                 \
        for (j = iy; j < iy + ny; j++)                                  \
        {                                                               \
          for (i = ix; i < ix + nx; i++)                                \
          {                                                             \
            i1 = INC_IDX(i1_start, (i - ix), (j - iy), (k - iz),        \
                         nx, ny, sx1, PV_jinc_1, PV_kinc_1);            \
            i2 = INC_IDX(i2_start, (i - ix), (j - iy), (k - iz),        \
                         nx, ny, sx2, PV_jinc_2, PV_kinc_2);            \
            body;                                                       \
          }                                                             \
        }                                                               \
      }                                                                 \
  }

#undef _BoxLoopI3
#define _BoxLoopI3(pragma, locals,                                      \
                   i, j, k,                                             \
                   ix, iy, iz, nx, ny, nz,                              \
                   i1, nx1, ny1, nz1, sx1, sy1, sz1,                    \
                   i2, nx2, ny2, nz2, sx2, sy2, sz2,                    \
                   i3, nx3, ny3, nz3, sx3, sy3, sz3,                    \
                   body)                                                \
  {                                                                     \
    int i1_start = i1;                                                  \
    int i2_start = i2;                                                  \
    int i3_start = i3;                                                  \
    DeclareInc(PV_jinc_1, PV_kinc_1, nx, ny, nz, nx1, ny1, nz1, sx1, sy1, sz1); \
    DeclareInc(PV_jinc_2, PV_kinc_2, nx, ny, nz, nx2, ny2, nz2, sx2, sy2, sz2); \
    DeclareInc(PV_jinc_3, PV_kinc_3, nx, ny, nz, nx3, ny3, nz3, sx3, sy3, sz3); \
    PRAGMA(omp pragma collapse(3) private(i, j, k, i1, i2, i3 locals))  \
      for (k = iz; k < iz + nz; k++)                                    \
      {                                                                 \
        for (j = iy; j < iy + ny; j++)                                  \
        {                                                               \
          for (i = ix; i < ix + nx; i++)                                \
          {                                                             \
            i1 = INC_IDX(i1_start, (i - ix), (j - iy), (k - iz),        \
                         nx, ny, sx1, PV_jinc_1, PV_kinc_1);            \
            i2 = INC_IDX(i2_start, (i - ix), (j - iy), (k - iz),        \
                         nx, ny, sx2, PV_jinc_2, PV_kinc_2);            \
            i3 = INC_IDX(i3_start, (i - ix), (j - iy), (k - iz),        \
                         nx, ny, sx3, PV_jinc_3, PV_kinc_3);            \
            body;                                                       \
          }                                                             \
        }                                                               \
      }                                                                 \
  }

/**************************************************************************
 * SIMD Variants
 **************************************************************************/
#define SIMD_BoxLoopI0(pragma, locals,                      \
                       i, j, k,                             \
                       ix, iy, iz,                          \
                       nx, ny, nz,                          \
                       body)                                \
  {                                                         \
    PRAGMA(omp pragma collapse(2) private(i, j, k locals))  \
      for (k = iz; k < iz + nz; k++)                        \
      {                                                     \
        for (j = iy; j < iy + ny; j++)                      \
        {                                                   \
          PRAGMA(omp simd)                                  \
            for (i = ix; i < ix + nx; i++)                  \
            {                                               \
              body;                                         \
            }                                               \
        }                                                   \
      }                                                     \
  }

#define SIMD_BoxLoopI1(pragma, locals,                                  \
                       i, j, k,                                         \
                       ix, iy, iz, nx, ny, nz,                          \
                       i1, nx1, ny1, nz1, sx1, sy1, sz1,                \
                       body)                                            \
  {                                                                     \
    int i1_start = i1;                                                  \
    DeclareInc(PV_jinc_1, PV_kinc_1, nx, ny, nz, nx1, ny1, nz1, sx1, sy1, sz1); \
    PRAGMA(omp pragma collapse(3) private(i, j, k, i1 locals))          \
      for (k = iz; k < iz + nz; k++)                                    \
      {                                                                 \
        for (j = iy; j < iy + ny; j++)                                  \
        {                                                               \
          for (i = ix; i < ix + nx; i++)                                \
          {                                                             \
            i1 = INC_IDX(i1_start, (i - ix), (j - iy), (k - iz),        \
                         nx, ny, sx1, PV_jinc_1, PV_kinc_1);            \
            body;                                                       \
          }                                                             \
        }                                                               \
      }                                                                 \
  }

#define SIMD_BoxLoopI2(pragma, locals,                                  \
                       i, j, k,                                         \
                       ix, iy, iz, nx, ny, nz,                          \
                       i1, nx1, ny1, nz1, sx1, sy1, sz1,                \
                       i2, nx2, ny2, nz2, sx2, sy2, sz2,                \
                       body)                                            \
  {                                                                     \
    int i1_start = i1;                                                  \
    int i2_start = i2;                                                  \
    DeclareInc(PV_jinc_1, PV_kinc_1, nx, ny, nz, nx1, ny1, nz1, sx1, sy1, sz1); \
    DeclareInc(PV_jinc_2, PV_kinc_2, nx, ny, nz, nx2, ny2, nz2, sx2, sy2, sz2); \
    PRAGMA(omp pragma collapse(3) private(i, j, k, i1, i2 locals))      \
      for (k = iz; k < iz + nz; k++)                                    \
      {                                                                 \
        for (j = iy; j < iy + ny; j++)                                  \
        {                                                               \
          for (i = ix; i < ix + nx; i++)                                \
          {                                                             \
            i1 = INC_IDX(i1_start, (i - ix), (j - iy), (k - iz),        \
                         nx, ny, sx1, PV_jinc_1, PV_kinc_1);            \
            i2 = INC_IDX(i2_start, (i - ix), (j - iy), (k - iz),        \
                         nx, ny, sx2, PV_jinc_2, PV_kinc_2);            \
            body;                                                       \
          }                                                             \
        }                                                               \
      }                                                                 \
  }

#define SIMD_BoxLoopI2(pragma, locals,                                  \
                       i, j, k,                                         \
                       ix, iy, iz, nx, ny, nz,                          \
                       i1, nx1, ny1, nz1, sx1, sy1, sz1,                \
                       i2, nx2, ny2, nz2, sx2, sy2, sz2,                \
                       body)                                            \
  {                                                                     \
    int i1_start = i1;                                                  \
    int i2_start = i2;                                                  \
    DeclareInc(PV_jinc_1, PV_kinc_1, nx, ny, nz, nx1, ny1, nz1, sx1, sy1, sz1); \
    DeclareInc(PV_jinc_2, PV_kinc_2, nx, ny, nz, nx2, ny2, nz2, sx2, sy2, sz2); \
    PRAGMA(omp pragma collapse(3) private(i, j, k, i1, i2 locals))      \
      for (k = iz; k < iz + nz; k++)                                    \
      {                                                                 \
        for (j = iy; j < iy + ny; j++)                                  \
        {                                                               \
          for (i = ix; i < ix + nx; i++)                                \
          {                                                             \
            i1 = INC_IDX(i1_start, (i - ix), (j - iy), (k - iz),        \
                         nx, ny, sx1, PV_jinc_1, PV_kinc_1);            \
            i2 = INC_IDX(i2_start, (i - ix), (j - iy), (k - iz),        \
                         nx, ny, sx2, PV_jinc_2, PV_kinc_2);            \
            body;                                                       \
          }                                                             \
        }                                                               \
      }                                                                 \
  }

#define SIMD_BoxLoopI3(pragma, locals,                                  \
                       i, j, k,                                         \
                       ix, iy, iz, nx, ny, nz,                          \
                       i1, nx1, ny1, nz1, sx1, sy1, sz1,                \
                       i2, nx2, ny2, nz2, sx2, sy2, sz2,                \
                       i3, nx3, ny3, nz3, sx3, sy3, sz3,                \
                       body)                                            \
  {                                                                     \
    int i1_start = i1;                                                  \
    int i2_start = i2;                                                  \
    int i3_start = i3;                                                  \
    DeclareInc(PV_jinc_1, PV_kinc_1, nx, ny, nz, nx1, ny1, nz1, sx1, sy1, sz1); \
    DeclareInc(PV_jinc_2, PV_kinc_2, nx, ny, nz, nx2, ny2, nz2, sx2, sy2, sz2); \
    DeclareInc(PV_jinc_3, PV_kinc_3, nx, ny, nz, nx3, ny3, nz3, sx3, sy3, sz3); \
    PRAGMA(omp pragma collapse(3) private(i, j, k, i1, i2, i3 locals))  \
      for (k = iz; k < iz + nz; k++)                                    \
      {                                                                 \
        for (j = iy; j < iy + ny; j++)                                  \
        {                                                               \
          for (i = ix; i < ix + nx; i++)                                \
          {                                                             \
            i1 = INC_IDX(i1_start, (i - ix), (j - iy), (k - iz),        \
                         nx, ny, sx1, PV_jinc_1, PV_kinc_1);            \
            i2 = INC_IDX(i2_start, (i - ix), (j - iy), (k - iz),        \
                         nx, ny, sx2, PV_jinc_2, PV_kinc_2);            \
            i3 = INC_IDX(i3_start, (i - ix), (j - iy), (k - iz),        \
                         nx, ny, sx3, PV_jinc_3, PV_kinc_3);            \
            body;                                                       \
          }                                                             \
        }                                                               \
      }                                                                 \
  }

#undef BoxLoopReduceI1
#define BoxLoopReduceI1(locals, sum,                                    \
                        i, j, k,                                        \
                        ix, iy, iz, nx, ny, nz,                         \
                        i1, nx1, ny1, nz1, sx1, sy1, sz1,               \
                        body)                                           \
  {                                                                     \
    int i1_start = i1;                                                  \
    DeclareInc(PV_jinc_1, PV_kinc_1, nx, ny, nz, nx1, ny1, nz1, sx1, sy1, sz1); \
    PRAGMA(omp parallel for reduction(+:sum) collapse(3) private(i, j, k, i1 locals)) \
      for (k = iz; k < iz + nz; k++)                                    \
      {                                                                 \
        for (j = iy; j < iy + ny; j++)                                  \
        {                                                               \
          for (i = ix; i < ix + nx; i++)                                \
          {                                                             \
            i1 = INC_IDX(i1_start, (i - ix), (j - iy), (k - iz),        \
                         nx, ny, sx1, PV_jinc_1, PV_kinc_1);            \
            body;                                                       \
          }                                                             \
        }                                                               \
      }                                                                 \
  }

#define __BoxLoopReduceI1(locals, sum,                                  \
                          i, j, k,                                      \
                          ix, iy, iz, nx, ny, nz,                       \
                          i1, nx1, ny1, nz1, sx1, sy1, sz1,             \
                          body)                                         \
  {                                                                     \
    int i1_start = i1;                                                  \
    DeclareInc(PV_jinc_1, PV_kinc_1, nx, ny, nz, nx1, ny1, nz1, sx1, sy1, sz1); \
    PRAGMA(omp for reduction(+:sum) collapse(3) private(i, j, k, i1 locals)) \
      for (k = iz; k < iz + nz; k++)                                    \
      {                                                                 \
        for (j = iy; j < iy + ny; j++)                                  \
        {                                                               \
          for (i = ix; i < ix + nx; i++)                                \
          {                                                             \
            i1 = INC_IDX(i1_start, (i - ix), (j - iy), (k - iz),        \
                         nx, ny, sx1, PV_jinc_1, PV_kinc_1);            \
            body;                                                       \
          }                                                             \
        }                                                               \
      }                                                                 \
  }

#undef BoxLoopReduceI2
#define BoxLoopReduceI2(locals, sum,                                    \
                        i, j, k,                                        \
                        ix, iy, iz, nx, ny, nz,                         \
                        i1, nx1, ny1, nz1, sx1, sy1, sz1,               \
                        i2, nx2, ny2, nz2, sx2, sy2, sz2,               \
                        body)                                           \
  {                                                                     \
    int i1_start = i1;                                                  \
    int i2_start = i2;                                                  \
    DeclareInc(PV_jinc_1, PV_kinc_1, nx, ny, nz, nx1, ny1, nz1, sx1, sy1, sz1); \
    DeclareInc(PV_jinc_2, PV_kinc_2, nx, ny, nz, nx2, ny2, nz2, sx2, sy2, sz2); \
    PRAGMA(omp parallel for reduction(+:sum) collapse(3) private(i, j, k, i1, i2 locals)) \
      for (k = iz; k < iz + nz; k++)                                    \
      {                                                                 \
        for (j = iy; j < iy + ny; j++)                                  \
        {                                                               \
          for (i = ix; i < ix + nx; i++)                                \
          {                                                             \
            i1 = INC_IDX(i1_start, (i - ix), (j - iy), (k - iz),        \
                         nx, ny, sx1, PV_jinc_1, PV_kinc_1);            \
            i2 = INC_IDX(i2_start, (i - ix), (j - iy), (k - iz),        \
                         nx, ny, sx2, PV_jinc_2, PV_kinc_2);            \
            body;                                                       \
          }                                                             \
        }                                                               \
      }                                                                 \
  }

#endif // HAVE_OMP


#endif // _PF_OMPLOOPS_H
