#ifdef HAVE_CUDA

#include "cuda.h"
#include "cuda_runtime.h"

/*--------------------------------------------------------------------------
 * CUDA error handling macro
 *--------------------------------------------------------------------------*/

#define CUDA_ERR( err ) (gpu_error( err, __FILE__, __LINE__ ))
static void gpu_error(cudaError_t err, const char *file, int line) {
	if (err != cudaSuccess) {
		printf("\n\n%s in %s at line %d\n", cudaGetErrorString(err), file, line);
		exit(1);
	}
}

/*--------------------------------------------------------------------------
 * Define unified memory allocation routines for CUDA
 *--------------------------------------------------------------------------*/

static inline void *tallocCUDA(size_t size)
{
   void *ptr = NULL;

   CUDA_ERR(cudaMallocManaged((void**)&ptr, size, cudaMemAttachGlobal));

   return ptr;
}

static inline void *ctallocCUDA(size_t size)
{
   void *ptr = NULL;

   CUDA_ERR(cudaMallocManaged((void**)&ptr, size, cudaMemAttachGlobal));
   CUDA_ERR(cudaMemset(ptr, 0, size));

   return ptr;
}
static inline void tfreeCUDA(void *ptr)
{
   CUDA_ERR(cudaFree(ptr));
}

/*--------------------------------------------------------------------------
 * Refefine allocation macros for CUDA
 *--------------------------------------------------------------------------*/

// Redefine amps.h definitions
#undef amps_TAlloc
#define amps_TAlloc(type, count) \
  ((count) ? (type*)tallocCUDA(sizeof(type) * (unsigned int)(count)) : NULL)

#undef amps_CTAlloc
#define amps_CTAlloc(type, count) \
  ((count) ? (type*)ctallocCUDA(sizeof(type) * (unsigned int)(count)) : NULL)

#undef amps_TFree
#define amps_TFree(ptr) if (ptr) tfreeCUDA(ptr); else {}

// Redefine general.h definitions
#undef talloc
#define talloc(type, count) \
  ((count) ? (type*)tallocCUDA(sizeof(type) * (unsigned int)(count)) : NULL)

#undef ctalloc
#define ctalloc(type, count) \
  ((count) ? (type*)ctallocCUDA(sizeof(type) * (unsigned int)(count)) : NULL)

#undef tfree
#define tfree(ptr) if (ptr) tfreeCUDA(ptr); else {}

#endif