/**********************************************************************
 *
 *  Please read the LICENSE file for the GNU Lesser General Public License.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License (as published
 *  by the Free Software Foundation) version 2.1 dated February 1999.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms
 *  and conditions of the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 *  USA
 ***********************************************************************/

/* @file
 * @brief Contains macros, functions, and structs for CUDA compute kernels.
 */

#ifndef PFCUDALOOPS_H
#define PFCUDALOOPS_H

/*--------------------------------------------------------------------------
 * Include CUDA headers
 *--------------------------------------------------------------------------*/
#include "pfcudaerr.h"
#include "pfcudamalloc.h"

extern "C++"{

#include <tuple>
#include "cub.cuh"

/*--------------------------------------------------------------------------
 * CUDA block size definitions
 *--------------------------------------------------------------------------*/
#define BLOCKSIZE_MAX 1024
#define BLOCKSIZE_X 32
#define BLOCKSIZE_Y 8
#define BLOCKSIZE_Z 4

/*--------------------------------------------------------------------------
 * CUDA lambda definition (visible for host and device functions)
 *--------------------------------------------------------------------------*/
#define GPU_LAMBDA [=] __host__ __device__

/*--------------------------------------------------------------------------
 * CUDA helper macros and functions
 *--------------------------------------------------------------------------*/
#define RAND48_SEED_0   (0x330e)
#define RAND48_SEED_1   (0xabcd)
#define RAND48_SEED_2   (0x1234)
#define RAND48_MULT_0   (0xe66d)
#define RAND48_MULT_1   (0xdeec)
#define RAND48_MULT_2   (0x0005)
#define RAND48_ADD      (0x000b)

 /** Helper struct for type comparison. @note Not for direct use! */
template <typename T>
struct function_traits
    : public function_traits<decltype(&T::operator())>
{};
// For generic types, directly use the result of the signature of its 'operator()'

 /** Helper struct for type comparison. @note Not for direct use! */
template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits<ReturnType(ClassType::*)(Args...) const>
// we specialize for pointers to member function
{
    enum { arity = sizeof...(Args) };
    // arity is the number of arguments.

    typedef ReturnType result_type;

    template <size_t i>
    struct arg
    {
        typedef typename std::tuple_element<i, std::tuple<Args...>>::type type;
        // the i-th argument is equivalent to the i-th tuple element of a tuple
        // composed of those arguments.
    };
};

/** Device-callable dorand48() function for CUDA compute kernels. */
__host__ __device__ __forceinline__ static void dev_dorand48(unsigned short xseed[3])
{
  unsigned long accu;
  
  unsigned short _rand48_mult[3] = {
    RAND48_MULT_0,
    RAND48_MULT_1,
    RAND48_MULT_2
  };
  unsigned short _rand48_add = RAND48_ADD;
  unsigned short temp[2];

  accu = (unsigned long)_rand48_mult[0] * (unsigned long)xseed[0] +
         (unsigned long)_rand48_add;
  temp[0] = (unsigned short)accu;               /* lower 16 bits */
  accu >>= sizeof(unsigned short) * 8;
  accu += (unsigned long)_rand48_mult[0] * (unsigned long)xseed[1] +
          (unsigned long)_rand48_mult[1] * (unsigned long)xseed[0];
  temp[1] = (unsigned short)accu;               /* middle 16 bits */
  accu >>= sizeof(unsigned short) * 8;
  accu += _rand48_mult[0] * xseed[2] + _rand48_mult[1] * xseed[1] + _rand48_mult[2] * xseed[0];
  xseed[0] = temp[0];
  xseed[1] = temp[1];
  xseed[2] = (unsigned short)accu;
}

/** Device-callable erand48() function for CUDA compute kernels. */
__host__ __device__ __forceinline__ static double dev_erand48(unsigned short xseed[3])
{
  dev_dorand48(xseed);
  return ldexp((double)xseed[0], -48) +
         ldexp((double)xseed[1], -32) +
         ldexp((double)xseed[2], -16);
}

/** Device-callable drand48() function for CUDA compute kernels. */
__host__ __device__ __forceinline__ static double dev_drand48(void)
{
  unsigned short _rand48_seed[3] = {
    RAND48_SEED_0,
    RAND48_SEED_1,
    RAND48_SEED_2
  };
  return dev_erand48(_rand48_seed);
}

/** Device-callable RPowerR() function for CUDA compute kernels. */
template <typename T>
__host__ __device__ __forceinline__ static T RPowerR(T base, T exponent)
{
  if (base <= 0.0)
    return(0.0);

  return((T)pow((double)base, (double)exponent));
}

/**
 * Thread-safe function to find minimum value in compute kernels.
 * The function definition depends whether called from host or device code.
 *
 * @note Not for direct use!
 *
 * @param address pointer to first value [IN], pointer to min value [OUT]
 * @param val2 second value [IN]
 */
__host__ __device__ __forceinline__ static void AtomicMin(int *address, int val2)
{
#ifdef __CUDA_ARCH__
    atomicMin(address, val2);
#else
    if(*address > val2) *address = val2;
#endif
}

/**
 * Thread-safe function to find maximum value in compute kernels.
 * The function definition depends whether called from host or device code.
 *
 * @note Not for direct use!
 *
 * @param address pointer to first value [IN], pointer to max value [OUT]
 * @param val2 second value [IN]
 */
__host__ __device__ __forceinline__ static void AtomicMax(int *address, int val2)
{
#ifdef __CUDA_ARCH__
    atomicMax(address, val2);
#else
    if(*address < val2) *address = val2;
#endif
}

/**
 * Thread-safe function to find minimum value in compute kernels.
 * The function definition depends whether called from host or device code.
 *
 * @note Not for direct use!
 *
 * @param address pointer to first value [IN], pointer to min value [OUT]
 * @param val2 second value [IN]
 */
__host__ __device__ __forceinline__ static void AtomicMin(double *address, double val2)
{
#ifdef __CUDA_ARCH__
    unsigned long long ret = __double_as_longlong(*address);
    while(val2 < __longlong_as_double(ret))
    {
        unsigned long long old = ret;
        if((ret = atomicCAS((unsigned long long *)address, old, __double_as_longlong(val2))) == old)
            break;
    }
    //return __longlong_as_double(ret);
#else
    if(*address > val2) *address = val2;
#endif
}

/**
 * Thread-safe function to find maximum value in compute kernels.
 * The function definition depends whether called from host or device code.
 *
 * @note Not for direct use!
 *
 * @param address pointer to first value [IN], pointer to max value [OUT]
 * @param val2 second value [IN]
 */
__host__ __device__ __forceinline__ static void AtomicMax(double *address, double val2)
{
#ifdef __CUDA_ARCH__
    unsigned long long ret = __double_as_longlong(*address);
    while(val2 > __longlong_as_double(ret))
    {
        unsigned long long old = ret;
        if((ret = atomicCAS((unsigned long long *)address, old, __double_as_longlong(val2))) == old)
            break;
    }
    //return __longlong_as_double(ret);
#else
    if(*address < val2) *address = val2;
#endif
}

/**
 * Thread-safe addition assignment for compute kernels.
 * The function definition depends whether called from host or device code.
 *
 * @note Not for direct use!
 *
 * @param array_loc original value [IN], sum result [OUT]
 * @param value value to be added [IN]
 */
template <typename T>
__host__ __device__ __forceinline__ static void AtomicAdd(T *array_loc, T value)
{
    //Define this function depending on whether it runs on GPU or CPU
#ifdef __CUDA_ARCH__
    atomicAdd(array_loc, value);
#else
    *array_loc += value;
#endif
}

 /** Helper struct for type comparison. @note Not for direct use! */
struct SkipParallelSync {const int dummy = 0;};
#undef SKIP_PARALLEL_SYNC
#define SKIP_PARALLEL_SYNC struct SkipParallelSync sync_struct; return sync_struct;

#undef PARALLEL_SYNC
#define PARALLEL_SYNC CUDA_ERR(cudaStreamSynchronize(0)); 

#undef PlusEquals
#define PlusEquals(a, b) AtomicAdd(&(a), b)

 /** Helper struct for type comparison. @note Not for direct use! */
template <typename T>
struct ReduceMaxType {T value;};
#undef ReduceMax
#define ReduceMax(a, b) struct ReduceMaxType<std::decay<decltype(a)>::type> reduce_struct {.value = b}; return reduce_struct;

 /** Helper struct for type comparison. @note Not for direct use! */
template <typename T>
struct ReduceMinType {T value;};
#undef ReduceMin
#define ReduceMin(a, b) struct ReduceMinType<std::decay<decltype(a)>::type> reduce_struct {.value = b}; return reduce_struct;

 /** Helper struct for type comparison. @note Not for direct use! */
template <typename T>
struct ReduceSumType {T value;};
#undef ReduceSum
#define ReduceSum(a, b) struct ReduceSumType<std::decay<decltype(a)>::type> reduce_struct {.value = b}; return reduce_struct;

/*--------------------------------------------------------------------------
 * CUDA loop kernels
 *--------------------------------------------------------------------------*/

/**
 * @brief CUDA compute kernel.
 *
 * @param loop_body A lambda function that evaluates the loop body [IN/OUT]
 * @param ix The loop starting index for the first dim [IN]
 * @param iy The loop starting index for the second dim [IN]
 * @param iz The loop starting index for the third dim [IN]
 * @param nx The size of the first dim [IN]
 * @param ny The size of the second dim [IN]
 * @param nz The size of the third dim [IN]
 */
template <typename LambdaBody>
__global__ static void 
__launch_bounds__(BLOCKSIZE_MAX)
BoxKernelI0(LambdaBody loop_body, 
    const int ix, const int iy, const int iz, const int nx, const int ny, const int nz)
{

    int i = ((blockIdx.x*blockDim.x)+threadIdx.x);
    if(i >= nx)return;
    int j = ((blockIdx.y*blockDim.y)+threadIdx.y);
    if(j >= ny)return;
    int k = ((blockIdx.z*blockDim.z)+threadIdx.z);
    if(k >= nz)return;

    i += ix;
    j += iy;
    k += iz;
    
    loop_body(i, j, k);
}


/**
 * @brief CUDA compute kernel.
 *
 * @param loop_init A lambda function that calculates a thread-local index [IN]
 * @param loop_body A lambda function that evaluates the loop body [IN/OUT]
 * @param ix The loop starting index for the first dim [IN]
 * @param iy The loop starting index for the second dim [IN]
 * @param iz The loop starting index for the third dim [IN]
 * @param nx The size of the first dim [IN]
 * @param ny The size of the second dim [IN]
 * @param nz The size of the third dim [IN]
 */
template <typename LambdaInit, typename LambdaBody>
__global__ static void 
__launch_bounds__(BLOCKSIZE_MAX)
BoxKernelI1(LambdaInit loop_init, LambdaBody loop_body, 
    const int ix, const int iy, const int iz, const int nx, const int ny, const int nz)
{

    int i = ((blockIdx.x*blockDim.x)+threadIdx.x);
    if(i >= nx)return;
    int j = ((blockIdx.y*blockDim.y)+threadIdx.y);
    if(j >= ny)return;
    int k = ((blockIdx.z*blockDim.z)+threadIdx.z);
    if(k >= nz)return;

    const int i1 = loop_init(i, j, k);

    i += ix;
    j += iy;
    k += iz;
    
    loop_body(i, j, k, i1);       
}

/**
 * @brief CUDA compute kernel.
 *
 * @param loop_init1 A lambda function that calculates a thread-local index [IN]
 * @param loop_init2 A lambda function that calculates a thread-local index [IN]
 * @param loop_body A lambda function that evaluates the loop body [IN/OUT]
 * @param ix The loop starting index for the first dim [IN]
 * @param iy The loop starting index for the second dim [IN]
 * @param iz The loop starting index for the third dim [IN]
 * @param nx The size of the first dim [IN]
 * @param ny The size of the second dim [IN]
 * @param nz The size of the third dim [IN]
 */
template <typename LambdaInit1, typename LambdaInit2, typename LambdaBody>
__global__ static void 
__launch_bounds__(BLOCKSIZE_MAX)
BoxKernelI2(LambdaInit1 loop_init1, LambdaInit2 loop_init2, LambdaBody loop_body, 
    const int ix, const int iy, const int iz, const int nx, const int ny, const int nz)
{

    int i = ((blockIdx.x*blockDim.x)+threadIdx.x);
    if(i >= nx)return;
    int j = ((blockIdx.y*blockDim.y)+threadIdx.y);
    if(j >= ny)return;
    int k = ((blockIdx.z*blockDim.z)+threadIdx.z);
    if(k >= nz)return;

    const int i1 = loop_init1(i, j, k);
    const int i2 = loop_init2(i, j, k);

    i += ix;
    j += iy;
    k += iz;
    
    loop_body(i, j, k, i1, i2);    
}

/**
 * @brief CUDA compute kernel.
 *
 * @param loop_init1 A lambda function that calculates a thread-local index [IN]
 * @param loop_init2 A lambda function that calculates a thread-local index [IN]
 * @param loop_init3 A lambda function that calculates a thread-local index [IN]
 * @param loop_body A lambda function that evaluates the loop body [IN/OUT]
 * @param ix The loop starting index for the first dim [IN]
 * @param iy The loop starting index for the second dim [IN]
 * @param iz The loop starting index for the third dim [IN]
 * @param nx The size of the first dim [IN]
 * @param ny The size of the second dim [IN]
 * @param nz The size of the third dim [IN]
 */
template <typename LambdaInit1, typename LambdaInit2, typename LambdaInit3, typename LambdaBody>
__global__ static void 
__launch_bounds__(BLOCKSIZE_MAX)
BoxKernelI3(LambdaInit1 loop_init1, LambdaInit2 loop_init2, LambdaInit3 loop_init3, 
    LambdaBody loop_body, const int ix, const int iy, const int iz, const int nx, const int ny, const int nz)
{

    int i = ((blockIdx.x*blockDim.x)+threadIdx.x);
    if(i >= nx)return;
    int j = ((blockIdx.y*blockDim.y)+threadIdx.y);
    if(j >= ny)return;
    int k = ((blockIdx.z*blockDim.z)+threadIdx.z);
    if(k >= nz)return;

    const int i1 = loop_init1(i, j, k);
    const int i2 = loop_init2(i, j, k);
    const int i3 = loop_init3(i, j, k);

    i += ix;
    j += iy;
    k += iz;
    
    loop_body(i, j, k, i1, i2, i3);    
}

/**
 * @brief CUDA compute kernel for parallel reductions.
 *
 * @param loop_init1 A lambda function that calculates a thread-local index [IN]
 * @param loop_init2 A lambda function that calculates a thread-local index [IN]
 * @param loop_fun A lambda function that evaluates the loop body [IN/OUT]
 * @param init_val The initial value of the reduction variable [IN]
 * @param rslt A pointer to the result/sum variable [OUT]
 * @param ix The loop starting index for the first dim [IN]
 * @param iy The loop starting index for the second dim [IN]
 * @param iz The loop starting index for the third dim [IN]
 * @param nx The size of the first dim [IN]
 * @param ny The size of the second dim [IN]
 * @param nz The size of the third dim [IN]
 */
template <typename ReduceOp, typename LambdaInit1, typename LambdaInit2, typename LambdaFun, typename T>
__global__ static void 
__launch_bounds__(BLOCKSIZE_MAX)
DotKernelI2(LambdaInit1 loop_init1, LambdaInit2 loop_init2, LambdaFun loop_fun, 
    const T init_val, T * __restrict__ rslt, const int ix, const int iy, const int iz, 
    const int nx, const int ny, const int nz)
{
    // Specialize BlockReduce for a 1D block of BLOCKSIZE_X * 1 * 1 threads on type T
#ifdef __CUDA_ARCH__
    typedef cub::BlockReduce<T, BLOCKSIZE_MAX, cub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY, 1, 1, __CUDA_ARCH__> BlockReduce;
#else
    typedef cub::BlockReduce<T, BLOCKSIZE_MAX, cub::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY, 1, 1> BlockReduce;
#endif

    __shared__ typename BlockReduce::TempStorage temp_storage;

    const int idx = ((blockIdx.x*blockDim.x)+threadIdx.x);

    const int i = idx % nx;
    const int j = (idx / nx) % ny;
    const int k = idx / (nx * ny);
    const int ntot = nx * ny * nz;

    T thread_data;
    
    // Initialize thread_data depending on reduction operation
    if(std::is_same<ReduceOp, struct ReduceSumType<T>>::value)
      thread_data = 0;
    else 
      thread_data = init_val;

    // Evaluate the loop body
    if (idx < ntot)
    {
      const int i1 = loop_init1(i, j, k);
      const int i2 = loop_init2(i, j, k);

      thread_data = loop_fun(i, j, k, i1, i2).value;
    }

    // Perform reductions
    if(std::is_same<ReduceOp, struct ReduceSumType<T>>::value)
    {
      // Compute the block-wide sum for thread0
      T aggregate = BlockReduce(temp_storage).Sum(thread_data);

      // Store aggregate
      if(threadIdx.x == 0) 
      {
        atomicAdd(rslt, aggregate);
      }
    }
    else if(std::is_same<ReduceOp, struct ReduceMaxType<T>>::value)
    {
      // Compute the block-wide sum for thread0
      T aggregate = BlockReduce(temp_storage).Reduce(thread_data, cub::Max());

      // Store aggregate
      if(threadIdx.x == 0) 
      {
        AtomicMax(rslt, aggregate);
      }

      // Write to global memory directly from all threads
      // if (idx < ntot)
      // {
      //   AtomicMax(rslt, thread_data);
      // }
    }
    else if(std::is_same<ReduceOp, struct ReduceMinType<T>>::value)
    {
      // Compute the block-wide sum for thread0
      T aggregate = BlockReduce(temp_storage).Reduce(thread_data, cub::Min());

      // Store aggregate
      if(threadIdx.x == 0) 
      {
        AtomicMin(rslt, aggregate);
      }
      
      // Write to global memory directly from all threads
      // if (idx < ntot)
      // {
      //   AtomicMin(rslt, thread_data);
      // }
    }
    else 
    {
      printf("ERROR at %s:%d: Invalid reduction identifier, likely a problem with a BoxLoopReduce body.", __FILE__, __LINE__);
    }
}

/*--------------------------------------------------------------------------
 * CUDA loop macro redefinitions
 *--------------------------------------------------------------------------*/

/**
 * A macro for finding the 3D CUDA grid and block dimensions.
 *
 * @note Not for direct use!
 *
 * @param grid grid dimensions [OUT]
 * @param block block dimensions [OUT]
 * @param nx The size of the first dim [IN]
 * @param ny The size of the second dim [IN]
 * @param nz The size of the third dim [IN]
 * @param dyn_blocksize Runtime adjustment of y- and z-blocksizes [IN]
 */
#define FindDims(grid, block, nx, ny, nz, dyn_blocksize)                            \
{                                                                                   \
  int blocksize_x = BLOCKSIZE_X;                                                    \
  int blocksize_y = BLOCKSIZE_Y;                                                    \
  int blocksize_z = BLOCKSIZE_Z;                                                    \
  grid = dim3(((nx - 1) + blocksize_x) / blocksize_x,                               \
      ((ny - 1) + blocksize_y) / blocksize_y,                                       \
      ((nz - 1) + blocksize_z) / blocksize_z);                                      \
  while(dyn_blocksize && (grid.x*grid.y*grid.z < 80)                                \
          && ((blocksize_y * blocksize_z) >= 4 ))                                   \
  {                                                                                 \
      if ( blocksize_z >= 2 )                                                       \
          blocksize_z /= 2;                                                         \
      else                                                                          \
          blocksize_y /= 2;                                                         \
      grid = dim3(((nx - 1) + blocksize_x) / blocksize_x,                           \
      ((ny - 1) + blocksize_y) / blocksize_y,                                       \
      ((nz - 1) + blocksize_z) / blocksize_z);                                      \
  }                                                                                 \
  block = dim3(blocksize_x, blocksize_y, blocksize_z);                              \
}  

#undef BoxLoopI1
 /** Loop redefinition for CUDA. */
#define BoxLoopI1(i, j, k,                                                          \
  ix, iy, iz, nx, ny, nz,                                                           \
  i1, nx1, ny1, nz1, sx1, sy1, sz1,                                                 \
  loop_body)                                                                        \
{                                                                                   \
  if(nx > 0 && ny > 0 && nz > 0)                                                    \
  {                                                                                 \
    DeclareInc(PV_jinc_1, PV_kinc_1, nx, ny, nz, nx1, ny1, nz1, sx1, sy1, sz1);     \
                                                                                    \
    dim3 block, grid;                                                               \
    FindDims(grid, block, nx, ny, nz, 1);                                           \
                                                                                    \
    auto lambda_init =                                                              \
        GPU_LAMBDA(const int i, const int j, const int k)                           \
        {                                                                           \
            return k * PV_kinc_1 + (k * ny + j) * PV_jinc_1                         \
            + (k * ny * nx + j * nx + i) * sx1 + i1;                                \
        };                                                                          \
    auto lambda_body =                                                              \
        GPU_LAMBDA(const int i, const int j, const int k, const int i1)             \
        loop_body;                                                                  \
                                                                                    \
    BoxKernelI1<<<grid, block>>>(                                                   \
        lambda_init, lambda_body, ix, iy, iz, nx, ny, nz);                          \
    CUDA_ERR(cudaPeekAtLastError());                                                \
                                                                                    \
    typedef function_traits<decltype(lambda_body)> traits;                          \
    if(!std::is_same<traits::result_type, struct SkipParallelSync>::value)          \
      CUDA_ERR(cudaStreamSynchronize(0));                                           \
  }                                                                                 \
  (void)i;(void)j;(void)k;                                                          \
}

#undef BoxLoopI2
 /** Loop redefinition for CUDA. */
#define BoxLoopI2(i, j, k,                                                          \
  ix, iy, iz, nx, ny, nz,                                                           \
  i1, nx1, ny1, nz1, sx1, sy1, sz1,                                                 \
  i2, nx2, ny2, nz2, sx2, sy2, sz2,                                                 \
  loop_body)                                                                        \
{                                                                                   \
  if(nx > 0 && ny > 0 && nz > 0)                                                    \
  {                                                                                 \
    DeclareInc(PV_jinc_1, PV_kinc_1, nx, ny, nz, nx1, ny1, nz1, sx1, sy1, sz1);     \
    DeclareInc(PV_jinc_2, PV_kinc_2, nx, ny, nz, nx2, ny2, nz2, sx2, sy2, sz2);     \
                                                                                    \
    dim3 block, grid;                                                               \
    FindDims(grid, block, nx, ny, nz, 1);                                           \
                                                                                    \
    auto lambda_init1 =                                                             \
        GPU_LAMBDA(const int i, const int j, const int k)                           \
        {                                                                           \
            return k * PV_kinc_1 + (k * ny + j) * PV_jinc_1                         \
            + (k * ny * nx + j * nx + i) * sx1 + i1;                                \
        };                                                                          \
    auto lambda_init2 =                                                             \
        GPU_LAMBDA(const int i, const int j, const int k)                           \
        {                                                                           \
            return k * PV_kinc_2 + (k * ny + j) * PV_jinc_2                         \
            + (k * ny * nx + j * nx + i) * sx2 + i2;                                \
        };                                                                          \
    auto lambda_body =                                                              \
        GPU_LAMBDA(const int i, const int j, const int k,                           \
        const int i1, const int i2)loop_body;                                       \
                                                                                    \
    BoxKernelI2<<<grid, block>>>(                                                   \
        lambda_init1, lambda_init2, lambda_body, ix, iy, iz, nx, ny, nz);           \
    CUDA_ERR(cudaPeekAtLastError());                                                \
                                                                                    \
    typedef function_traits<decltype(lambda_body)> traits;                          \
    if(!std::is_same<traits::result_type, struct SkipParallelSync>::value)          \
      CUDA_ERR(cudaStreamSynchronize(0));                                           \
  }                                                                                 \
  (void)i;(void)j;(void)k;                                                          \
}

#undef BoxLoopI3
 /** Loop redefinition for CUDA. */
#define BoxLoopI3(i, j, k,                                                          \
  ix, iy, iz, nx, ny, nz,                                                           \
  i1, nx1, ny1, nz1, sx1, sy1, sz1,                                                 \
  i2, nx2, ny2, nz2, sx2, sy2, sz2,                                                 \
  i3, nx3, ny3, nz3, sx3, sy3, sz3,                                                 \
  loop_body)                                                                        \
{                                                                                   \
  if(nx > 0 && ny > 0 && nz > 0)                                                    \
  {                                                                                 \
    DeclareInc(PV_jinc_1, PV_kinc_1, nx, ny, nz, nx1, ny1, nz1, sx1, sy1, sz1);     \
    DeclareInc(PV_jinc_2, PV_kinc_2, nx, ny, nz, nx2, ny2, nz2, sx2, sy2, sz2);     \
    DeclareInc(PV_jinc_3, PV_kinc_3, nx, ny, nz, nx3, ny3, nz3, sx3, sy3, sz3);     \
                                                                                    \
    dim3 block, grid;                                                               \
    FindDims(grid, block, nx, ny, nz, 1);                                           \
                                                                                    \
    auto lambda_init1 =                                                             \
        GPU_LAMBDA(const int i, const int j, const int k)                           \
        {                                                                           \
            return k * PV_kinc_1 + (k * ny + j) * PV_jinc_1                         \
            + (k * ny * nx + j * nx + i) * sx1 + i1;                                \
        };                                                                          \
    auto lambda_init2 =                                                             \
        GPU_LAMBDA(const int i, const int j, const int k)                           \
        {                                                                           \
            return k * PV_kinc_2 + (k * ny + j) * PV_jinc_2                         \
            + (k * ny * nx + j * nx + i) * sx2 + i2;                                \
        };                                                                          \
    auto lambda_init3 =                                                             \
        GPU_LAMBDA(const int i, const int j, const int k)                           \
        {                                                                           \
            return k * PV_kinc_3 + (k * ny + j) * PV_jinc_3                         \
            + (k * ny * nx + j * nx + i) * sx3 + i3;                                \
        };                                                                          \
    auto lambda_body =                                                              \
        GPU_LAMBDA(const int i, const int j, const int k,                           \
        const int i1, const int i2, const int i3)loop_body;                         \
                                                                                    \
    BoxKernelI3<<<grid, block>>>(lambda_init1, lambda_init2, lambda_init3,          \
        lambda_body, ix, iy, iz, nx, ny, nz);                                       \
                                                                                    \
    typedef function_traits<decltype(lambda_body)> traits;                          \
    if(!std::is_same<traits::result_type, struct SkipParallelSync>::value)          \
      CUDA_ERR(cudaStreamSynchronize(0));                                           \
  }                                                                                 \
  (void)i;(void)j;(void)k;                                                          \
}

#undef BoxLoopReduceI1
 /** Loop redefinition for CUDA. */
#define BoxLoopReduceI1(rslt, i, j, k,                                              \
  ix, iy, iz, nx, ny, nz,                                                           \
  i1, nx1, ny1, nz1, sx1, sy1, sz1,                                                 \
  loop_body)                                                                        \
{                                                                                   \
  if(nx > 0 && ny > 0 && nz > 0)                                                    \
  {                                                                                 \
    DeclareInc(PV_jinc_1, PV_kinc_1, nx, ny, nz, nx1, ny1, nz1, sx1, sy1, sz1);     \
                                                                                    \
    int block = BLOCKSIZE_MAX;                                                      \
    int grid = ((nx * ny * nz - 1) + block) / block;                                \
                                                                                    \
    auto lambda_init =                                                              \
        GPU_LAMBDA(const int i, const int j, const int k)                           \
        {                                                                           \
            return k * PV_kinc_1 + (k * ny + j) * PV_jinc_1                         \
            + (k * ny * nx + j * nx + i) * sx1 + i1;                                \
        };                                                                          \
    auto &ref_rslt = rslt;                                                          \
    auto lambda_body =                                                              \
        GPU_LAMBDA(const int i, const int j, const int k,                           \
                   const int i1, const int i2)                                      \
        {                                                                           \
            auto rslt = ref_rslt;                                                   \
            loop_body;                                                              \
        };                                                                          \
                                                                                    \
    decltype(rslt) *ptr_rslt = (decltype(rslt)*)talloc_cuda(sizeof(decltype(rslt)));\
    MemPrefetchDeviceToHost(ptr_rslt, sizeof(decltype(rslt)), 0);                   \
    *ptr_rslt = rslt;                                                               \
    MemPrefetchHostToDevice(ptr_rslt, sizeof(decltype(rslt)), 0);                   \
                                                                                    \
    typedef function_traits<decltype(lambda_body)> traits;                          \
    DotKernelI2<traits::result_type><<<grid, block>>>(lambda_init, lambda_init,     \
      lambda_body, rslt, ptr_rslt, ix, iy, iz, nx, ny, nz);                         \
    CUDA_ERR(cudaPeekAtLastError());                                                \
    CUDA_ERR(cudaStreamSynchronize(0));                                             \
                                                                                    \
    MemPrefetchDeviceToHost(ptr_rslt, sizeof(decltype(rslt)), 0);                   \
    rslt = *ptr_rslt;                                                               \
    tfree_cuda(ptr_rslt);                                                           \
  }                                                                                 \
  (void)i;(void)j;(void)k;                                                          \
}

#undef BoxLoopReduceI2
 /** Loop redefinition for CUDA. */
#define BoxLoopReduceI2(rslt, i, j, k,                                              \
  ix, iy, iz, nx, ny, nz,                                                           \
  i1, nx1, ny1, nz1, sx1, sy1, sz1,                                                 \
  i2, nx2, ny2, nz2, sx2, sy2, sz2,                                                 \
  loop_body)                                                                        \
{                                                                                   \
  if(nx > 0 && ny > 0 && nz > 0)                                                    \
  {                                                                                 \
    DeclareInc(PV_jinc_1, PV_kinc_1, nx, ny, nz, nx1, ny1, nz1, sx1, sy1, sz1);     \
    DeclareInc(PV_jinc_2, PV_kinc_2, nx, ny, nz, nx2, ny2, nz2, sx2, sy2, sz2);     \
                                                                                    \
    int block = BLOCKSIZE_MAX;                                                      \
    int grid = ((nx * ny * nz - 1) + block) / block;                                \
                                                                                    \
    auto lambda_init1 =                                                             \
        GPU_LAMBDA(const int i, const int j, const int k)                           \
        {                                                                           \
            return k * PV_kinc_1 + (k * ny + j) * PV_jinc_1                         \
            + (k * ny * nx + j * nx + i) * sx1 + i1;                                \
        };                                                                          \
    auto lambda_init2 =                                                             \
        GPU_LAMBDA(const int i, const int j, const int k)                           \
        {                                                                           \
            return k * PV_kinc_2 + (k * ny + j) * PV_jinc_2                         \
            + (k * ny * nx + j * nx + i) * sx2 + i2;                                \
        };                                                                          \
    auto &ref_rslt = rslt;                                                          \
    auto lambda_body =                                                              \
        GPU_LAMBDA(const int i, const int j, const int k,                           \
                   const int i1, const int i2)                                      \
        {                                                                           \
            auto rslt = ref_rslt;                                                   \
            loop_body;                                                              \
        };                                                                          \
                                                                                    \
    decltype(rslt) *ptr_rslt = (decltype(rslt)*)talloc_cuda(sizeof(decltype(rslt)));\
    MemPrefetchDeviceToHost(ptr_rslt, sizeof(decltype(rslt)), 0);                   \
    *ptr_rslt = rslt;                                                               \
    MemPrefetchHostToDevice(ptr_rslt, sizeof(decltype(rslt)), 0);                   \
                                                                                    \
    typedef function_traits<decltype(lambda_body)> traits;                          \
    DotKernelI2<traits::result_type><<<grid, block>>>(lambda_init1, lambda_init2,   \
      lambda_body, rslt, ptr_rslt, ix, iy, iz, nx, ny, nz);                         \
    CUDA_ERR(cudaPeekAtLastError());                                                \
    CUDA_ERR(cudaStreamSynchronize(0));                                             \
                                                                                    \
    MemPrefetchDeviceToHost(ptr_rslt, sizeof(decltype(rslt)), 0);                   \
    rslt = *ptr_rslt;                                                               \
    tfree_cuda(ptr_rslt);                                                           \
  }                                                                                 \
  (void)i;(void)j;(void)k;                                                          \
}

#undef GrGeomInLoopBoxes
 /** Loop redefinition for CUDA. */
#define GrGeomInLoopBoxes(i, j, k,                                                  \
  grgeom, ix, iy, iz, nx, ny, nz, loop_body)                                        \
{                                                                                   \
  BoxArray* boxes = GrGeomSolidInteriorBoxes(grgeom);                               \
  for (int PV_box = 0; PV_box < BoxArraySize(boxes); PV_box++)                      \
  {                                                                                 \
    Box box = BoxArrayGetBox(boxes, PV_box);                                        \
    /* find octree and region intersection */                                       \
    int PV_ixl = pfmax(ix, box.lo[0]);                                              \
    int PV_iyl = pfmax(iy, box.lo[1]);                                              \
    int PV_izl = pfmax(iz, box.lo[2]);                                              \
    int PV_ixu = pfmin((ix + nx - 1), box.up[0]);                                   \
    int PV_iyu = pfmin((iy + ny - 1), box.up[1]);                                   \
    int PV_izu = pfmin((iz + nz - 1), box.up[2]);                                   \
                                                                                    \
    if(PV_ixl <= PV_ixu && PV_iyl <= PV_iyu && PV_izl <= PV_izu)                    \
    {                                                                               \
      int PV_diff_x = PV_ixu - PV_ixl;                                              \
      int PV_diff_y = PV_iyu - PV_iyl;                                              \
      int PV_diff_z = PV_izu - PV_izl;                                              \
                                                                                    \
      dim3 block, grid;                                                             \
      FindDims(grid, block, PV_diff_x + 1, PV_diff_y + 1, PV_diff_z + 1, 1);        \
                                                                                    \
      int nx = PV_diff_x + 1;                                                       \
      int ny = PV_diff_y + 1;                                                       \
      int nz = PV_diff_z + 1;                                                       \
                                                                                    \
      BoxKernelI0<<<grid, block>>>(                                                 \
          GPU_LAMBDA(const int i, const int j, const int k)loop_body,               \
          PV_ixl, PV_iyl, PV_izl, nx, ny, nz);                                      \
      CUDA_ERR(cudaPeekAtLastError());                                              \
      CUDA_ERR(cudaStreamSynchronize(0));                                           \
    }                                                                               \
  }                                                                                 \
  (void)i;(void)j;(void)k;                                                          \
}

#undef GrGeomSurfLoopBoxes
 /** Loop redefinition for CUDA. */
#define GrGeomSurfLoopBoxes(i, j, k, fdir, grgeom,                                  \
  ix, iy, iz, nx, ny, nz, loop_body)                                                \
{                                                                                   \
  int PV_fdir[3];                                                                   \
  fdir = PV_fdir;                                                                   \
  for (int PV_f = 0; PV_f < GrGeomOctreeNumFaces; PV_f++)                           \
  {                                                                                 \
    switch (PV_f)                                                                   \
    {                                                                               \
    case GrGeomOctreeFaceL:                                                         \
        fdir[0] = -1; fdir[1] = 0; fdir[2] = 0;                                     \
        break;                                                                      \
    case GrGeomOctreeFaceR:                                                         \
        fdir[0] = 1; fdir[1] = 0; fdir[2] = 0;                                      \
        break;                                                                      \
    case GrGeomOctreeFaceD:                                                         \
        fdir[0] = 0; fdir[1] = -1; fdir[2] = 0;                                     \
        break;                                                                      \
    case GrGeomOctreeFaceU:                                                         \
        fdir[0] = 0; fdir[1] = 1; fdir[2] = 0;                                      \
        break;                                                                      \
    case GrGeomOctreeFaceB:                                                         \
        fdir[0] = 0; fdir[1] = 0; fdir[2] = -1;                                     \
        break;                                                                      \
    case GrGeomOctreeFaceF:                                                         \
        fdir[0] = 0; fdir[1] = 0; fdir[2] = 1;                                      \
        break;                                                                      \
    default:                                                                        \
        fdir[0] = -9999; fdir[1] = -9999; fdir[2] = -99999;                         \
        break;                                                                      \
    }                                                                               \
                                                                                    \
    BoxArray* boxes = GrGeomSolidSurfaceBoxes(grgeom, PV_f);                        \
    for (int PV_box = 0; PV_box < BoxArraySize(boxes); PV_box++)                    \
    {                                                                               \
      Box box = BoxArrayGetBox(boxes, PV_box);                                      \
      /* find octree and region intersection */                                     \
      int PV_ixl = pfmax(ix, box.lo[0]);                                            \
      int PV_iyl = pfmax(iy, box.lo[1]);                                            \
      int PV_izl = pfmax(iz, box.lo[2]);                                            \
      int PV_ixu = pfmin((ix + nx - 1), box.up[0]);                                 \
      int PV_iyu = pfmin((iy + ny - 1), box.up[1]);                                 \
      int PV_izu = pfmin((iz + nz - 1), box.up[2]);                                 \
                                                                                    \
      if(PV_ixl <= PV_ixu && PV_iyl <= PV_iyu && PV_izl <= PV_izu)                  \
      {                                                                             \
        int PV_diff_x = PV_ixu - PV_ixl;                                            \
        int PV_diff_y = PV_iyu - PV_iyl;                                            \
        int PV_diff_z = PV_izu - PV_izl;                                            \
                                                                                    \
        dim3 block, grid;                                                           \
        FindDims(grid, block, PV_diff_x + 1, PV_diff_y + 1, PV_diff_z + 1, 1);      \
                                                                                    \
        int nx = PV_diff_x + 1;                                                     \
        int ny = PV_diff_y + 1;                                                     \
        int nz = PV_diff_z + 1;                                                     \
                                                                                    \
        const int fdir_capt0 = fdir[0];                                             \
        const int fdir_capt1 = fdir[1];                                             \
        const int fdir_capt2 = fdir[2];                                             \
                                                                                    \
        auto lambda_body =                                                          \
          GPU_LAMBDA(const int i, const int j, const int k)                         \
          {                                                                         \
            const int fdir[3] = {fdir_capt0, fdir_capt1, fdir_capt2};               \
            loop_body;                                                              \
            (void)fdir;                                                             \
          };                                                                        \
                                                                                    \
        BoxKernelI0<<<grid, block>>>(lambda_body,                                   \
            PV_ixl, PV_iyl, PV_izl, nx, ny, nz);                                    \
        CUDA_ERR(cudaPeekAtLastError());                                            \
        CUDA_ERR(cudaStreamSynchronize(0));                                         \
      }                                                                             \
    }                                                                               \
  }                                                                                 \
  (void)i;(void)j;(void)k;                                                          \
}

#undef GrGeomPatchLoopBoxes
 /** Loop redefinition for CUDA. */
#define GrGeomPatchLoopBoxes(i, j, k, fdir, grgeom, patch_num,                      \
  ix, iy, iz, nx, ny, nz, loop_body)                                                \
{                                                                                   \
  int PV_fdir[3];                                                                   \
  fdir = PV_fdir;                                                                   \
  for (int PV_f = 0; PV_f < GrGeomOctreeNumFaces; PV_f++)                           \
  {                                                                                 \
    switch (PV_f)                                                                   \
    {                                                                               \
      case GrGeomOctreeFaceL:                                                       \
        fdir[0] = -1; fdir[1] = 0; fdir[2] = 0;                                     \
        break;                                                                      \
      case GrGeomOctreeFaceR:                                                       \
        fdir[0] = 1; fdir[1] = 0; fdir[2] = 0;                                      \
        break;                                                                      \
      case GrGeomOctreeFaceD:                                                       \
        fdir[0] = 0; fdir[1] = -1; fdir[2] = 0;                                     \
        break;                                                                      \
      case GrGeomOctreeFaceU:                                                       \
        fdir[0] = 0; fdir[1] = 1; fdir[2] = 0;                                      \
        break;                                                                      \
      case GrGeomOctreeFaceB:                                                       \
        fdir[0] = 0; fdir[1] = 0; fdir[2] = -1;                                     \
        break;                                                                      \
      case GrGeomOctreeFaceF:                                                       \
        fdir[0] = 0; fdir[1] = 0; fdir[2] = 1;                                      \
        break;                                                                      \
      default:                                                                      \
        fdir[0] = -9999; fdir[1] = -9999; fdir[2] = -99999;                         \
        break;                                                                      \
    }                                                                               \
                                                                                    \
    int n_prev = 0;                                                                 \
    BoxArray* boxes = GrGeomSolidPatchBoxes(grgeom, patch_num, PV_f);               \
    for (int PV_box = 0; PV_box < BoxArraySize(boxes); PV_box++)                    \
    {                                                                               \
      Box box = BoxArrayGetBox(boxes, PV_box);                                      \
      /* find octree and region intersection */                                     \
      int PV_ixl = pfmax(ix, box.lo[0]);                                            \
      int PV_iyl = pfmax(iy, box.lo[1]);                                            \
      int PV_izl = pfmax(iz, box.lo[2]);                                            \
      int PV_ixu = pfmin((ix + nx - 1), box.up[0]);                                 \
      int PV_iyu = pfmin((iy + ny - 1), box.up[1]);                                 \
      int PV_izu = pfmin((iz + nz - 1), box.up[2]);                                 \
                                                                                    \
      if(PV_ixl <= PV_ixu && PV_iyl <= PV_iyu && PV_izl <= PV_izu)                  \
      {                                                                             \
        int PV_diff_x = PV_ixu - PV_ixl;                                            \
        int PV_diff_y = PV_iyu - PV_iyl;                                            \
        int PV_diff_z = PV_izu - PV_izl;                                            \
                                                                                    \
        dim3 block, grid;                                                           \
        FindDims(grid, block, PV_diff_x + 1, PV_diff_y + 1, PV_diff_z + 1, 1);      \
                                                                                    \
        int nx = PV_diff_x + 1;                                                     \
        int ny = PV_diff_y + 1;                                                     \
        int nz = PV_diff_z + 1;                                                     \
                                                                                    \
        auto lambda_init =                                                          \
            GPU_LAMBDA(const int i, const int j, const int k)                       \
            {                                                                       \
                return n_prev + k * ny * nx + j * nx + i;                           \
            };                                                                      \
        n_prev += nz * ny * nx;                                                     \
                                                                                    \
        const int fdir_capt0 = fdir[0];                                             \
        const int fdir_capt1 = fdir[1];                                             \
        const int fdir_capt2 = fdir[2];                                             \
                                                                                    \
        auto lambda_body =                                                          \
             GPU_LAMBDA(const int i, const int j, const int k, int ival)            \
             {                                                                      \
               const int fdir[3] = {fdir_capt0, fdir_capt1, fdir_capt2};            \
               loop_body;                                                           \
               (void)fdir;                                                          \
             };                                                                     \
                                                                                    \
        BoxKernelI1<<<grid, block>>>(                                               \
            lambda_init, lambda_body, PV_ixl, PV_iyl, PV_izl, nx, ny, nz);          \
        CUDA_ERR(cudaPeekAtLastError());                                            \
        CUDA_ERR(cudaStreamSynchronize(0));                                         \
      }                                                                             \
    }                                                                               \
  }                                                                                 \
  (void)i;(void)j;(void)k;                                                          \
}

#undef GrGeomPatchLoopBoxesNoFdir
 /** Loop redefinition for CUDA. */
#define GrGeomPatchLoopBoxesNoFdir(i, j, k, grgeom, patch_num,                      \
  ix, iy, iz, nx, ny, nz, locals, setup,                                            \
  f_left, f_right, f_down, f_up, f_back, f_front, finalize)                         \
{                                                                                   \
  for (int PV_f = 0; PV_f < GrGeomOctreeNumFaces; PV_f++)                           \
  {                                                                                 \
    int n_prev = 0;                                                                 \
    BoxArray* boxes = GrGeomSolidPatchBoxes(grgeom, patch_num, PV_f);               \
    for (int PV_box = 0; PV_box < BoxArraySize(boxes); PV_box++)                    \
    {                                                                               \
      Box box = BoxArrayGetBox(boxes, PV_box);                                      \
      /* find octree and region intersection */                                     \
      int PV_ixl = pfmax(ix, box.lo[0]);                                            \
      int PV_iyl = pfmax(iy, box.lo[1]);                                            \
      int PV_izl = pfmax(iz, box.lo[2]);                                            \
      int PV_ixu = pfmin((ix + nx - 1), box.up[0]);                                 \
      int PV_iyu = pfmin((iy + ny - 1), box.up[1]);                                 \
      int PV_izu = pfmin((iz + nz - 1), box.up[2]);                                 \
                                                                                    \
      if(PV_ixl <= PV_ixu && PV_iyl <= PV_iyu && PV_izl <= PV_izu)                  \
      {                                                                             \
        int PV_diff_x = PV_ixu - PV_ixl;                                            \
        int PV_diff_y = PV_iyu - PV_iyl;                                            \
        int PV_diff_z = PV_izu - PV_izl;                                            \
                                                                                    \
        dim3 block, grid;                                                           \
        FindDims(grid, block, PV_diff_x + 1, PV_diff_y + 1, PV_diff_z + 1, 1);      \
                                                                                    \
        int nx = PV_diff_x + 1;                                                     \
        int ny = PV_diff_y + 1;                                                     \
        int nz = PV_diff_z + 1;                                                     \
                                                                                    \
        auto lambda_init =                                                          \
            GPU_LAMBDA(const int i, const int j, const int k)                       \
            {                                                                       \
                return n_prev + k * ny * nx + j * nx + i;                           \
            };                                                                      \
        n_prev += nz * ny * nx;                                                     \
                                                                                    \
        auto lambda_body =                                                          \
             GPU_LAMBDA(const int i, const int j, const int k, int ival)            \
             {                                                                      \
                UNPACK(locals);                                                     \
                setup;                                                              \
                switch(PV_f)                                                        \
                {                                                                   \
                  f_left;                                                           \
                  f_right;                                                          \
                  f_down;                                                           \
                  f_up;                                                             \
                  f_back;                                                           \
                  f_front;                                                          \
                }                                                                   \
                finalize;                                                           \
             };                                                                     \
                                                                                    \
        BoxKernelI1<<<grid, block>>>(                                               \
            lambda_init, lambda_body, PV_ixl, PV_iyl, PV_izl, nx, ny, nz);          \
        CUDA_ERR(cudaPeekAtLastError());                                            \
        CUDA_ERR(cudaStreamSynchronize(0));                                         \
      }                                                                             \
    }                                                                               \
  }                                                                                 \
  (void)i;(void)j;(void)k;                                                          \
}

#undef GrGeomOctreeExteriorNodeLoop
 /** Loop redefinition for CUDA. */
#define GrGeomOctreeExteriorNodeLoop(i, j, k, node, octree, level,                  \
  ix, iy, iz, nx, ny, nz, val_test, loop_body)                                      \
{                                                                                   \
  int PV_i, PV_j, PV_k, PV_l;                                                       \
  int PV_ixl, PV_iyl, PV_izl, PV_ixu, PV_iyu, PV_izu;                               \
                                                                                    \
  PV_i = i;                                                                         \
  PV_j = j;                                                                         \
  PV_k = k;                                                                         \
                                                                                    \
  GrGeomOctreeExteriorLoop(PV_i, PV_j, PV_k, PV_l, node, octree, level, val_test,   \
  {                                                                                 \
    if ((PV_i >= ix) && (PV_i < (ix + nx)) &&                                       \
        (PV_j >= iy) && (PV_j < (iy + ny)) &&                                       \
        (PV_k >= iz) && (PV_k < (iz + nz)))                                         \
    {                                                                               \
      i = PV_i;                                                                     \
      j = PV_j;                                                                     \
      k = PV_k;                                                                     \
      loop_body;                                                                    \
    }                                                                               \
  },                                                                                \
  {                                                                                 \
    /* find octree and region intersection */                                       \
    PV_ixl = pfmax(ix, PV_i);                                                       \
    PV_iyl = pfmax(iy, PV_j);                                                       \
    PV_izl = pfmax(iz, PV_k);                                                       \
    PV_ixu = pfmin((ix + nx), (PV_i + (int)PV_inc));                                \
    PV_iyu = pfmin((iy + ny), (PV_j + (int)PV_inc));                                \
    PV_izu = pfmin((iz + nz), (PV_k + (int)PV_inc));                                \
                                                                                    \
    if(PV_ixl < PV_ixu && PV_iyl < PV_iyu && PV_izl < PV_izu)                       \
    {                                                                               \
      const int PV_diff_x = PV_ixu - PV_ixl;                                        \
      const int PV_diff_y = PV_iyu - PV_iyl;                                        \
      const int PV_diff_z = PV_izu - PV_izl;                                        \
                                                                                    \
      dim3 block;                                                                   \
      dim3 grid;                                                                    \
      FindDims(grid, block, PV_diff_x, PV_diff_y, PV_diff_z, 1);                    \
                                                                                    \
      auto lambda_body = GPU_LAMBDA(const int i, const int j, const int k)loop_body;\
                                                                                    \
      (BoxKernelI0<<<grid, block>>>(lambda_body,                                    \
          PV_ixl, PV_iyl, PV_izl, PV_diff_x, PV_diff_y, PV_diff_z));                \
      CUDA_ERR(cudaPeekAtLastError());                                              \
      CUDA_ERR(cudaStreamSynchronize(0));                                           \
    }                                                                               \
    i = PV_ixu;                                                                     \
    j = PV_iyu;                                                                     \
    k = PV_izu;                                                                     \
  })                                                                                \
  (void)i;(void)j;(void)k;                                                          \
}

#undef GrGeomOutLoop
 /** Loop redefinition for CUDA. */
#define GrGeomOutLoop(i, j, k, grgeom, r,                                           \
  ix, iy, iz, nx, ny, nz, body)                                                     \
{                                                                                   \
  if(nx > 0 && ny > 0 && nz > 0)                                                    \
  {                                                                                 \
    if(!GrGeomSolidOutflag(grgeom))                                                 \
    {                                                                               \
      GrGeomOctree  *PV_node;                                                       \
      double PV_ref = pow(2.0, r);                                                  \
      unsigned outflag_size = sizeof(int) * (unsigned)(nz * ny * nx);               \
      int *outflag = (int*)ctalloc_cuda(outflag_size);                              \
                                                                                    \
      i = GrGeomSolidOctreeIX(grgeom) * (int)PV_ref;                                \
      j = GrGeomSolidOctreeIY(grgeom) * (int)PV_ref;                                \
      k = GrGeomSolidOctreeIZ(grgeom) * (int)PV_ref;                                \
      GrGeomOctreeExteriorNodeLoop(i, j, k, PV_node,                                \
                                   GrGeomSolidData(grgeom),                         \
                                   GrGeomSolidOctreeBGLevel(grgeom) + r,            \
                                   ix, iy, iz, nx, ny, nz,                          \
                                   TRUE,                                            \
      {                                                                             \
        body;                                                                       \
        outflag[(k - iz) * ny * nx + (j - iy) * nx + (i - ix)] = 1;                 \
      });                                                                           \
      GrGeomSolidOutflag(grgeom) = outflag;                                         \
    }                                                                               \
    else                                                                            \
    {                                                                               \
      dim3 block, grid;                                                             \
      FindDims(grid, block, nx, ny, nz, 1);                                         \
                                                                                    \
      int *outflag = GrGeomSolidOutflag(grgeom);                                    \
      auto lambda_body =                                                            \
        GPU_LAMBDA(const int i, const int j, const int k)                           \
        {                                                                           \
          if(outflag[(k - iz) * ny * nx + (j - iy) * nx + (i - ix)] == 1)           \
          {                                                                         \
            body;                                                                   \
          }                                                                         \
        };                                                                          \
                                                                                    \
      BoxKernelI0<<<grid, block>>>(                                                 \
          lambda_body, ix, iy, iz, nx, ny, nz);                                     \
      CUDA_ERR(cudaPeekAtLastError());                                              \
      CUDA_ERR(cudaStreamSynchronize(0));                                           \
    }                                                                               \
  }                                                                                 \
  (void)i;(void)j;(void)k;                                                          \
}

}
#endif // PFCUDALOOPS_H
