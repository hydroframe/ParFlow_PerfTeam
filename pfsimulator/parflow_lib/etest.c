typedef struct {
  SubgridArray *subgrids;
  GrGeomSolid *gr_domain;
  int num_patches;
  int *patch_indexes;
  int *bc_types;
  double ***values;
} BCStruct;
typedef struct {
  double x, y, z;
} GeomVertex;
typedef struct {
  int v0, v1, v2;
} GeomTriangle;
typedef struct {
  GeomVertex **vertices;
  int nV;
  int num_ptrs_to;
} GeomVertexArray;
typedef struct {
  GeomVertexArray *vertex_array;
  GeomTriangle **triangles;
  int nT;
} GeomTIN;
typedef struct {
  GeomTIN *surface;
  int **patches;
  int num_patches;
  int *num_patch_triangles;
} GeomTSolid;
typedef struct {
  void *data;
  int type;
  NameArray patches;
} GeomSolid;
typedef struct grgeom_octree {
  unsigned char flags;
  unsigned char faces;
  struct grgeom_octree *parent;
  struct grgeom_octree **children;
} GrGeomOctree;

typedef unsigned char __u_char;
typedef unsigned short int __u_short;
typedef unsigned int __u_int;
typedef unsigned long int __u_long;
typedef signed char __int8_t;
typedef unsigned char __uint8_t;
typedef signed short int __int16_t;
typedef unsigned short int __uint16_t;
typedef signed int __int32_t;
typedef unsigned int __uint32_t;
typedef signed long int __int64_t;
typedef unsigned long int __uint64_t;
typedef long int __quad_t;
typedef unsigned long int __u_quad_t;
typedef long int __intmax_t;
typedef unsigned long int __uintmax_t;
typedef unsigned long int __dev_t;
typedef unsigned int __uid_t;
typedef unsigned int __gid_t;
typedef unsigned long int __ino_t;
typedef unsigned long int __ino64_t;
typedef unsigned int __mode_t;
typedef unsigned long int __nlink_t;
typedef long int __off_t;
typedef long int __off64_t;
typedef int __pid_t;
typedef struct { int __val[2]; } __fsid_t;
typedef long int __clock_t;
typedef unsigned long int __rlim_t;
typedef unsigned long int __rlim64_t;
typedef unsigned int __id_t;
typedef long int __time_t;
typedef unsigned int __useconds_t;
typedef long int __suseconds_t;
typedef int __daddr_t;
typedef int __key_t;
typedef int __clockid_t;
typedef void * __timer_t;
typedef long int __blksize_t;
typedef long int __blkcnt_t;
typedef long int __blkcnt64_t;
typedef unsigned long int __fsblkcnt_t;
typedef unsigned long int __fsblkcnt64_t;
typedef unsigned long int __fsfilcnt_t;
typedef unsigned long int __fsfilcnt64_t;
typedef long int __fsword_t;
typedef long int __ssize_t;
typedef long int __syscall_slong_t;
typedef unsigned long int __syscall_ulong_t;
typedef __off64_t __loff_t;
typedef char *__caddr_t;
typedef long int __intptr_t;
typedef unsigned int __socklen_t;
typedef int __sig_atomic_t;
typedef __u_char u_char;
typedef __u_short u_short;
typedef __u_int u_int;
typedef __u_long u_long;
typedef __quad_t quad_t;
typedef __u_quad_t u_quad_t;
typedef __fsid_t fsid_t;
typedef __loff_t loff_t;
typedef __ino_t ino_t;
typedef __dev_t dev_t;
typedef __gid_t gid_t;
typedef __mode_t mode_t;
typedef __nlink_t nlink_t;
typedef __uid_t uid_t;
typedef __off_t off_t;
typedef __pid_t pid_t;
typedef __id_t id_t;
typedef __ssize_t ssize_t;
typedef __daddr_t daddr_t;
typedef __caddr_t caddr_t;
typedef __key_t key_t;
typedef __clock_t clock_t;
typedef __clockid_t clockid_t;
typedef __time_t time_t;
typedef __timer_t timer_t;
typedef long unsigned int size_t;
typedef unsigned long int ulong;
typedef unsigned short int ushort;
typedef unsigned int uint;
typedef __int8_t int8_t;
typedef __int16_t int16_t;
typedef __int32_t int32_t;
typedef __int64_t int64_t;
typedef unsigned int u_int8_t __attribute__ ((__mode__ (__QI__)));
typedef unsigned int u_int16_t __attribute__ ((__mode__ (__HI__)));
typedef unsigned int u_int32_t __attribute__ ((__mode__ (__SI__)));
typedef unsigned int u_int64_t __attribute__ ((__mode__ (__DI__)));
typedef int register_t __attribute__ ((__mode__ (__word__)));
static __inline unsigned int
__bswap_32 (unsigned int __bsx)
{
  return __builtin_bswap32 (__bsx);
}
static __inline __uint64_t
__bswap_64 (__uint64_t __bsx)
{
  return __builtin_bswap64 (__bsx);
}
static __inline __uint16_t
__uint16_identity (__uint16_t __x)
{
  return __x;
}
static __inline __uint32_t
__uint32_identity (__uint32_t __x)
{
  return __x;
}
static __inline __uint64_t
__uint64_identity (__uint64_t __x)
{
  return __x;
}
typedef struct
{
  unsigned long int __val[(1024 / (8 * sizeof (unsigned long int)))];
} __sigset_t;
typedef __sigset_t sigset_t;
struct timeval
{
  __time_t tv_sec;
  __suseconds_t tv_usec;
};
struct timespec
{
  __time_t tv_sec;
  __syscall_slong_t tv_nsec;
};
typedef __suseconds_t suseconds_t;
typedef long int __fd_mask;
typedef struct
  {
    __fd_mask __fds_bits[1024 / (8 * (int) sizeof (__fd_mask))];
  } fd_set;
typedef __fd_mask fd_mask;

extern int select (int __nfds, fd_set *__restrict __readfds,
     fd_set *__restrict __writefds,
     fd_set *__restrict __exceptfds,
     struct timeval *__restrict __timeout);
extern int pselect (int __nfds, fd_set *__restrict __readfds,
      fd_set *__restrict __writefds,
      fd_set *__restrict __exceptfds,
      const struct timespec *__restrict __timeout,
      const __sigset_t *__restrict __sigmask);


extern unsigned int gnu_dev_major (__dev_t __dev) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern unsigned int gnu_dev_minor (__dev_t __dev) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern __dev_t gnu_dev_makedev (unsigned int __major, unsigned int __minor) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

typedef __blksize_t blksize_t;
typedef __blkcnt_t blkcnt_t;
typedef __fsblkcnt_t fsblkcnt_t;
typedef __fsfilcnt_t fsfilcnt_t;
struct __pthread_rwlock_arch_t
{
  unsigned int __readers;
  unsigned int __writers;
  unsigned int __wrphase_futex;
  unsigned int __writers_futex;
  unsigned int __pad3;
  unsigned int __pad4;
  int __cur_writer;
  int __shared;
  signed char __rwelision;
  unsigned char __pad1[7];
  unsigned long int __pad2;
  unsigned int __flags;
};
typedef struct __pthread_internal_list
{
  struct __pthread_internal_list *__prev;
  struct __pthread_internal_list *__next;
} __pthread_list_t;
struct __pthread_mutex_s
{
  int __lock ;
  unsigned int __count;
  int __owner;
  unsigned int __nusers;
  int __kind;
 
  short __spins; short __elision;
  __pthread_list_t __list;
 
};
struct __pthread_cond_s
{
  __extension__ union
  {
    __extension__ unsigned long long int __wseq;
    struct
    {
      unsigned int __low;
      unsigned int __high;
    } __wseq32;
  };
  __extension__ union
  {
    __extension__ unsigned long long int __g1_start;
    struct
    {
      unsigned int __low;
      unsigned int __high;
    } __g1_start32;
  };
  unsigned int __g_refs[2] ;
  unsigned int __g_size[2];
  unsigned int __g1_orig_size;
  unsigned int __wrefs;
  unsigned int __g_signals[2];
};
typedef unsigned long int pthread_t;
typedef union
{
  char __size[4];
  int __align;
} pthread_mutexattr_t;
typedef union
{
  char __size[4];
  int __align;
} pthread_condattr_t;
typedef unsigned int pthread_key_t;
typedef int pthread_once_t;
union pthread_attr_t
{
  char __size[56];
  long int __align;
};
typedef union pthread_attr_t pthread_attr_t;
typedef union
{
  struct __pthread_mutex_s __data;
  char __size[40];
  long int __align;
} pthread_mutex_t;
typedef union
{
  struct __pthread_cond_s __data;
  char __size[48];
  __extension__ long long int __align;
} pthread_cond_t;
typedef union
{
  struct __pthread_rwlock_arch_t __data;
  char __size[56];
  long int __align;
} pthread_rwlock_t;
typedef union
{
  char __size[8];
  long int __align;
} pthread_rwlockattr_t;
typedef volatile int pthread_spinlock_t;
typedef union
{
  char __size[32];
  long int __align;
} pthread_barrier_t;
typedef union
{
  char __size[4];
  int __align;
} pthread_barrierattr_t;


typedef __sig_atomic_t sig_atomic_t;
union sigval
{
  int sival_int;
  void *sival_ptr;
};
typedef union sigval __sigval_t;
typedef struct
  {
    int si_signo;
    int si_errno;
    int si_code;
    int __pad0;
    union
      {
 int _pad[((128 / sizeof (int)) - 4)];
 struct
   {
     __pid_t si_pid;
     __uid_t si_uid;
   } _kill;
 struct
   {
     int si_tid;
     int si_overrun;
     __sigval_t si_sigval;
   } _timer;
 struct
   {
     __pid_t si_pid;
     __uid_t si_uid;
     __sigval_t si_sigval;
   } _rt;
 struct
   {
     __pid_t si_pid;
     __uid_t si_uid;
     int si_status;
     __clock_t si_utime;
     __clock_t si_stime;
   } _sigchld;
 struct
   {
     void *si_addr;
    
     short int si_addr_lsb;
     union
       {
  struct
    {
      void *_lower;
      void *_upper;
    } _addr_bnd;
  __uint32_t _pkey;
       } _bounds;
   } _sigfault;
 struct
   {
     long int si_band;
     int si_fd;
   } _sigpoll;
 struct
   {
     void *_call_addr;
     int _syscall;
     unsigned int _arch;
   } _sigsys;
      } _sifields;
  } siginfo_t ;
enum
{
  SI_ASYNCNL = -60,
  SI_TKILL = -6,
  SI_SIGIO,
  SI_ASYNCIO,
  SI_MESGQ,
  SI_TIMER,
  SI_QUEUE,
  SI_USER,
  SI_KERNEL = 0x80
};
enum
{
  ILL_ILLOPC = 1,
  ILL_ILLOPN,
  ILL_ILLADR,
  ILL_ILLTRP,
  ILL_PRVOPC,
  ILL_PRVREG,
  ILL_COPROC,
  ILL_BADSTK
};
enum
{
  FPE_INTDIV = 1,
  FPE_INTOVF,
  FPE_FLTDIV,
  FPE_FLTOVF,
  FPE_FLTUND,
  FPE_FLTRES,
  FPE_FLTINV,
  FPE_FLTSUB
};
enum
{
  SEGV_MAPERR = 1,
  SEGV_ACCERR,
  SEGV_BNDERR,
  SEGV_PKUERR
};
enum
{
  BUS_ADRALN = 1,
  BUS_ADRERR,
  BUS_OBJERR,
  BUS_MCEERR_AR,
  BUS_MCEERR_AO
};
enum
{
  CLD_EXITED = 1,
  CLD_KILLED,
  CLD_DUMPED,
  CLD_TRAPPED,
  CLD_STOPPED,
  CLD_CONTINUED
};
enum
{
  POLL_IN = 1,
  POLL_OUT,
  POLL_MSG,
  POLL_ERR,
  POLL_PRI,
  POLL_HUP
};
typedef __sigval_t sigval_t;
typedef struct sigevent
  {
    __sigval_t sigev_value;
    int sigev_signo;
    int sigev_notify;
    union
      {
 int _pad[((64 / sizeof (int)) - 4)];
 __pid_t _tid;
 struct
   {
     void (*_function) (__sigval_t);
     pthread_attr_t *_attribute;
   } _sigev_thread;
      } _sigev_un;
  } sigevent_t;
enum
{
  SIGEV_SIGNAL = 0,
  SIGEV_NONE,
  SIGEV_THREAD,
  SIGEV_THREAD_ID = 4
};
typedef void (*__sighandler_t) (int);
extern __sighandler_t __sysv_signal (int __sig, __sighandler_t __handler)
     __attribute__ ((__nothrow__ , __leaf__));
extern __sighandler_t signal (int __sig, __sighandler_t __handler)
     __attribute__ ((__nothrow__ , __leaf__));
extern int kill (__pid_t __pid, int __sig) __attribute__ ((__nothrow__ , __leaf__));
extern int killpg (__pid_t __pgrp, int __sig) __attribute__ ((__nothrow__ , __leaf__));
extern int raise (int __sig) __attribute__ ((__nothrow__ , __leaf__));
extern __sighandler_t ssignal (int __sig, __sighandler_t __handler)
     __attribute__ ((__nothrow__ , __leaf__));
extern int gsignal (int __sig) __attribute__ ((__nothrow__ , __leaf__));
extern void psignal (int __sig, const char *__s);
extern void psiginfo (const siginfo_t *__pinfo, const char *__s);
extern int sigblock (int __mask) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__deprecated__));
extern int sigsetmask (int __mask) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__deprecated__));
extern int siggetmask (void) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__deprecated__));
typedef __sighandler_t sig_t;
extern int sigemptyset (sigset_t *__set) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));
extern int sigfillset (sigset_t *__set) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));
extern int sigaddset (sigset_t *__set, int __signo) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));
extern int sigdelset (sigset_t *__set, int __signo) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));
extern int sigismember (const sigset_t *__set, int __signo)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));
struct sigaction
  {
    union
      {
 __sighandler_t sa_handler;
 void (*sa_sigaction) (int, siginfo_t *, void *);
      }
    __sigaction_handler;
    __sigset_t sa_mask;
    int sa_flags;
    void (*sa_restorer) (void);
  };
extern int sigprocmask (int __how, const sigset_t *__restrict __set,
   sigset_t *__restrict __oset) __attribute__ ((__nothrow__ , __leaf__));
extern int sigsuspend (const sigset_t *__set) __attribute__ ((__nonnull__ (1)));
extern int sigaction (int __sig, const struct sigaction *__restrict __act,
        struct sigaction *__restrict __oact) __attribute__ ((__nothrow__ , __leaf__));
extern int sigpending (sigset_t *__set) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));
extern int sigwait (const sigset_t *__restrict __set, int *__restrict __sig)
     __attribute__ ((__nonnull__ (1, 2)));
extern int sigwaitinfo (const sigset_t *__restrict __set,
   siginfo_t *__restrict __info) __attribute__ ((__nonnull__ (1)));
extern int sigtimedwait (const sigset_t *__restrict __set,
    siginfo_t *__restrict __info,
    const struct timespec *__restrict __timeout)
     __attribute__ ((__nonnull__ (1)));
extern int sigqueue (__pid_t __pid, int __sig, const union sigval __val)
     __attribute__ ((__nothrow__ , __leaf__));
extern const char *const _sys_siglist[(64 + 1)];
extern const char *const sys_siglist[(64 + 1)];
struct _fpx_sw_bytes
{
  __uint32_t magic1;
  __uint32_t extended_size;
  __uint64_t xstate_bv;
  __uint32_t xstate_size;
  __uint32_t __glibc_reserved1[7];
};
struct _fpreg
{
  unsigned short significand[4];
  unsigned short exponent;
};
struct _fpxreg
{
  unsigned short significand[4];
  unsigned short exponent;
  unsigned short __glibc_reserved1[3];
};
struct _xmmreg
{
  __uint32_t element[4];
};
struct _fpstate
{
  __uint16_t cwd;
  __uint16_t swd;
  __uint16_t ftw;
  __uint16_t fop;
  __uint64_t rip;
  __uint64_t rdp;
  __uint32_t mxcsr;
  __uint32_t mxcr_mask;
  struct _fpxreg _st[8];
  struct _xmmreg _xmm[16];
  __uint32_t __glibc_reserved1[24];
};
struct sigcontext
{
  __uint64_t r8;
  __uint64_t r9;
  __uint64_t r10;
  __uint64_t r11;
  __uint64_t r12;
  __uint64_t r13;
  __uint64_t r14;
  __uint64_t r15;
  __uint64_t rdi;
  __uint64_t rsi;
  __uint64_t rbp;
  __uint64_t rbx;
  __uint64_t rdx;
  __uint64_t rax;
  __uint64_t rcx;
  __uint64_t rsp;
  __uint64_t rip;
  __uint64_t eflags;
  unsigned short cs;
  unsigned short gs;
  unsigned short fs;
  unsigned short __pad0;
  __uint64_t err;
  __uint64_t trapno;
  __uint64_t oldmask;
  __uint64_t cr2;
  __extension__ union
    {
      struct _fpstate * fpstate;
      __uint64_t __fpstate_word;
    };
  __uint64_t __reserved1 [8];
};
struct _xsave_hdr
{
  __uint64_t xstate_bv;
  __uint64_t __glibc_reserved1[2];
  __uint64_t __glibc_reserved2[5];
};
struct _ymmh_state
{
  __uint32_t ymmh_space[64];
};
struct _xstate
{
  struct _fpstate fpstate;
  struct _xsave_hdr xstate_hdr;
  struct _ymmh_state ymmh;
};
extern int sigreturn (struct sigcontext *__scp) __attribute__ ((__nothrow__ , __leaf__));
typedef struct
  {
    void *ss_sp;
    int ss_flags;
    size_t ss_size;
  } stack_t;
__extension__ typedef long long int greg_t;
typedef greg_t gregset_t[23];
struct _libc_fpxreg
{
  unsigned short int significand[4];
  unsigned short int exponent;
  unsigned short int __glibc_reserved1[3];
};
struct _libc_xmmreg
{
  __uint32_t element[4];
};
struct _libc_fpstate
{
  __uint16_t cwd;
  __uint16_t swd;
  __uint16_t ftw;
  __uint16_t fop;
  __uint64_t rip;
  __uint64_t rdp;
  __uint32_t mxcsr;
  __uint32_t mxcr_mask;
  struct _libc_fpxreg _st[8];
  struct _libc_xmmreg _xmm[16];
  __uint32_t __glibc_reserved1[24];
};
typedef struct _libc_fpstate *fpregset_t;
typedef struct
  {
    gregset_t gregs;
    fpregset_t fpregs;
    __extension__ unsigned long long __reserved1 [8];
} mcontext_t;
typedef struct ucontext_t
  {
    unsigned long int uc_flags;
    struct ucontext_t *uc_link;
    stack_t uc_stack;
    mcontext_t uc_mcontext;
    sigset_t uc_sigmask;
    struct _libc_fpstate __fpregs_mem;
  } ucontext_t;
extern int siginterrupt (int __sig, int __interrupt) __attribute__ ((__nothrow__ , __leaf__));
enum
{
  SS_ONSTACK = 1,
  SS_DISABLE
};
extern int sigaltstack (const stack_t *__restrict __ss,
   stack_t *__restrict __oss) __attribute__ ((__nothrow__ , __leaf__));
struct sigstack
  {
    void *ss_sp;
    int ss_onstack;
  };
extern int sigstack (struct sigstack *__ss, struct sigstack *__oss)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__deprecated__));
extern int pthread_sigmask (int __how,
       const __sigset_t *__restrict __newmask,
       __sigset_t *__restrict __oldmask)__attribute__ ((__nothrow__ , __leaf__));
extern int pthread_kill (pthread_t __threadid, int __signo) __attribute__ ((__nothrow__ , __leaf__));
extern int __libc_current_sigrtmin (void) __attribute__ ((__nothrow__ , __leaf__));
extern int __libc_current_sigrtmax (void) __attribute__ ((__nothrow__ , __leaf__));


struct _IO_FILE;
typedef struct _IO_FILE __FILE;
struct _IO_FILE;
typedef struct _IO_FILE FILE;
typedef struct
{
  int __count;
  union
  {
    unsigned int __wch;
    char __wchb[4];
  } __value;
} __mbstate_t;
typedef struct
{
  __off_t __pos;
  __mbstate_t __state;
} _G_fpos_t;
typedef struct
{
  __off64_t __pos;
  __mbstate_t __state;
} _G_fpos64_t;
typedef __builtin_va_list __gnuc_va_list;
struct _IO_jump_t; struct _IO_FILE;
typedef void _IO_lock_t;
struct _IO_marker {
  struct _IO_marker *_next;
  struct _IO_FILE *_sbuf;
  int _pos;
};
enum __codecvt_result
{
  __codecvt_ok,
  __codecvt_partial,
  __codecvt_error,
  __codecvt_noconv
};
struct _IO_FILE {
  int _flags;
  char* _IO_read_ptr;
  char* _IO_read_end;
  char* _IO_read_base;
  char* _IO_write_base;
  char* _IO_write_ptr;
  char* _IO_write_end;
  char* _IO_buf_base;
  char* _IO_buf_end;
  char *_IO_save_base;
  char *_IO_backup_base;
  char *_IO_save_end;
  struct _IO_marker *_markers;
  struct _IO_FILE *_chain;
  int _fileno;
  int _flags2;
  __off_t _old_offset;
  unsigned short _cur_column;
  signed char _vtable_offset;
  char _shortbuf[1];
  _IO_lock_t *_lock;
  __off64_t _offset;
  void *__pad1;
  void *__pad2;
  void *__pad3;
  void *__pad4;
  size_t __pad5;
  int _mode;
  char _unused2[15 * sizeof (int) - 4 * sizeof (void *) - sizeof (size_t)];
};
typedef struct _IO_FILE _IO_FILE;
struct _IO_FILE_plus;
extern struct _IO_FILE_plus _IO_2_1_stdin_;
extern struct _IO_FILE_plus _IO_2_1_stdout_;
extern struct _IO_FILE_plus _IO_2_1_stderr_;
typedef __ssize_t __io_read_fn (void *__cookie, char *__buf, size_t __nbytes);
typedef __ssize_t __io_write_fn (void *__cookie, const char *__buf,
     size_t __n);
typedef int __io_seek_fn (void *__cookie, __off64_t *__pos, int __w);
typedef int __io_close_fn (void *__cookie);
extern int __underflow (_IO_FILE *);
extern int __uflow (_IO_FILE *);
extern int __overflow (_IO_FILE *, int);
extern int _IO_getc (_IO_FILE *__fp);
extern int _IO_putc (int __c, _IO_FILE *__fp);
extern int _IO_feof (_IO_FILE *__fp) __attribute__ ((__nothrow__ , __leaf__));
extern int _IO_ferror (_IO_FILE *__fp) __attribute__ ((__nothrow__ , __leaf__));
extern int _IO_peekc_locked (_IO_FILE *__fp);
extern void _IO_flockfile (_IO_FILE *) __attribute__ ((__nothrow__ , __leaf__));
extern void _IO_funlockfile (_IO_FILE *) __attribute__ ((__nothrow__ , __leaf__));
extern int _IO_ftrylockfile (_IO_FILE *) __attribute__ ((__nothrow__ , __leaf__));
extern int _IO_vfscanf (_IO_FILE * __restrict, const char * __restrict,
   __gnuc_va_list, int *__restrict);
extern int _IO_vfprintf (_IO_FILE *__restrict, const char *__restrict,
    __gnuc_va_list);
extern __ssize_t _IO_padn (_IO_FILE *, int, __ssize_t);
extern size_t _IO_sgetn (_IO_FILE *, void *, size_t);
extern __off64_t _IO_seekoff (_IO_FILE *, __off64_t, int, int);
extern __off64_t _IO_seekpos (_IO_FILE *, __off64_t, int);
extern void _IO_free_backup_area (_IO_FILE *) __attribute__ ((__nothrow__ , __leaf__));
typedef __gnuc_va_list va_list;
typedef _G_fpos_t fpos_t;
extern struct _IO_FILE *stdin;
extern struct _IO_FILE *stdout;
extern struct _IO_FILE *stderr;
extern int remove (const char *__filename) __attribute__ ((__nothrow__ , __leaf__));
extern int rename (const char *__old, const char *__new) __attribute__ ((__nothrow__ , __leaf__));
extern int renameat (int __oldfd, const char *__old, int __newfd,
       const char *__new) __attribute__ ((__nothrow__ , __leaf__));
extern FILE *tmpfile (void) ;
extern char *tmpnam (char *__s) __attribute__ ((__nothrow__ , __leaf__)) ;
extern char *tmpnam_r (char *__s) __attribute__ ((__nothrow__ , __leaf__)) ;
extern char *tempnam (const char *__dir, const char *__pfx)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__malloc__)) ;
extern int fclose (FILE *__stream);
extern int fflush (FILE *__stream);
extern int fflush_unlocked (FILE *__stream);
extern FILE *fopen (const char *__restrict __filename,
      const char *__restrict __modes) ;
extern FILE *freopen (const char *__restrict __filename,
        const char *__restrict __modes,
        FILE *__restrict __stream) ;
extern FILE *fdopen (int __fd, const char *__modes) __attribute__ ((__nothrow__ , __leaf__)) ;
extern FILE *fmemopen (void *__s, size_t __len, const char *__modes)
  __attribute__ ((__nothrow__ , __leaf__)) ;
extern FILE *open_memstream (char **__bufloc, size_t *__sizeloc) __attribute__ ((__nothrow__ , __leaf__)) ;
extern void setbuf (FILE *__restrict __stream, char *__restrict __buf) __attribute__ ((__nothrow__ , __leaf__));
extern int setvbuf (FILE *__restrict __stream, char *__restrict __buf,
      int __modes, size_t __n) __attribute__ ((__nothrow__ , __leaf__));
extern void setbuffer (FILE *__restrict __stream, char *__restrict __buf,
         size_t __size) __attribute__ ((__nothrow__ , __leaf__));
extern void setlinebuf (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__));
extern int fprintf (FILE *__restrict __stream,
      const char *__restrict __format, ...);
extern int printf (const char *__restrict __format, ...);
extern int sprintf (char *__restrict __s,
      const char *__restrict __format, ...) __attribute__ ((__nothrow__));
extern int vfprintf (FILE *__restrict __s, const char *__restrict __format,
       __gnuc_va_list __arg);
extern int vprintf (const char *__restrict __format, __gnuc_va_list __arg);
extern int vsprintf (char *__restrict __s, const char *__restrict __format,
       __gnuc_va_list __arg) __attribute__ ((__nothrow__));
extern int snprintf (char *__restrict __s, size_t __maxlen,
       const char *__restrict __format, ...)
     __attribute__ ((__nothrow__)) __attribute__ ((__format__ (__printf__, 3, 4)));
extern int vsnprintf (char *__restrict __s, size_t __maxlen,
        const char *__restrict __format, __gnuc_va_list __arg)
     __attribute__ ((__nothrow__)) __attribute__ ((__format__ (__printf__, 3, 0)));
extern int vdprintf (int __fd, const char *__restrict __fmt,
       __gnuc_va_list __arg)
     __attribute__ ((__format__ (__printf__, 2, 0)));
extern int dprintf (int __fd, const char *__restrict __fmt, ...)
     __attribute__ ((__format__ (__printf__, 2, 3)));
extern int fscanf (FILE *__restrict __stream,
     const char *__restrict __format, ...) ;
extern int scanf (const char *__restrict __format, ...) ;
extern int sscanf (const char *__restrict __s,
     const char *__restrict __format, ...) __attribute__ ((__nothrow__ , __leaf__));
extern int fscanf (FILE *__restrict __stream, const char *__restrict __format, ...) __asm__ ("" "__isoc99_fscanf") ;
extern int scanf (const char *__restrict __format, ...) __asm__ ("" "__isoc99_scanf") ;
extern int sscanf (const char *__restrict __s, const char *__restrict __format, ...) __asm__ ("" "__isoc99_sscanf") __attribute__ ((__nothrow__ , __leaf__));
extern int vfscanf (FILE *__restrict __s, const char *__restrict __format,
      __gnuc_va_list __arg)
     __attribute__ ((__format__ (__scanf__, 2, 0))) ;
extern int vscanf (const char *__restrict __format, __gnuc_va_list __arg)
     __attribute__ ((__format__ (__scanf__, 1, 0))) ;
extern int vsscanf (const char *__restrict __s,
      const char *__restrict __format, __gnuc_va_list __arg)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__format__ (__scanf__, 2, 0)));
extern int vfscanf (FILE *__restrict __s, const char *__restrict __format, __gnuc_va_list __arg) __asm__ ("" "__isoc99_vfscanf")
     __attribute__ ((__format__ (__scanf__, 2, 0))) ;
extern int vscanf (const char *__restrict __format, __gnuc_va_list __arg) __asm__ ("" "__isoc99_vscanf")
     __attribute__ ((__format__ (__scanf__, 1, 0))) ;
extern int vsscanf (const char *__restrict __s, const char *__restrict __format, __gnuc_va_list __arg) __asm__ ("" "__isoc99_vsscanf") __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__format__ (__scanf__, 2, 0)));
extern int fgetc (FILE *__stream);
extern int getc (FILE *__stream);
extern int getchar (void);
extern int getc_unlocked (FILE *__stream);
extern int getchar_unlocked (void);
extern int fgetc_unlocked (FILE *__stream);
extern int fputc (int __c, FILE *__stream);
extern int putc (int __c, FILE *__stream);
extern int putchar (int __c);
extern int fputc_unlocked (int __c, FILE *__stream);
extern int putc_unlocked (int __c, FILE *__stream);
extern int putchar_unlocked (int __c);
extern int getw (FILE *__stream);
extern int putw (int __w, FILE *__stream);
extern char *fgets (char *__restrict __s, int __n, FILE *__restrict __stream)
     ;
extern __ssize_t __getdelim (char **__restrict __lineptr,
          size_t *__restrict __n, int __delimiter,
          FILE *__restrict __stream) ;
extern __ssize_t getdelim (char **__restrict __lineptr,
        size_t *__restrict __n, int __delimiter,
        FILE *__restrict __stream) ;
extern __ssize_t getline (char **__restrict __lineptr,
       size_t *__restrict __n,
       FILE *__restrict __stream) ;
extern int fputs (const char *__restrict __s, FILE *__restrict __stream);
extern int puts (const char *__s);
extern int ungetc (int __c, FILE *__stream);
extern size_t fread (void *__restrict __ptr, size_t __size,
       size_t __n, FILE *__restrict __stream) ;
extern size_t fwrite (const void *__restrict __ptr, size_t __size,
        size_t __n, FILE *__restrict __s);
extern size_t fread_unlocked (void *__restrict __ptr, size_t __size,
         size_t __n, FILE *__restrict __stream) ;
extern size_t fwrite_unlocked (const void *__restrict __ptr, size_t __size,
          size_t __n, FILE *__restrict __stream);
extern int fseek (FILE *__stream, long int __off, int __whence);
extern long int ftell (FILE *__stream) ;
extern void rewind (FILE *__stream);
extern int fseeko (FILE *__stream, __off_t __off, int __whence);
extern __off_t ftello (FILE *__stream) ;
extern int fgetpos (FILE *__restrict __stream, fpos_t *__restrict __pos);
extern int fsetpos (FILE *__stream, const fpos_t *__pos);
extern void clearerr (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__));
extern int feof (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;
extern int ferror (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;
extern void clearerr_unlocked (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__));
extern int feof_unlocked (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;
extern int ferror_unlocked (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;
extern void perror (const char *__s);
extern int sys_nerr;
extern const char *const sys_errlist[];
extern int fileno (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;
extern int fileno_unlocked (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;
extern FILE *popen (const char *__command, const char *__modes) ;
extern int pclose (FILE *__stream);
extern char *ctermid (char *__s) __attribute__ ((__nothrow__ , __leaf__));
extern void flockfile (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__));
extern int ftrylockfile (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;
extern void funlockfile (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__));

typedef struct _HBT_element {
  void *obj;
  short balance;
  struct _HBT_element *left;
  struct _HBT_element *right;
} HBT_element;
typedef struct _HBT {
  unsigned int height;
  unsigned int num;
  HBT_element *root;
  int (*compare)(void *, void *);
  void (*free)(void *);
  void (*printf)(FILE *file, void *);
  int (*scanf)(FILE *file, void **);
  int malloc_flag;
} HBT;
typedef struct _IDB_Entry {
  char *key;
  char *value;
  char used;
} IDB_Entry;
typedef HBT IDB;
typedef struct NameArray__ {
  int num;
  char **names;
  char *tok_string;
  char *string;
} NameArrayStruct;
typedef NameArrayStruct *NameArray;
typedef struct {
  double X, Y, Z;
  double DX, DY, DZ;
  int IX, IY, IZ;
  int NX, NY, NZ;
} Background;
typedef struct {
  int num_send_invoices;
  int *send_ranks;
  amps_Invoice *send_invoices;
  int num_recv_invoices;
  int *recv_ranks;
  amps_Invoice *recv_invoices;
  amps_Package package;
  int *loop_array;
} CommPkg;
typedef amps_Handle CommHandle;
typedef struct {
  int ix, iy, iz;
  int nx, ny, nz;
  int sx, sy, sz;
  int rx, ry, rz;
  int level;
  int process;
} Subregion;
typedef struct {
  Subregion **subregions;
  int size;
} SubregionArray;
typedef struct {
  SubregionArray **subregion_arrays;
  int size;
} Region;
typedef struct {
  Region *send_region;
  Region *recv_region;
  Region *ind_region;
  Region *dep_region;
} ComputePkg;
typedef Subregion Subgrid;
typedef SubregionArray SubgridArray;
typedef struct {
  SubgridArray *subgrids;
  SubgridArray *all_subgrids;
  int size;
  ComputePkg **compute_pkgs;
  Subgrid *background;
} Grid;
enum matrix_type {
  matrix_cell_centered,
  matrix_non_samrai
};
typedef int StencilElt[3];
typedef struct {
  StencilElt *shape;
  int size;
} Stencil;
typedef struct {
  double* data;
  int allocated;
  int* data_index;
  int data_size;
  Subregion* data_space;
} Submatrix;
typedef struct {
  Submatrix **submatrices;
  Grid *grid;
  SubregionArray *range;
  SubregionArray *data_space;
  Stencil *stencil;
  int *data_stencil;
  int data_stencil_size;
  int symmetric;
  int size;
  CommPkg *comm_pkg;
  enum matrix_type type;
} Matrix;
enum vector_type {
  vector_cell_centered,
  vector_cell_centered_2D,
  vector_side_centered_x,
  vector_side_centered_y,
  vector_side_centered_z,
  vector_clm_topsoil,
  vector_met,
  vector_non_samrai
};
typedef struct {
  double *data;
  int allocated;
  Subgrid *data_space;
  int data_size;
} Subvector;
typedef struct _Vector {
  Subvector **subvectors;
  int data_size;
  Grid *grid;
  SubgridArray *data_space;
  int size;
  CommPkg *comm_pkg[10];
  enum vector_type type;
} Vector;
typedef Vector *N_Vector;
typedef struct _VectorUpdateCommHandle {
  Vector *vector;
  CommHandle *comm_handle;
} VectorUpdateCommHandle;
typedef struct {
  void (*call)();
  void (*init_instance_xtra)();
  void (*free_instance_xtra)();
  void (*new_public_xtra)();
  void (*free_public_xtra)();
  int (*sizeof_temp_data)();
  void *instance_xtra;
  void *public_xtra;
} PFModule;
amps_ThreadLocalDcl(extern PFModule *, global_ptr_this_pf_module);
typedef struct {
  int nc;
  double *x, *y, *z;
  double *v;
} RFCondData;
typedef struct {
  double mean, sigma;
  double lambdaX, lambdaY, lambdaZ;
  int lognormal;
} Statistics;
typedef struct {
  char *data;
  int *data_index;
  Subgrid *data_space;
  int nc;
} Subcharvector;
typedef struct _CharVector {
  Subcharvector **subcharvectors;
  char *data;
  int data_size;
  Grid *grid;
  SubgridArray *data_space;
  int nc;
  int size;
  CommPkg *comm_pkg[10];
} CharVector;
typedef struct {
  Grid *grid;
  double **e;
  double *c;
  Vector *pressure;
  Vector *perm;
  Vector *phi;
  double **Ktensor;
  double *Kscale;
  CharVector *cellType;
  double beta_perm;
  double beta_pore;
  double beta_fracture;
  double beta_fluid;
  double viscosity;
  double density;
  Vector *pwork;
  double *bforce;
  double t;
  double start;
  double stop;
  double step;
  double dump;
  double tscale;
  double cfl;
  int comp_compress_flag;
} Lattice;
enum ParflowGridType {
      invalid_grid_type,
      flow_3D_grid_type,
      surface_2D_grid_type,
      met_2D_grid_type,
      vector_clm_topsoil_grid_type
};
typedef struct _Globals {
  char run_name[256];
  char in_file_name[256];
  char out_file_name[256];
  int logging_level;
  int num_procs;
  int num_procs_x;
  int num_procs_y;
  int num_procs_z;
  int p;
  int q;
  int r;
  Background *background;
  Grid *user_grid;
  int max_ref_level;
  NameArray geom_names;
  GeomSolid **geometries;
  NameArray phase_names;
  NameArray contaminant_names;
  NameArray cycle_names;
  int num_cycles;
  NameArray *interval_names;
  int *interval_divisions;
  int **intervals;
  int *repeat_counts;
  Grid *grid3d;
  Grid *grid2d;
  int use_clustering;
} Globals;
amps_ThreadLocalDcl(extern Globals *, globals_ptr);
amps_ThreadLocalDcl(extern IDB *, input_database);
typedef struct {
  int number_of_cycles;
  int *interval_divisions;
  int **intervals;
  int *repeat_counts;
  int *cycle_lengths;
} TimeCycleData;
typedef struct {
  int num_f_points, num_h_points;
  double *f_points, *f_values;
  double *h_points, *h_values;
} EvalStruct;
typedef struct {
  int number;
  char *name;
  double x_lower, y_lower, z_lower;
  double x_upper, y_upper, z_upper;
  double diameter;
  Subgrid *subgrid;
  double size;
  int action;
  int method;
  int cycle_number;
  double average_permeability_x;
  double average_permeability_y;
  double average_permeability_z;
} WellDataPhysical;
typedef struct {
  double *phase_values;
  double *saturation_values;
  double *delta_saturation_ptrs;
  double *contaminant_values;
  double *delta_contaminant_ptrs;
  double *contaminant_fractions;
} WellDataValue;
typedef struct {
  double *delta_phases;
  double *phase_stats;
  double *delta_saturations;
  double *saturation_stats;
  double *delta_contaminants;
  double *contaminant_stats;
} WellDataStat;
typedef struct {
  int num_phases;
  int num_contaminants;
  int num_wells;
  int num_press_wells;
  WellDataPhysical **press_well_physicals;
  WellDataValue ***press_well_values;
  WellDataStat **press_well_stats;
  int num_flux_wells;
  WellDataPhysical **flux_well_physicals;
  WellDataValue ***flux_well_values;
  WellDataStat **flux_well_stats;
  TimeCycleData *time_cycle_data;
} WellData;
typedef struct { int reference_solid; int reference_patch; double value; double *value_at_interfaces; } BCPressureTypeDirEquilRefPatch; typedef struct { double xlower; double ylower; double xupper; double yupper; int num_points; double *points; double *values; double *value_at_interfaces; } BCPressureTypeDirEquilPLinear; typedef struct { double value; } BCPressureTypeFluxConst; typedef struct { double value; } BCPressureTypeFluxVolumetric; typedef struct { char *filename; } BCPressureTypePressureFile; typedef struct { char *filename; } BCPressureTypeFluxFile; typedef struct { int function_type; } BCPressureTypeExactSolution; typedef struct { double value; } BCPressureTypeOverlandFlow; typedef struct { char *filename; } BCPressureTypeOverlandFlowPFB; typedef struct { double value; } BCPressureTypeSeepageFace; typedef struct { double value; } BCPressureTypeOverlandKinematic; typedef struct { double value; } BCPressureTypeOverlandDiffusive;
typedef struct {
  int num_phases;
  int num_patches;
  int *types;
  int *patch_indexes;
  int *cycle_numbers;
  int *bc_types;
  void ***values;
  TimeCycleData *time_cycle_data;
} BCPressureData;
typedef struct {
  PFModule *geometries;
  PFModule *domain;
  int num_phases;
  int num_contaminants;
  double base_time_unit;
  int start_count;
  double start_time;
  double stop_time;
  double dump_interval;
  int dump_interval_execution_time_limit;
  int dump_at_end;
  PFModule *select_time_step;
  double gravity;
  double *phase_viscosity;
  double *contaminant_degradation;
  PFModule *phase_density;
  PFModule *permeability;
  PFModule *porosity;
  PFModule *retardation;
  PFModule *phase_mobility;
  PFModule *phase_rel_perm;
  PFModule *phase_source;
  PFModule *specific_storage;
  PFModule *FBx;
  PFModule *FBy;
  PFModule *FBz;
  PFModule *capillary_pressure;
  PFModule *saturation;
  PFModule *bc_internal;
  PFModule *bc_pressure;
  PFModule *bc_pressure_package;
  PFModule *bc_phase_saturation;
  PFModule *ic_phase_concen;
  PFModule *ic_phase_satur;
  PFModule *ic_phase_pressure;
  PFModule *l2_error_norm;
  PFModule *constitutive;
  PFModule *well_package;
  PFModule *x_slope;
  PFModule *y_slope;
  PFModule *mann;
  PFModule *overlandflow_eval;
  PFModule *overlandflow_eval_diff;
  PFModule *overlandflow_eval_kin;
  PFModule *dz_mult;
  PFModule *real_space_z;
} Problem;
typedef struct {
  int num_solids;
  GeomSolid **solids;
  GrGeomSolid **gr_solids;
  GeomSolid *domain;
  GrGeomSolid *gr_domain;
  Vector *index_of_domain_top;
  Vector *permeability_x;
  Vector *permeability_y;
  Vector *permeability_z;
  Vector *porosity;
  Vector *specific_storage;
  Vector *FBx;
  Vector *FBy;
  Vector *FBz;
  WellData *well_data;
  BCPressureData *bc_pressure_data;
  Vector *x_slope;
  Vector *y_slope;
  Vector *mann;
  Vector *x_sslope;
  Vector *y_sslope;
  Vector *dz_mult;
  Vector *rsz;
} ProblemData;
typedef struct {
  PFModule *nl_function_eval;
  PFModule *richards_jacobian_eval;
  PFModule *precond;
  PFModule *bc_pressure;
  ProblemData *problem_data;
  Matrix *jacobian_matrix;
  Matrix *jacobian_matrix_C;
  Matrix *jacobian_matrix_E;
  Matrix *jacobian_matrix_F;
  Vector *old_density;
  Vector *old_saturation;
  Vector *old_pressure;
  Vector *density;
  Vector *saturation;
  double dt;
  double time;
  double *outflow;
  Vector *evap_trans;
  Vector *ovrl_bc_flx;
  Vector *x_velocity;
  Vector *y_velocity;
  Vector *z_velocity;
} State;
typedef void (*AdvectionConcentrationInvoke) (ProblemData *problem_data, int phase, int concentration, Vector *old_concentration, Vector *new_concentration, Vector *x_velocity, Vector *y_velocity, Vector *z_velocity, Vector *solid_mass_factor, double time, double deltat, int order);
typedef PFModule *(*AdvectionConcentrationInitInstanceXtraType) (Problem *problem, Grid *grid, double *temp_data);
void Godunov(ProblemData *problem_data, int phase, int concentration, Vector *old_concentration, Vector *new_concentration, Vector *x_velocity, Vector *y_velocity, Vector *z_velocity, Vector *solid_mass_factor, double time, double deltat, int order);
PFModule *GodunovInitInstanceXtra(Problem *problem, Grid *grid, double *temp_data);
void GodunovFreeInstanceXtra(void);
PFModule *GodunovNewPublicXtra(void);
void GodunovFreePublicXtra(void);
int GodunovSizeOfTempData(void);
void Axpy(double alpha, Vector *x, Vector *y);
Background *ReadBackground(void);
void FreeBackground(Background *background);
void SetBackgroundBounds(Background *background, Grid *grid);
void LBInitializeBC(Lattice *lattice, Problem *problem, ProblemData *problem_data);
BCPressureData *NewBCPressureData(void);
void FreeBCPressureData(BCPressureData *bc_pressure_data);
void PrintBCPressureData(BCPressureData *bc_pressure_data);
typedef void (*BCPressurePackageInvoke) (ProblemData *problem_data);
typedef PFModule *(*BCPressurePackageInitInstanceXtraInvoke) (Problem *problem);
typedef PFModule *(*BCPressurePackageNewPublicXtraInvoke) (int num_phases);
void BCPressurePackage(ProblemData *problem_data);
PFModule *BCPressurePackageInitInstanceXtra(Problem *problem);
void BCPressurePackageFreeInstanceXtra(void);
PFModule *BCPressurePackageNewPublicXtra(int num_phases);
void BCPressurePackageFreePublicXtra(void);
int BCPressurePackageSizeOfTempData(void);
double **CalcElevations(GeomSolid *geom_solid, int ref_patch, SubgridArray *subgrids, ProblemData *problem_data);
typedef void (*LinearSolverInvoke) (Vector *x, Vector *b, double tol, int zero);
typedef PFModule *(*LinearSolverInitInstanceXtraInvoke) (Problem *problem, Grid *grid, ProblemData *problem_data, Matrix *A, double *temp_data);
typedef PFModule *(*LinearSolverNewPublicXtraInvoke) (char *name);
void CGHS(Vector *x, Vector *b, double tol, int zero);
PFModule *CGHSInitInstanceXtra(Problem *problem, Grid *grid, ProblemData *problem_data, Matrix *A, double *temp_data);
void CGHSFreeInstanceXtra(void);
PFModule *CGHSNewPublicXtra(char *name);
void CGHSFreePublicXtra(void);
int CGHSSizeOfTempData(void);
CommPkg *NewCharVectorUpdatePkg(CharVector *charvector, int update_mode);
CommHandle *InitCharVectorUpdate(CharVector *charvector, int update_mode);
void FinalizeCharVectorUpdate(CommHandle *handle);
CharVector *NewTempCharVector(Grid *grid, int nc, int num_ghost);
void SetTempCharVectorData(CharVector *charvector, char *data);
CharVector *NewCharVector(Grid *grid, int nc, int num_ghost);
void FreeTempCharVector(CharVector *charvector);
void FreeCharVector(CharVector *charvector);
void InitCharVector(CharVector *v, char value);
void InitCharVectorAll(CharVector *v, char value);
void InitCharVectorInc(CharVector *v, char value, int inc);
typedef void (*ChebyshevInvoke) (Vector *x, Vector *b, double tol, int zero, double ia, double ib, int num_iter);
typedef PFModule *(*ChebyshevInitInstanceXtraInvoke) (Problem *problem, Grid *grid, ProblemData *problem_data, Matrix *A, double *temp_data);
typedef PFModule *(*ChebyshevNewPublicXtraInvoke) (char *name);
void Chebyshev(Vector *x, Vector *b, double tol, int zero, double ia, double ib, int num_iter);
PFModule *ChebyshevInitInstanceXtra(Problem *problem, Grid *grid, ProblemData *problem_data, Matrix *A, double *temp_data);
void ChebyshevFreeInstanceXtra(void);
PFModule *ChebyshevNewPublicXtra(char *name);
void ChebyshevFreePublicXtra(void);
int ChebyshevSizeOfTempData(void);
void ProjectRegion(Region *region, int sx, int sy, int sz, int ix, int iy, int iz);
Region *ProjectRBPoint(Region *region, int rb [4 ][3 ]);
void CreateComputePkgs(Grid *grid);
void FreeComputePkgs(Grid *grid);
int NewCommPkgInfo(Subregion *data_sr, Subregion *comm_sr, int index, int num_vars, int *loop_array);
CommPkg *NewCommPkg(Region *send_region, Region *recv_region, SubregionArray *data_space, int num_vars, double *data);
void FreeCommPkg(CommPkg *pkg);
CommHandle *InitCommunication(CommPkg *comm_pkg);
void FinalizeCommunication(CommHandle *handle);
ComputePkg *NewComputePkg(Region *send_reg, Region *recv_reg, Region *dep_reg, Region *ind_reg);
void FreeComputePkg(ComputePkg *compute_pkg);
double ComputePhaseMaximum(double phase_u_max, double dx, double phase_v_max, double dy, double phase_w_max, double dz);
double ComputeTotalMaximum(Problem *problem, EvalStruct *eval_struct, double s_lower, double s_upper, double total_u_max, double dx, double total_v_max, double dy, double total_w_max, double beta_max, double dz);
double ComputeTotalConcen(GrGeomSolid *gr_domain, Grid *grid, Vector *substance);
typedef void (*ConstantRFInvoke) (GeomSolid *geounit, GrGeomSolid *gr_geounit, Vector *field, RFCondData *cdata);
typedef PFModule *(*ConstantRFInitInstanceXtraInvoke) (Grid *grid, double *temp_data);
void ConstantRF(GeomSolid *geounit, GrGeomSolid *gr_geounit, Vector *field, RFCondData *cdata);
PFModule *ConstantRFInitInstanceXtra(Grid *grid, double *temp_data);
void ConstantRFFreeInstanceXtra(void);
PFModule *ConstantRFNewPublicXtra(char *geom_name);
void ConstantRFFreePublicXtra(void);
int ConstantRFSizeOfTempData(void);
typedef void (*PorosityFieldInvoke) (GeomSolid *geounit, GrGeomSolid *gr_geounit, Vector *field);
typedef PFModule *(*PorosityFieldInitInstanceXtraInvoke) (Grid *grid, double *temp_data);
typedef PFModule *(*PorosityFieldNewPublicXtraInvoke) (char *geom_name);
void ConstantPorosity(GeomSolid *geounit, GrGeomSolid *gr_geounit, Vector *field);
PFModule *ConstantPorosityInitInstanceXtra(Grid *grid, double *temp_data);
void ConstantPorosityFreeInstanceXtra(void);
PFModule *ConstantPorosityNewPublicXtra(char *geom_name);
void ConstantPorosityFreePublicXtra(void);
int ConstantPorositySizeOfTempData(void);
void Copy(Vector *x, Vector *y);
SubgridArray *GetGridSubgrids(SubgridArray *all_subgrids);
Grid *CreateGrid(Grid *user_grid);
typedef void (*DiagScaleInvoke) (Vector *x, Matrix *A, Vector *b, Vector *d);
void DiagScale(Vector *x, Matrix *A, Vector *b, Vector *d);
void DiffuseLB(Lattice *lattice, Problem *problem, int max_iterations, char *file_prefix);
void LatticeFlowInit(Lattice *lattice, Problem *problem);
double MaxVectorValue(Vector *field);
double MaxVectorDividend(Vector *field1, Vector *field2);
typedef void (*DiscretizePressureInvoke) (Matrix **ptr_to_A, Vector **ptr_to_f, ProblemData *problem_data, double time, Vector *total_mobility_x, Vector *total_mobility_y, Vector *total_mobility_z, Vector **phase_saturations);
typedef PFModule *(*DiscretizePressureInitInstanceXtraInvoke) (Problem *problem, Grid *grid, double *temp_data);
void DiscretizePressure(Matrix **ptr_to_A, Vector **ptr_to_f, ProblemData *problem_data, double time, Vector *total_mobility_x, Vector *total_mobility_y, Vector *total_mobility_z, Vector **phase_saturations);
PFModule *DiscretizePressureInitInstanceXtra(Problem *problem, Grid *grid, double *temp_data);
void DiscretizePressureFreeInstanceXtra(void);
PFModule *DiscretizePressureNewPublicXtra(void);
void DiscretizePressureFreePublicXtra(void);
int DiscretizePressureSizeOfTempData(void);
SubgridArray *DistributeUserGrid(Grid *user_grid);
int dpofa_(double *a, int *lda, int *n, int *info);
double ddot_(int *n, double *dx, int *incx, double *dy, int *incy);
int dposl_(double *a, int *lda, int *n, double *b);
int daxpy_(int *n, double *da, double *dx, int *incx, double *dy, int *incy);
int gauinv_(double *p, double *xp, int *ierr);
char *malloc_chk(int size, char *file, int line);
char *calloc_chk(int count, int elt_size, char *file, int line);
int Exp2(int p);
void printMemoryInfo(FILE *log_file);
void recordMemoryInfo();
void printMaxMemory(FILE *log_file);
GeomTSolid *GeomNewTSolid(GeomTIN *surface, int **patches, int num_patches, int *num_patch_triangles);
void GeomFreeTSolid(GeomTSolid *solid);
int GeomReadTSolids(GeomTSolid ***solids_data_ptr, char *geom_input_name);
GeomTSolid *GeomTSolidFromBox(double xl, double yl, double zl, double xu, double yu, double zu);
GeomVertexArray *GeomNewVertexArray(GeomVertex **vertices, int nV);
void GeomFreeVertexArray(GeomVertexArray *vertex_array);
GeomTIN *GeomNewTIN(GeomVertexArray *vertex_array, GeomTriangle **triangles, int nT);
void GeomFreeTIN(GeomTIN *surface);
GeomSolid *GeomNewSolid(void *data, int type);
void GeomFreeSolid(GeomSolid *solid);
int GeomReadSolids(GeomSolid ***solids_ptr, char *geom_input_name, int type);
GeomSolid *GeomSolidFromBox(double xl, double yl, double zl, double xu, double yu, double zu, int type);
void IntersectLineWithTriangle(unsigned int line_direction, double coord_0, double coord_1, double v0_x, double v0_y, double v0_z, double v1_x, double v1_y, double v1_z, double v2_x, double v2_y, double v2_z, int *intersects, double *point, int *normal_component);
void NewGlobals(char *run_name);
void FreeGlobals(void);
void LogGlobals(void);
ListMember *NewListMember(double value, int normal_component, int triangle_id);
void FreeListMember(ListMember *member);
void ListInsert(ListMember **head, ListMember *member);
int ListDelete(ListMember **head, ListMember *member);
ListMember *ListSearch(ListMember *head, double value, int normal_component, int triangle_id);
ListMember *ListValueSearch(ListMember *head, double value);
ListMember *ListValueNormalComponentSearch(ListMember *head, double value, int normal_component);
ListMember *ListTriangleIDSearch(ListMember *head, int triangle_id);
void ListFree(ListMember **head);
int ListLength(ListMember *head);
void ListPrint(ListMember *head);
int GrGeomCheckOctree(GrGeomOctree *grgeom_octree);
void GrGeomFixOctree(GrGeomOctree *grgeom_octree, GrGeomOctree **patch_octrees, int num_patches, int level, int num_indices);
GrGeomOctree *GrGeomNewOctree(void);
void GrGeomNewOctreeChildren(GrGeomOctree *grgeom_octree);
void GrGeomFreeOctree(GrGeomOctree *grgeom_octree);
GrGeomOctree *GrGeomOctreeFind(int *new_level, GrGeomOctree *grgeom_octree_root, int ix, int iy, int iz, int level);
GrGeomOctree *GrGeomOctreeAddCell(GrGeomOctree *grgeom_octree_root, unsigned int cell, int ix, int iy, int iz, int level);
GrGeomOctree *GrGeomOctreeAddFace(GrGeomOctree *grgeom_octree_root, int line_direction, int cell_index0, int cell_index1, int face_index, int extent_lower, int extent_upper, int level, int normal_in_direction);
void GrGeomOctreeFromTIN(GrGeomOctree **solid_octree_ptr, GrGeomOctree ***patch_octrees_ptr, GeomTIN *solid, int **patches, int num_patches, int *num_patch_triangles, GrGeomExtentArray *extent_array, double xlower, double ylower, double zlower, double xupper, double yupper, double zupper, int min_level, int max_level);
void GrGeomOctreeFromInd(GrGeomOctree **solid_octree_ptr, Vector *indicator_field, int indicator, double xlower, double ylower, double zlower, double xupper, double yupper, double zupper, int octree_bg_level, int octree_ix, int octree_iy, int octree_iz);
void GrGeomPrintOctreeStruc(amps_File file, GrGeomOctree *grgeom_octree);
int GrGeomPrintOctreeLevel(amps_File file, GrGeomOctree *grgeom_octree, int level, int current_level);
void GrGeomPrintOctree(char *filename, GrGeomOctree *grgeom_octree_root);
void GrGeomPrintOctreeCells(char *filename, GrGeomOctree *octree, int last_level);
void GrGeomOctreeFree(GrGeomOctree *grgeom_octree_root);
int GrGeomGetOctreeInfo(double *xlp, double *ylp, double *zlp, double *xup, double *yup, double *zup, int *ixp, int *iyp, int *izp);
GrGeomExtentArray *GrGeomNewExtentArray(GrGeomExtents *extents, int size);
void GrGeomFreeExtentArray(GrGeomExtentArray *extent_array);
GrGeomExtentArray *GrGeomCreateExtentArray(SubgridArray *subgrids, int xl_ghost, int xu_ghost, int yl_ghost, int yu_ghost, int zl_ghost, int zu_ghost);
GrGeomSolid *GrGeomNewSolid(GrGeomOctree *data, GrGeomOctree **patches, int num_patches, int octree_bg_level, int octree_ix, int octree_iy, int octree_iz);
void GrGeomFreeSolid(GrGeomSolid *solid);
void GrGeomSolidFromInd(GrGeomSolid **solid_ptr, Vector *indicator_field, int indicator);
void GrGeomSolidFromGeom(GrGeomSolid **solid_ptr, GeomSolid *geom_solid, GrGeomExtentArray *extent_array);
Grid *NewGrid(SubgridArray *subgrids, SubgridArray *all_subgrids);
void FreeGrid(Grid *grid);
int ProjectSubgrid(Subgrid *subgrid, int sx, int sy, int sz, int ix, int iy, int iz);
Subgrid *ConvertToSubgrid(Subregion *subregion);
Subgrid *ExtractSubgrid(int rx, int ry, int rz, Subgrid *subgrid);
Subgrid *IntersectSubgrids(Subgrid *subgrid1, Subgrid *subgrid2);
SubgridArray *SubtractSubgrids(Subgrid *subgrid1, Subgrid *subgrid2);
SubgridArray *UnionSubgridArray(SubgridArray *subgrids);
HBT *HBT_new(
             int (*compare_method)(void *, void *),
             void (*free_method)(void *),
             void (*printf_method)(FILE *, void *),
             int (*scanf_method)(FILE *, void **),
             int malloc_flag);
HBT_element *_new_HBT_element(HBT *tree, void *object, int sizeof_obj);
void _free_HBT_element(HBT *tree, HBT_element *el);
void _HBT_free(HBT *tree, HBT_element *subtree);
void HBT_free(HBT *tree);
void *HBT_lookup(HBT *tree, void *obj);
void *HBT_replace(HBT *tree, void *obj, int sizeof_obj);
int HBT_insert(HBT *tree, void *obj, int sizeof_obj);
void *HBT_delete(HBT *tree, void *obj);
void *HBT_successor(HBT *tree, void *obj);
void HBT_printf(FILE *file, HBT *tree);
void HBT_scanf(FILE *file, HBT *tree);
double InfinityNorm(Vector *x);
double InnerProd(Vector *x, Vector *y);
void InputPorosity(GeomSolid *geounit, GrGeomSolid *gr_geounit, Vector *field);
PFModule *InputPorosityInitInstanceXtra(Grid *grid, double *temp_data);
void InputPorosityFreeInstanceXtra(void);
PFModule *InputPorosityNewPublicXtra(char *geom_name);
void InputPorosityFreePublicXtra(void);
int InputPorositySizeOfTempData(void);
void InputRF(GeomSolid *geounit, GrGeomSolid *gr_geounit, Vector *field, RFCondData *cdata);
PFModule *InputRFInitInstanceXtra(Grid *grid, double *temp_data);
void InputRFFreeInstanceXtra(void);
PFModule *InputRFNewPublicXtra(char *geom_name);
void InputRFFreePublicXtra(void);
int InputRFSizeOfTempData(void);
void IDB_Print(FILE *file, void *entry);
int IDB_Compare(void *a, void *b);
void IDB_Free(void *a);
IDB_Entry *IDB_NewEntry(char *key, char *value);
IDB *IDB_NewDB(char *filename);
void IDB_FreeDB(IDB *database);
void IDB_PrintUsage(FILE *file, IDB *database);
char *IDB_GetString(IDB *database, const char *key);
char *IDB_GetStringDefault(IDB *database, const char *key, char *default_value);
double IDB_GetDoubleDefault(IDB *database, const char *key, double default_value);
double IDB_GetDouble(IDB *database, const char *key);
int IDB_GetIntDefault(IDB *database, const char *key, int default_value);
int IDB_GetInt(IDB *database, const char *key);
NameArray NA_NewNameArray(char *string);
int NA_AppendToArray(NameArray name_array, char *string);
void NA_FreeNameArray(NameArray name_array);
int NA_NameToIndex(NameArray name_array, char *name);
char *NA_IndexToName(NameArray name_array, int index);
int NA_Sizeof(NameArray name_array);
void InputError(const char *format, const char *s1, const char *s2);
typedef int (*NonlinSolverInvoke) (Vector *pressure, Vector *density, Vector *old_density, Vector *saturation, Vector *old_saturation, double t, double dt, ProblemData *problem_data, Vector *old_pressure, Vector *evap_trans, Vector *ovrl_bc_flx, Vector *x_velocity, Vector *y_velocity, Vector *z_velocity);
typedef PFModule *(*NonlinSolverInitInstanceXtraInvoke) (Problem *problem, Grid *grid, ProblemData *problem_data, double *temp_data);
int KINSolInitPC(int neq, N_Vector pressure, N_Vector uscale, N_Vector fval, N_Vector fscale, N_Vector vtemp1, N_Vector vtemp2, void *nl_function, double uround, long int *nfePtr, void *current_state);
int KINSolCallPC(int neq, N_Vector pressure, N_Vector uscale, N_Vector fval, N_Vector fscale, N_Vector vtem, N_Vector ftem, void *nl_function, double uround, long int *nfePtr, void *current_state);
void PrintFinalStats(FILE *out_file, long int *integer_outputs_now, long int *integer_outputs_total);
int KinsolNonlinSolver(Vector *pressure, Vector *density, Vector *old_density, Vector *saturation, Vector *old_saturation, double t, double dt, ProblemData *problem_data, Vector *old_pressure, Vector *evap_trans, Vector *ovrl_bc_flx, Vector *x_velocity, Vector *y_velocity, Vector *z_velocity);
PFModule *KinsolNonlinSolverInitInstanceXtra(Problem *problem, Grid *grid, ProblemData *problem_data, double *temp_data);
void KinsolNonlinSolverFreeInstanceXtra(void);
PFModule *KinsolNonlinSolverNewPublicXtra(void);
void KinsolNonlinSolverFreePublicXtra(void);
int KinsolNonlinSolverSizeOfTempData(void);
typedef void (*KinsolPCInvoke) (Vector *rhs);
typedef PFModule * (*KinsolPCInitInstanceXtraInvoke) (Problem *problem, Grid *grid, ProblemData *problem_data, double *temp_data, Vector *pressure, Vector *old_pressure, Vector *saturation, Vector *density, double dt, double time);
typedef PFModule *(*KinsolPCNewPublicXtraInvoke) (char *name, char *pc_name);
void KinsolPC(Vector *rhs);
PFModule *KinsolPCInitInstanceXtra(Problem *problem, Grid *grid, ProblemData *problem_data, double *temp_data, Vector *pressure, Vector *old_pressure, Vector *saturation, Vector *density, double dt, double time);
void KinsolPCFreeInstanceXtra(void);
PFModule *KinsolPCNewPublicXtra(char *name, char *pc_name);
void KinsolPCFreePublicXtra(void);
int KinsolPCSizeOfTempData(void);
typedef void (*L2ErrorNormInvoke) (double time, Vector *pressure, ProblemData *problem_data, double *l2_error_norm);
void L2ErrorNorm(double time, Vector *pressure, ProblemData *problem_data, double *l2_error_norm);
PFModule *L2ErrorNormInitInstanceXtra(void);
void L2ErrorNormFreeInstanceXtra(void);
PFModule *L2ErrorNormNewPublicXtra(void);
void L2ErrorNormFreePublicXtra(void);
int L2ErrorNormSizeOfTempData(void);
void LineProc(double *Z, double phi, double theta, double dzeta, int izeta, int nzeta, double Kmax, double dK);
void NewLogging(void);
void FreeLogging(void);
FILE *OpenLogFile(char *module_name);
int CloseLogFile(FILE *log_file);
void PrintVersionInfo(FILE *log_file);
typedef void (*MatrixDiagScaleInvoke) (Vector *x, Matrix *A, Vector *b, int flag);
typedef PFModule *(*MatrixDiagScaleInitInstanceXtraInvoke) (Grid *grid);
typedef PFModule *(*MatrixDiagScaleNewPublicXtraInvoke) (char *name);
void MatDiagScale(Vector *x, Matrix *A, Vector *b, int flag);
PFModule *MatDiagScaleInitInstanceXtra(Grid *grid);
void MatDiagScaleFreeInstanceXtra(void);
PFModule *MatDiagScaleNewPublicXtra(char *name);
void MatDiagScaleFreePublicXtra(void);
int MatDiagScaleSizeOfTempData(void);
Stencil *NewStencil(int shape [][3 ], int sz);
CommPkg *NewMatrixUpdatePkg(Matrix *matrix, Stencil *ghost);
CommHandle *InitMatrixUpdate(Matrix *matrix);
void FinalizeMatrixUpdate(CommHandle *handle);
Matrix *NewMatrix(Grid *grid, SubregionArray *range, Stencil *stencil, int symmetry, Stencil *ghost);
Matrix *NewMatrixType(Grid *grid, SubregionArray *range, Stencil *stencil, int symmetry, Stencil *ghost, enum matrix_type type);
void FreeStencil(Stencil *stencil);
void FreeMatrix(Matrix *matrix);
void InitMatrix(Matrix *A, double value);
void Matvec(double alpha, Matrix *A, Vector *x, double beta, Vector *y);
void MatvecSubMat(void * current_state,
                  double alpha,
                  Matrix *JB,
                  Matrix *JC,
                  Vector *x,
                  double beta,
                  Vector *y);
void MatvecJacF(
                           ProblemData *problem_data,
                           double alpha,
                           Matrix * JF,
                           Vector * x,
                           double beta,
                           Vector * y);
void MatvecJacE(
                           ProblemData *problem_data,
                           double alpha,
                           Matrix * JE,
                           Vector * x,
                           double beta,
                           Vector * y);
double MaxFieldValue(Vector *field, Vector *phi, int dir);
double MaxPhaseFieldValue(Vector *x_velocity, Vector *y_velocity, Vector *z_velocity, Vector *phi);
double MaxTotalFieldValue(Problem *problem, EvalStruct *eval_struct, Vector *saturation, Vector *x_velocity, Vector *y_velocity, Vector *z_velocity, Vector *beta, Vector *phi);
typedef void (*PrecondInvoke) (Vector *x, Vector *b, double tol, int zero);
typedef PFModule * (*PrecondInitInstanceXtraInvoke) (Problem *problem, Grid *grid, ProblemData *problem_data, Matrix *A, Matrix *B, double *temp_data);
typedef PFModule *(*PrecondNewPublicXtra) (char *name);
void MGSemi(Vector *x, Vector *b, double tol, int zero);
void SetupCoarseOps(Matrix **A_l, Matrix **P_l, int num_levels, SubregionArray **f_sra_l, SubregionArray **c_sra_l);
PFModule *MGSemiInitInstanceXtra(Problem *problem, Grid *grid, ProblemData *problem_data, Matrix *A, double *temp_data);
void MGSemiFreeInstanceXtra(void);
PFModule *MGSemiNewPublicXtra(char *name);
void MGSemiFreePublicXtra(void);
int MGSemiSizeOfTempData(void);
void MGSemiProlong(Matrix *A_f, Vector *e_f, Vector *e_c, Matrix *P, SubregionArray *f_sr_array, SubregionArray *c_sr_array, ComputePkg *compute_pkg, CommPkg *e_f_comm_pkg);
ComputePkg *NewMGSemiProlongComputePkg(Grid *grid, Stencil *stencil, int sx, int sy, int sz, int c_index, int f_index);
void MGSemiRestrict(Matrix *A_f, Vector *r_f, Vector *r_c, Matrix *P, SubregionArray *f_sr_array, SubregionArray *c_sr_array, ComputePkg *compute_pkg, CommPkg *r_f_comm_pkg);
ComputePkg *NewMGSemiRestrictComputePkg(Grid *grid, Stencil *stencil, int sx, int sy, int sz, int c_index, int f_index);
void SetPf2KinsolData(Grid *grid, int num_ghost);
N_Vector N_VNew(int N, void *machEnv);
void N_VPrint(N_Vector x);
void FreeTempVector(Vector *vector);
void NewEndpts(double *alpha, double *beta, double *pp, int *size_ptr, int n, double *a_ptr, double *b_ptr, double *cond_ptr, double ereps);
typedef void (*NlFunctionEvalInvoke) (Vector *pressure, Vector *fval, ProblemData *problem_data, Vector *saturation, Vector *old_saturation, Vector *density, Vector *old_density, double dt, double time, Vector *old_pressure, Vector *evap_trans, Vector *ovrl_bc_flx, Vector *x_velocity, Vector *y_velocity, Vector *z_velocity);
typedef PFModule *(*NlFunctionEvalInitInstanceXtraInvoke) (Problem *problem, Grid *grid, double *temp_data);
void KINSolFunctionEval(int size, N_Vector pressure, N_Vector fval, void *current_state);
void NlFunctionEval(Vector *pressure, Vector *fval, ProblemData *problem_data, Vector *saturation, Vector *old_saturation, Vector *density, Vector *old_density, double dt, double time, Vector *old_pressure, Vector *evap_trans, Vector *ovrl_bc_flx, Vector *x_velocity, Vector *y_velocity, Vector *z_velocity);
PFModule *NlFunctionEvalInitInstanceXtra(Problem *problem, Grid *grid, double *temp_data);
void NlFunctionEvalFreeInstanceXtra(void);
PFModule *NlFunctionEvalNewPublicXtra(char *name);
void NlFunctionEvalFreePublicXtra(void);
int NlFunctionEvalSizeOfTempData(void);
void NoDiagScale(Vector *x, Matrix *A, Vector *b, int flag);
PFModule *NoDiagScaleInitInstanceXtra(Grid *grid);
void NoDiagScaleFreeInstanceXtra(void);
PFModule *NoDiagScaleNewPublicXtra(char *name);
void NoDiagScaleFreePublicXtra(void);
int NoDiagScaleSizeOfTempData(void);
int main(int argc, char *argv []);
void PCG(Vector *x, Vector *b, double tol, int zero);
PFModule *PCGInitInstanceXtra(Problem *problem, Grid *grid, ProblemData *problem_data, Matrix *A, Matrix *C, double *temp_data);
void PCGFreeInstanceXtra(void);
PFModule *PCGNewPublicXtra(char *name);
void PCGFreePublicXtra(void);
int PCGSizeOfTempData(void);
typedef void (*PermeabilityFaceInvoke) (Vector *zperm, Vector *permeability);
typedef PFModule *(*PermeabilityFaceInitInstanceXtraInvoke) (Grid *z_grid);
void PermeabilityFace(Vector *zperm, Vector *permeability);
PFModule *PermeabilityFaceInitInstanceXtra(Grid *z_grid);
void PermeabilityFaceFreeInstanceXtra(void);
PFModule *PermeabilityFaceNewPublicXtra(void);
void PermeabilityFaceFreePublicXtra(void);
int PermeabilityFaceSizeOfTempData(void);
void PerturbSystem(Lattice *lattice, Problem *problem);
PFModule *NewPFModule(void *call, void *init_instance_xtra, void *free_instance_xtra, void *new_public_xtra, void *free_public_xtra, void *sizeof_temp_data, void *instance_xtra, void *public_xtra);
PFModule *DupPFModule(PFModule *pf_module);
void FreePFModule(PFModule *pf_module);
void PFMG(Vector *soln, Vector *rhs, double tol, int zero);
PFModule *PFMGInitInstanceXtra(Problem *problem, Grid *grid, ProblemData *problem_data, Matrix *pf_Bmat, Matrix *pf_Cmat, double *temp_data);
void PFMGFreeInstanceXtra(void);
PFModule *PFMGNewPublicXtra(char *name);
void PFMGFreePublicXtra(void);
int PFMGSizeOfTempData(void);
void PFMGOctree(Vector *soln, Vector *rhs, double tol, int zero);
PFModule *PFMGOctreeInitInstanceXtra(
                                      Problem * problem,
                                      Grid * grid,
                                      ProblemData *problem_data,
                                      Matrix * pf_Bmat,
                                      Matrix * pf_Cmat,
                                      double * temp_data);
void PFMGOctreeFreeInstanceXtra(void);
PFModule *PFMGOctreeNewPublicXtra(char *name);
void PFMGOctreeFreePublicXtra(void);
int PFMGOctreeSizeOfTempData(void);
void SMG(Vector *soln, Vector *rhs, double tol, int zero);
PFModule *SMGInitInstanceXtra(Problem *problem, Grid *grid, ProblemData *problem_data, Matrix *pf_matrix, double *temp_data);
void SMGFreeInstanceXtra(void);
PFModule *SMGNewPublicXtra(char *name);
void SMGFreePublicXtra(void);
int SMGSizeOfTempData(void);
void PField(Grid *grid, GeomSolid *geounit, GrGeomSolid *gr_geounit, Vector *field, RFCondData *cdata, Statistics *stats);
void PGSRF(GeomSolid *geounit, GrGeomSolid *gr_geounit, Vector *field, RFCondData *cdata);
PFModule *PGSRFInitInstanceXtra(Grid *grid, double *temp_data);
void PGSRFFreeInstanceXtra(void);
PFModule *PGSRFNewPublicXtra(char *geom_name);
void PGSRFFreePublicXtra(void);
int PGSRFSizeOfTempData(void);
typedef void (*PhaseVelocityFaceInvoke) (Vector *xvel, Vector *yvel, Vector *zvel, ProblemData *problem_data, Vector *pressure, Vector **saturations, int phase);
typedef PFModule *(*PhaseVelocityFaceInitInstanceXtraInvoke) (Problem *problem, Grid *grid, Grid *x_grid, Grid *y_grid, Grid *z_grid, double *temp_data);
void PhaseVelocityFace(Vector *xvel, Vector *yvel, Vector *zvel, ProblemData *problem_data, Vector *pressure, Vector **saturations, int phase);
PFModule *PhaseVelocityFaceInitInstanceXtra(Problem *problem, Grid *grid, Grid *x_grid, Grid *y_grid, Grid *z_grid, double *temp_data);
void PhaseVelocityFaceFreeInstanceXtra(void);
PFModule *PhaseVelocityFaceNewPublicXtra(void);
void PhaseVelocityFaceFreePublicXtra(void);
int PhaseVelocityFaceSizeOfTempData(void);
void PPCG(Vector *x, Vector *b, double tol, int zero);
PFModule *PPCGInitInstanceXtra(Problem *problem, Grid *grid, ProblemData *problem_data, Matrix *A, double *temp_data);
void PPCGFreeInstanceXtra(void);
PFModule *PPCGNewPublicXtra(char *name);
void PPCGFreePublicXtra(void);
int PPCGSizeOfTempData(void);
void PrintGrid(char *filename, Grid *grid);
void PrintSubmatrixAll(amps_File file, Submatrix *submatrix, Stencil *stencil);
void PrintMatrixAll(char *filename, Matrix *A);
void PrintSubmatrix(amps_File file, Submatrix *submatrix, Subregion *subregion, Stencil *stencil);
void PrintMatrix(char *filename, Matrix *A);
void PrintSortMatrix(char *filename, Matrix *A, int all);
void PrintSubvectorAll(amps_File file, Subvector *subvector);
void PrintVectorAll(char *filename, Vector *v);
void PrintSubvector(amps_File file, Subvector *subvector, Subgrid *subgrid);
void PrintVector(char *filename, Vector *v);
Problem *NewProblem(int solver);
void FreeProblem(Problem *problem, int solver);
ProblemData *NewProblemData(Grid *grid, Grid *grid2d);
void FreeProblemData(ProblemData *problem_data);
BCStruct *NewBCStruct(SubgridArray *subgrids, GrGeomSolid *gr_domain, int num_patches, int *patch_indexes, int *bc_types, double ***values);
void FreeBCStruct(BCStruct *bc_struct);
typedef void (*BCInternalInvoke) (Problem *problem, ProblemData *problem_data, Matrix *A, Vector *f, double time);
void BCInternal(Problem *problem, ProblemData *problem_data, Matrix *A, Vector *f, double time);
PFModule *BCInternalInitInstanceXtra(void);
void BCInternalFreeInstanceXtra(void);
PFModule *BCInternalNewPublicXtra(void);
void BCInternalFreePublicXtra(void);
int BCInternalSizeOfTempData(void);
typedef void (*BCPhaseSaturationInvoke) (Vector *saturation, int phase, GrGeomSolid *gr_domain);
typedef PFModule *(*BCPhaseSaturationNewPublicXtraInvoke) (int num_phases);
void BCPhaseSaturation(Vector *saturation, int phase, GrGeomSolid *gr_domain);
PFModule *BCPhaseSaturationInitInstanceXtra(void);
void BCPhaseSaturationFreeInstanceXtra(void);
PFModule *BCPhaseSaturationNewPublicXtra(int num_phases);
void BCPhaseSaturationFreePublicXtra(void);
int BCPhaseSaturationSizeOfTempData(void);
typedef BCStruct *(*BCPressureInvoke) (ProblemData *problem_data, Grid *grid, GrGeomSolid *gr_domain, double time);
typedef PFModule *(*BCPressureInitInstanceXtraInvoke) (Problem *problem);
typedef PFModule *(*BCPressureNewPublicXtraInvoke) (int num_phases);
BCStruct *BCPressure(ProblemData *problem_data, Grid *grid, GrGeomSolid *gr_domain, double time);
PFModule *BCPressureInitInstanceXtra(Problem *problem);
void BCPressureFreeInstanceXtra(void);
PFModule *BCPressureNewPublicXtra(int num_phases);
void BCPressureFreePublicXtra(void);
int BCPressureSizeOfTempData(void);
typedef void (*CapillaryPressureInvoke) (Vector *capillary_pressure, int phase_i, int phase_j, ProblemData *problem_data, Vector *phase_saturation);
typedef PFModule *(*CapillaryPressureNewPublicXtraInvoke) (int num_phases);
void CapillaryPressure(Vector *capillary_pressure, int phase_i, int phase_j, ProblemData *problem_data, Vector *phase_saturation);
PFModule *CapillaryPressureInitInstanceXtra(void);
void CapillaryPressureFreeInstanceXtra(void);
PFModule *CapillaryPressureNewPublicXtra(int num_phases);
void CapillaryPressureFreePublicXtra(void);
int CapillaryPressureSizeOfTempData(void);
typedef void (*DomainInvoke) (ProblemData *problem_data);
typedef PFModule *(*DomainInitInstanceXtraInvoke) (Grid *grid);
void Domain(ProblemData *problem_data);
PFModule *DomainInitInstanceXtra(Grid *grid);
void DomainFreeInstanceXtra(void);
PFModule *DomainNewPublicXtra(void);
void DomainFreePublicXtra(void);
int DomainSizeOfTempData(void);
EvalStruct *NewEvalStruct(Problem *problem);
void FreeEvalStruct(EvalStruct *eval_struct);
typedef void (*GeometriesInvoke) (ProblemData *problem_data);
typedef PFModule *(*GeometriesInitInstanceXtraInvoke) (Grid *grid);
void Geometries(ProblemData *problem_data);
PFModule *GeometriesInitInstanceXtra(Grid *grid);
void GeometriesFreeInstanceXtra(void);
PFModule *GeometriesNewPublicXtra(void);
void GeometriesFreePublicXtra(void);
int GeometriesSizeOfTempData(void);
typedef void (*ICPhaseConcenInvoke) (Vector *ic_phase_concen, int phase, int contaminant, ProblemData *problem_data);
typedef PFModule *(*ICPhaseConcenNewPublicXtraInvoke) (int num_phases, int num_contaminants);
void ICPhaseConcen(Vector *ic_phase_concen, int phase, int contaminant, ProblemData *problem_data);
PFModule *ICPhaseConcenInitInstanceXtra(void);
void ICPhaseConcenFreeInstanceXtra(void);
PFModule *ICPhaseConcenNewPublicXtra(int num_phases, int num_contaminants);
void ICPhaseConcenFreePublicXtra(void);
int ICPhaseConcenSizeOfTempData(void);
typedef void (*ICPhasePressureInvoke) (Vector *ic_pressure, Vector *mask, ProblemData *problem_data, Problem *problem);
typedef PFModule *(*ICPhasePressureInitInstanceXtraInvoke) (Problem *problem, Grid *grid, double *temp_data);
void ICPhasePressure(Vector *ic_pressure, Vector *mask, ProblemData *problem_data, Problem *problem);
PFModule *ICPhasePressureInitInstanceXtra(Problem *problem, Grid *grid, double *temp_data);
void ICPhasePressureFreeInstanceXtra(void);
PFModule *ICPhasePressureNewPublicXtra(void);
void ICPhasePressureFreePublicXtra(void);
int ICPhasePressureSizeOfTempData(void);
typedef void (*ManningsInvoke) (ProblemData *problem_data, Vector *mann, Vector *dummy);
typedef PFModule *(*ManningsInitInstanceXtraInvoke) (Grid *grid3d, Grid *grid2d);
void Mannings(ProblemData *problem_data, Vector *mann, Vector *dummy);
PFModule *ManningsInitInstanceXtra(Grid *grid3d, Grid *grid2d);
void ManningsFreeInstanceXtra(void);
PFModule *ManningsNewPublicXtra(void);
void ManningsFreePublicXtra(void);
int ManningsSizeOfTempData(void);
typedef void (*SpecStorageInvoke) (ProblemData *problem_data, Vector *specific_storage);
void SpecStorage(ProblemData *problem_data, Vector *specific_storage);
PFModule *SpecStorageInitInstanceXtra(void);
void SpecStorageFreeInstanceXtra(void);
PFModule *SpecStorageNewPublicXtra(void);
void SpecStorageFreePublicXtra(void);
int SpecStorageSizeOfTempData(void);
typedef void (*dzScaleInvoke) (ProblemData *problem_data, Vector *dz_mult);
void dzScale(ProblemData *problem_data, Vector *dz_mult);
PFModule *dzScaleInitInstanceXtra(void);
void dzScaleFreeInstanceXtra(void);
PFModule *dzScaleNewPublicXtra(void);
void dzScaleFreePublicXtra(void);
int dzScaleSizeOfTempData(void);
typedef void (*FBxInvoke) (ProblemData *problem_data, Vector *FBx);
void FBx(ProblemData *problem_data, Vector *FBx);
PFModule *FBxInitInstanceXtra(void);
void FBxFreeInstanceXtra(void);
PFModule *FBxNewPublicXtra(void);
void FBxFreePublicXtra(void);
int FBxSizeOfTempData(void);
typedef void (*FByInvoke) (ProblemData *problem_data, Vector *FBy);
void FBy(ProblemData *problem_data, Vector *FBy);
PFModule *FByInitInstanceXtra(void);
void FByFreeInstanceXtra(void);
PFModule *FByNewPublicXtra(void);
void FByFreePublicXtra(void);
int FBySizeOfTempData(void);
typedef void (*FBzInvoke) (ProblemData *problem_data, Vector *FBz);
void FBz(ProblemData *problem_data, Vector *FBz);
PFModule *FBzInitInstanceXtra(void);
void FBzFreeInstanceXtra(void);
PFModule *FBzNewPublicXtra(void);
void FBzFreePublicXtra(void);
int FBzSizeOfTempData(void);
typedef void (*realSpaceZInvoke) (ProblemData *problem_data, Vector *rsz);
void realSpaceZ(ProblemData *problem_data, Vector *rsz);
PFModule *realSpaceZInitInstanceXtra(void);
void realSpaceZFreeInstanceXtra(void);
PFModule *realSpaceZNewPublicXtra(void);
void realSpaceZFreePublicXtra(void);
int realSpaceZSizeOfTempData(void);
typedef void (*OverlandFlowEvalInvoke) (Grid * grid,
                                        int sg,
                                        BCStruct * bc_struct,
                                        int ipatch,
                                        ProblemData *problem_data,
                                        Vector * pressure,
                                        Vector * old_pressure,
                                        double * ke_v,
                                        double * kw_v,
                                        double * kn_v,
                                        double * ks_v,
                                        double * qx_v,
                                        double * qy_v,
                                        int fcn);
void OverlandFlowEval(Grid * grid,
                      int sg,
                      BCStruct * bc_struct,
                      int ipatch,
                      ProblemData *problem_data,
                      Vector * pressure,
                      Vector * old_pressure,
                      double * ke_v,
                      double * kw_v,
                      double * kn_v,
                      double * ks_v,
                      double * qx_v,
                      double * qy_v,
                      int fcn);
PFModule *OverlandFlowEvalInitInstanceXtra(void);
void OverlandFlowEvalFreeInstanceXtra(void);
PFModule *OverlandFlowEvalNewPublicXtra(void);
void OverlandFlowEvalFreePublicXtra(void);
int OverlandFlowEvalSizeOfTempData(void);
typedef void (*OverlandFlowEvalDiffInvoke) (Grid * grid,
                                            int sg,
                                            BCStruct * bc_struct,
                                            int ipatch,
                                            ProblemData *problem_data,
                                            Vector * pressure,
         Vector * old_pressure,
                                            double * ke_v,
                                            double * kw_v,
                                            double * kn_v,
                                            double * ks_v,
                                            double * ke_vns,
                                            double * kw_vns,
                                            double * kn_vns,
                                            double * ks_vns,
                                            double * qx_v,
                                            double * qy_v,
                                            int fcn);
void OverlandFlowEvalDiff(Grid * grid,
                          int sg,
                          BCStruct * bc_struct,
                          int ipatch,
                          ProblemData *problem_data,
                          Vector * pressure,
     Vector * old_pressure,
                          double * ke_v,
                          double * kw_v,
                          double * kn_v,
                          double * ks_v,
                          double * ke_vns,
                          double * kw_vns,
                          double * kn_vns,
                          double * ks_vns,
                          double * qx_v,
                          double * qy_v,
                          int fcn);
PFModule *OverlandFlowEvalDiffInitInstanceXtra(void);
void OverlandFlowEvalDiffFreeInstanceXtra(void);
PFModule *OverlandFlowEvalDiffNewPublicXtra(void);
void OverlandFlowEvalDiffFreePublicXtra(void);
int OverlandFlowEvalDiffSizeOfTempData(void);
typedef void (*OverlandFlowEvalKinInvoke) (Grid * grid,
                                            int sg,
                                            BCStruct * bc_struct,
                                            int ipatch,
                                            ProblemData *problem_data,
                                            Vector * pressure,
                                            double * ke_v,
                                            double * kw_v,
                                            double * kn_v,
                                            double * ks_v,
                                            double * ke_vns,
                                            double * kw_vns,
                                            double * kn_vns,
                                            double * ks_vns,
                                            double * qx_v,
                                            double * qy_v,
                                            int fcn);
void OverlandFlowEvalKin(Grid * grid,
                          int sg,
                          BCStruct * bc_struct,
                          int ipatch,
                          ProblemData *problem_data,
                          Vector * pressure,
                          double * ke_v,
                          double * kw_v,
                          double * kn_v,
                          double * ks_v,
                          double * ke_vns,
                          double * kw_vns,
                          double * kn_vns,
                          double * ks_vns,
                          double * qx_v,
                          double * qy_v,
                          int fcn);
PFModule *OverlandFlowEvalKinInitInstanceXtra(void);
void OverlandFlowEvalKinFreeInstanceXtra(void);
PFModule *OverlandFlowEvalKinNewPublicXtra(void);
void OverlandFlowEvalKinFreePublicXtra(void);
int OverlandFlowEvalKinSizeOfTempData(void);
typedef void (*ICPhaseSaturInvoke) (Vector *ic_phase_satur, int phase, ProblemData *problem_data);
typedef PFModule *(*ICPhaseSaturNewPublicXtraInvoke) (int num_phases);
void ICPhaseSatur(Vector *ic_phase_satur, int phase, ProblemData *problem_data);
PFModule *ICPhaseSaturInitInstanceXtra(void);
void ICPhaseSaturFreeInstanceXtra(void);
PFModule *ICPhaseSaturNewPublicXtra(int num_phases);
void ICPhaseSaturFreePublicXtra(void);
int ICPhaseSaturSizeOfTempData(void);
typedef void (*PhaseDensityInvoke) (int phase, Vector *phase_pressure, Vector *density_v, double *pressure_d, double *density_d, int fcn);
typedef PFModule *(*PhaseDensityNewPublicXtraInvoke) (int num_phases);
void PhaseDensity(int phase, Vector *phase_pressure, Vector *density_v, double *pressure_d, double *density_d, int fcn);
PFModule *PhaseDensityInitInstanceXtra(void);
void PhaseDensityFreeInstanceXtra(void);
PFModule *PhaseDensityNewPublicXtra(int num_phases);
void PhaseDensityFreePublicXtra(void);
int PhaseDensitySizeOfTempData(void);
typedef void (*PhaseMobilityInvoke) (Vector *phase_mobility_x, Vector *phase_mobility_y, Vector *phase_mobility_z, Vector *perm_x, Vector *perm_y, Vector *perm_z, int phase, Vector *phase_saturation, double phase_viscosity);
typedef PFModule *(*PhaseMobilityNewPublicXtraInvoke) (int num_phases);
void PhaseMobility(Vector *phase_mobility_x, Vector *phase_mobility_y, Vector *phase_mobility_z, Vector *perm_x, Vector *perm_y, Vector *perm_z, int phase, Vector *phase_saturation, double phase_viscosity);
PFModule *PhaseMobilityInitInstanceXtra(void);
void PhaseMobilityFreeInstanceXtra(void);
PFModule *PhaseMobilityNewPublicXtra(int num_phases);
void PhaseMobilityFreePublicXtra(void);
int PhaseMobilitySizeOfTempData(void);
typedef void (*PhaseRelPermInvoke) (Vector *phase_rel_perm, Vector *phase_pressure, Vector *phase_density, double gravity, ProblemData *problem_data, int fcn);
typedef PFModule *(*PhaseRelPermInitInstanceXtraInvoke) (Grid *grid, double *temp_data);
void PhaseRelPerm(Vector *phase_rel_perm, Vector *phase_pressure, Vector *phase_density, double gravity, ProblemData *problem_data, int fcn);
PFModule *PhaseRelPermInitInstanceXtra(Grid *grid, double *temp_data);
void PhaseRelPermFreeInstanceXtra(void);
PFModule *PhaseRelPermNewPublicXtra(void);
void PhaseRelPermFreePublicXtra(void);
int PhaseRelPermSizeOfTempData(void);
typedef void (*PhaseSourceInvoke) (Vector *phase_source, int phase, Problem *problem, ProblemData *problem_data, double time);
typedef PFModule *(*PhaseSourceNewPublicXtraInvoke) (int num_phases);
void PhaseSource(Vector *phase_source, int phase, Problem *problem, ProblemData *problem_data, double time);
PFModule *PhaseSourceInitInstanceXtra(void);
void PhaseSourceFreeInstanceXtra(void);
PFModule *PhaseSourceNewPublicXtra(int num_phases);
void PhaseSourceFreePublicXtra(void);
int PhaseSourceSizeOfTempData(void);
typedef void (*PorosityInvoke) (ProblemData *problem_data, Vector *porosity, int num_geounits, GeomSolid **geounits, GrGeomSolid **gr_geounits);
typedef PFModule *(*PorosityInitInstanceXtraInvoke) (Grid *grid, double *temp_data);
void Porosity(ProblemData *problem_data, Vector *porosity, int num_geounits, GeomSolid **geounits, GrGeomSolid **gr_geounits);
PFModule *PorosityInitInstanceXtra(Grid *grid, double *temp_data);
void PorosityFreeInstanceXtra(void);
PFModule *PorosityNewPublicXtra(void);
void PorosityFreePublicXtra(void);
int PorositySizeOfTempData(void);
typedef void (*RetardationInvoke) (Vector *solidmassfactor, int contaminant, ProblemData *problem_data);
typedef PFModule *(*RetardationInitInstanceXtraInvoke) (double *temp_data);
typedef PFModule *(*RetardationNewPublicXtraInvoke) (int num_contaminants);
void Retardation(Vector *solidmassfactor, int contaminant, ProblemData *problem_data);
PFModule *RetardationInitInstanceXtra(double *temp_data);
void RetardationFreeInstanceXtra(void);
PFModule *RetardationNewPublicXtra(int num_contaminants);
void RetardationFreePublicXtra(void);
int RetardationSizeOfTempData(void);
typedef void (*RichardsBCInternalInvoke) (Problem *problem, ProblemData *problem_data, Vector *f, Matrix *A, double time, Vector *pressure, int fcn);
void RichardsBCInternal(Problem *problem, ProblemData *problem_data, Vector *f, Matrix *A, double time, Vector *pressure, int fcn);
PFModule *RichardsBCInternalInitInstanceXtra(void);
void RichardsBCInternalFreeInstanceXtra(void);
PFModule *RichardsBCInternalNewPublicXtra(void);
void RichardsBCInternalFreePublicXtra(void);
int RichardsBCInternalSizeOfTempData(void);
typedef void (*SaturationInvoke) (Vector *phase_saturation, Vector *phase_pressure, Vector *phase_density, double gravity, ProblemData *problem_data, int fcn);
typedef PFModule *(*SaturationInitInstanceXtraInvoke) (Grid *grid, double *temp_data);
void Saturation(Vector *phase_saturation, Vector *phase_pressure, Vector *phase_density, double gravity, ProblemData *problem_data, int fcn);
PFModule *SaturationInitInstanceXtra(Grid *grid, double *temp_data);
void SaturationFreeInstanceXtra(void);
PFModule *SaturationNewPublicXtra(void);
void SaturationFreePublicXtra(void);
int SaturationSizeOfTempData(void);
typedef void (*SaturationConstitutiveInvoke) (Vector **phase_saturations);
typedef PFModule *(*SaturationConstitutiveInitInstanceXtraInvoke) (Grid *grid);
typedef PFModule *(*SaturationConstitutiveNewPublicXtraInvoke) (int num_phases);
void SaturationConstitutive(Vector **phase_saturations);
PFModule *SaturationConstitutiveInitInstanceXtra(Grid *grid);
void SaturationConstitutiveFreeInstanceXtra(void);
PFModule *SaturationConstitutiveNewPublicXtra(int num_phases);
void SaturationConstitutiveFreePublicXtra(void);
int SaturationConstitutiveSizeOfTempData(void);
typedef void (*SlopeInvoke) (ProblemData *problem_data, Vector *x_sl, Vector *dummy);
typedef PFModule *(*SlopeInitInstanceXtraInvoke) (Grid *grid3d, Grid *grid2d);
void XSlope(ProblemData *problem_data, Vector *x_sl, Vector *dummy);
PFModule *XSlopeInitInstanceXtra(Grid *grid3d, Grid *grid2d);
void XSlopeFreeInstanceXtra(void);
PFModule *XSlopeNewPublicXtra(void);
void XSlopeFreePublicXtra(void);
int XSlopeSizeOfTempData(void);
void YSlope(ProblemData *problem_data, Vector *y_slope, Vector *dummy);
PFModule *YSlopeInitInstanceXtra(Grid *grid3d, Grid *grid2d);
void YSlopeFreeInstanceXtra(void);
PFModule *YSlopeNewPublicXtra(void);
void YSlopeFreePublicXtra(void);
int YSlopeSizeOfTempData(void);
void SeedRand(int seed);
double Rand(void);
int ratqr_(int *n, double *eps1, double *d, double *e, double *e2, int *m, double *w, int *ind, double *bd, int *type, int *idef, int *ierr);
double epslon_(double *x);
void RedBlackGSPoint(Vector *x, Vector *b, double tol, int zero);
PFModule *RedBlackGSPointInitInstanceXtra(Problem *problem, Grid *grid, ProblemData *problem_data, Matrix *A, double *temp_data);
void RedBlackGSPointFreeInstanceXtra(void);
PFModule *RedBlackGSPointNewPublicXtra(char *name);
void RedBlackGSPointFreePublicXtra(void);
int RedBlackGSPointSizeOfTempData(void);
void ReadPFBinary_Subvector(amps_File file, Subvector *subvector, Subgrid *subgrid);
void ReadPFBinary(char *filename, Vector *v);
void ComputeRegFromStencil(Region **dep_reg_ptr, Region **ind_reg_ptr, SubregionArray *cr_array, Region *send_reg, Region *recv_reg, Stencil *stencil);
SubgridArray *GetGridNeighbors(SubgridArray *subgrids, SubgridArray *all_subgrids, Stencil *stencil);
void CommRegFromStencil(Region **send_region_ptr, Region **recv_region_ptr, Grid *grid, Stencil *stencil);
Subregion *NewSubregion(int ix, int iy, int iz, int nx, int ny, int nz, int sx, int sy, int sz, int rx, int ry, int rz, int process);
SubregionArray *NewSubregionArray(void);
Region *NewRegion(int size);
void FreeSubregion(Subregion *subregion);
void FreeSubregionArray(SubregionArray *subregion_array);
void FreeRegion(Region *region);
Subregion *DuplicateSubregion(Subregion *subregion);
SubregionArray *DuplicateSubregionArray(SubregionArray *subregion_array);
Region *DuplicateRegion(Region *region);
void AppendSubregion(Subregion *subregion, SubregionArray *sr_array);
void DeleteSubregion(SubregionArray *sr_array, int index);
void AppendSubregionArray(SubregionArray *sr_array_0, SubregionArray *sr_array_1);
typedef void (*RichardsJacobianEvalInvoke) (Vector *pressure, Vector *old_pressure, Matrix **ptr_to_J, Matrix **ptr_to_JC, Vector *saturation, Vector *density, ProblemData *problem_data, double dt, double time, int symm_part);
typedef PFModule *(*RichardsJacobianEvalInitInstanceXtraInvoke) (Problem *problem, Grid *grid, ProblemData *problem_data, double *temp_data, int symmetric_jac);
typedef PFModule *(*RichardsJacobianEvalNewPublicXtraInvoke) (char *name);
int KINSolMatVec(void *current_state, N_Vector x, N_Vector y, int *recompute, N_Vector pressure);
void RichardsJacobianEval(Vector *pressure, Vector *old_pressure, Matrix **ptr_to_J, Matrix **ptr_to_JC, Vector *saturation, Vector *density, ProblemData *problem_data, double dt, double time, int symm_part);
PFModule *RichardsJacobianEvalInitInstanceXtra(Problem *problem, Grid *grid, ProblemData *problem_data, double *temp_data, int symmetric_jac);
void RichardsJacobianEvalFreeInstanceXtra(void);
PFModule *RichardsJacobianEvalNewPublicXtra(char *name);
void RichardsJacobianEvalFreePublicXtra(void);
int RichardsJacobianEvalSizeOfTempData(void);
typedef void (*AdvectionSaturationInvoke) (ProblemData *problem_data, int phase, Vector *old_saturation, Vector *new_saturation, Vector *x_velocity, Vector *y_velocity, Vector *z_velocity, Vector *z_permeability, Vector *solid_mass_factor, double *viscosity, double *density, double gravity, double time, double deltat, int order);
typedef PFModule *(*AdvectionSaturationInitInstanceXtraInvoke) (Problem *problem, Grid *grid, double *temp_data);
void SatGodunov(ProblemData *problem_data, int phase, Vector *old_saturation, Vector *new_saturation, Vector *x_velocity, Vector *y_velocity, Vector *z_velocity, Vector *z_permeability, Vector *solid_mass_factor, double *viscosity, double *density, double gravity, double time, double deltat, int order);
PFModule *SatGodunovInitInstanceXtra(Problem *problem, Grid *grid, double *temp_data);
void SatGodunovFreeInstanceXtra(void);
PFModule *SatGodunovNewPublicXtra(void);
void SatGodunovFreePublicXtra(void);
int SatGodunovSizeOfTempData(void);
void Scale(double alpha, Vector *y);
typedef void (*SelectTimeStepInvoke) (double *dt, char *dt_info, double time, Problem *problem, ProblemData *problem_data);
void SelectTimeStep(double *dt, char *dt_info, double time, Problem *problem, ProblemData *problem_data);
PFModule *SelectTimeStepInitInstanceXtra(void);
void SelectTimeStepFreeInstanceXtra(void);
PFModule *SelectTimeStepNewPublicXtra(void);
void SelectTimeStepFreePublicXtra(void);
int SelectTimeStepSizeOfTempData(void);
PFModule *WRFSelectTimeStepNewPublicXtra(
                                          double initial_step,
                                          double growth_factor,
                                          double max_step,
                                          double min_step);
void WRFSelectTimeStepFreePublicXtra();
PFModule *WRFSelectTimeStepInitInstanceXtra();
PFModule *WRFSelectTimeStepNewPublicXtra(
                                          double initial_step,
                                          double growth_factor,
                                          double max_step,
                                          double min_step);
void WRFSelectTimeStepFreePublicXtra();
typedef void (*SetProblemDataInvoke) (ProblemData *problem_data);
typedef PFModule *(*SetProblemDataInitInstanceXtraInvoke) (Problem *problem, Grid *grid, Grid *grid2d, double *temp_data);
void SetProblemData(ProblemData *problem_data);
PFModule *SetProblemDataInitInstanceXtra(Problem *problem, Grid *grid, Grid *grid2d, double *temp_data);
void SetProblemDataFreeInstanceXtra(void);
PFModule *SetProblemDataNewPublicXtra(void);
void SetProblemDataFreePublicXtra(void);
int SetProblemDataSizeOfTempData(void);
double **SimShear(double **shear_min_ptr, double **shear_max_ptr, GeomSolid *geom_solid, SubgridArray *subgrids, int type);
void Solve(void);
void NewSolver(void);
void FreeSolver(void);
typedef void (*SolverInvoke)(void);
typedef PFModule *(*SolverImpesNewPublicXtraInvoke) (char *name);
void SolverImpes(void);
PFModule *SolverImpesInitInstanceXtra(void);
void SolverImpesFreeInstanceXtra(void);
PFModule *SolverImpesNewPublicXtra(char *name);
void SolverImpesFreePublicXtra(void);
int SolverImpesSizeOfTempData(void);
void SolverDiffusion(void);
PFModule *SolverDiffusionInitInstanceXtra(void);
void SolverDiffusionFreeInstanceXtra(void);
PFModule *SolverDiffusionNewPublicXtra(char *name);
void SolverDiffusionFreePublicXtra(void);
int SolverDiffusionSizeOfTempData(void);
typedef PFModule *(*SolverNewPublicXtraInvoke) (char *name);
void SolverRichards(void);
PFModule *SolverRichardsInitInstanceXtra(void);
void SolverRichardsFreeInstanceXtra(void);
PFModule *SolverRichardsNewPublicXtra(char *name);
void SolverRichardsFreePublicXtra(void);
int SolverRichardsSizeOfTempData(void);
ProblemData *GetProblemDataRichards(PFModule *this_module);
Problem *GetProblemRichards(PFModule *this_module);
PFModule *GetICPhasePressureRichards(PFModule *this_module);
void AdvanceRichards(PFModule *this_module,
                     double start_time,
                     double stop_time,
                     PFModule *time_step_control,
                     Vector * evap_trans,
                     Vector ** pressure_out,
                     Vector ** porosity_out,
                     Vector ** saturation_out
                     );
void SetupRichards(PFModule *this_module);
typedef void (*SubsrfSimInvoke) (ProblemData *problem_data, Vector *perm_x, Vector *perm_y, Vector *perm_z, int num_geounits, GeomSolid **geounits, GrGeomSolid **gr_geounits);
typedef PFModule *(*SubsrfSimInitInstanceXtraInvoke) (Grid *grid, double *temp_data);
void SubsrfSim(ProblemData *problem_data, Vector *perm_x, Vector *perm_y, Vector *perm_z, int num_geounits, GeomSolid **geounits, GrGeomSolid **gr_geounits);
PFModule *SubsrfSimInitInstanceXtra(Grid *grid, double *temp_data);
void SubsrfSimFreeInstanceXtra(void);
PFModule *SubsrfSimNewPublicXtra(void);
void SubsrfSimFreePublicXtra(void);
int SubsrfSimSizeOfTempData(void);
TimeCycleData *NewTimeCycleData(int number_of_cycles, int *number_of_intervals);
void FreeTimeCycleData(TimeCycleData *time_cycle_data);
void PrintTimeCycleData(TimeCycleData *time_cycle_data);
int TimeCycleDataComputeIntervalNumber(Problem *problem, double time, TimeCycleData *time_cycle_data, int cycle_number);
double TimeCycleDataComputeNextTransition(Problem *problem, double time, TimeCycleData *time_cycle_data);
void ReadGlobalTimeCycleData(void);
void FreeGlobalTimeCycleData(void);
typedef void (*TotalVelocityFaceInvoke) (Vector *xvel, Vector *yvel, Vector *zvel, ProblemData *problem_data, Vector *total_mobility_x, Vector *total_mobility_y, Vector *total_mobility_z, Vector *pressure, Vector **saturations);
typedef PFModule *(*TotalVelocityFaceInitInstanceXtraInvoke) (Problem *problem, Grid *grid, Grid *x_grid, Grid *y_grid, Grid *z_grid, double *temp_data);
void TotalVelocityFace(Vector *xvel, Vector *yvel, Vector *zvel, ProblemData *problem_data, Vector *total_mobility_x, Vector *total_mobility_y, Vector *total_mobility_z, Vector *pressure, Vector **saturations);
PFModule *TotalVelocityFaceInitInstanceXtra(Problem *problem, Grid *grid, Grid *x_grid, Grid *y_grid, Grid *z_grid, double *temp_data);
void TotalVelocityFaceFreeInstanceXtra(void);
PFModule *TotalVelocityFaceNewPublicXtra(void);
void TotalVelocityFaceFreePublicXtra(void);
int TotalVelocityFaceSizeOfTempData(void);
void Turn(Vector *field, void *vxtra);
typedef void (*TurnInvoke) (Vector *field, void *vxtra);
int InitTurn(void);
void *NewTurn(char *geom_name);
void FreeTurn(void *xtra);
typedef void (*KFieldSimulatorInvoke) (GeomSolid *geounit, GrGeomSolid *gr_geounit, Vector *field, RFCondData *cdata);
typedef PFModule *(*KFieldSimulatorInitInstanceXtraInvoke) (Grid *grid, double *temp_data);
typedef PFModule *(*KFieldSimulatorNewPublicXtra) (char *geom_name);
void TurningBandsRF(GeomSolid *geounit, GrGeomSolid *gr_geounit, Vector *field, RFCondData *cdata);
PFModule *TurningBandsRFInitInstanceXtra(Grid *grid, double *temp_data);
void TurningBandsRFFreeInstanceXtra(void);
PFModule *TurningBandsRFNewPublicXtra(char *geom_name);
void TurningBandsRFFreePublicXtra(void);
int TurningBandsRFSizeOfTempData(void);
Subgrid *ReadUserSubgrid(void);
Grid *ReadUserGrid(void);
void FreeUserGrid(Grid *user_grid);
CommPkg *NewVectorCommPkg(Vector *vector, ComputePkg *compute_pkg);
VectorUpdateCommHandle *InitVectorUpdate(
                                          Vector *vector,
                                          int update_mode);
void FinalizeVectorUpdate(
                                  VectorUpdateCommHandle *handle);
Vector *NewVector(
                   Grid *grid,
                   int nc,
                   int num_ghost);
Vector *NewVectorType(
                       Grid * grid,
                       int nc,
                       int num_ghost,
                       enum vector_type type);
void FreeVector(Vector *vector);
void InitVector(Vector *v, double value);
void InitVectorAll(Vector *v, double value);
void InitVectorInc(Vector *v, double value, double inc);
void InitVectorRandom(Vector *v, long seed);
void PFVLinearSum(double a, Vector *x, double b, Vector *y, Vector *z);
void PFVConstInit(double c, Vector *z);
void PFVProd(Vector *x, Vector *y, Vector *z);
void PFVDiv(Vector *x, Vector *y, Vector *z);
void PFVScale(double c, Vector *x, Vector *z);
void PFVAbs(Vector *x, Vector *z);
void PFVInv(Vector *x, Vector *z);
void PFVAddConst(Vector *x, double b, Vector *z);
double PFVDotProd(Vector *x, Vector *y);
double PFVMaxNorm(Vector *x);
double PFVWrmsNorm(Vector *x, Vector *w);
double PFVWL2Norm(Vector *x, Vector *w);
double PFVL1Norm(Vector *x);
double PFVMin(Vector *x);
double PFVMax(Vector *x);
int PFVConstrProdPos(Vector *c, Vector *x);
void PFVCompare(double c, Vector *x, Vector *z);
int PFVInvTest(Vector *x, Vector *z);
void PFVCopy(Vector *x, Vector *y);
void PFVSum(Vector *x, Vector *y, Vector *z);
void PFVDiff(Vector *x, Vector *y, Vector *z);
void PFVNeg(Vector *x, Vector *z);
void PFVScaleSum(double c, Vector *x, Vector *y, Vector *z);
void PFVScaleDiff(double c, Vector *x, Vector *y, Vector *z);
void PFVLin1(double a, Vector *x, Vector *y, Vector *z);
void PFVLin2(double a, Vector *x, Vector *y, Vector *z);
void PFVAxpy(double a, Vector *x, Vector *y);
void PFVScaleBy(double a, Vector *x);
void PFVLayerCopy(int a, int b, Vector *x, Vector *y);
void WJacobi(Vector *x, Vector *b, double tol, int zero);
PFModule *WJacobiInitInstanceXtra(Problem *problem, Grid *grid, ProblemData *problem_data, Matrix *A, double *temp_data);
void WJacobiFreeInstanceXtra(void);
PFModule *WJacobiNewPublicXtra(char *name);
void WJacobiFreePublicXtra(void);
int WJacobiSizeOfTempData(void);
WellData *NewWellData(void);
void FreeWellData(WellData *well_data);
void PrintWellData(WellData *well_data, unsigned int print_mask);
void WriteWells(char *file_prefix, Problem *problem, WellData *well_data, double time, int write_header);
typedef void (*WellPackageInvoke) (ProblemData *problem_data);
typedef PFModule *(*WellPackageNewPublicXtraInvoke) (int num_phases, int num_contaminants);
void WellPackage(ProblemData *problem_data);
PFModule *WellPackageInitInstanceXtra(void);
void WellPackageFreeInstanceXtra(void);
PFModule *WellPackageNewPublicXtra(int num_phases, int num_contaminants);
void WellPackageFreePublicXtra(void);
int WellPackageSizeOfTempData(void);
void LBWells(Lattice *lattice, Problem *problem, ProblemData *problem_data);
long SizeofPFBinarySubvector(Subvector *subvector, Subgrid *subgrid);
void WritePFBinary_Subvector(amps_File file, Subvector *subvector, Subgrid *subgrid);
void WritePFBinary(char *file_prefix, char *file_suffix, Vector *v);
long SizeofPFSBinarySubvector(Subvector *subvector, Subgrid *subgrid, double drop_tolerance);
void WritePFSBinary_Subvector(amps_File file, Subvector *subvector, Subgrid *subgrid, double drop_tolerance);
void WritePFSBinary(char *file_prefix, char *file_suffix, Vector *v, double drop_tolerance);
void WriteSilo(char * file_prefix,
               char * file_type,
               char * file_suffix,
               Vector *v,
               double time,
               int step,
               char * variable_name);
void WriteSiloInit(char *file_prefix);
void pf_mk_dir(char *filename);
void WriteSiloPMPIO(char * file_prefix,
                        char * file_type,
                        char * file_suffix,
                        Vector *v,
                        double time,
                        int step,
                        char * variable_name);
void WriteSiloPMPIOInit(char *file_prefix);
void wrfparflowinit_();
void wrfparflowadvance_(double *current_time,
                        double *dt,
                        float * wrf_flux,
                        float * wrf_pressure,
                        float * wrf_porosity,
                        float * wrf_saturation,
                        int * num_soil_layers,
                        int * ghost_size_i_lower,
                        int * ghost_size_j_lower,
                        int * ghost_size_i_upper,
                        int * ghost_size_j_upper);
void WRF2PF(float * wrf_array,
            int wrf_depth,
            int ghost_size_i_lower,
            int ghost_size_j_lower,
            int ghost_size_i_upper,
            int ghost_size_j_upper,
            Vector *pf_vector,
            Vector *top);
void PF2WRF(Vector *pf_vector,
            float * wrf_array,
            int wrf_depth,
            int ghost_size_i_lower,
            int ghost_size_j_lower,
            int ghost_size_i_upper,
            int ghost_size_j_upper,
            Vector *top);
void ComputeTop(Problem * problem,
                ProblemData *problem_data
                );
int CheckTime(Problem *problem, char *key, double time);
void EvapTransSum(ProblemData *problem_data, double dt, Vector *evap_trans_sum, Vector *evap_trans);
void OverlandSum(ProblemData *problem_data,
                 Vector * pressure,
                 double dt,
                 Vector * overland_sum);
Grid *ReadProcessGrid();

void advect_(double *s, double *sn,
            double *uedge, double *vedge, double *wedge, double *phi,
            double *slx, double *sly, double *slz,
            int *lo, int *hi, int *dlo, int *dhi, double *hx, double *dt, int *fstord,
            double *sbot, double *stop, double *sbotp,
            double *sfrt, double *sbck,
            double *sleft, double *sright, double *sfluxz,
            double *dxscr, double *dyscr, double *dzscr, double *dzfrm);
void sadvect_(double *s, double *sn,
             double *uedge, double *vedge, double *wedge, double *betaedge, double *phi,
             double *viscosity, double *density, double *gravity,
             double *slx, double *sly, double *slz,
             int *lohi, int *dlohi, double *hx, double *dt,
             double *sbot, double *stop, double *sbotp,
             double *sfrt, double *sbck,
             double *sleft, double *sright, double *sfluxz,
             double *dxscr, double *dyscr, double *dzscr, double *dzfrm);
void clm_lsm_(double *pressure_data, double *saturation_data, double *evap_trans_data, double *mask, double *porosity_data,
             double *dz_mult_data, int *istep, double *dt, double *t, double *start_time,
             double *dx, double *dy, double *dz, int *ix, int *iy, int *nx, int *ny, int *nz,
             int *nx_f, int *ny_f, int *nz_f, int *nz_rz, int *ip, int *p, int *q, int *r, int *gnx, int *gny, int *rank,
             double *sw_data, double *lw_data, double *prcp_data, double *tas_data, double *u_data, double *v_data, double *patm_data, double *qatm_data,
             double *lai_data, double *sai_data, double *z0m_data, double *displa_data,
             double *eflx_lh_tot_data, double *eflx_lwrad_out_data, double *eflx_sh_tot_data, double *eflx_soil_grnd_data, double *qflx_eval_tot_data,
             double *qflx_evap_grnd_data, double *qflx_evap_soi_data, double *qflx_evap_veg_data, double *qflx_tran_veg_data,
             double *qflx_infl_data, double *swe_out_data, double *t_grnd_data, double *t_soil_data, int *clm_dump_interval, int *clm_1d_out,
             int *clm_forc_veg, char *clm_file_dir, int *clm_file_dir_length, int *clm_bin_out_dir, int *write_CLM_binary, int *clm_beta_function,
             int *clm_veg_function, double *clm_veg_wilting, double *clm_veg_fieldc, double *clm_res_sat,
             int *clm_irr_type, int *clm_irr_cycle, double *clm_irr_rate, double *clm_irr_start, double *clm_irr_stop,
             double *clm_irr_threshold, double *qirr, double *qirr_inst, double *iflag, int *clm_irr_thresholdtype, int *soi_z,
             int *clm_next, int *clm_write_logs, int *clm_last_rst, int *clm_daily_rst, int *clm_nlevsoi, int *clm_nlevlak);


typedef float float_t;
typedef double double_t;
extern int __fpclassify (double __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));
extern int __signbit (double __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));
extern int __isinf (double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern int __finite (double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern int __isnan (double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern int __iseqsig (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__));
extern int __issignaling (double __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));
extern double acos (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __acos (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double asin (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __asin (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double atan (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __atan (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double atan2 (double __y, double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __atan2 (double __y, double __x) __attribute__ ((__nothrow__ , __leaf__));
 extern double cos (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __cos (double __x) __attribute__ ((__nothrow__ , __leaf__));
 extern double sin (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __sin (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double tan (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __tan (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double cosh (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __cosh (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double sinh (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __sinh (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double tanh (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __tanh (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double acosh (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __acosh (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double asinh (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __asinh (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double atanh (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __atanh (double __x) __attribute__ ((__nothrow__ , __leaf__));
 extern double exp (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __exp (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double frexp (double __x, int *__exponent) __attribute__ ((__nothrow__ , __leaf__)); extern double __frexp (double __x, int *__exponent) __attribute__ ((__nothrow__ , __leaf__));
extern double ldexp (double __x, int __exponent) __attribute__ ((__nothrow__ , __leaf__)); extern double __ldexp (double __x, int __exponent) __attribute__ ((__nothrow__ , __leaf__));
 extern double log (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __log (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double log10 (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __log10 (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double modf (double __x, double *__iptr) __attribute__ ((__nothrow__ , __leaf__)); extern double __modf (double __x, double *__iptr) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2)));
extern double expm1 (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __expm1 (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double log1p (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __log1p (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double logb (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __logb (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double exp2 (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __exp2 (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double log2 (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __log2 (double __x) __attribute__ ((__nothrow__ , __leaf__));
 extern double pow (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)); extern double __pow (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__));
extern double sqrt (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __sqrt (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double hypot (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)); extern double __hypot (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__));
extern double cbrt (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __cbrt (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double ceil (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __ceil (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern double fabs (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __fabs (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern double floor (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __floor (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern double fmod (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)); extern double __fmod (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__));
extern int isinf (double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern int finite (double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern double drem (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)); extern double __drem (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__));
extern double significand (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __significand (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double copysign (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __copysign (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern double nan (const char *__tagb) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __nan (const char *__tagb) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern int isnan (double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern double j0 (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __j0 (double) __attribute__ ((__nothrow__ , __leaf__));
extern double j1 (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __j1 (double) __attribute__ ((__nothrow__ , __leaf__));
extern double jn (int, double) __attribute__ ((__nothrow__ , __leaf__)); extern double __jn (int, double) __attribute__ ((__nothrow__ , __leaf__));
extern double y0 (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __y0 (double) __attribute__ ((__nothrow__ , __leaf__));
extern double y1 (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __y1 (double) __attribute__ ((__nothrow__ , __leaf__));
extern double yn (int, double) __attribute__ ((__nothrow__ , __leaf__)); extern double __yn (int, double) __attribute__ ((__nothrow__ , __leaf__));
extern double erf (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __erf (double) __attribute__ ((__nothrow__ , __leaf__));
extern double erfc (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __erfc (double) __attribute__ ((__nothrow__ , __leaf__));
extern double lgamma (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __lgamma (double) __attribute__ ((__nothrow__ , __leaf__));
extern double tgamma (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __tgamma (double) __attribute__ ((__nothrow__ , __leaf__));
extern double gamma (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __gamma (double) __attribute__ ((__nothrow__ , __leaf__));
extern double lgamma_r (double, int *__signgamp) __attribute__ ((__nothrow__ , __leaf__)); extern double __lgamma_r (double, int *__signgamp) __attribute__ ((__nothrow__ , __leaf__));
extern double rint (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __rint (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double nextafter (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)); extern double __nextafter (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__));
extern double nexttoward (double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern double __nexttoward (double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));
extern double remainder (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)); extern double __remainder (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__));
extern double scalbn (double __x, int __n) __attribute__ ((__nothrow__ , __leaf__)); extern double __scalbn (double __x, int __n) __attribute__ ((__nothrow__ , __leaf__));
extern int ilogb (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern int __ilogb (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double scalbln (double __x, long int __n) __attribute__ ((__nothrow__ , __leaf__)); extern double __scalbln (double __x, long int __n) __attribute__ ((__nothrow__ , __leaf__));
extern double nearbyint (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __nearbyint (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double round (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __round (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern double trunc (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __trunc (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern double remquo (double __x, double __y, int *__quo) __attribute__ ((__nothrow__ , __leaf__)); extern double __remquo (double __x, double __y, int *__quo) __attribute__ ((__nothrow__ , __leaf__));
extern long int lrint (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long int __lrint (double __x) __attribute__ ((__nothrow__ , __leaf__));
__extension__
extern long long int llrint (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long long int __llrint (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long int lround (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long int __lround (double __x) __attribute__ ((__nothrow__ , __leaf__));
__extension__
extern long long int llround (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long long int __llround (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern double fdim (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)); extern double __fdim (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__));
extern double fmax (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __fmax (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern double fmin (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __fmin (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern double fma (double __x, double __y, double __z) __attribute__ ((__nothrow__ , __leaf__)); extern double __fma (double __x, double __y, double __z) __attribute__ ((__nothrow__ , __leaf__));
extern double scalb (double __x, double __n) __attribute__ ((__nothrow__ , __leaf__)); extern double __scalb (double __x, double __n) __attribute__ ((__nothrow__ , __leaf__));
extern int __fpclassifyf (float __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));
extern int __signbitf (float __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));
extern int __isinff (float __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern int __finitef (float __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern int __isnanf (float __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern int __iseqsigf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__));
extern int __issignalingf (float __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));
extern float acosf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __acosf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float asinf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __asinf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float atanf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __atanf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float atan2f (float __y, float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __atan2f (float __y, float __x) __attribute__ ((__nothrow__ , __leaf__));
 extern float cosf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __cosf (float __x) __attribute__ ((__nothrow__ , __leaf__));
 extern float sinf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __sinf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float tanf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __tanf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float coshf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __coshf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float sinhf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __sinhf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float tanhf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __tanhf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float acoshf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __acoshf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float asinhf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __asinhf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float atanhf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __atanhf (float __x) __attribute__ ((__nothrow__ , __leaf__));
 extern float expf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __expf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float frexpf (float __x, int *__exponent) __attribute__ ((__nothrow__ , __leaf__)); extern float __frexpf (float __x, int *__exponent) __attribute__ ((__nothrow__ , __leaf__));
extern float ldexpf (float __x, int __exponent) __attribute__ ((__nothrow__ , __leaf__)); extern float __ldexpf (float __x, int __exponent) __attribute__ ((__nothrow__ , __leaf__));
 extern float logf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __logf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float log10f (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __log10f (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float modff (float __x, float *__iptr) __attribute__ ((__nothrow__ , __leaf__)); extern float __modff (float __x, float *__iptr) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2)));
extern float expm1f (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __expm1f (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float log1pf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __log1pf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float logbf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __logbf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float exp2f (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __exp2f (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float log2f (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __log2f (float __x) __attribute__ ((__nothrow__ , __leaf__));
 extern float powf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)); extern float __powf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__));
extern float sqrtf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __sqrtf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float hypotf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)); extern float __hypotf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__));
extern float cbrtf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __cbrtf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float ceilf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __ceilf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern float fabsf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __fabsf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern float floorf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __floorf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern float fmodf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)); extern float __fmodf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__));
extern int isinff (float __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern int finitef (float __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern float dremf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)); extern float __dremf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__));
extern float significandf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __significandf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float copysignf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __copysignf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern float nanf (const char *__tagb) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __nanf (const char *__tagb) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern int isnanf (float __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern float j0f (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __j0f (float) __attribute__ ((__nothrow__ , __leaf__));
extern float j1f (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __j1f (float) __attribute__ ((__nothrow__ , __leaf__));
extern float jnf (int, float) __attribute__ ((__nothrow__ , __leaf__)); extern float __jnf (int, float) __attribute__ ((__nothrow__ , __leaf__));
extern float y0f (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __y0f (float) __attribute__ ((__nothrow__ , __leaf__));
extern float y1f (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __y1f (float) __attribute__ ((__nothrow__ , __leaf__));
extern float ynf (int, float) __attribute__ ((__nothrow__ , __leaf__)); extern float __ynf (int, float) __attribute__ ((__nothrow__ , __leaf__));
extern float erff (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __erff (float) __attribute__ ((__nothrow__ , __leaf__));
extern float erfcf (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __erfcf (float) __attribute__ ((__nothrow__ , __leaf__));
extern float lgammaf (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __lgammaf (float) __attribute__ ((__nothrow__ , __leaf__));
extern float tgammaf (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __tgammaf (float) __attribute__ ((__nothrow__ , __leaf__));
extern float gammaf (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __gammaf (float) __attribute__ ((__nothrow__ , __leaf__));
extern float lgammaf_r (float, int *__signgamp) __attribute__ ((__nothrow__ , __leaf__)); extern float __lgammaf_r (float, int *__signgamp) __attribute__ ((__nothrow__ , __leaf__));
extern float rintf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __rintf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float nextafterf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)); extern float __nextafterf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__));
extern float nexttowardf (float __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern float __nexttowardf (float __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));
extern float remainderf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)); extern float __remainderf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__));
extern float scalbnf (float __x, int __n) __attribute__ ((__nothrow__ , __leaf__)); extern float __scalbnf (float __x, int __n) __attribute__ ((__nothrow__ , __leaf__));
extern int ilogbf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern int __ilogbf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float scalblnf (float __x, long int __n) __attribute__ ((__nothrow__ , __leaf__)); extern float __scalblnf (float __x, long int __n) __attribute__ ((__nothrow__ , __leaf__));
extern float nearbyintf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __nearbyintf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float roundf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __roundf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern float truncf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __truncf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern float remquof (float __x, float __y, int *__quo) __attribute__ ((__nothrow__ , __leaf__)); extern float __remquof (float __x, float __y, int *__quo) __attribute__ ((__nothrow__ , __leaf__));
extern long int lrintf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern long int __lrintf (float __x) __attribute__ ((__nothrow__ , __leaf__));
__extension__
extern long long int llrintf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern long long int __llrintf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern long int lroundf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern long int __lroundf (float __x) __attribute__ ((__nothrow__ , __leaf__));
__extension__
extern long long int llroundf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern long long int __llroundf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern float fdimf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)); extern float __fdimf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__));
extern float fmaxf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __fmaxf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern float fminf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __fminf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern float fmaf (float __x, float __y, float __z) __attribute__ ((__nothrow__ , __leaf__)); extern float __fmaf (float __x, float __y, float __z) __attribute__ ((__nothrow__ , __leaf__));
extern float scalbf (float __x, float __n) __attribute__ ((__nothrow__ , __leaf__)); extern float __scalbf (float __x, float __n) __attribute__ ((__nothrow__ , __leaf__));
extern int __fpclassifyl (long double __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));
extern int __signbitl (long double __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));
extern int __isinfl (long double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern int __finitel (long double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern int __isnanl (long double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern int __iseqsigl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));
extern int __issignalingl (long double __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));
extern long double acosl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __acosl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double asinl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __asinl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double atanl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __atanl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double atan2l (long double __y, long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __atan2l (long double __y, long double __x) __attribute__ ((__nothrow__ , __leaf__));
 extern long double cosl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __cosl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
 extern long double sinl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __sinl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double tanl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __tanl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double coshl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __coshl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double sinhl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __sinhl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double tanhl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __tanhl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double acoshl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __acoshl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double asinhl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __asinhl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double atanhl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __atanhl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
 extern long double expl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __expl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double frexpl (long double __x, int *__exponent) __attribute__ ((__nothrow__ , __leaf__)); extern long double __frexpl (long double __x, int *__exponent) __attribute__ ((__nothrow__ , __leaf__));
extern long double ldexpl (long double __x, int __exponent) __attribute__ ((__nothrow__ , __leaf__)); extern long double __ldexpl (long double __x, int __exponent) __attribute__ ((__nothrow__ , __leaf__));
 extern long double logl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __logl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double log10l (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __log10l (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double modfl (long double __x, long double *__iptr) __attribute__ ((__nothrow__ , __leaf__)); extern long double __modfl (long double __x, long double *__iptr) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2)));
extern long double expm1l (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __expm1l (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double log1pl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __log1pl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double logbl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __logbl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double exp2l (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __exp2l (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double log2l (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __log2l (long double __x) __attribute__ ((__nothrow__ , __leaf__));
 extern long double powl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double __powl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));
extern long double sqrtl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __sqrtl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double hypotl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double __hypotl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));
extern long double cbrtl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __cbrtl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double ceill (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __ceill (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern long double fabsl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __fabsl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern long double floorl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __floorl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern long double fmodl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double __fmodl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));
extern int isinfl (long double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern int finitel (long double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern long double dreml (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double __dreml (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));
extern long double significandl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __significandl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double copysignl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __copysignl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern long double nanl (const char *__tagb) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __nanl (const char *__tagb) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern int isnanl (long double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern long double j0l (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __j0l (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double j1l (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __j1l (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double jnl (int, long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __jnl (int, long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double y0l (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __y0l (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double y1l (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __y1l (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double ynl (int, long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __ynl (int, long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double erfl (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __erfl (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double erfcl (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __erfcl (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double lgammal (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __lgammal (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double tgammal (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __tgammal (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double gammal (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __gammal (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double lgammal_r (long double, int *__signgamp) __attribute__ ((__nothrow__ , __leaf__)); extern long double __lgammal_r (long double, int *__signgamp) __attribute__ ((__nothrow__ , __leaf__));
extern long double rintl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __rintl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double nextafterl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double __nextafterl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));
extern long double nexttowardl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double __nexttowardl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));
extern long double remainderl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double __remainderl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));
extern long double scalbnl (long double __x, int __n) __attribute__ ((__nothrow__ , __leaf__)); extern long double __scalbnl (long double __x, int __n) __attribute__ ((__nothrow__ , __leaf__));
extern int ilogbl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern int __ilogbl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double scalblnl (long double __x, long int __n) __attribute__ ((__nothrow__ , __leaf__)); extern long double __scalblnl (long double __x, long int __n) __attribute__ ((__nothrow__ , __leaf__));
extern long double nearbyintl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __nearbyintl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double roundl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __roundl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern long double truncl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __truncl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern long double remquol (long double __x, long double __y, int *__quo) __attribute__ ((__nothrow__ , __leaf__)); extern long double __remquol (long double __x, long double __y, int *__quo) __attribute__ ((__nothrow__ , __leaf__));
extern long int lrintl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long int __lrintl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
__extension__
extern long long int llrintl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long long int __llrintl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long int lroundl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long int __lroundl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
__extension__
extern long long int llroundl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long long int __llroundl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long double fdiml (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double __fdiml (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));
extern long double fmaxl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __fmaxl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern long double fminl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __fminl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern long double fmal (long double __x, long double __y, long double __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double __fmal (long double __x, long double __y, long double __z) __attribute__ ((__nothrow__ , __leaf__));
extern long double scalbl (long double __x, long double __n) __attribute__ ((__nothrow__ , __leaf__)); extern long double __scalbl (long double __x, long double __n) __attribute__ ((__nothrow__ , __leaf__));
extern int __fpclassifyf128 (_Float128 __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));
extern int __signbitf128 (_Float128 __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));
extern int __isinff128 (_Float128 __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern int __finitef128 (_Float128 __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern int __isnanf128 (_Float128 __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern int __iseqsigf128 (_Float128 __x, _Float128 __y) __attribute__ ((__nothrow__ , __leaf__));
extern int __issignalingf128 (_Float128 __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));
extern int signgam;
enum
  {
    FP_NAN =
      0,
    FP_INFINITE =
      1,
    FP_ZERO =
      2,
    FP_SUBNORMAL =
      3,
    FP_NORMAL =
      4
  };

typedef struct list_member {
  double value;
  int normal_component;
  int triangle_id;
  struct list_member *next_list_member;
} ListMember;
typedef int Point[3];
typedef struct
{
  Point lo;
  Point up;
} Box;
typedef struct _BoxListElement
{
  Box box;
  struct _BoxListElement* next;
  struct _BoxListElement* prev;
} BoxListElement;
typedef struct _BoxList
{
  BoxListElement* head;
  BoxListElement* tail;
  unsigned int size;
} BoxList;
typedef struct _BoxArray
{
  Box* boxes;
  unsigned int size;
} BoxArray;
void PointCopy(Point dst, Point src);
int BoxSize(Box *box);
void BoxNumberCells(Box* box, Point* number_cells);
void BoxClear(Box *box);
void BoxSet(Box *box, Point lo, Point up);
void BoxCopy(Box *dst, Box *src);
void BoxPrint(Box* box);
BoxList* NewBoxList(void);
void FreeBoxList(BoxList *box_list);
int BoxListSize(BoxList *box_list);
int BoxListIsEmpty(BoxList *box_list);
Box* BoxListFront(BoxList *box_list);
void BoxListAppend(BoxList* box_list, Box* box);
void BoxListConcatenate(BoxList *box_list, BoxList *concatenate_list);
void BoxListClearItems(BoxList* box_list);
void BoxListPrint(BoxList* box_list);
BoxArray* NewBoxArray(BoxList *box_list);
void FreeBoxArray(BoxArray* box_array);
typedef int GrGeomExtents[6];
typedef struct {
  GrGeomExtents *extents;
  int size;
} GrGeomExtentArray;
typedef struct {
  GrGeomOctree *data;
  GrGeomOctree **patches;
  int num_patches;
  int octree_bg_level;
  int octree_ix, octree_iy, octree_iz;
  BoxArray* interior_boxes;
  BoxArray* surface_boxes[6];
  BoxArray** patch_boxes[6];
} GrGeomSolid;
for (ipatch = 0; ipatch < ((bc_struct)->num_patches); ipatch++) { bc_patches_values = ((bc_struct)->values[ipatch][is]); switch(((bc_struct)->bc_types[ipatch])) { case 0: { { GrGeomSolid *PV_gr_domain = ((bc_struct)->gr_domain); int PV_patch_index = ((bc_struct)->patch_indexes[ipatch]); Subgrid *PV_subgrid = ((Subgrid*)((((bc_struct)->subgrids))->subregions[(is)])); int PV_r = ((PV_subgrid)->rx); int PV_ix = ((PV_subgrid)->ix); int PV_iy = ((PV_subgrid)->iy); int PV_iz = ((PV_subgrid)->iz); int PV_nx = ((PV_subgrid)->nx); int PV_ny = ((PV_subgrid)->ny); int PV_nz = ((PV_subgrid)->nz); ival = 0; { if(PV_r == 0 && ((PV_gr_domain)->patch_boxes[(6 -1)][(PV_patch_index)])) { { int PV_ixl, PV_iyl, PV_izl, PV_ixu, PV_iyu, PV_izu; int *PV_visiting = ((void *)0); for(int PV_f=0; PV_f < 6; PV_f++) { BoxArray* boxes = ((PV_gr_domain)->patch_boxes[(PV_f)][(PV_patch_index)]); for(int PV_box = 0; PV_box < (boxes -> size); PV_box++) { Box box = (boxes -> boxes[PV_box]); PV_ixl = (((PV_ix) < (box.lo[0])) ? (box.lo[0]) : (PV_ix)); PV_iyl = (((PV_iy) < (box.lo[1])) ? (box.lo[1]) : (PV_iy)); PV_izl = (((PV_iz) < (box.lo[2])) ? (box.lo[2]) : (PV_iz)); PV_ixu = ((((PV_ix + PV_nx - 1)) < (box.up[0])) ? ((PV_ix + PV_nx - 1)) : (box.up[0])); PV_iyu = ((((PV_iy + PV_ny - 1)) < (box.up[1])) ? ((PV_iy + PV_ny - 1)) : (box.up[1])); PV_izu = ((((PV_iz + PV_nz - 1)) < (box.up[2])) ? ((PV_iz + PV_nz - 1)) : (box.up[2])); for(k = PV_izl; k <= PV_izu; k++) for(j =PV_iyl; j <= PV_iyu; j++) for(i = PV_ixl; i <= PV_ixu; i++) { { ip = 1; }; switch ((PV_f)) { ; case 0: { { op = wp; }; break; } }; { { ip = 2; }; ival++; }; } } } }; } else { GrGeomOctree *PV_node; double PV_ref = pow(2.0, PV_r); i = ((PV_gr_domain)->octree_ix) * (int)PV_ref; j = ((PV_gr_domain)->octree_iy) * (int)PV_ref; k = ((PV_gr_domain)->octree_iz) * (int)PV_ref; { int PV_f; { int PV_i, PV_j, PV_k, PV_l; int PV_ixl, PV_iyl, PV_izl, PV_ixu, PV_iyu, PV_izu; PV_i = i; PV_j = j; PV_k = k; { int PV_level = ((PV_gr_domain)->octree_bg_level) + PV_r; unsigned int PV_inc; int *PV_visiting; int PV_visit_child; PV_node = ((PV_gr_domain)->patches[(PV_patch_index)]); PV_l = 0; PV_inc = 1 << PV_level; PV_visiting = ((PV_level + 2) ? (int*)calloc((unsigned int)(PV_level + 2), (unsigned int)sizeof(int)) : ((void *)0)); PV_visiting++; PV_visiting[0] = 0; while (PV_l >= 0) { if (PV_l == PV_level) { if ((((((PV_node)->flags) & ((unsigned char)0x04)) == ((unsigned char)0x04))) && (1)) { if ((PV_i >= PV_ix) && (PV_i < (PV_ix + PV_nx)) && (PV_j >= PV_iy) && (PV_j < (PV_iy + PV_ny)) && (PV_k >= PV_iz) && (PV_k < (PV_iz + PV_nz))) { i = PV_i; j = PV_j; k = PV_k; { for (PV_f = 0; PV_f < 6; PV_f++) if (((((PV_node)->faces) & (((unsigned char)0x01) << PV_f)) != ((unsigned char)0x00))) { { ip = 1; }; switch ((PV_f)) { ; case 0: { { op = wp; }; break; } }; { { ip = 2; }; ival++; }; } }; } }; PV_visit_child = 0; } else if (((((PV_node)->flags) & ((unsigned char)0x10)) == ((unsigned char)0x10))) { if (((((PV_node)->flags) & ((unsigned char)0x04)) == ((unsigned char)0x04)) && (1)) { PV_ixl = (((PV_ix) < (PV_i)) ? (PV_i) : (PV_ix)); PV_iyl = (((PV_iy) < (PV_j)) ? (PV_j) : (PV_iy)); PV_izl = (((PV_iz) < (PV_k)) ? (PV_k) : (PV_iz)); PV_ixu = ((((PV_ix + PV_nx)) < ((PV_i + (int)PV_inc))) ? ((PV_ix + PV_nx)) : ((PV_i + (int)PV_inc))); PV_iyu = ((((PV_iy + PV_ny)) < ((PV_j + (int)PV_inc))) ? ((PV_iy + PV_ny)) : ((PV_j + (int)PV_inc))); PV_izu = ((((PV_iz + PV_nz)) < ((PV_k + (int)PV_inc))) ? ((PV_iz + PV_nz)) : ((PV_k + (int)PV_inc))); for (k = PV_izl; k < PV_izu; k++) for (j = PV_iyl; j < PV_iyu; j++) for (i = PV_ixl; i < PV_ixu; i++) { { for (PV_f = 0; PV_f < 6; PV_f++) if (((((PV_node)->faces) & (((unsigned char)0x01) << PV_f)) != ((unsigned char)0x00))) { { ip = 1; }; switch ((PV_f)) { ; case 0: { { op = wp; }; break; } }; { { ip = 2; }; ival++; }; } }; } }; PV_visit_child = 0; } else if (!((((PV_node)->flags) & ((unsigned char)0x04)) == ((unsigned char)0x04))) { PV_visit_child = 0; } else if (PV_visiting[PV_l] < 8) PV_visit_child = 1; else PV_visit_child = 0; if (PV_visit_child) { PV_node = ((PV_node)->children[PV_visiting[PV_l]]); PV_inc = PV_inc >> 1; PV_i += (int)(PV_inc) * ((PV_visiting[PV_l] & 1) ? 1 : 0); PV_j += (int)(PV_inc) * ((PV_visiting[PV_l] & 2) ? 1 : 0); PV_k += (int)(PV_inc) * ((PV_visiting[PV_l] & 4) ? 1 : 0); PV_l++; PV_visiting[PV_l] = 0; } else { PV_l--; PV_i -= (int)(PV_inc) * ((PV_visiting[PV_l] & 1) ? 1 : 0); PV_j -= (int)(PV_inc) * ((PV_visiting[PV_l] & 2) ? 1 : 0); PV_k -= (int)(PV_inc) * ((PV_visiting[PV_l] & 4) ? 1 : 0); PV_inc = PV_inc << 1; PV_node = ((PV_node)->parent); PV_visiting[PV_l]++; } } if (PV_visiting - 1) free(PV_visiting - 1); else {}; } }; }; } }; }; break; } case 1: { { GrGeomSolid *PV_gr_domain = ((bc_struct)->gr_domain); int PV_patch_index = ((bc_struct)->patch_indexes[ipatch]); Subgrid *PV_subgrid = ((Subgrid*)((((bc_struct)->subgrids))->subregions[(is)])); int PV_r = ((PV_subgrid)->rx); int PV_ix = ((PV_subgrid)->ix); int PV_iy = ((PV_subgrid)->iy); int PV_iz = ((PV_subgrid)->iz); int PV_nx = ((PV_subgrid)->nx); int PV_ny = ((PV_subgrid)->ny); int PV_nz = ((PV_subgrid)->nz); ival = 0; { if(PV_r == 0 && ((PV_gr_domain)->patch_boxes[(6 -1)][(PV_patch_index)])) { { int PV_ixl, PV_iyl, PV_izl, PV_ixu, PV_iyu, PV_izu; int *PV_visiting = ((void *)0); for(int PV_f=0; PV_f < 6; PV_f++) { BoxArray* boxes = ((PV_gr_domain)->patch_boxes[(PV_f)][(PV_patch_index)]); for(int PV_box = 0; PV_box < (boxes -> size); PV_box++) { Box box = (boxes -> boxes[PV_box]); PV_ixl = (((PV_ix) < (box.lo[0])) ? (box.lo[0]) : (PV_ix)); PV_iyl = (((PV_iy) < (box.lo[1])) ? (box.lo[1]) : (PV_iy)); PV_izl = (((PV_iz) < (box.lo[2])) ? (box.lo[2]) : (PV_iz)); PV_ixu = ((((PV_ix + PV_nx - 1)) < (box.up[0])) ? ((PV_ix + PV_nx - 1)) : (box.up[0])); PV_iyu = ((((PV_iy + PV_ny - 1)) < (box.up[1])) ? ((PV_iy + PV_ny - 1)) : (box.up[1])); PV_izu = ((((PV_iz + PV_nz - 1)) < (box.up[2])) ? ((PV_iz + PV_nz - 1)) : (box.up[2])); for(k = PV_izl; k <= PV_izu; k++) for(j =PV_iyl; j <= PV_iyu; j++) for(i = PV_ixl; i <= PV_ixu; i++) { { ip = 3; }; switch ((PV_f)) { ; case 1: { { op = ep; }; break; } }; { { ip = 4; }; ival++; }; } } } }; } else { GrGeomOctree *PV_node; double PV_ref = pow(2.0, PV_r); i = ((PV_gr_domain)->octree_ix) * (int)PV_ref; j = ((PV_gr_domain)->octree_iy) * (int)PV_ref; k = ((PV_gr_domain)->octree_iz) * (int)PV_ref; { int PV_f; { int PV_i, PV_j, PV_k, PV_l; int PV_ixl, PV_iyl, PV_izl, PV_ixu, PV_iyu, PV_izu; PV_i = i; PV_j = j; PV_k = k; { int PV_level = ((PV_gr_domain)->octree_bg_level) + PV_r; unsigned int PV_inc; int *PV_visiting; int PV_visit_child; PV_node = ((PV_gr_domain)->patches[(PV_patch_index)]); PV_l = 0; PV_inc = 1 << PV_level; PV_visiting = ((PV_level + 2) ? (int*)calloc((unsigned int)(PV_level + 2), (unsigned int)sizeof(int)) : ((void *)0)); PV_visiting++; PV_visiting[0] = 0; while (PV_l >= 0) { if (PV_l == PV_level) { if ((((((PV_node)->flags) & ((unsigned char)0x04)) == ((unsigned char)0x04))) && (1)) { if ((PV_i >= PV_ix) && (PV_i < (PV_ix + PV_nx)) && (PV_j >= PV_iy) && (PV_j < (PV_iy + PV_ny)) && (PV_k >= PV_iz) && (PV_k < (PV_iz + PV_nz))) { i = PV_i; j = PV_j; k = PV_k; { for (PV_f = 0; PV_f < 6; PV_f++) if (((((PV_node)->faces) & (((unsigned char)0x01) << PV_f)) != ((unsigned char)0x00))) { { ip = 3; }; switch ((PV_f)) { ; case 1: { { op = ep; }; break; } }; { { ip = 4; }; ival++; }; } }; } }; PV_visit_child = 0; } else if (((((PV_node)->flags) & ((unsigned char)0x10)) == ((unsigned char)0x10))) { if (((((PV_node)->flags) & ((unsigned char)0x04)) == ((unsigned char)0x04)) && (1)) { PV_ixl = (((PV_ix) < (PV_i)) ? (PV_i) : (PV_ix)); PV_iyl = (((PV_iy) < (PV_j)) ? (PV_j) : (PV_iy)); PV_izl = (((PV_iz) < (PV_k)) ? (PV_k) : (PV_iz)); PV_ixu = ((((PV_ix + PV_nx)) < ((PV_i + (int)PV_inc))) ? ((PV_ix + PV_nx)) : ((PV_i + (int)PV_inc))); PV_iyu = ((((PV_iy + PV_ny)) < ((PV_j + (int)PV_inc))) ? ((PV_iy + PV_ny)) : ((PV_j + (int)PV_inc))); PV_izu = ((((PV_iz + PV_nz)) < ((PV_k + (int)PV_inc))) ? ((PV_iz + PV_nz)) : ((PV_k + (int)PV_inc))); for (k = PV_izl; k < PV_izu; k++) for (j = PV_iyl; j < PV_iyu; j++) for (i = PV_ixl; i < PV_ixu; i++) { { for (PV_f = 0; PV_f < 6; PV_f++) if (((((PV_node)->faces) & (((unsigned char)0x01) << PV_f)) != ((unsigned char)0x00))) { { ip = 3; }; switch ((PV_f)) { ; case 1: { { op = ep; }; break; } }; { { ip = 4; }; ival++; }; } }; } }; PV_visit_child = 0; } else if (!((((PV_node)->flags) & ((unsigned char)0x04)) == ((unsigned char)0x04))) { PV_visit_child = 0; } else if (PV_visiting[PV_l] < 8) PV_visit_child = 1; else PV_visit_child = 0; if (PV_visit_child) { PV_node = ((PV_node)->children[PV_visiting[PV_l]]); PV_inc = PV_inc >> 1; PV_i += (int)(PV_inc) * ((PV_visiting[PV_l] & 1) ? 1 : 0); PV_j += (int)(PV_inc) * ((PV_visiting[PV_l] & 2) ? 1 : 0); PV_k += (int)(PV_inc) * ((PV_visiting[PV_l] & 4) ? 1 : 0); PV_l++; PV_visiting[PV_l] = 0; } else { PV_l--; PV_i -= (int)(PV_inc) * ((PV_visiting[PV_l] & 1) ? 1 : 0); PV_j -= (int)(PV_inc) * ((PV_visiting[PV_l] & 2) ? 1 : 0); PV_k -= (int)(PV_inc) * ((PV_visiting[PV_l] & 4) ? 1 : 0); PV_inc = PV_inc << 1; PV_node = ((PV_node)->parent); PV_visiting[PV_l]++; } } if (PV_visiting - 1) free(PV_visiting - 1); else {}; } }; }; } }; }; break; }; } };
