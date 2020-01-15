#ifndef _BC_LOOP_UTIL_H
#define _BC_LOOP_UTIL_H

#define Left GrGeomOctreeFaceL
#define Right GrGeomOctreeFaceR
#define Down GrGeomOctreeFaceD
#define Up GrGeomOctreeFaceU
#define Back GrGeomOctreeFaceB
#define Front GrGeomOctreeFaceF

#define AllBCTypes default
#define FACE(a, b) XCASE(a, b)
#define PROLOGUE(x) x
#define EPILOGUE(x) x
#define NO_PROLOGUE {}
#define NO_EPILOGUE {}

#define PRECONDITION(x) x
#define POSTCONDITION(x) x
#define NO_PRECOND {}
#define NO_POSTCOND {}

#define LOOP_VARS(...) __VA_ARGS__

#define OctreeFacePhysics(key, ...) XSWITCH(key, __VA_ARGS__)

#define ApplyBCPatch(_case, loop_vars, precond, postcond, prologue, epilogue, ...) \
  case _case:                                                           \
  {                                                                     \
    precond;                                                            \
    BCStructPatchLoopNoFdir(loop_vars, prologue, epilogue, __VA_ARGS__); \
    postcond;                                                           \
    break;                                                              \
  }

#define ApplyAllPatches(loop_vars, ...)               \
  {                                                   \
    BCStructPatchLoopNoFdir(loop_vars, __VA_ARGS__);  \
  }

#define AnyOf(_case, ...) \
  _case:                  \
  EXPAND_CASES(COUNT_VARARGS(__VA_ARGS__))(__VA_ARGS__)


#define WITH_BODY(_case, body) \
  case _case:                  \
  {                            \
    body;                      \
    break;                     \
  }

#define ONLY_CASE(_case) case _case:

#define CASE_BODY_SELECTION(arg1, arg2, arg3, ...) arg3
#define XCASE(...) CASE_BODY_SELECTION(__VA_ARGS__, WITH_BODY, ONLY_CASE)(__VA_ARGS__)

#define XSWITCH(key, ...)                                 \
  switch (key)                                            \
  {                                                       \
    EXPAND_FACES(COUNT_VARARGS(__VA_ARGS__))(__VA_ARGS__) \
                                                          \
    default:                                              \
    continue;                                             \
  }

#define EVAL_COUNT_VARARGS(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, N, ...) N
#define COUNT_VARARGS(...)                                              \
  EVAL_COUNT_VARARGS("not_evaluated", ##__VA_ARGS__, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)

#define EXPAND_INDIRECT(a, n) a ## n
#define EXPAND_FACES(n) EXPAND_INDIRECT(EXPAND_FACES_, n)
#define EXPAND_CASES(n) EXPAND_INDIRECT(EXPAND_CASES_, n)

/*
  Used with the AnyOf macro.  EX: AnyOf(DirichletBC, OverlandBC)
  Unlikely we'll need more than 7.  If there are extras, the EXPAND_CASES_1 will
  spit out the extra symbols and give a compiler error.
*/
#define EXPAND_CASES_1(a, ...) case a __VA_ARGS__
#define EXPAND_CASES_2(a, ...) case a: EXPAND_CASES_1(__VA_ARGS__)
#define EXPAND_CASES_3(a, ...) case a: EXPAND_CASES_2(__VA_ARGS__)
#define EXPAND_CASES_4(a, ...) case a: EXPAND_CASES_3(__VA_ARGS__)
#define EXPAND_CASES_5(a, ...) case a: EXPAND_CASES_4(__VA_ARGS__)
#define EXPAND_CASES_6(a, ...) case a: EXPAND_CASES_5(__VA_ARGS__)
#define EXPAND_CASES_7(a, ...) case a: EXPAND_CASES_6(__VA_ARGS__)

/* @MCB: We only have 6 directions, so this should be fine for now */
#define EXPAND_FACES_0(...) {};
#define EXPAND_FACES_1(a, ...) a
#define EXPAND_FACES_2(a, ...) a EXPAND_FACES_1(__VA_ARGS__)
#define EXPAND_FACES_3(a, ...) a EXPAND_FACES_2(__VA_ARGS__)
#define EXPAND_FACES_4(a, ...) a EXPAND_FACES_3(__VA_ARGS__)
#define EXPAND_FACES_5(a, ...) a EXPAND_FACES_4(__VA_ARGS__)
#define EXPAND_FACES_6(a, ...) a EXPAND_FACES_5(__VA_ARGS__)
#define EXPAND_FACES_7(a, ...) a EXPAND_FACES_6(__VA_ARGS__)

#endif // _BC_LOOP_UTIL_H
