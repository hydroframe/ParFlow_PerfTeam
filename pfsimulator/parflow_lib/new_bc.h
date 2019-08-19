#ifndef _NEW_BC_H
#define _NEW_BC_H

#define Left GrGeomOctreeFaceL
#define Right GrGeomOctreeFaceR
#define Down GrGeomOctreeFaceD
#define Up GrGeomOctreeFaceU
#define Back GrGeomOctreeFaceB
#define Front GrGeomOctreeFaceF

#define FACE(a,b) XCASE(a,b)
#define PROLOGUE(x) x
#define EPILOGUE(x) x
#define NO_PROLOGUE PROLOGUE({})
#define NO_EPILOGUE EPILOGUE({})
#define PRECONDITION(x) x
#define POSTCONDITION(x) x

#define AnyOf(_case, ...)                       \
  _case:                                        \
  EXPAND_CASES(COUNT_VARARGS(__VA_ARGS__))(__VA_ARGS__)

#define ApplyBCPatch(_case, ...)                \
  case _case:                                   \
  {                                             \
    EXPAND_FACE_PHYSICS(__VA_ARGS__);           \
    break;                                      \
  }

#define ApplyBCPatchSubtypes(_case, _subcase, equations)  \
  case _case:                                             \
  {                                                       \
    switch ((_subcase))                                   \
    {                                                     \
      equations;                                          \
    }                                                     \
  }

#define ApplyBCPatch_WithPrecondition(_case, precondition, ...) \
  case _case:                                                   \
  {                                                             \
    precondition;                                               \
    EXPAND_FACE_PHYSICS(__VA_ARGS__);                           \
    break;                                                      \
  }

#define ApplyBCPatch_WithPostcondition(_case, postcondition, ...) \
  case _case:                                                   \
  {                                                             \
    EXPAND_FACE_PHYSICS(__VA_ARGS__);                           \
    postcondition;                                              \
    break;                                                      \
  }

#define OctreeFacePhysics(key, ...)             \
  XSWITCH(key, __VA_ARGS__)


#define Do_BCContrib(i, j, k, ival, bc_struct, ipatch, is, bc_patch_values, body) \
  ForBCStructNumPatches(bc_struct, ipatch)                                     \
  {                                                                     \
    bc_patch_values = BCStructPatchValues(bc_struct, ipatch, is);       \
    {                                                                   \
      switch(BCStructBCType(bc_struct, ipatch))                         \
      {                                                                 \
        body;                                                           \
      }                                                                 \
    }                                                                   \
  }

#define XSWITCH(key, ...)                       \
  switch ((key))                                \
  {                                                         \
    EXPAND_FACES(COUNT_VARARGS(__VA_ARGS__))(__VA_ARGS__)  \
  }

#define EVAL_COUNT_VARARGS(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, N, ...) N
#define COUNT_VARARGS(...)                      \
  EVAL_COUNT_VARARGS("not_evaluated", ##__VA_ARGS__, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)

#define EXPAND_INDIRECT(a, n) a ## n
#define EXPAND_FACES(n) EXPAND_INDIRECT(EXPAND_FACES_, n)
#define EXPAND_CASES(n) EXPAND_INDIRECT(EXPAND_CASES_, n)

#define EXPAND_CASES_1(a, ...) case a
#define EXPAND_CASES_2(a, ...) case a: EXPAND_CASES_1(__VA_ARGS__)
#define EXPAND_CASES_3(a, ...) case a: EXPAND_CASES_2(__VA_ARGS__)
#define EXPAND_CASES_4(a, ...) case a: EXPAND_CASES_3(__VA_ARGS__)
#define EXPAND_CASES_5(a, ...) case a: EXPAND_CASES_4(__VA_ARGS__)
#define EXPAND_CASES_6(a, ...) case a: EXPAND_CASES_5(__VA_ARGS__)
#define EXPAND_CASES_7(a, ...) case a: EXPAND_CASES_6(__VA_ARGS__)

#define EXPAND_FACES_0(...) {};
#define EXPAND_FACES_1(a, ...) a
#define EXPAND_FACES_2(a, ...) a EXPAND_FACES_1(__VA_ARGS__)
#define EXPAND_FACES_3(a, ...) a EXPAND_FACES_2(__VA_ARGS__)
#define EXPAND_FACES_4(a, ...) a EXPAND_FACES_3(__VA_ARGS__)
#define EXPAND_FACES_5(a, ...) a EXPAND_FACES_4(__VA_ARGS__)
#define EXPAND_FACES_6(a, ...) a EXPAND_FACES_5(__VA_ARGS__)

#define EXPAND_FACE_PHYSICS(prologue, epilogue, ...)  \
  BCStructPatchLoopNoFdir(i, j, k, ival, bc_struct, ipatch, is, prologue, epilogue, __VA_ARGS__);

#define WITH_BODY(_case, body);                 \
  case _case:                                   \
  {                                             \
    body;                                       \
    break;                                      \
  }

#define ONLY_CASE(_case)                        \
  case _case:                                   \

#define CASE_BODY_SELECTION(arg1, arg2, arg3, ...) arg3
#define XCASE(...)                              \
  CASE_BODY_SELECTION(__VA_ARGS__, WITH_BODY, ONLY_CASE)(__VA_ARGS__)

#endif // _NEW_BC_H
