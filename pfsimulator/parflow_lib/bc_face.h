#ifndef _BC_FACE_H
#define _BC_FACE_H

/*********************************************
 *
 * Utility header file for constructing boundary condition face physics selection.
 * For use in conjunction with BCStructLoopX.
 *
 *********************************************/



/*********************************************
 *
 * Safe to call macros for defining BC Patch Face Physics
 *
 *********************************************/

#define Left GrGeomOctreeFaceL
#define Right GrGeomOctreeFaceR
#define Down GrGeomOctreeFaceD
#define Up GrGeomOctreeFaceU
#define Back GrGeomOctreeFaceB
#define Front GrGeomOctreeFaceF
#define AllBCTypes default
#define FACE(a,b) XCASE(a,b)
#define PROLOGUE(x) x
#define EPILOGUE(x) x
#define NO_PROLOGUE PROLOGUE({})
#define NO_EPILOGUE EPILOGUE({})
#define PRECONDITION(x) x
#define POSTCONDITION(x) x

/* Apply same physics for a set of boundary conditions */
#define AnyOf(_case, ...)                               \
  _case:                                                \
  EXPAND_CASES(COUNT_VARARGS(__VA_ARGS__))(__VA_ARGS__)

/* Apply physics for a boundary condition */
#define ApplyBCPatch(_case, ...)                \
  case _case:                                   \
  {                                             \
    EXPAND_FACE_PHYSICS(__VA_ARGS__);           \
    break;                                      \
  }

/* NOTE: This will be removed once additional OverlandFlow BC conditions are added */
#define ApplyBCPatchSubtypes(_case, _subcase, equations)  \
  case _case:                                             \
  {                                                       \
    switch ((_subcase))                                   \
    {                                                     \
      equations;                                          \
    }                                                     \
  }

/* Apply physics for a boundary condition with an action before the loop body */
/* EX: Applying kinematic or diffusive function call in nl_function_eval for
 * the OverlandBC case */
#define ApplyBCPatch_WithPrecondition(_case, precondition, ...) \
  case _case:                                                   \
  {                                                             \
    precondition;                                               \
    EXPAND_FACE_PHYSICS(__VA_ARGS__);                           \
    break;                                                      \
  }

/* Apply physics for a boundary condition with an action after the loop body */
/* EX: Applying kinematic or diffusive function call in richards_jacobian_eval for
   the OverlandBC case */
#define ApplyBCPatch_WithPostcondition(_case, postcondition, ...) \
  case _case:                                                     \
  {                                                               \
    EXPAND_FACE_PHYSICS(__VA_ARGS__);                             \
    postcondition;                                                \
    break;                                                        \
  }

/* To be called from GrGeomOctreeFaceLoopX */
#define OctreeFacePhysics(key, ...)              \
  XSWITCH(key, __VA_ARGS__)

/* Used in RichardsJacobianEval and NlFunctionEval to add
 *    Boundary Condition patch value contributions */
#define Do_BCContrib(i, j, k, ival, bc_struct, ipatch, is,          \
                              bc_patch_values, body)                \
  ForBCStructNumPatches(ipatch, bc_struct)                          \
  {                                                                 \
    bc_patch_values = BCStructPatchValues(bc_struct, ipatch, is);   \
    {                                                               \
      switch(BCStructBCType(bc_struct, ipatch))                     \
      {                                                             \
        body;                                                       \
      }                                                             \
    }                                                               \
  }

/*********************************************
 *
 * Do not call any of the below macros directly!
 *
 *********************************************/

/* Case statement with body */
#define WITH_BODY(_case, body)                  \
  case _case:                                   \
  {                                             \
    body;                                       \
    break;                                      \
  }

/* Case statement for fallthrough */
#define ONLY_CASE(_case)                        \
  case _case:

/* Variadic macros to determine what type of case to build */
#define CASE_BODY_SELECTION(arg1, arg2, arg3, ...) arg3
#define XCASE(...)                              \
  CASE_BODY_SELECTION(__VA_ARGS__, WITH_BODY, ONLY_CASE)(__VA_ARGS__)

/* Expand to switch on different boundary condition types */
#define XSWITCH(key, ...)                       \
  switch (key)                                  \
  {                                                       \
    EXPAND_FACES(COUNT_VARARGS(__VA_ARGS__))(__VA_ARGS__) \
  }

// NOTE: This macro will break if too many arguments are passed in
#define EVAL_COUNT_VARARGS(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, N, ...) N
#define COUNT_VARARGS(...) \
  EVAL_COUNT_VARARGS("not_evaluated", ##__VA_ARGS__, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)

/* Macro expansion indirection to get us the right EXPAND_FACES_N macro */
#define EXPAND_INDIRECTION(a, n) a ## n
#define EXPAND_FACES(n) EXPAND_INDIRECTION(EXPAND_FACES_, n)
#define EXPAND_CASES(n) EXPAND_INDIRECTION(EXPAND_CASES_, n)

#define EXPAND_CASES_1(a, ...) case a /* Drop last colon to prevent case a: : */
#define EXPAND_CASES_2(a, ...) case a: EXPAND_CASES_1(__VA_ARGS__)
#define EXPAND_CASES_3(a, ...) case a: EXPAND_CASES_2(__VA_ARGS__)
#define EXPAND_CASES_4(a, ...) case a: EXPAND_CASES_3(__VA_ARGS__)
#define EXPAND_CASES_5(a, ...) case a: EXPAND_CASES_4(__VA_ARGS__)
#define EXPAND_CASES_6(a, ...) case a: EXPAND_CASES_5(__VA_ARGS__)


/* This is a bit silly, might be a proper preprocessor way to generate these */
#define EXPAND_FACES_0(...) {};
#define EXPAND_FACES_1(a, ...) a
#define EXPAND_FACES_2(a, ...) a EXPAND_FACES_1(__VA_ARGS__)
#define EXPAND_FACES_3(a, ...) a EXPAND_FACES_2(__VA_ARGS__)
#define EXPAND_FACES_4(a, ...) a EXPAND_FACES_3(__VA_ARGS__)
#define EXPAND_FACES_5(a, ...) a EXPAND_FACES_4(__VA_ARGS__)
#define EXPAND_FACES_6(a, ...) a EXPAND_FACES_5(__VA_ARGS__)

/* Expand the boundary condition patch equations into the loop */
/* TODO: nl_function_eval has a case for Overland that uses a different loop than this.
 *       Will need to consider a different design?
 */
#define EXPAND_FACE_PHYSICS(prologue, epilogue, ...)  \
  BCStructPatchLoopX(i, j, k, ival, bc_struct, ipatch, is, prologue, epilogue, __VA_ARGS__);




/******************************
 *
 * Example of what the macro can expand to
 *
 ******************************/
#if 0

#define ApplyBCPatchExample                     \
  Do_ApplyBCPatch_Example(i, j, k, ival, bc_struct, ipatch, is, bc_patch_values, \
  {                                                                     \
    ApplyBCPatch(DirichletBC,                                           \
                 PROLOGUE({                                             \
                     ip = SubvectorEltIndex(p_sub, i, j, k);            \
                   }),                                                  \
                 EPILOGUE({                                             \
                     cp[im] = op[im];                                   \
                     op[im] = 0.0;                                      \
                   }),                                                  \
                 FACE(Left,  { /*Left Physics */ }),                    \
                 FACE(Right, { /*Right Physics */ }),                   \
                 FACE(Down,  { /*Down Physics */ }),                    \
                 FACE(Up,    { /*Up Physics */ }),                      \
                 FACE(Back,  { /*Back Physics */ }),                    \
                 FACE(Front, { /*Front Physics */ })                    \
      );                                                                \
  });

#define ApplyBCPatchTemplate                                            \
  ForBCStructNumPatches(ipatch, bc_struct)                              \
  {                                                                     \
    switch (BCStructPatchType(bc_struct, ipatch))                       \
    {                                                                   \
      /* ApplyBCPatch(DirichletBC, ...) */                              \
      case DirichletBC:                                                 \
      {                                                                 \
        BCStructPatchLoopX(i, j, k, ival, bc_struct, ipatch, is,        \
        {                                                               \
          /* PV_f and node come from GrGeomOctreeFaceLoopX expanded from BCStructLoopX */ \
          for (PV_f = 0; PV_f < GrGeomOctreeNumFaces; PV_f++)           \
            if (GrGeomOctreeHasFace(node, PV_f))                        \
            {                                                           \
              {                                                         \
                /* Prologue */                                          \
                ip = SubvectorEltIndex(p_sub, i, j, k);                 \
              }                                                         \
              switch (PV_f)                                             \
              {                                                         \
                case Left:                                              \
                {                                                       \
                  /* Left Physics */                                    \
                  break;                                                \
                }                                                       \
                case Right:                                             \
                {                                                       \
                  /* Right Physics */                                   \
                  break;                                                \
                }                                                       \
                case Down:                                              \
                {                                                       \
                  /* Down Physics */                                    \
                  break;                                                \
                }                                                       \
                case Up:                                                \
                {                                                       \
                  /* Up Physics */                                      \
                  break;                                                \
                }                                                       \
                case Back:                                              \
                {                                                       \
                  /* Back Physics */                                    \
                  break;                                                \
                }                                                       \
                case Front:                                             \
                {                                                       \
                  /* Front Physics */                                   \
                  break;                                                \
                }                                                       \
              }                                                         \
              {                                                         \
                /* Epilogue */                                          \
                cp[im] += op[im];                                       \
                op[im] = 0.0;                                           \
              }                                                         \
            }                                                           \
        });                                                             \
      }                                                                 \
    }                                                                   \
  }

#endif

#endif // _BC_FACE_H
