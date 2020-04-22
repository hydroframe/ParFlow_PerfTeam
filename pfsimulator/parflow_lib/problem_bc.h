/*BHEADER*********************************************************************
 *
 *  Copyright (c) 1995-2009, Lawrence Livermore National Security,
 *  LLC. Produced at the Lawrence Livermore National Laboratory. Written
 *  by the Parflow Team (see the CONTRIBUTORS file)
 *  <parflow@lists.llnl.gov> CODE-OCEC-08-103. All rights reserved.
 *
 *  This file is part of Parflow. For details, see
 *  http://www.llnl.gov/casc/parflow
 *
 *  Please read the COPYRIGHT file or Our Notice and the LICENSE file
 *  for the GNU Lesser General Public License.
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
 **********************************************************************EHEADER*/
/*****************************************************************************
*
*****************************************************************************/

#ifndef _PROBLEM_BC_HEADER
#define _PROBLEM_BC_HEADER

#define DirichletBC  0
#define FluxBC       1
#define OverlandBC   2   //sk
#define SeepageFaceBC   3   //rmm

/* @MCB: Additional overlandflow cases per LEC */
#define OverlandKinematicBC 4
#define OverlandDiffusiveBC 5

/*----------------------------------------------------------------
 * BCStruct structure
 *----------------------------------------------------------------*/

typedef struct {
  SubgridArray    *subgrids;   /* subgrids that BC data is defined on */

  GrGeomSolid     *gr_domain;

  int num_patches;
  int             *patch_indexes;  /* num_patches patch indexes */
  int             *bc_types;          /* num_patches BC types */

  double        ***values;   /* num_patches x num_subgrids data arrays */
} BCStruct;


/*--------------------------------------------------------------------------
 * Accessor macros:
 *--------------------------------------------------------------------------*/

#define BCStructSubgrids(bc_struct)           ((bc_struct)->subgrids)
#define BCStructNumSubgrids(bc_struct) \
  SubgridArrayNumSubgrids(BCStructSubgrids(bc_struct))
#define BCStructGrDomain(bc_struct)           ((bc_struct)->gr_domain)
#define BCStructNumPatches(bc_struct)         ((bc_struct)->num_patches)
#define BCStructPatchIndexes(bc_struct)       ((bc_struct)->patch_indexes)
#define BCStructBCTypes(bc_struct)            ((bc_struct)->bc_types)
#define BCStructValues(bc_struct)             ((bc_struct)->values)
#define BCStructSubgrid(bc_struct, i) \
  SubgridArraySubgrid(BCStructSubgrids(bc_struct), i)
#define BCStructPatchIndex(bc_struct, p)      ((bc_struct)->patch_indexes[p])
#define BCStructBCType(bc_struct, p)          ((bc_struct)->bc_types[p])
#define BCStructPatchValues(bc_struct, p, s)  ((bc_struct)->values[p][s])

/**
 * @brief Iterate over number of Boundary Condition patches
 *
 * @param ipatch Iterator variable to use
 * @param bc_struct BCStruct struct to read patch count from
 */
#define ForBCStructNumPatches(ipatch, bc_struct)										\
  for(ipatch = 0; ipatch < BCStructNumPatches(bc_struct); ipatch++)

/*--------------------------------------------------------------------------
 * Looping macro:
 *--------------------------------------------------------------------------*/

#define BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, body) \
  {                                                                         \
    GrGeomSolid  *PV_gr_domain = BCStructGrDomain(bc_struct);               \
    int PV_patch_index = BCStructPatchIndex(bc_struct, ipatch);             \
    Subgrid      *PV_subgrid = BCStructSubgrid(bc_struct, is);              \
                                                                            \
    int PV_r = SubgridRX(PV_subgrid);                                       \
    int PV_ix = SubgridIX(PV_subgrid);                                      \
    int PV_iy = SubgridIY(PV_subgrid);                                      \
    int PV_iz = SubgridIZ(PV_subgrid);                                      \
    int PV_nx = SubgridNX(PV_subgrid);                                      \
    int PV_ny = SubgridNY(PV_subgrid);                                      \
    int PV_nz = SubgridNZ(PV_subgrid);                                      \
                                                                            \
    ival = 0;                                                               \
    GrGeomPatchLoop(i, j, k, fdir, PV_gr_domain, PV_patch_index,            \
                    PV_r, PV_ix, PV_iy, PV_iz, PV_nx, PV_ny, PV_nz,         \
    {                                                                       \
      body;                                                                 \
      ival++;                                                               \
    });                                                                     \
  }

#define BCStructPatchLoopOvrlnd(i, j, k, fdir, ival, bc_struct, ipatch, is, body) \
  {                                                                               \
    GrGeomSolid  *PV_gr_domain = BCStructGrDomain(bc_struct);                     \
    int PV_patch_index = BCStructPatchIndex(bc_struct, ipatch);                   \
    Subgrid      *PV_subgrid = BCStructSubgrid(bc_struct, is);                    \
                                                                                  \
    int PV_r = SubgridRX(PV_subgrid);                                             \
    int PV_ix = SubgridIX(PV_subgrid) - 1;                                        \
    int PV_iy = SubgridIY(PV_subgrid) - 1;                                        \
    int PV_iz = SubgridIZ(PV_subgrid) - 1;                                        \
    int PV_nx = SubgridNX(PV_subgrid) + 2;                                        \
    int PV_ny = SubgridNY(PV_subgrid) + 2;                                        \
    int PV_nz = SubgridNZ(PV_subgrid) + 2;                                        \
                                                                                  \
    ival = 0;                                                                     \
    GrGeomPatchLoop(i, j, k, fdir, PV_gr_domain, PV_patch_index,                  \
                    PV_r, PV_ix, PV_iy, PV_iz, PV_nx, PV_ny, PV_nz,               \
    {                                                                             \
      body;                                                                       \
      ival++;                                                                     \
    });                                                                           \
  }

/**
 * @brief Variant of BCStructPatchLoop that doesn't utilize an fdir variable, used in ForPatchCells loops.
 *
 * @note Do not call directly! Not intended for use code.
 */
#define BCStructPatchLoopNoFdir(i, j, k, ival, bc_struct, ipatch, is,		\
                                locals, setup,                          \
                                f_left, f_right,                        \
                                f_down, f_up,                           \
                                f_back, f_front,                        \
                                finalize)                               \
  {                                                                     \
    GrGeomSolid  *PV_gr_domain = BCStructGrDomain(bc_struct);           \
    int PV_patch_index = BCStructPatchIndex(bc_struct, ipatch);         \
    Subgrid      *PV_subgrid = BCStructSubgrid(bc_struct, is);          \
                                                                        \
    int PV_r = SubgridRX(PV_subgrid);                                   \
    int PV_ix = SubgridIX(PV_subgrid);                                  \
    int PV_iy = SubgridIY(PV_subgrid);                                  \
    int PV_iz = SubgridIZ(PV_subgrid);                                  \
    int PV_nx = SubgridNX(PV_subgrid);                                  \
    int PV_ny = SubgridNY(PV_subgrid);                                  \
    int PV_nz = SubgridNZ(PV_subgrid);                                  \
                                                                        \
    ival = 0;                                                           \
    GrGeomPatchLoopNoFdir(i, j, k, PV_gr_domain, PV_patch_index,        \
                          PV_r, PV_ix, PV_iy, PV_iz, PV_nx, PV_ny, PV_nz, \
                          locals, setup,                                \
                          f_left, f_right,                              \
                          f_down, f_up,                                 \
                          f_back, f_front,                              \
    {                                                                   \
      finalize;                                                         \
      ival++;                                                           \
    });                                                                 \
  }

#define BCStructPatchLoopOvrlndNoFdir(i, j, k, ival, bc_struct, ipatch, is, \
                                      locals, setup,                    \
                                      f_left, f_right,                  \
                                      f_down, f_up,                     \
                                      f_back, f_front,                  \
                                      finalize)                         \
  {                                                                     \
    GrGeomSolid  *PV_gr_domain = BCStructGrDomain(bc_struct);           \
    int PV_patch_index = BCStructPatchIndex(bc_struct, ipatch);         \
    Subgrid      *PV_subgrid = BCStructSubgrid(bc_struct, is);          \
                                                                        \
    int PV_r = SubgridRX(PV_subgrid);                                   \
    int PV_ix = SubgridIX(PV_subgrid) - 1;                              \
    int PV_iy = SubgridIY(PV_subgrid) - 1;                              \
    int PV_iz = SubgridIZ(PV_subgrid) - 1;                              \
    int PV_nx = SubgridNX(PV_subgrid) + 2;                              \
    int PV_ny = SubgridNY(PV_subgrid) + 2;                              \
    int PV_nz = SubgridNZ(PV_subgrid) + 2;                              \
                                                                        \
    ival = 0;                                                           \
    GrGeomPatchLoopNoFdir(i, j, k, PV_gr_domain, PV_patch_index,        \
                          PV_r, PV_ix, PV_iy, PV_iz, PV_nx, PV_ny, PV_nz, \
                          locals, setup,                                \
                          f_left, f_right,                              \
                          f_down, f_up,                                 \
                          f_back, f_front,                              \
    {                                                                   \
      finalize;                                                         \
      ival++;                                                           \
    });                                                                 \
  }

/*--------------------------------------------------------------------------
 * ForPatch loops and macros
 *--------------------------------------------------------------------------*/

/**
 * @brief For use when a BCLoop should execute no matter what the patch type is
 */
#define ALL -1

/**
 * @brief For use when a statement body is unnecessary in ForPatchCellsPerFace.
 */
#define DoNothing

/**
 * @name Face directions
 * @brief Declare which face a given FACE body should apply to
 *
 * These are aliases of the internal GrGeomOctreeFace definitions
 * @{
 */
#define LeftFace GrGeomOctreeFaceL
#define RightFace GrGeomOctreeFaceR
#define DownFace GrGeomOctreeFaceD
#define UpFace GrGeomOctreeFaceU
#define BackFace GrGeomOctreeFaceB
#define FrontFace GrGeomOctreeFaceF
/** @} */

/**
 * @brief Unconditionally executes body at the beginning of each boundary cell iteration
 */
#define CellSetup(body) { body; };

/**
 * @brief Unconditionally executes body at the end of each boundary cell iteration
 */
#define CellFinalize(body) { body; };

/**
 * @brief Unconditionally executes body before all boundary cell iterations
 */
#define BeforeAllCells(body) { body; };

/**
 * @brief Unconditionally executes body after all boundary cell iterations
 */
#define AfterAllCells(body) { body; }

/**
 * @brief Used to pass loop variables (ex: i, j, k, ival, etc.) through to inner loop macros
 */
#define LoopVars(...) __VA_ARGS__

/**
 * @brief Packs arbitrary number of statements into paranthesis to pass through to other macros.
 *
 * Used in ForPatchCellsPerFace loop to declare local variables.  Provides architecture portability and scope safety.
 * Expand the packed arguments using the UNPACK macro at the appropriate location in code.
 */
#define Locals(...) (__VA_ARGS__)

/**
 * @brief For use when no locals variables need to be declared in a loop.
 */
#define NoLocals ()
#define UNPACK(locals) _UNPACK locals
#define _UNPACK(...) __VA_ARGS__ ;

/**
 * @brief For use with ForPatchCellsPerFace loop, manages conditional branching internally
 *
 * Expands to a case statement containing the statement body.
 *
 * @param[in] fdir Face direction to execute the body on (e.g. FaceLeft)
 * @param[in] body Arbitrary statement body block to execute
 */
#define FACE(fdir, body)      \
  case fdir:                  \
  {                           \
    body;                     \
    break;                    \
  }

#if 0
/* @MCB: Template for copy+pasting for new BC loops.
 Body blocks intentionally contain ~~ as a statement to cause syntax errors if left alone.
 This is to ensure each cell body statement has been set and not left as a copy+paste.
 */
ForPatchCellsPerFace(InsertBCTypeHere,
                     BeforeAllCells(DoNothing),
                     LoopVars(i, j, k, ival, bc_struct, ipatch, is),
                     NoLocals,
                     CellSetup({ ~~ }),
                     FACE(LeftFace, { ~~ }),
                     FACE(RightFace, { ~~ }),
                     FACE(DownFace, { ~~ }),
                     FACE(UpFace, { ~~ }),
                     FACE(BackFace, { ~~ }),
                     FACE(FrontFace, { ~~ }),
                     CellFinalize({ ~~ }),
                     AfterAllCells(DoNothing)
  );
#endif

/**
 * @brief Iterates over the cells of a boundary condition patch with conditional branching on each face direction
 *
 *
 *
 * @param bctype Boundary condition type this loop should apply computations to.  (e.g. OverlandBC)
 * @param before_loop See BeforeAllCells macro.
 * @param loopvars Variables to pass through to inner looping macro.  See LoopVars macro.
 * @param locals Locally defined variables for use in loop.  See Locals macro.
 * @param setup See CellSetup macro
 * @param f_left Statement body to execute on Left cell face.  See FACE macro.
 * @param f_right Statement body to execute on Right cell face.  See FACE macro.
 * @param f_down Statement body to execute on Down cell face.  See FACE macro.
 * @param f_up Statement body to execute on Up cell face.  See FACE macro.
 * @param f_back Statement body to execute on Back cell face.  See FACE macro.
 * @param f_front Statement body to execute on Front cell face.  See FACE macro.
 * @param finalize See CellFinalize macro
 * @param after_loop See AfterAllCells macro
 */
#define ForPatchCellsPerFace(bctype,													\
                             before_loop,											\
														 loopvars, locals,								\
                             setup,                           \
                             f_left, f_right,                 \
                             f_down, f_up,                    \
                             f_back, f_front,                 \
                             finalize,                        \
                             after_loop)                      \
  {                                                           \
    if ( ((bctype) == ALL) ||                                 \
         ((bctype) == _GetCurrentPatch(loopvars)))            \
    {                                                         \
      before_loop;                                            \
      BCStructPatchLoopNoFdir(loopvars,                       \
                              locals, setup,                  \
                              f_left, f_right,                \
                              f_down, f_up,                   \
                              f_back, f_front,                \
                              finalize);                      \
      after_loop;                                             \
    }                                                         \
  }

/**
 * @brief Variation of ForPatchCellsPerFace that extends loop bounds to include ghost cells
 * See ForPatchCellsPerFace for further documentation
 */
#define ForPatchCellsPerFaceWithGhost(bctype, \
                                      before_loop, loopvars,          \
																			locals,													\
                                      setup,                          \
                                      f_left, f_right,                \
                                      f_down, f_up,                   \
                                      f_back, f_front,                \
                                      finalize,                       \
                                      after_loop)                     \
  {                                                                   \
    if ( ((bctype) == ALL) ||                                         \
         ((bctype) == _GetCurrentPatch(loopvars)))                    \
    {                                                                 \
      before_loop;                                                    \
      BCStructPatchLoopOvrlndNoFdir(loopvars,                         \
																		locals,														\
                                    setup,                            \
                                    f_left, f_right,                  \
                                    f_down, f_up,                     \
                                    f_back, f_front,                  \
                                    finalize);                        \
      after_loop;                                                     \
    }                                                                 \
  }

/**
 * @brief Iterates over cells of a boundary condition patch without any conditional branching
 *
 * @note Any local variable definitions should be made explicitly in the body of the loop
 *
 * @param[in,out] i X index
 * @param[in,out] j Y index
 * @param[in,out] k Z index
 * @param[in,out] ival Index for patch value array
 * @param bc_struct BCStruct to use in internal looping structures
 * @param ipatch Current patch number
 * @param is Current subgrid number
 * @param body Statement body to execute
*/
#define ForEachPatchCell(i, j, k, ival, bc_struct, ipatch, is, body)		\
  BCStructPatchLoopNoFdir(i, j, k, ival, bc_struct, ipatch, is,					\
													NoLocals,																			\
													body,																					\
                          DoNothing, DoNothing,													\
                          DoNothing, DoNothing,													\
                          DoNothing, DoNothing,													\
                          DoNothing);

#endif
