#ifndef _PF_OMP_BCLOOPS_H
#define _PF_OMP_BCLOOPS_H

/* Used to calculate ival */
#define CALC_IVAL(diff, a, b) ((diff) * (a) + (a) + (b))

#define _GetBCDomain(i, j, k, ival, bc_struct, ipatch, is)\
  BCStructGrDomain(bc_struct)

#define _GetBCIndex(i, j, k, ival, bc_struct, ipatch, is)  \
  BCStructPatchIndex(bc_struct, ipatch)

#define _GetBCSubgrid(i, j, k, ival, bc_struct, ipatch, is)  \
  BCStructSubgrid(bc_struct, is)

#define _SetIvalZero(i, j, k, ival, bc_struct, ipatch, is)  \
  ival = 0

#define _GetIJK(i, j, k, ival, bc_struct, ipatch, is) i, j, k

#undef ForPatchCellsPerFace
#define ForPatchCellsPerFace(bctype, pragma, locals,                    \
                             before_loop,                               \
                             loopvars,                                  \
                             setup,                                     \
                             f_left, f_right,                           \
                             f_down, f_up,                              \
                             f_back, f_front,                           \
                             finalize,                                  \
                             after_loop)                                \
  {                                                                     \
    if ( ((bctype) == ALL) ||                                           \
         ((bctype) == _GetCurrentPatch(loopvars)) )                     \
    {                                                                   \
      GrGeomSolid  *PV_gr_domain = _GetBCDomain(loopvars);              \
      int PV_patch_index = _GetBCIndex(loopvars);                       \
      Subgrid      *PV_subgrid = _GetBCSubgrid(loopvars);               \
                                                                        \
      int PV_r = SubgridRX(PV_subgrid);                                 \
      int PV_ix = SubgridIX(PV_subgrid);                                \
      int PV_iy = SubgridIY(PV_subgrid);                                \
      int PV_iz = SubgridIZ(PV_subgrid);                                \
      int PV_nx = SubgridNX(PV_subgrid);                                \
      int PV_ny = SubgridNY(PV_subgrid);                                \
      int PV_nz = SubgridNZ(PV_subgrid);                                \
                                                                        \
      _SetIvalZero(loopvars);                                           \
      if (PV_r != 0 ||                                                  \
          !GrGeomSolidPatchBoxes(PV_gr_domain, PV_patch_index, GrGeomOctreeNumFaces -1)) \
      {                                                                 \
        if(!amps_Rank(amps_CommWorld))                                  \
        {                                                               \
          amps_Printf("Use of OpenMP requires box clustering! Aborting from:\n%s:%d\n", \
                      __FILE__, __LINE__);                              \
        }                                                               \
        exit(1);                                                        \
      }                                                                 \
                                                                        \
      pragma ## _GrGeomPatchLoopBoxesNoFdir(locals,                     \
                                            loopvars,                   \
                                            PV_gr_domain, PV_patch_index, \
                                            PV_r, PV_ix, PV_iy, PV_iz,  \
                                            PV_nx, PV_ny, PV_nz,        \
                                            before_loop,                \
                                            setup,                      \
                                            f_left, f_right,            \
                                            f_down, f_up,               \
                                            f_back, f_front,            \
                                            finalize,                   \
                                            after_loop);                \
    }                                                                   \
  }

#undef ForPatchCellsPerFaceWithGhost
#define ForPatchCellsPerFaceWithGhost(bctype, pragma, locals,                    \
                                      before_loop,                      \
                                      loopvars,                         \
                                      setup,                            \
                                      f_left, f_right,                  \
                                      f_down, f_up,                     \
                                      f_back, f_front,                  \
                                      finalize,                         \
                                      after_loop)                       \
  {                                                                     \
    if ( ((bctype) == ALL) ||                                           \
         ((bctype) == _GetCurrentPatch(loopvars)) )                     \
    {                                                                   \
      GrGeomSolid  *PV_gr_domain = _GetBCDomain(loopvars);              \
      int PV_patch_index = _GetBCIndex(loopvars);                       \
      Subgrid      *PV_subgrid = _GetBCSubgrid(loopvars);               \
                                                                        \
      int PV_r = SubgridRX(PV_subgrid);                                 \
      int PV_ix = SubgridIX(PV_subgrid) - 1;                            \
      int PV_iy = SubgridIY(PV_subgrid) - 1;                            \
      int PV_iz = SubgridIZ(PV_subgrid) - 1;                            \
      int PV_nx = SubgridNX(PV_subgrid) + 2;                            \
      int PV_ny = SubgridNY(PV_subgrid) + 2;                            \
      int PV_nz = SubgridNZ(PV_subgrid) + 2;                            \
                                                                        \
      _SetIvalZero(loopvars);                                           \
      if (PV_r != 0 ||                                                  \
          !GrGeomSolidPatchBoxes(PV_gr_domain, PV_patch_index, GrGeomOctreeNumFaces -1)) \
      {                                                                 \
        if(!amps_Rank(amps_CommWorld))                                  \
        {                                                               \
          amps_Printf("Use of OpenMP requires box clustering! Aborting from:\n%s:%d\n", \
                      __FILE__, __LINE__);                              \
        }                                                               \
        exit(1);                                                        \
      }                                                                 \
                                                                        \
      pragma ## _GrGeomPatchLoopBoxesNoFdir(locals,                     \
                                            loopvars,                   \
                                            PV_gr_domain, PV_patch_index, \
                                            PV_r, PV_ix, PV_iy, PV_iz,  \
                                            PV_nx, PV_ny, PV_nz,        \
                                            before_loop,                \
                                            setup,                      \
                                            f_left, f_right,            \
                                            f_down, f_up,               \
                                            f_back, f_front,            \
                                            finalize,                   \
                                            after_loop);                \
    }                                                                   \
  }

#undef ForEachPatchCell
#define ForEachPatchCell(pragma, locals, i, j, k, ival, bc_struct, ipatch, is, body) \
  ForPatchCellsPerFace(ALL, pragma, locals, DoNothing,\
                       LoopVars(i, j, k, ival, bc_struct, ipatch, is), \
                       body,                                           \
                       DoNothing, DoNothing,                           \
                       DoNothing, DoNothing,                           \
                       DoNothing, DoNothing,                           \
                       DoNothing,                                      \
                       DoNothing);

#define NewParallel_GrGeomPatchLoopBoxesNoFdir(locals, i, j, k, ival, bc_struct, ipatch, is, \
                                               PV_gr_domain, PV_patch_index, \
                                               PV_r, PV_ix, PV_iy, PV_iz, \
                                               PV_nx, PV_ny, PV_nz,     \
                                               before_loop,             \
                                               setup,                   \
                                               f_left, f_right,         \
                                               f_down, f_up,            \
                                               f_back, f_front,         \
                                               finalize,                \
                                               after_loop)              \
  PRAGMA(omp parallel firstprivate(i, j, k, ival locals))               \
  {                                                                     \
    before_loop;                                                        \
    MainBody_GrGeomPatchLoopBoxesNoFdir(for,                            \
                                        i, j, k,                        \
                                        PV_gr_domain, PV_patch_index,   \
                                        PV_r, PV_ix, PV_iy, PV_iz,      \
                                        PV_nx, PV_ny, PV_nz,            \
                                        setup,                          \
                                        f_left, f_right,                \
                                        f_down, f_up,                   \
                                        f_back, f_front,                \
                                        finalize);                      \
    after_loop;                                                         \
  }

#define InParallel_GrGeomPatchLoopBoxesNoFdir(locals, i, j, k, ival, bc_struct, ipatch, is, \
                                              PV_gr_domain, PV_patch_index, \
                                              PV_r, PV_ix, PV_iy, PV_iz, \
                                              PV_nx, PV_ny, PV_nz,      \
                                              before_loop,              \
                                              setup,                    \
                                              f_left, f_right,          \
                                              f_down, f_up,             \
                                              f_back, f_front,          \
                                              finalize,                 \
                                              after_loop)               \
  {                                                                     \
    before_loop;                                                        \
    MainBody_GrGeomPatchLoopBoxesNoFdir(for, i, j, k,                   \
                                        PV_gr_domain, PV_patch_index,   \
                                        PV_r, PV_ix, PV_iy, PV_iz,      \
                                        PV_nx, PV_ny, PV_nz,            \
                                        setup,                          \
                                        f_left, f_right,                \
                                        f_down, f_up,                   \
                                        f_back, f_front,                \
                                        finalize);                      \
    after_loop;                                                         \
  }

#define NoWait_GrGeomPatchLoopBoxesNoFdir(locals, i, j, k, ival, bc_struct, ipatch, is, \
                                          PV_gr_domain, PV_patch_index, \
                                          PV_r, PV_ix, PV_iy, PV_iz,    \
                                          PV_nx, PV_ny, PV_nz,          \
                                          before_loop,                  \
                                          setup,                        \
                                          f_left, f_right,              \
                                          f_down, f_up,                 \
                                          f_back, f_front,              \
                                          finalize,                     \
                                          after_loop)                   \
  {                                                                     \
    before_loop;                                                        \
    MainBody_GrGeomPatchLoopBoxesNoFdir(for nowait, i, j, k,            \
                                        PV_gr_domain, PV_patch_index,   \
                                        PV_r, PV_ix, PV_iy, PV_iz,      \
                                        PV_nx, PV_ny, PV_nz,            \
                                        setup,                          \
                                        f_left, f_right,                \
                                        f_down, f_up,                   \
                                        f_back, f_front,                \
                                        finalize);                      \
    after_loop;                                                         \
  }


#define MainBody_GrGeomPatchLoopBoxesNoFdir(pragma_clause,              \
                                            i, j, k,                    \
                                            grgeom, patch_num,          \
                                            r, ix, iy, iz,              \
                                            nx, ny, nz,                 \
                                            setup,                      \
                                            f_left, f_right,            \
                                            f_down, f_up,               \
                                            f_back, f_front,            \
                                            finalize)                   \
  {                                                                     \
    int PV_ixl, PV_iyl, PV_izl, PV_ixu, PV_iyu, PV_izu;                 \
    int *PV_visiting = NULL;                                            \
    for (int PV_f = 0; PV_f < GrGeomOctreeNumFaces; PV_f++)             \
    {                                                                   \
      BoxArray* boxes = GrGeomSolidPatchBoxes(grgeom, patch_num, PV_f); \
      for (int PV_box = 0; PV_box < BoxArraySize(boxes); PV_box++)      \
      {                                                                 \
        Box box = BoxArrayGetBox(boxes, PV_box);                        \
        /* find octree and region intersection */                       \
        PV_ixl = pfmax(ix, box.lo[0]);                                  \
        PV_iyl = pfmax(iy, box.lo[1]);                                  \
        PV_izl = pfmax(iz, box.lo[2]);                                  \
        PV_ixu = pfmin((ix + nx - 1), box.up[0]);                       \
        PV_iyu = pfmin((iy + ny - 1), box.up[1]);                       \
        PV_izu = pfmin((iz + nz - 1), box.up[2]);                       \
                                                                        \
        int PV_diff_x = PV_ixu - PV_ixl;                                \
        int PV_diff_y = PV_iyu - PV_iyl;                                \
        int PV_diff_z = PV_izu - PV_izl;                                \
        int x_scale = !!PV_diff_x;                                      \
        int y_scale = !!PV_diff_y;                                      \
        int z_scale = !!PV_diff_z;                                      \
        if (PV_diff_x * PV_diff_y * PV_diff_z != 0) {                   \
          fprintf(stderr, "ERROR: Diff not 0 at %s %d\n", __FILE__, __LINE__); \
          exit(-1);                                                     \
        }                                                               \
        PRAGMA(omp pragma_clause collapse(3) private(i, j, k, ival))    \
          for (k = PV_izl; k <= PV_izu; k++)                            \
          {                                                             \
            for (j = PV_iyl; j <= PV_iyu; j++)                          \
            {                                                           \
              for (i = PV_ixl; i <= PV_ixu; i++)                        \
              {                                                         \
                int PV_tmp_i = i - PV_ixl;                              \
                int PV_tmp_j = j - PV_iyl;                              \
                int PV_tmp_k = k - PV_izl;                              \
                if (!z_scale) {                                         \
                  ival = CALC_IVAL(PV_diff_x, PV_tmp_j, PV_tmp_i);      \
                } else if (!y_scale) {                                  \
                  ival = CALC_IVAL(PV_diff_x, PV_tmp_k, PV_tmp_i);      \
                } else {                                                \
                  ival = CALC_IVAL(PV_diff_y, PV_tmp_k, PV_tmp_j);      \
                }                                                       \
                setup;                                                  \
                                                                        \
                switch(PV_f)                                            \
                {                                                       \
                  f_left;                                               \
                  f_right;                                              \
                  f_down;                                               \
                  f_up;                                                 \
                  f_back;                                               \
                  f_front;                                              \
                }                                                       \
                                                                        \
                finalize;                                               \
              }                                                         \
            }                                                           \
          }                                                             \
      }                                                                 \
    }                                                                   \
  }

#undef _BCStructPatchLoop
#define _BCStructPatchLoop(pragma, locals,                              \
                           i, j, k, fdir, ival, bc_struct,              \
                           ipatch, is, body)                            \
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
    if (PV_r != 0 ||                                                    \
        !GrGeomSolidPatchBoxes(PV_gr_domain, PV_patch_index, GrGeomOctreeNumFaces -1)) \
    {                                                                   \
      if(!amps_Rank(amps_CommWorld))                                    \
      {                                                                 \
        amps_Printf("Use of OpenMP requires box clustering! Aborting from:\n%s:%d\n", __FILE__, __LINE__); \
      }                                                                 \
      exit(1);                                                          \
    }                                                                   \
    pragma ## _GrGeomPatchLoopBoxes(locals, i, j, k,                    \
                                    fdir, PV_gr_domain, PV_patch_index, \
                                    PV_r, PV_ix, PV_iy, PV_iz,          \
                                    PV_nx, PV_ny, PV_nz, body);         \
  }

#undef _BCStructPatchLoopOvrlnd
#define _BCStructPatchLoopOvrlnd(pragma, locals,                        \
                                 i, j, k, fdir, ival, bc_struct,        \
                                 ipatch, is, body)                      \
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
    if (PV_r != 0 ||                                                    \
        !GrGeomSolidPatchBoxes(PV_gr_domain, PV_patch_index, GrGeomOctreeNumFaces -1)) \
    {                                                                   \
      if(!amps_Rank(amps_CommWorld))                                    \
      {                                                                 \
        amps_Printf("Use of OpenMP requires box clustering! Aborting from:\n%s:%d\n", __FILE__, __LINE__); \
      }                                                                 \
      exit(1);                                                          \
    }                                                                   \
    pragma ## _GrGeomPatchLoopBoxes(locals, i, j, k,                    \
                                    fdir, PV_gr_domain, PV_patch_index, \
                                    PV_r, PV_ix, PV_iy, PV_iz,          \
                                    PV_nx, PV_ny, PV_nz, body);         \
  }

#define MainBody_GrGeomPatchLoopBoxes(pragma_clause, locals,            \
                                      i, j, k, fdir,                    \
                                      grgeom, patch_num,                \
                                      r, ix, iy, iz,                    \
                                      nx, ny, nz, body)                 \
  {                                                                     \
    int PV_ixl, PV_iyl, PV_izl, PV_ixu, PV_iyu, PV_izu;                 \
    int *PV_visiting = NULL;                                            \
    int PV_fdir[3];                                                     \
    fdir = PV_fdir;                                                     \
    for (int PV_f = 0; PV_f < GrGeomOctreeNumFaces; PV_f++)             \
    {                                                                   \
      switch (PV_f)                                                     \
      {                                                                 \
        case GrGeomOctreeFaceL:                                         \
          fdir[0] = -1; fdir[1] = 0; fdir[2] = 0;                       \
          break;                                                        \
        case GrGeomOctreeFaceR:                                         \
          fdir[0] = 1; fdir[1] = 0; fdir[2] = 0;                        \
          break;                                                        \
        case GrGeomOctreeFaceD:                                         \
          fdir[0] = 0; fdir[1] = -1; fdir[2] = 0;                       \
          break;                                                        \
        case GrGeomOctreeFaceU:                                         \
          fdir[0] = 0; fdir[1] = 1; fdir[2] = 0;                        \
          break;                                                        \
        case GrGeomOctreeFaceB:                                         \
          fdir[0] = 0; fdir[1] = 0; fdir[2] = -1;                       \
          break;                                                        \
        case GrGeomOctreeFaceF:                                         \
          fdir[0] = 0; fdir[1] = 0; fdir[2] = 1;                        \
          break;                                                        \
        default:                                                        \
          fdir[0] = -9999; fdir[1] = -9999; fdir[2] = -99999;           \
          break;                                                        \
      }                                                                 \
                                                                        \
      BoxArray* boxes = GrGeomSolidPatchBoxes(grgeom, patch_num, PV_f); \
      for (int PV_box = 0; PV_box < BoxArraySize(boxes); PV_box++)      \
      {                                                                 \
        Box box = BoxArrayGetBox(boxes, PV_box);                        \
        /* find octree and region intersection */                       \
        PV_ixl = pfmax(ix, box.lo[0]);                                  \
        PV_iyl = pfmax(iy, box.lo[1]);                                  \
        PV_izl = pfmax(iz, box.lo[2]);                                  \
        PV_ixu = pfmin((ix + nx - 1), box.up[0]);                       \
        PV_iyu = pfmin((iy + ny - 1), box.up[1]);                       \
        PV_izu = pfmin((iz + nz - 1), box.up[2]);                       \
                                                                        \
        int PV_diff_x = PV_ixu - PV_ixl;                                \
        int PV_diff_y = PV_iyu - PV_iyl;                                \
        int PV_diff_z = PV_izu - PV_izl;                                \
        int x_scale = !!PV_diff_x;                                      \
        int y_scale = !!PV_diff_y;                                      \
        int z_scale = !!PV_diff_z;                                      \
        if (PV_diff_x * PV_diff_y * PV_diff_z != 0) {                   \
          fprintf(stderr, "ERROR: Diff not 0 at %s %d\n", __FILE__, __LINE__); \
          exit(-1);                                                     \
        }                                                               \
        PRAGMA(omp pragma_clause collapse(3) private(i, j, k, ival locals)) \
          for (k = PV_izl; k <= PV_izu; k++)                            \
          {                                                             \
            for (j = PV_iyl; j <= PV_iyu; j++)                          \
            {                                                           \
              for (i = PV_ixl; i <= PV_ixu; i++)                        \
              {                                                         \
                int PV_tmp_i = i - PV_ixl;                              \
                int PV_tmp_j = j - PV_iyl;                              \
                int PV_tmp_k = k - PV_izl;                              \
                if (!z_scale) {                                         \
                  ival = CALC_IVAL(PV_diff_x, PV_tmp_j, PV_tmp_i);      \
                } else if (!y_scale) {                                  \
                  ival = CALC_IVAL(PV_diff_x, PV_tmp_k, PV_tmp_i);      \
                } else {                                                \
                  ival = CALC_IVAL(PV_diff_y, PV_tmp_k, PV_tmp_j);      \
                }                                                       \
                body;                                                   \
              }                                                         \
            }                                                           \
          }                                                             \
      }                                                                 \
    }                                                                   \
  }

/*------------------------------------------------------------------------
 * BCPatchLoop Redefinitions
 *------------------------------------------------------------------------*/
#define NewParallel_GrGeomPatchLoopBoxes(locals,                        \
                                         i, j, k, fdir,                 \
                                         grgeom, patch_num,             \
                                         r, ix, iy, iz,                 \
                                         nx, ny, nz, body)              \
  PRAGMA(omp parallel private(fdir))                                    \
  {                                                                     \
    MainBody_GrGeomPatchLoopBoxes(for, locals,                          \
                                  i, j, k, fdir,                        \
                                  grgeom, patch_num,                    \
                                  r, ix, iy, iz,                        \
                                  nx, ny, nz, body);                    \
  }

#define InParallel_GrGeomPatchLoopBoxes(...) \
  MainBody_GrGeomPatchLoopBoxes(for, __VA_ARGS__)

#define NoWait_GrGeomPatchLoopBoxes(...)  \
  MainBody_GrGeomPatchLoopBoxes(for nowait, __VA_ARGS__)


#endif // _PF_OMP_BCLOOPS_H
