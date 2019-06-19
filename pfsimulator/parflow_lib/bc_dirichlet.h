#ifndef _BC_DIRICHLET_H
#define _BC_DIRICHLET_H

#define Do_AddDirichletBC_Pressure(i, j, k, is, ival, ipatch,           \
                                   bc_struct, bc_patch_values, body)    \
  ForBCStructNumPatches(ipatch, bc_struct)                              \
  {                                                                     \
    bc_patch_values = BCStructPatchValues(bc_struct, ipatch, is);       \
    {                                                                   \
      body;                                                             \
    }                                                                   \
  }                                                                     \



#endif // _BC_DIRICHLET_H
