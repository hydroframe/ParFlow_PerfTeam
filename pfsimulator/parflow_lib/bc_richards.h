#ifndef _BC_RICHARDS_H
#define _BC_RICHARDS_H

#define Do_Richards_SymmContrib(i, j, k, ival, bc_struct, ipatch, is, body) \
  ForBCStructNumPatches(ipatch, bc_struct)                              \
  {                                                                     \
    switch(BCStructBCType(bc_struct, ipatch))                           \
    {                                                                   \
      body;                                                             \
    }                                                                   \
  }

#define Do_Richards_BuildJCMatrix(i, j, k, ival, bc_struct, ipatch, is, body) \
  ForBCStructNumPatches(ipatch, bc_struct)                              \
  {                                                                     \
    switch(BCStructBCType(bc_struct, ipatch))                           \
    {                                                                   \
      body;                                                             \
    }                                                                   \
  }

#define IfSurfaceNode(k1, k, body)              \
  if (((k) >= 0) && ((k1) == (k)))              \
  {                                             \
    body;                                       \
  }

#endif // _BC_RICHARDS_H
