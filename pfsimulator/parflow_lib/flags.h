#ifndef _FLAGS_H
#define _FLAGS_H

typedef struct {
  unsigned char type_of; /* Internal, BC, etc */
  unsigned char ipatch_types[6]; /* BC Type for each ipatch */
  unsigned char ipatch_faces[6]; /* Faces for each ipatch */
  unsigned short ipatch_ivals[6][6]; /* Ivals for each ipatch */
} Flag;

#define XYZStride(i, j, k, x, y, z) ((i) + (x) * ((j) + (y) * (k)))
#define FlagStride(flags, i, j, k, x, y, z) flags[ XYZStride(i, j, k, x, y, z) ]
#define GetFlag(flags, i, j, k, x, y, z) (Flag)FlagStride(flags, i, j, k, x, y, z);

#define Intr 0x01

#define DirBC 0x02
#define FluBC 0x04
#define OvlBC 0x08
#define SeeBC 0x10
#define OvlKBC 0x20
#define OvlDBC 0x40
#define AnyBC 0xFE

#define LeftByte 0x01
#define RightByte 0x02
#define DownByte 0x04
#define UpByte 0x08
#define BackByte 0x10
#define FrontByte 0x20
#define AnyFace 0xFF

#define SetTypeof(flag, byte) ((flag.type_of) |= byte)
#define SetInterior(flag) SetTypeof(flag, Intr)
#define SetBCType(flag, bc) SetTypeof(flag, bc)
#define SetPatchType(flag, ipatch, bc) ((flag.ipatch_types[ipatch]) = bc)
#define SetPatchFace(flag, ipatch, face) ((flag.ipatch_faces[ipatch]) |= (0x01 << face))
#define SetPatchIval(flag, ipatch, face, ival) ((flag.ipatch_ivals[ipatch][face]) = (unsigned short)ival)

#define GetFlagPatchType(flag, ipatch) (flag.ipatch_types[ipatch])
#define FlagPatchHasFace(flag, ipatch, face) (flag.ipatch_faces[ipatch] & (0x01 << face))

#define FlagInterior(flag) ((flag.type_of) & Intr)
#define FlagAnyBC(flag) ((flag.type_of) & AnyBC)
#define FlagIsBC(flag, bc) ((flag.type_of) & (bc))
#define FlagHasPatch(flag, ipatch) ((flag.ipatch_faces[ipatch]) & AnyFace)
#define FlagHasFace(flag, ipatch, face) ((flag.ipatch_faces[ipatch]) & (0x01 << face))
#define FlagGetIval(flag, ipatch, face) ((flag.ipatch_ivals[ipatch][face]))

#define SetInteriorFlags(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz, flags) \
  {                                                                     \
    GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,         \
    {                                                                   \
      Flag flag = GetFlag(flags, i, j, k, nx, ny, nz);                  \
      SetInterior(flag);                                                \
      FlagStride(flags, i, j, k, nx, ny, nz) = flag;                    \
    });                                                                 \
  }

#define SetBCFlags(ipatch, bc_struct, body)                              \
  {                                                                     \
    for (ipatch = 0; ipatch < BCStructNumPatches(bc_struct); ipatch++)  \
    {                                                                   \
      switch(BCStructBCType(bc_struct, ipatch))                         \
      {                                                                 \
        body;                                                           \
      }                                                                 \
    }                                                                   \
  }

#define SetBCFlag(bctype, type, i, j, k, nx, ny, nz, fdir, ival, bc_struct, ipatch, is, flags) \
  case bctype:                                                          \
  {                                                                     \
    BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,       \
    {                                                                   \
      Flag flag = GetFlag(flags, i, j, k, nx, ny, nz);                  \
      SetBCType(flag, type);                                            \
      SetPatchType(flag, ipatch, type);                                 \
      SetPatchFace(flag, ipatch, PV_f);                                 \
      SetPatchIval(flag, ipatch, PV_f, ival);                           \
      FlagStride(flags, i, j, k, nx, ny, nz) = flag;                    \
    });                                                                 \
    break;                                                              \
  }

#define PROLOGUE(x) { x }
#define EPILOGUE(x) { x }
#define NO_PROLOGUE {}
#define NO_EPILOGUE {}

#define Left 0
#define Right 1
#define Down 2
#define Up 3
#define Back 4
#define Front 5

#define FACE(_case, body)                         \
  case _case:                                     \
  {                                               \
    body;                                         \
    break;                                        \
  }

#define ApplyBCFlag(flag, type, ipatch, prologue, epilogue, ...)  \
  {                                                               \
    if (FlagIsBC(flag, type) && FlagHasPatch(flag, ipatch)) {     \
      for (int idx = 0; idx < 6; idx++) {                         \
        if (FlagHasFace(flag, ipatch, idx)) {                     \
          prologue;                                               \
          switch(idx) {                                           \
            __VA_ARGS__;                                          \
          }                                                       \
          epilogue;                                               \
        }                                                         \
      }                                                           \
    }                                                             \
  }

#endif // _FLAGS_H
