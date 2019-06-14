/* LS_fft_mex.c
 *
 * compile with "mex -largeArrayDims -lmwlapack"
 *
 * As the FFT trick grid may be large, here we enforce using
 *   size_t:    for (large-sized) grid indexing
 *   ptrdiff_t: to address negative offsets from grid center
 * instead of mwSize and mwSignedIndex
 *
 * some acronyms used:
 * _ctor: constructor
 * _dtor: destructor
 * _B:    a Block of a block matrix
 * _P:    Pointer
 * _U:    Upper matrix
 * _L:    Lower matrix
 * _t:    Temporary (when appended to a variable)
 * _i:    Iterator
 * _R:    Real
 * _I:    Imag
 */
#if !defined(_WIN32)
#define cposv cposv_
/* #define chesv chesv_ */
#endif

#include <mex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lapack.h>
#include <stddef.h> /* common macros: size_t, ptrdiff_t, etc. */

/* def of GRID
 * GIM GCO GNO are defined w/ ptrdiff_t as they are used in signed index
 * calculation in LS_init()
 */
struct GRID
{
  float      Tik;  /* (1,), Tik regularizor coeff */
  size_t     nC;   /* (1,), number of coils */
  size_t     nDim; /* (1,), number of dimensions (1|2|3) */
  size_t     nVtx; /* (1,), number of lin-interp ngbs, (2|4|8) */
  mwSize    *nd;   /* (nDim,), grid size along each dim */
  mxArray   *G;    /* (nC*(nC+1)/2,), Grid cells */
  float     *GIU;  /* (nDim,), _Inv _Unit, normalizing shifts for grids */
  ptrdiff_t *GIM;  /* (nDim,), _Index _Multiplier */
  ptrdiff_t  GCO;  /* (1,),    _Center _Offset */
  ptrdiff_t *GNO;  /* (nVtx,), _Neighbors _Offsets */
};

struct GRID *Grid_ctor(
  const size_t   nDim,
  const size_t   nC,
  const float    Tik,
  const float   *GIU,
  const mxArray *G)
{
  struct GRID *grid = mxCalloc(1, sizeof(struct GRID));

  size_t nVtx_t = 2, i_t;
  for (i_t = 1; i_t < nDim; i_t++) {nVtx_t *= 2;}
  const size_t nVtx = nVtx_t;

  grid->Tik  = Tik;
  grid->nC   = nC;
  grid->nDim = nDim;
  grid->nVtx = nVtx;

  grid->nd = mxCalloc(nDim, sizeof(size_t));
  memcpy(grid->nd, mxGetDimensions(mxGetCell(G,0)), nDim*sizeof(mwSize));
  grid->G  = mxDuplicateArray(G);

  grid->GIU = mxCalloc(nDim, sizeof(float));
  memcpy(grid->GIU, GIU, nDim*sizeof(float));
  grid->GIM = mxCalloc(nDim, sizeof(ptrdiff_t));
  grid->GCO = 0;
  grid->GNO = mxCalloc(nVtx, sizeof(ptrdiff_t));

  const mwSize *nd = grid->nd;
  /* setup GCO GNO GIM */
  ptrdiff_t *GCO, *GNO, *GIM;
  GCO = &(grid->GCO);
  GNO =   grid->GNO + nVtx-1;
  GIM =   grid->GIM + nDim;

  switch(nDim)
  {
    case 3:
      /* ->GCO = ->GCO + floor(pz/2.0)*px*py */
      *GCO  += ( (*(nd+2))/2 ) * (*(nd+1)) * (*nd);
      *--GIM = *(nd+1) * (*nd);
      *GNO-- = *GIM + *nd + 1;
      *GNO-- = *GIM + *nd;
      *GNO-- = *GIM + 1;
      *GNO-- = *GIM; /* no break */
    case 2:
      *GCO  += ( (*(nd+1))/2 ) * (*nd);
      *--GIM = *nd;
      *GNO-- = *nd + 1;
      *GNO-- = *nd; /* no break */
    case 1:
      *GCO  += ( (*nd)/2+!!(*nd) ) - 1; /* -1 convert to 0-based indexing */
      *--GIM = 1;
      *GNO-- = 1;
      *GNO   = 0; /* break only in case 1 */
      break;
    default:
      mexErrMsgTxt("Only support nDim == 1 or 2 or 3");
  } /* [0, 1(, nx, nx+1(, ny*nx, ny*nx+1, ny*nx+nx, ny*nx+nx+1))] */

  return grid;
}

void Grid_dtor(struct GRID *grid)
{
  mxFree(grid->nd);

  mxDestroyArray(grid->G);

  mxFree(grid->GIU);
  mxFree(grid->GIM);
  mxFree(grid->GNO);

  mxFree(grid);
}

/* def of LS
 * NhN_GI NhC_GI are defined w/ size_t, as grid may be large.
 * Their calculation involves signed addition
 */
struct LS
{
  size_t  nNgb;

  float  *NhN;       /* (2*nC*nNgb, nC*nNgb), 2 for complex */
  float  *NhC;       /* (2*nC*nNgb, nC), 2 for complex */
  /* _Block _Ptr, pointers indexing the blocks of NhN & NhC */
  float **NhN_BPL; /* L->(nC*(nC+1)/2,): L-triangle w/ Diag */
  float **NhN_BPU; /* U->(nC*(nC-1)/2,): U-triangle */
  float **NhC_BPL;
  float **NhC_BPU;
  /* _Interpolation _Coefficients (_Reversed ptr) */
  float  *IC;        /* (nVtx,) container for interpolation coeffs */
  float  *GD;        /* (nDim,) container for floored _Grid _Distance */
};

struct LS *LS_ctor(const struct GRID *grid)
{
  struct LS *ls = mxCalloc(1, sizeof(struct LS));
  ls->nNgb = 0;
  ls->NhN = NULL;
  ls->NhC = NULL;

  const size_t nC = grid->nC;
  ls->NhN_BPL = mxCalloc((nC*(nC+1))/2, sizeof(float*));
  ls->NhN_BPU = mxCalloc((nC*(nC-1))/2, sizeof(float*));
  ls->NhC_BPL = mxCalloc((nC*(nC+1))/2, sizeof(float*));
  ls->NhC_BPU = mxCalloc((nC*(nC-1))/2, sizeof(float*));

  ls->IC = mxCalloc(grid->nVtx, sizeof(float));
  ls->GD = mxCalloc(grid->nDim, sizeof(float));

  return ls;
}

/* prepare block indices */
void LS_BP_prep(
  const size_t  nNgb,
  const size_t  nC,
  struct LS    *ls)
{
  /* setup ls->NhN_BP & ls->NhC_BP, 2 for Fortran */
  const size_t nNgb2     = 2*nNgb;
  const size_t nNgbSq2   = nNgb*nNgb2;
  const size_t nCnNgb2   = nC*nNgb2;
  const size_t nCnNgbSq2 = nC*nNgbSq2;

  size_t  CoilRef_i, CoilQry_i;
  float **NhN_BPL_t = ls->NhN_BPL;
  float **NhN_BPU_t = ls->NhN_BPU;
  float **NhC_BPL_t = ls->NhC_BPL;
  float **NhC_BPU_t = ls->NhC_BPU;
  float  *NhN_Lt, *NhC_Lt, *NhN_Ut, *NhC_Ut; /* _Lower & _Upper triangle _tmp */

  for (CoilRef_i = 0; CoilRef_i < nC; CoilRef_i++)
  {
    /* _Upper triangle _temp current value is not used, just to update offset */
    *NhN_BPL_t++ = (NhN_Ut = NhN_Lt = ls->NhN + CoilRef_i*(nCnNgbSq2+nNgb2));
    *NhC_BPL_t++ = (NhC_Ut = NhC_Lt = ls->NhC + CoilRef_i*(nCnNgb2+nNgb2));

    for (CoilQry_i = CoilRef_i+1; CoilQry_i < nC; CoilQry_i++)
    {
      *NhN_BPL_t++ = (NhN_Lt += nNgb2);
      *NhC_BPL_t++ = (NhC_Lt += nNgb2);

      *NhN_BPU_t++ = (NhN_Ut += nCnNgbSq2);
      *NhC_BPU_t++ = (NhC_Ut += nCnNgb2);
    }
  }
}

void LS_init(
  const size_t       nNgb,
  const struct GRID *grid,
        struct LS   *ls)
{
  const size_t nC = grid->nC;

  ls->nNgb = nNgb;
  ls->NhN  = mxCalloc(2*(nC*nNgb)*(nC*nNgb), sizeof(float)); /* 2 for cmplx */
  ls->NhC  = mxCalloc(2*(nC*nNgb)*nC,        sizeof(float));

  LS_BP_prep(nNgb, nC, ls);
}

/* This only frees some struct members */
/* Double free for this struct is allowed, as it has member being persistent */
void LS_free(struct LS *ls)
{
  if(NULL != ls->NhN)
  {
    mxFree(ls->NhN);
    ls->NhN = NULL;
  }
  if(NULL != ls->NhC)
  {
    mxFree(ls->NhC);
    ls->NhC = NULL;
  }
}

/* Double destruct is not allowed tho */
void LS_dtor(struct LS *ls)
{
  /* LS_free(ls); */ /* non-persistent var are taken cared by matlab */

  mxFree(ls->NhN_BPL);
  mxFree(ls->NhN_BPU);
  mxFree(ls->NhC_BPL);
  mxFree(ls->NhC_BPU);

  mxFree(ls->IC);
  mxFree(ls->GD);

  mxFree(ls);
}

static struct GRID *grid;
static struct LS   *ls;
static int          constructed = 0;

/* make-persistent & free persistent*/
void makePersistent(void)
{
  mexMakeMemoryPersistent(grid);

  mexMakeMemoryPersistent(grid->nd);

  mexMakeArrayPersistent(grid->G);

  mexMakeMemoryPersistent(grid->GIU);
  mexMakeMemoryPersistent(grid->GIM);
  mexMakeMemoryPersistent(grid->GNO);

  mexMakeMemoryPersistent(ls);

  mexMakeMemoryPersistent(ls->NhN_BPL);
  mexMakeMemoryPersistent(ls->NhN_BPU);
  mexMakeMemoryPersistent(ls->NhC_BPL);
  mexMakeMemoryPersistent(ls->NhC_BPU);

  mexMakeMemoryPersistent(ls->IC);
  mexMakeMemoryPersistent(ls->GD);
}

void freePersistent(void)
{
  Grid_dtor(grid);
  LS_dtor(ls);
}

/*
 * core part
 */

void IC_calc(
  const size_t  nVtx, /* (1,), # of vertices used for interpolation 2/4/8 */
        float  *GD,   /* (nDim,), shortest Grid Distance of each _dim */
        float  *IC)   /* (nVtx,), storing the coefficients*/
{
  switch(nVtx)
  {
    case 2:
      *IC++ =                         (1-*GD);
      *IC   =                         (  *GD);
      break;
    case 4:
      *IC++ =             (1-*(GD+1))*(1-*GD);
      *IC++ =             (1-*(GD+1))*(  *GD);
      *IC++ =             (  *(GD+1))*(1-*GD);
      *IC   =             (  *(GD+1))*(  *GD);
      break;
    case 8:
      *IC++ = (1-*(GD+2))*(1-*(GD+1))*(1-*GD);
      *IC++ = (1-*(GD+2))*(1-*(GD+1))*(  *GD);
      *IC++ = (1-*(GD+2))*(  *(GD+1))*(1-*GD);
      *IC++ = (1-*(GD+2))*(  *(GD+1))*(  *GD);
      *IC++ = (  *(GD+2))*(1-*(GD+1))*(1-*GD);
      *IC++ = (  *(GD+2))*(1-*(GD+1))*(  *GD);
      *IC++ = (  *(GD+2))*(  *(GD+1))*(1-*GD);
      *IC   = (  *(GD+2))*(  *(GD+1))*(  *GD);
      break;
    default:
      mexErrMsgTxt("Only support nVtx == 2 or 4 or 8");
  }
}

void NhN_BP_NextRef( /* update NhN block matrices ptr's offsets */
  const mwSize  nC,
  const size_t  nNgb,
        float **NhN_BPL,
        float **NhN_BPU)
{
  /* 2 for Fortran */
  const size_t inc_L = 2*((nC-1)*nNgb);
  const size_t dec_U = 2*(nC*nNgb*nNgb - 1);
  size_t i;
  const size_t nNhN_Lt = (nC*(nC+1))/2;
  const size_t nNhN_Ut = (nC*(nC-1))/2;

  for (i = 0; i < nNhN_Lt; i++)
  {
    *NhN_BPL += inc_L;
    NhN_BPL++;
  }
  for (i = 0; i < nNhN_Ut; i++ )
  {
    *NhN_BPU -= dec_U; /* go backwards */
    NhN_BPU++;
  }
}

float doInterp(
  const float     *Gdata,
  const size_t     GI0,
  const ptrdiff_t *GNO,
  const float     *IC,
  const size_t     nVtx,
  const int        isRev)
{
  float res = 0;
  size_t i;
  if (isRev) { for(i=0;i<nVtx;i++) { res += *(IC++) * Gdata[GI0 - *(GNO++)]; } }
  else       { for(i=0;i<nVtx;i++) { res += *(IC++) * Gdata[GI0 + *(GNO++)]; } }
  return res;
}

void update_GDGI(
  const size_t     nDim,
  const float     *GIU,
  const ptrdiff_t *GIM,
  const float     *shQry,
  const float     *shRef,
        float     *GD,
        ptrdiff_t *GI)
{
  float  diff;
  size_t i = 0;
  *GI = 0;
  if (NULL == shRef)
  {
    for(i = 0; i < nDim; i++)
    {
      diff  = (GIU[i]) * (         - shQry[i]);
      GD[i] = diff - floorf(diff);
      *GI  += (GIM[i]) * (ptrdiff_t)floorf(diff);
    }
  }
  else
  {
    for(i = 0; i < nDim; i++)
    {
      diff  = (GIU[i]) * (shRef[i] - shQry[i]);
      GD[i] = diff - floorf(diff);
      *GI  += (GIM[i]) * (ptrdiff_t)floorf(diff);
    }
  }
}

void LS_NhN_form(
  const size_t       nNgb,
  const float       *sh,
  const struct GRID *grid,
        struct LS   *ls)
{
  const float     *shRef   = sh;
  const float     *shQry_t = NULL; /* needs reset for each Ref */
  /* Para from grid */
  const float     *GIU = grid->GIU;
  const ptrdiff_t *GIM = grid->GIM;
  const ptrdiff_t *GNO = grid->GNO;
  const mxArray   *G   = grid->G;
  const float     *G_R, *G_I;
  const float      Tik = grid->Tik;
  const size_t nC = grid->nC, nDim = grid->nDim, nVtx = grid->nVtx;

  /* interp containers: coefficients, grid_dist, grid_ind */
  float *IC = ls->IC, *GD = ls->GD;
  ptrdiff_t GI = 0;  /* _Grid _Ind to the 1st neighbor for interp */
  size_t    GI0 = 0; /* _Grid _Ind, starting at 0 */

  /* BPL & BPU iteraters */
  size_t   inc = 2*nC*nNgb - 1; /* 2 for cplx, -1 shift from imag to real */
  size_t   NhN_BPL_i, NhN_BPU_i;
  float  **NhN_BPL_t = ls->NhN_BPL, **NhN_BPU_t = ls->NhN_BPU;

  /* for-loop vars */
  size_t   NgbRef_i, NgbQry_i, CoilRef_i, CoilQry_i;
  float    R, I;

  /* mexPrintf("GIO: %d\n", grid->GCO); */
  /* mexPrintf("nNgb: %d\n", nNgb); */
  for(NgbRef_i = 0; NgbRef_i < nNgb; NgbRef_i++)
  {
    shQry_t = sh;     /* reset to the 1st shift  */
    for(NgbQry_i = 0; NgbQry_i < nNgb; NgbQry_i++)
    {
      update_GDGI(nDim, GIU, GIM, shQry_t, shRef, GD, &GI);
      GI0 = grid->GCO + GI; /* center offset */

      /* update coefficients & do interpolation */
      IC_calc(nVtx, GD, IC);
      NhN_BPL_i = NhN_BPU_i = 0;

      for(CoilRef_i = 0; CoilRef_i < nC; CoilRef_i++)
      {
        G_R = (const float*)    mxGetData(mxGetCell(G, NhN_BPL_i));
        G_I = (const float*)mxGetImagData(mxGetCell(G, NhN_BPL_i)) ;

        *(NhN_BPL_t[NhN_BPL_i]++) = (R=doInterp(G_R, GI0, GNO, IC, nVtx, 0));
        *(NhN_BPL_t[NhN_BPL_i]++) = (I=doInterp(G_I, GI0, GNO, IC, nVtx, 0));

        NhN_BPL_i++;

        for(CoilQry_i = CoilRef_i+1; CoilQry_i < nC; CoilQry_i++)
        {
          G_R = (const float*)    mxGetData(mxGetCell(G, NhN_BPL_i));
          G_I = (const float*)mxGetImagData(mxGetCell(G, NhN_BPL_i));

          *(NhN_BPL_t[NhN_BPL_i]++) = (R=doInterp(G_R, GI0, GNO, IC, nVtx, 0));
          *(NhN_BPL_t[NhN_BPL_i]++) = (I=doInterp(G_I, GI0, GNO, IC, nVtx, 0));

          *(NhN_BPU_t[NhN_BPU_i]++) = R;
          *(NhN_BPU_t[NhN_BPU_i]) = -I; /* pointing in-row next element */
          NhN_BPU_t[NhN_BPU_i] += inc;

          NhN_BPL_i++;
          NhN_BPU_i++;
        }
      }
      shQry_t += nDim; /* next query shift */
    }
    /* update BP to next ref column in NhN */
    NhN_BP_NextRef(nC, nNgb, NhN_BPL_t, NhN_BPU_t);

    shRef += nDim;     /* next reference shift */
  }
  /* apply the Tik reg */
  size_t Tik_i = 0, Tik_inc = 2*(nC*nNgb+1), nTik = nC*nNgb;
  float  *NhN_t = ls->NhN, nNgbTik = nNgb*Tik;
  if (Tik)
  {
    for(Tik_i = 0; Tik_i < nTik; Tik_i++)
    {
      *NhN_t += nNgbTik;
      NhN_t  += Tik_inc;
    }
  }
}

void LS_NhC_form(
  const size_t       nNgb,
  const float       *sh,
  const struct GRID *grid,
        struct LS   *ls)
{
  const float     *shQry_t = sh;
  /* Paras from grid */
  const float     *GIU = grid->GIU;
  const ptrdiff_t *GIM = grid->GIM;
  const ptrdiff_t *GNO = grid->GNO;
  const mxArray   *G = grid->G;
  const float     *G_R, *G_I;
  const size_t nC = grid->nC, nDim = grid->nDim, nVtx = grid->nVtx;

  /* interp containers: coefficients, grid_dist, grid_ind */
  float *IC = ls->IC, *GD = ls->GD;
  ptrdiff_t GI = 0;  /* _Grid _Ind to the 1st neighbor for interpolation */
  size_t    GI0_L = 0, GI0_U = 0; /* _Grid _Ind, starting at 0 */

  /* BPL & BPU iteraters */
  size_t   NhC_BPL_i, NhC_BPU_i;
  float  **NhC_BPL_t = ls->NhC_BPL, **NhC_BPU_t = ls->NhC_BPU;

  /* for-loop vars */
  size_t  CtrQry_i = 0, CoilRef_i = 0, CoilQry_i = 0;

  for(CtrQry_i = 0; CtrQry_i < nNgb; CtrQry_i++)
  {
    update_GDGI(nDim, GIU, GIM, shQry_t, NULL, GD, &GI);

    GI0_L = grid->GCO + GI; /* nd-array center offset */
    GI0_U = grid->GCO - GI;
    /* update coefficients & do interpolation */
    IC_calc(nVtx, GD, IC);
    NhC_BPL_i = NhC_BPU_i = 0;
    for(CoilRef_i = 0; CoilRef_i < nC; CoilRef_i++)
    {
      G_R = (const float*)    mxGetData(mxGetCell(G, NhC_BPL_i));
      G_I = (const float*)mxGetImagData(mxGetCell(G, NhC_BPL_i));

      *(NhC_BPL_t[NhC_BPL_i]++) = doInterp(G_R, GI0_L, GNO, IC, nVtx, 0);
      *(NhC_BPL_t[NhC_BPL_i]++) = doInterp(G_I, GI0_L, GNO, IC, nVtx, 0);

      NhC_BPL_i++;

      for(CoilQry_i = CoilRef_i+1; CoilQry_i < nC; CoilQry_i++)
      {
        G_R = (const float*)    mxGetData(mxGetCell(G, NhC_BPL_i));
        G_I = (const float*)mxGetImagData(mxGetCell(G, NhC_BPL_i));

        *(NhC_BPL_t[NhC_BPL_i]++) =  doInterp(G_R, GI0_L, GNO, IC, nVtx, 0);
        *(NhC_BPL_t[NhC_BPL_i]++) =  doInterp(G_I, GI0_L, GNO, IC, nVtx, 0);

        *(NhC_BPU_t[NhC_BPU_i]++) =  doInterp(G_R, GI0_U, GNO, IC, nVtx, 1);
        *(NhC_BPU_t[NhC_BPU_i]++) = -doInterp(G_I, GI0_U, GNO, IC, nVtx, 1);

        NhC_BPL_i++;
        NhC_BPU_i++;
      }
    }
    shQry_t += nDim; /* next query shift */
  }
}

/*
 * Convert Fortran complex storage to MATLAB real and imaginary parts.
 * X = fort2mat(Z,ldz,m,n) copies Z to X, producing a complex mxArray
 * with mxGetM(X) = m and mxGetN(X) = n.
 */
mxArray* fort2mat(
  const size_t    m,
  const size_t    n,
  const ptrdiff_t ldz,
        float    *Z)
{
  size_t i,j;
  ptrdiff_t incz;
  float *xr,*xi,*zp;
  mxArray *X;

  X = mxCreateNumericMatrix(m,n, mxSINGLE_CLASS, mxCOMPLEX);
  xr = (float *)mxGetData(X);
  xi = (float *)mxGetImagData(X);
  zp = Z;
  incz = 2*(ldz-m);

  for (j = 0; j < n; j++)
  {
    for (i = 0; i < m; i++)
    {
      *xr++ = *zp++;
      *xi++ = *zp++;
    }
    zp += incz;
  }

  return(X);
}

/*
 * Function [coeff] = LS_fft_mex(sh, nC, Tik, F_c, dgrid, iofst, m2, coilSub)
 * INPUTS:
 *  - sh (nDim, nNgb) single, shifts, NOTICE the dimensions
 *  - ->Tik (1,), Tiknov regularizer
 *  - ->nC  (1,), number of coils
 *  - F_c (->nC*(->nC+1)/2,) cell, spectrum of conj correlated ACS
 *  - idgrid (nDim,) single, inv grid size for normalizing shifts to calc coefs
 *
 */
void mexFunction(
  int nlhs,       mxArray *plhs[],
  int nrhs, const mxArray *prhs[])
{
  if (!mxIsSingle(prhs[0])) {mexErrMsgTxt("Require single shifts input");}
  const size_t nDim = mxGetM(prhs[0]);
  const size_t nNgb = mxGetN(prhs[0]);

  /* handling INPUTS */
  if( !mxIsSingle(prhs[0]) ) { mexErrMsgTxt( "sh: single precision" ); }
  const float *sh = (float*)mxGetData(prhs[0]);

  size_t nC;
  float *GIU;
  float  Tik = (float) mxGetScalar(prhs[2]);
  if (!constructed) /* dealing with the remaining prhs */
  {
    nC  = (size_t)mxGetScalar(prhs[1]);

    if( !mxIsSingle(prhs[4]) ) { mexErrMsgTxt( "dgrid: single precision" ); }
    GIU = (float *)mxGetData(prhs[4]);

    /* Initialize the two essential structs */
    grid = Grid_ctor(nDim, nC, Tik, GIU, prhs[3]); /* persistent static var */
    ls   = LS_ctor(grid);

    makePersistent();
    mexAtExit(freePersistent); /* remember to call clear mex after all is done */

    constructed = 1;
  }
  else { grid->Tik = Tik; }

  LS_init(nNgb, grid, ls); /* update ls by the input shifts */

  LS_NhN_form(nNgb, sh, grid, ls);
  LS_NhC_form(nNgb, sh, grid, ls);

  float *NhN = ls->NhN;
  float *NhC = ls->NhC;

  size_t N = grid->nC*nNgb, NRHS = grid->nC;

  if (nlhs > 1) {plhs[1] = fort2mat(N, N, N, NhN);}
  if (nlhs > 2) {plhs[2] = fort2mat(N, NRHS, N, NhC);}

  /* Call LAPACK function */
  char uplo = 'U';
  mwSignedIndex info;
  cposv(&uplo, &N, &NRHS, NhN, &N, NhC, &N, &info);
  if (info != 0) {mexPrintf("LS_fft_mex: ill-conditioned, try a larger Tik\n");}

  /* LAPACK always do in-place operation, so NhC is now coeffs */
  plhs[0] = fort2mat(N, NRHS, N, NhC);

  /* clean up */
  LS_free(ls);
}

