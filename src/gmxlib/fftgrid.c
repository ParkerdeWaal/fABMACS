#include "fftgrid.h"
#include "smalloc.h"
#include "futil.h"

#define TOL -1

t_fftgrid *mk_fftgrid(int nx,int ny,int nz)
{
  t_fftgrid *grid;
  
  snew(grid,1);
  grid->nx   = nx;
  grid->ny   = ny;
  grid->nz   = nz;
  grid->nxyz = nx*ny*nz;
  
#ifdef USE_SGI_FFT
  if (nx != nz) {
    fprintf(stderr,"You can't use SGI optimized FFT routines when the number of grid points is not the same in X and Z directions. Sorry\n");
    exit(1);
  }
  grid->la1  = nx+2;
  grid->la2  = ny;
  grid->nptr = nz*grid->la2*grid->la1;
#else
  grid->la1  = ny;
  grid->la2  = nz;
  grid->nptr = nx*ny*nz;
#endif
  grid->la12 = grid->la1*grid->la2;

  snew(grid->ptr,grid->nptr);
  
  return grid;
}

void gmxfft3D(FILE *fp,bool bVerbose,t_fftgrid *grid,int dir)
{
  static bool bFirst=TRUE;
  
#ifdef USE_SGI_FFT
  static real *coeff;

  if (bVerbose) {
    if (dir == FFTW_FORWARD)
      fprintf(stderr,"Doing forward 3D-FFT (%d)\n",dir);
    else
      fprintf(stderr,"Doing backward 3D-FFT (%d)\n",dir);
  }
  
  if (bFirst) {
    fprintf(fp,"Going to use SGI optimized FFT routines.\n");
    
#ifdef DOUBLE
    coeff  = dzfft3dui(grid->nx,grid->ny,grid->nz,NULL);
#else
    coeff  = scfft3dui(grid->nx,grid->ny,grid->nz,NULL);
#endif
    bFirst = FALSE;
  }
  if (dir == FFTW_FORWARD) {
#ifdef DOUBLE
    dzfft3du(dir,grid->nx,grid->ny,grid->nz,
	     grid->ptr,grid->la1,grid->la2,coeff);
#else
    scfft3du(dir,grid->nx,grid->ny,grid->nz,
	     grid->ptr,grid->la1,grid->la2,coeff);
#endif
  }
  else if (dir == FFTW_BACKWARD) {
#ifdef DOUBLE
    zdfft3du(dir,grid->nx,grid->ny,grid->nz,
	     grid->ptr,grid->la1,grid->la2,coeff);
#else
    csfft3du(dir,grid->nx,grid->ny,grid->nz,
	      grid->ptr,grid->la1,grid->la2,coeff);
#endif
  }
#else
  static fftwnd_plan forward_plan,backward_plan;
  
  if (bFirst) {
    fprintf(fp,"Using the FFTW library (Fastest Fourier Transform in the West)\n");
    forward_plan  = fftw3d_create_plan(grid->nx,grid->ny,grid->nz,
				       FFTW_FORWARD,
				       FFTW_ESTIMATE | FFTW_IN_PLACE);
    backward_plan = fftw3d_create_plan(grid->nx,grid->ny,grid->nz,
				       FFTW_BACKWARD,
				       FFTW_ESTIMATE | FFTW_IN_PLACE);
    bFirst        = FALSE;
  }
  if (dir == FFTW_FORWARD)
    fftwnd(forward_plan, 1,(FFTW_COMPLEX *)grid->ptr,1,0,NULL,0,0);
  else if (dir == FFTW_BACKWARD)
    fftwnd(backward_plan,1,(FFTW_COMPLEX *)grid->ptr,1,0,NULL,0,0);
  else
    fatal_error(0,"Invalid direction for FFT: %d",dir);
#endif
}

void clear_fftgrid(t_fftgrid *grid)
{
  int      i,ngrid;
  t_fft_tp *ptr;
  
  ngrid = grid->nptr;
  ptr   = grid->ptr;
  
  for (i=0; (i<ngrid); i++) {
#ifdef USE_SGI_FFT
    ptr[i] = 0;
#else
    ptr[i].re = ptr[i].im = 0;
#endif
  }
}

void unpack_fftgrid(t_fftgrid *grid,int *nx,int *ny,int *nz,
		    int *la1,int *la2,int *la12,t_fft_tp **ptr)
{
  *nx  = grid->nx;
  *ny  = grid->ny;
  *nz  = grid->nz;
  *la1 = grid->la1;
  *la2 = grid->la2;
  *la12= grid->la12;
  *ptr = grid->ptr;
}

void print_fftgrid(FILE *out,char *title,t_fftgrid *grid,real factor,char *pdb,
		   rvec box,bool bReal)
{
  static char *pdbformat="%-6s%5d  %-4.4s%3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n";
  FILE     *fp;
  int      i,ix,iy,iz;
  real     fac=50.0,value;
  rvec     boxfac;
  int      nx,ny,nz,la1,la2,la12;
  t_fft_tp *ptr,g;
  
  if (pdb)
    fp = ffopen(pdb,"w");
  else
    fp = out;
  if (!fp)
    return;

  unpack_fftgrid(grid,&nx,&ny,&nz,&la1,&la2,&la12,&ptr);
    
  boxfac[XX] = fac*box[XX]/nx;
  boxfac[YY] = fac*box[YY]/ny;
  boxfac[ZZ] = fac*box[ZZ]/nz;
  
  if (pdb)
    fprintf(fp,"REMARK ");
  
  fprintf(fp,"Printing all non-zero %s elements of %s\n",
	  bReal ? "Real" : "Imaginary",title);
  for(i=ix=0; (ix<nx); ix++)
    for(iy=0; (iy<ny); iy++)
      for(iz=0; (iz<nz); iz++,i++) {
	g = ptr[INDEX(ix,iy,iz)];
	if (pdb) {
#ifdef USE_SGI_FFT
	  value = g;
#else
	  value = bReal ? g.re : g.im;
#endif
	  if (fabs(value) > TOL)
	    fprintf(fp,pdbformat,"ATOM",i,"H","H",' ',
		    i,ix*boxfac[XX],iy*boxfac[YY],iz*boxfac[ZZ],
		    1.0,factor*value);
	} 
	else {
#ifdef USE_SGI_FFT
	  if (fabs(g) > TOL)
	    fprintf(fp,"%s[%2d][%2d][%2d] = %12.5e\n",
		    title,ix,iy,iz,g*factor);
#else
	  if ((fabs(g.re) > TOL) || (fabs(g.im) > TOL))
	    fprintf(fp,"%s[%2d][%2d][%2d] = %12.5e + i %12.5e%s\n",
		    title,ix,iy,iz,g.re*factor,g.im*factor,
		    (g.im != 0) ? " XXX" : "");
#endif
	}
      }
  fflush(fp);
}
