/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
 * University of Groningen, The Netherlands
 *
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */
 
#include <string.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "fatal.h"
#include "txtdump.h"
#include "names.h"
#include "futil.h"
#include "trnio.h"
#include "gmxfio.h"

#define BUFSIZE		128
#define GROMACS_MAGIC   1993

static int nFloatSize(t_trnheader *sh)
{
  int nflsize;
  
  if (sh->box_size)
    nflsize = sh->box_size/(DIM*DIM);
  else if (sh->x_size)
    nflsize = sh->x_size/(sh->natoms*DIM);
  else if (sh->v_size)
    nflsize = sh->v_size/(sh->natoms*DIM);
  else if (sh->f_size)
    nflsize = sh->f_size/(sh->natoms*DIM);
  else 
    fatal_error(0,"Can not determine precision of trn file, quit!");
  
  if (((nflsize != sizeof(float)) && (nflsize != sizeof(double))))
    fatal_error(0,"Float size %d. Maybe different CPU?",nflsize);
      
  return nflsize;
}

static bool do_trnheader(int fp,bool bRead,t_trnheader *sh)
{
  static int magic=GROMACS_MAGIC;
  static char *version = "GMX_trn_file";
  static bool bFirst=TRUE;
  char buf[256];
  
  fio_select(fp);
  if (!do_int(magic))
    return FALSE;
  
  if (bRead) {
    do_string(buf);
    if (bFirst) {
      fprintf(stderr,"trn version: %s\n",buf);
      bFirst = FALSE;
    }
  }
  else
    do_string(version);
  do_int(sh->ir_size);
  do_int(sh->e_size);
  do_int(sh->box_size);
  do_int(sh->vir_size);
  do_int(sh->pres_size);
  do_int(sh->top_size); 
  do_int(sh->sym_size); 
  do_int(sh->x_size); 
  do_int(sh->v_size); 
  do_int(sh->f_size); 
  
  fio_setprecision(fp,(nFloatSize(sh) == sizeof(double)));
  
  do_int(sh->natoms); 
  do_int(sh->step); 
  do_int(sh->nre); 
  do_real(sh->t); 
  do_real(sh->lambda); 
  
  return TRUE;
}

void pr_trnheader(FILE *fp,int indent,char *title,t_trnheader *sh)
{
  if (sh) {
    indent=pr_title(fp,indent,title);
    (void) pr_indent(fp,indent);
    (void) fprintf(fp,"box_size    = %d\n",sh->box_size);
    (void) pr_indent(fp,indent);
    (void) fprintf(fp,"x_size      = %d\n",sh->x_size);
    (void) pr_indent(fp,indent);
    (void) fprintf(fp,"v_size      = %d\n",sh->v_size);
    (void) pr_indent(fp,indent);
    (void) fprintf(fp,"f_size      = %d\n",sh->f_size);
    
    (void) pr_indent(fp,indent);
    (void) fprintf(fp,"natoms      = %d\n",sh->natoms);
    (void) pr_indent(fp,indent);
    (void) fprintf(fp,"step        = %d\n",sh->step);
    (void) pr_indent(fp,indent);
    (void) fprintf(fp,"t           = %e\n",sh->t);
    (void) pr_indent(fp,indent);
    (void) fprintf(fp,"lambda      = %e\n",sh->lambda);
  }
}

static bool do_htrn(int fp,bool bRead,t_trnheader *sh,
		    rvec *box,rvec *x,rvec *v,rvec *f)
{
  matrix pv;
  
  if (sh->box_size != 0) ndo_rvec(box,DIM);
  if (sh->vir_size != 0) ndo_rvec(pv,DIM);
  if (sh->pres_size!= 0) ndo_rvec(pv,DIM);
  if (sh->x_size   != 0) ndo_rvec(x,sh->natoms);
  if (sh->v_size   != 0) ndo_rvec(v,sh->natoms);
  if (sh->f_size   != 0) ndo_rvec(f,sh->natoms);

  return TRUE;
}

static bool do_trn(int fp,bool bRead,int *step,real *t,real *lambda,
		   rvec *box,int *natoms,rvec *x,rvec *v,rvec *f)
{
  t_trnheader *sh;
  bool bResult;
  
  snew(sh,1);
  if (!bRead) {
    sh->box_size=(box)?sizeof(matrix):0;
    sh->x_size=((x)?(*natoms*sizeof(x[0])):0);
    sh->v_size=((v)?(*natoms*sizeof(v[0])):0);
    sh->f_size=((f)?(*natoms*sizeof(f[0])):0);
    sh->natoms = *natoms;
    sh->step   = *step;
    sh->nre    = 0;
    sh->t      = *t;
    sh->lambda = *lambda;
  }
  if (!do_trnheader(fp,bRead,sh))
    return FALSE;
  if (bRead) {
    *natoms = sh->natoms;
    *step   = sh->step;
    *t      = sh->t;
    *lambda = sh->lambda;
    if (sh->ir_size)
      fatal_error(0,"inputrec in trn file");
    if (sh->e_size)
      fatal_error(0,"energies in trn file");
    if (sh->top_size)
      fatal_error(0,"topology in trn file");
    if (sh->sym_size)
      fatal_error(0,"symbol table in trn file");
  }
  bResult = do_htrn(fp,bRead,sh,box,x,v,f);

  sfree(sh);
  
  return bResult;
}

/************************************************************
 *
 *  The following routines are the exported ones
 *
 ************************************************************/
 
void read_trnheader(char *fn,t_trnheader *trn)
{
  int  fp;
  
  fp = open_trn(fn,"r");
  if (!do_trnheader(fp,TRUE,trn))
    fatal_error(0,"Empty file %s",fn);
  close_trn(fp);
}

bool fread_trnheader(int fp,t_trnheader *trn)
{
  return do_trnheader(fp,TRUE,trn);
}

void write_trn(char *fn,int step,real t,real lambda,
	       rvec *box,int natoms,rvec *x,rvec *v,rvec *f)
{
  int fp;
  
  fp = open_trn(fn,"w");
  do_trn(fp,FALSE,&step,&t,&lambda,box,&natoms,x,v,f);
  close_trn(fp);
}

void read_trn(char *fn,int *step,real *t,real *lambda,
	      rvec *box,int *natoms,rvec *x,rvec *v,rvec *f)
{
  int fp;
  
  fp = open_trn(fn,"r");
  (void) do_trn(fp,TRUE,step,t,lambda,box,natoms,x,v,f);
  close_trn(fp);
}

void fwrite_trn(int fp,int step,real t,real lambda,
		rvec *box,int natoms,rvec *x,rvec *v,rvec *f)
{
  (void) do_trn(fp,FALSE,&step,&t,&lambda,box,&natoms,x,v,f);
}


bool fread_trn(int fp,int *step,real *t,real *lambda,
	       rvec *box,int *natoms,rvec *x,rvec *v,rvec *f)
{
  return do_trn(fp,TRUE,step,t,lambda,box,natoms,x,v,f);
}

bool fread_htrn(int fp,t_trnheader *trn,rvec *box,rvec *x,rvec *v,rvec *f)
{
  return do_htrn(fp,TRUE,trn,box,x,v,f);
}

int open_trn(char *fn,char *mode)
{
  char *m;

  if (mode[0]=='r')
    m="rb";
  else
    m="wb";

  return fio_open(fn,m);
}

void close_trn(int fp)
{
  fio_close(fp);
}
