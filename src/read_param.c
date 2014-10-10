#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <dirent.h>

#include "common.h"

static int myrank_;

static int read_parameter_file(char const *fname,Parameters* const param);
static void bcast_string(char* string,int len);
static void bcast_array_double(double** parray,int* len);

//int read_parameters(const char filename[], Parameters* const param)
int read_parameters(char *filename)
{
  //Only process 0 reads params then communicates
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank_);
  if(myrank_==0) {
    int ret=read_parameter_file(filename,&Param);
    if(ret!=0)
      msg_abort(1001,"Error: Unable to read parameter file: %s\n",filename);
  }

  // Share parameters with other nodes
  MPI_Bcast(&Param,sizeof(Parameters),MPI_BYTE,0,MPI_COMM_WORLD);

  bcast_string(Param.snap_prefix,256);
  bcast_string(Param.snap_dir,256);
  bcast_string(Param.init_prefix,256);
  bcast_string(Param.init_dir,256);

  bcast_array_double(&(Param.aout),&(Param.n_aout));

  return 0;
}

static int linecount(FILE *f)
{
  //////
  // Counts #lines from file
  int i0=0;
  char ch[1000];
  while((fgets(ch,sizeof(ch),f))!=NULL) {
    i0++;
  }
  return i0;
}

static int read_parameter_file(char const *fname,Parameters *param)
{
  FILE *fi;
  int n_lin,ii;
  char fname_times[256];
  
  //Read parameters from file
  fi=fopen(fname,"r");
  if(fi==NULL)
    msg_abort(1000, "Error: couldn't open file %s\n",fname);

  n_lin=linecount(fi);
  rewind(fi);
  for(ii=0;ii<n_lin;ii++) {
    char s0[512],s1[64],s2[256];
    if(fgets(s0,sizeof(s0),fi)==NULL)
      msg_abort(1001,"Error reading line %d, file %s\n",ii+1,fname);
    if((s0[0]=='#')||(s0[0]=='\n')) continue;
    int sr=sscanf(s0,"%s %s",s1,s2);
    if(sr!=2)
      msg_abort(1002,"Error reading line %d, file %s\n",ii+1,fname);
    
    if(!strcmp(s1,"n_grid="))
      param->n_grid=atoi(s2);
    else if(!strcmp(s1,"ntimestep="))
      param->ntimestep=atoi(s2);
    else if(!strcmp(s1,"times_file="))
      sprintf(fname_times,"%s",s2);
    else if(!strcmp(s1,"np_alloc_factor="))
      param->np_alloc_factor=atof(s2);
    else if(!strcmp(s1,"loglevel="))
      param->loglevel=atoi(s2);
    else if(!strcmp(s1,"snap_dir="))
      sprintf(param->snap_dir,"%s",s2);
    else if(!strcmp(s1,"snap_prefix="))
      sprintf(param->snap_prefix,"%s",s2);
    else if(!strcmp(s1,"init_dir="))
      sprintf(param->init_dir,"%s",s2);
    else if(!strcmp(s1,"init_prefix="))
      sprintf(param->init_prefix,"%s",s2);
    else if(!strcmp(s1,"loglevel="))
      param->loglevel=atoi(s2);
    else
      msg_printf(info,"Unknown parameter %s\n",s1);
  }
  fclose(fi);


  //Read from times file
  fi=fopen(fname_times,"r");
  if(fi==NULL)
    msg_abort(1003, "Error: couldn't open file %s\n",fname_times);
  param->n_aout=linecount(fi);
  rewind(fi);
  param->aout=(double *)malloc(param->n_aout*sizeof(double));
  for(ii=0;ii<param->n_aout;ii++) {
    char s0[64];
    double z;
    if(fgets(s0,sizeof(s0),fi)==NULL)
      msg_abort(1004,"Error reading line %d, file %s\n",ii+1,fname_times);
    z=atof(s0);
    if(z<0 || z>100)
      msg_abort(1005,"Wrong output redshift %lf\n",z);
    param->aout[ii]=1.0/(1+z);
  }
  fclose(fi);

  double a_previous=param->aout[0];
  for(ii=1;ii<param->n_aout;ii++) {
    if(param->aout[ii]<=a_previous)
      msg_abort(1005,"Please provide redshifts in descending order\n");
    a_previous=param->aout[ii];
  }


  //Read IC header and check #files
  int nf=0;
  DIR *dir;
  struct dirent *ent;
  int prefix_len=strlen(param->init_prefix);
  GadgetHeader head=read_header_multiple(param->init_dir,param->init_prefix);

  print_header(head);
  param->n_part=(unsigned long long)(head.np_total[1])+
    ((unsigned long long)(head.np_total_highword[1]) << 32);
  param->boxsize=head.boxsize;
  param->omega_m=head.omega0;
  if(fabs(1-head.omega0-head.omega_lambda)>0.001)
    msg_printf(info,"WARNING: non-flat cosmology!!\n");
  param->h=head.hubble_param;
  param->n_files_init=head.num_files;
  param->a_init=1./(1+head.redshift);

  dir=opendir(param->init_dir);
  if(dir==NULL) {
    fprintf(stderr,"wrong directory %s\n",param->init_dir);
    exit(1);
  }
  while((ent=readdir(dir))!=NULL) {
    if(!strncmp(ent->d_name,param->init_prefix,prefix_len))
      nf++;
  }
  closedir(dir);
  if(nf!=param->n_files_init) {
    fprintf(stderr,"Wrong #files %d != %d\n",nf,param->n_files_init);
    exit(1);
  }


  return 0;
}

static void bcast_string(char* pstring,int len)
{
  const int n=len;

  const int ret2=MPI_Bcast(pstring,n,MPI_CHAR,0,MPI_COMM_WORLD);
  assert(ret2==MPI_SUCCESS);
}

static void bcast_array_double(double** parray,int* len)
{
  const int ret1= MPI_Bcast(len,1,MPI_INT,0,MPI_COMM_WORLD);
    assert(ret1==MPI_SUCCESS);

  const int n=*len;

  if(n==0) {
    *parray=0;
    return;
  }

  if(myrank_!=0) {
    *parray=malloc(sizeof(double)*n);
  }
  assert(*parray);

  const int ret2=MPI_Bcast(*parray,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
  assert(ret2 == MPI_SUCCESS);
}
