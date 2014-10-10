#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <dirent.h>

#include "common.h"

static void bad_block(void)
{
  msg_abort(123,"Bad block!\n");
}

static void my_fread(void *ptr,size_t size,size_t count,FILE *stream)
{
  size_t stat=fread(ptr,size,count,stream);
  if(stat!=count) {
    msg_abort(123,"Error freading\n");
  }
}

static FILE *my_fopen(char *fname,char *mode)
{
  FILE *f=fopen(fname,mode);
  if(f==NULL) {
    msg_abort(123,"Couldn't open file %s\n",fname);
  }

  return f;
}

void print_header(GadgetHeader header)
{
  int i;

  printf("Snapshot header\n");
  printf("  #part. in this file: ");
  for(i=0;i<6;i++) printf("%d, ",header.np[i]);
  printf("\n");

  printf("  #part. in total:     ");
  for(i=0;i<6;i++) {
    unsigned long long npa=(unsigned long long)(header.np_total[i]);
    unsigned long long npb=((unsigned long long)(header.np_total_highword[i]) << 32);
    unsigned long long np=npa+npb;
    printf("%llu, ",np);
  }
  printf("\n");

  printf("  #part. mass (10^10 M_sun/h) :");
  for(i=0;i<6;i++) printf("%.2lE, ",header.mass[i]);
  printf("\n");

  printf("  Redshift : %.2lE, Scale factor: %.2lE\n",
	 header.redshift,header.time);

  printf("  O_M = %.2lf, O_L = %.2lf, h = %.2lf\n",
	 header.omega0,header.omega_lambda,header.hubble_param);

  printf("  L_box = %.2lE Mpc/h\n",header.boxsize);

  printf("  Total number of files: %d\n",header.num_files);
}

GadgetHeader read_header_multiple(char *dirname,char *prefix)
{
  int i;
  struct dirent *ent;
  DIR *dir;
  char **fname_array;
  int n_files=0;
  int prefix_len=strlen(prefix);

  dir=opendir(dirname);
  //Find number of files
  if(dir==NULL) {
    msg_abort(123,"wrong directory %s\n",dirname);
  }
  while((ent=readdir(dir))!=NULL) {
    if(!strncmp(ent->d_name,prefix,prefix_len)) {
      printf("%s\n",ent->d_name);
      n_files++;
    }
  }
  closedir(dir);
  printf("Found %d files\n",n_files);

  //Collect filenames
  fname_array=(char **)malloc(n_files*sizeof(char *));
  for(i=0;i<n_files;i++)
    fname_array[i]=(char *)malloc(256*sizeof(char));

  dir=opendir(dirname);
  n_files=0;
  while((ent=readdir(dir))!=NULL) {
    if(!strncmp(ent->d_name,prefix,prefix_len)) {
      sprintf(fname_array[n_files],"%s/%s",dirname,ent->d_name);
      n_files++;
    }
  }
  closedir(dir);

  //Check first file's header
  int blklen1;
  GadgetHeader header;
  FILE *fp=my_fopen(fname_array[0],"r");
  my_fread(&blklen1,sizeof(int),1,fp);
  if(blklen1!=sizeof(GadgetHeader)) bad_block();
  my_fread(&header,sizeof(GadgetHeader),1,fp);
  my_fread(&blklen1,sizeof(int),1,fp);
  if(blklen1!=sizeof(GadgetHeader)) bad_block();
  fclose(fp);

  if(header.num_files!=n_files) {
    fprintf(stderr,"Wrong number of files! %d != %d\n",
	    header.num_files,n_files);
  }

  return header;
}

void read_initial_condition(Particles *particles)
{
  int i;
  struct dirent *ent;
  DIR *dir;
  char **fname_array;
  int n_files=0;
  int prefix_len=strlen(Param.init_prefix);
  double l_box=Param.boxsize;
  double i_l_box=1./l_box;
  int nnode=comm_nnode();
  int this_node=comm_this_node();

  dir=opendir(Param.init_dir);
  //Find number of files
  if(dir==NULL) {
    msg_abort(123,"wrong directory %s\n",Param.init_dir);
  }
  while((ent=readdir(dir))!=NULL) {
    if(!strncmp(ent->d_name,Param.init_prefix,prefix_len))
      n_files++;
  }
  closedir(dir);
  msg_printf(info,"Found %d files\n",n_files);

  //Collect filenames
  fname_array=(char **)malloc(n_files*sizeof(char *));
  for(i=0;i<n_files;i++)
    fname_array[i]=(char *)malloc(256*sizeof(char));

  dir=opendir(Param.init_dir);
  n_files=0;
  while((ent=readdir(dir))!=NULL) {
    if(!strncmp(ent->d_name,Param.init_prefix,prefix_len)) {
      sprintf(fname_array[n_files],"%s/%s",Param.init_dir,ent->d_name);
      n_files++;
    }
  }
  closedir(dir);

  //Check first file's header
  int blklen1,blklen2;
  unsigned long long np_tot,np_saved,np_read;
  GadgetHeader header;
  FILE *fp=my_fopen(fname_array[0],"r");
  my_fread(&blklen1,sizeof(int),1,fp);
  if(blklen1!=sizeof(GadgetHeader)) bad_block();
  my_fread(&header,sizeof(GadgetHeader),1,fp);
  my_fread(&blklen1,sizeof(int),1,fp);
  if(blklen1!=sizeof(GadgetHeader)) bad_block();
  fclose(fp);

  if(header.num_files!=n_files || Param.n_files_init) {
    msg_abort(123,"Wrong number of files! %d != %d != %d\n",
	      header.num_files,n_files,Param.n_files_init);
  }
  
  //Compute total number of particles
  np_tot=(unsigned long long)(header.np_total[1])+
    ((unsigned long long)(header.np_total_highword[1]) << 32);
  Particle *p=particles->p;
  
  //Iterate through files
  np_saved=0;
  np_read=0;
  for(i=0;i<n_files;i++) {
    int np_here;
    int j;

    fp=my_fopen(fname_array[i],"r");
    my_fread(&blklen1,sizeof(int),1,fp);
    if(blklen1!=sizeof(GadgetHeader)) bad_block();
    my_fread(&header,sizeof(GadgetHeader),1,fp);
    my_fread(&blklen1,sizeof(int),1,fp);
    if(blklen1!=sizeof(GadgetHeader)) bad_block();
    
    np_here=header.np[1];
    if(np_read+np_here>np_tot) {
      fprintf(stderr,"There's something wrong with the number of particles!\n");
      exit(1);
    }

    my_fread(&blklen1,sizeof(int),1,fp);
    for(j=0;j<np_here;j++) {
      int inode;
      my_fread(p,sizeof(Particle),1,fp);
      inode=(int)(p->x[0]*i_l_box*nnode);
      if(inode<0||inode>=nnode)
	msg_abort(123,"Wrong particle position %lE\n",p->x[0]);
      if(inode==this_node) {
	if(np_saved>=particles->np_allocated)
	  msg_abort(123,"Too many particles in file, enlarge alloc_factor\n");
	p++;
	p->v[0]=0;
	p->v[1]=0;
	p->v[2]=0;
	np_saved++;
      }
    }
    my_fread(&blklen2,sizeof(int),1,fp);
    if(blklen1!=blklen2) bad_block();

    fclose(fp);
    np_read+=np_here;
  }
  if(np_read!=np_tot) {
    fprintf(stderr,"There's something wrong with the number of particles!\n");
    exit(1);
  }

  particles->np_local=np_saved;
  particles->a_x=Param.a_init;
  particles->a_v=Param.a_init;

  unsigned long long np_global;
  MPI_Reduce(&np_saved,&np_global,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
  if(this_node == 0 && np_global!=np_tot)
    msg_abort(123,"Not all particles were saved %llu != %llu\n",np_global,np_tot);
}

void write_snapshot(const char filebase[],Snapshot const * const snapshot)
{
  char filename[256];
  sprintf(filename,"%s.%d",filebase,comm_this_node());

  FILE* fp=fopen(filename,"w");
  if(fp==0)
    msg_abort(9000,"Error: Unable to write to file: %s\n",filename);

  ParticleMinimum* const p=snapshot->p;
  const int np=snapshot->np_local;
  const double boxsize=Param.boxsize;
  const double omega_m=Param.omega_m;
  const double h=Param.h;

  msg_printf(normal,"Longid is used for GADGET snapshot. %d-byte.\n", 
	     sizeof(unsigned long long));

  long long np_send=np,np_total;
  MPI_Reduce(&np_send,&np_total,1,MPI_LONG_LONG, MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Bcast(&np_total,1,MPI_LONG_LONG,0,MPI_COMM_WORLD);

  if(np_total!=Param.n_part) {
    msg_abort(123,"Wrong number of particles %ll != %llu\n",
	      np_total,Param.n_part);
  }

  GadgetHeader header; assert(sizeof(GadgetHeader)==256);
  memset(&header,0,sizeof(GadgetHeader));

  const double rho_crit=27.7455;
  const double m=omega_m*rho_crit*pow(boxsize,3.0)/np_total;
  
  header.np[1]=np;
  header.mass[1]=m;
  header.time=snapshot->a;
  header.redshift=1.0/header.time-1;
  header.np_total[1]=(unsigned int)np_total;
  header.np_total_highword[1]=(unsigned int)(np_total >> 32);
  header.num_files=comm_nnode();
  header.boxsize=boxsize;
  header.omega0=omega_m;
  header.omega_lambda=1.0-omega_m;
  header.hubble_param=h;


  int blklen=sizeof(GadgetHeader);
  fwrite(&blklen,sizeof(blklen),1,fp);
  fwrite(&header,sizeof(GadgetHeader),1,fp);
  fwrite(&blklen,sizeof(blklen),1,fp);

  // position
  blklen=np*sizeof(float)*3;
  fwrite(&blklen,sizeof(blklen),1,fp);
  for(int i=0;i<np;i++) {
    fwrite(p[i].x,sizeof(float),3,fp);
  }
  fwrite(&blklen,sizeof(blklen),1,fp);

  // velocity
  const float vfac=1.0/sqrt(snapshot->a); // Gadget convention

  fwrite(&blklen,sizeof(blklen),1,fp);
  for(int i=0;i<np;i++) {
    float vout[]={vfac*p[i].v[0],vfac*p[i].v[1],vfac*p[i].v[2]};
    fwrite(vout,sizeof(float),3,fp);
  }
  fwrite(&blklen,sizeof(blklen),1,fp);

  // id
  blklen=np*sizeof(unsigned long long);
  fwrite(&blklen,sizeof(blklen),1,fp);
  for(int i=0;i<np;i++) {
    unsigned long long id_out=p[i].id;
    fwrite(&id_out,sizeof(unsigned long long),1,fp); 
  }
  fwrite(&blklen,sizeof(blklen),1,fp);

  fclose(fp);

  msg_printf(normal, "snapshot %s written\n", filebase);
}
