
Conversation opened. 1 unread message.

Skip to content
Using 연세대학교 Mail with screen readers
Pause

  More 
1 of 62
 
test_func.c 
Inbox
x 

남재욱
Attachments5:37 PM (28 minutes ago)
to me 

Images are not displayed. Display images below - Always display images from elqhkf@naver.com
   Translate message
Turn off for: Korean
Seq function화 한거하고 MPI 파일이에요 

2 Attachments
	
Click here to Reply or Forward
Using 3.72 GB
Manage
Program Policies
Powered by Google
Last account activity: 0 minutes ago
Open in 1 other location  Details


#include <petscksp.h>
#include <stdio.h>
#include <math.h>

int PETSCINIT(int argc,char **argv);
int MAT_GENER(int nl,Mat A);
int KSP_GENER(KSP ksp,Mat A);
int LI_SOLV(PetscInt npc,KSP ksp,Vec b,Vec x,long double *ig);
int DESTROY(KSP ksp,Mat A,Vec b,Vec x,Vec e);

int main(int argc,char **argv)
{
  KSP                ksp;
  Mat                A;
  Vec                x,b,e;
  PetscScalar        vp;
  PetscInt           ip,np,nps,npc;
  PetscErrorCode     ierr;
  int i,l,nl,ns,nc;
  long double *ig,*rhs;

  ierr=PETSC_INIT(argc,argv);
  nl=5; ns=pow(nl,2); nc=pow(nl,3);
  np=nl; nps=pow(np,2); npc=pow(np,3);

  ig=(long double*)malloc((nc-1)*sizeof(long double));
  rhs=(long double*)malloc((nc-1)*sizeof(long double));
  for (i=0; i<nc-1; i++) {
    ig[i]=0.0; rhs[i]=0.01;
  }
  
  MatCreate(PETSC_COMM_SELF,&A);
  MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,npc-1,npc-1);
  MatSetType(A,MATSEQAIJ);
  MatSeqAIJSetPreallocation(A,7,0);
  MatSetFromOptions(A);
  MatSetUp(A);
  ierr=MAT_GENER(nl,A);

  KSPCreate(PETSC_COMM_SELF,&ksp);
  ierr=KSP_GENER(ksp,A);

  VecCreate(PETSC_COMM_WORLD,&b);
  VecSetSizes(b,PETSC_DECIDE,npc-1);
  VecSetFromOptions(b);
  VecSetUp(b);
  VecDuplicate(b,&x);
  VecDuplicate(b,&e);

  for (ip=0; ip<npc-1; ip++) {
    i=ip; vp    = rhs[i];
    VecSetValues(b,1,&ip,&vp,INSERT_VALUES);
  }

  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  ierr=LI_SOLV(npc,ksp,b,x,ig);

  ierr=DESTROY(ksp,A,b,x,e);

  free(ig); free(rhs);

  return ierr;
}

int PETSC_INIT(int argc,char **argv)
{
  char help[] = "Solving Poisson Eqn for Pressure by function form()\n";
  PetscInitialize(&argc,&argv,0,help);

  return 0;
}

int MAT_GENER(int nl,Mat A)
{
  PetscScalar        vp;
  PetscInt           ip,jp,kp,sp,np,nps,npc;

  np=nl; nps=pow(np,2); npc=pow(np,3);

  for (ip=nps; ip<npc-1; ip++) {
    vp    = 1; jp=ip-nps;
    MatSetValues(A,1,&ip,1,&jp,&vp,INSERT_VALUES);
    MatSetValues(A,1,&jp,1,&ip,&vp,INSERT_VALUES);
  }
  for (ip=np; ip<nps-1; ip++) {
    vp    = 1; jp=ip-np;
    MatSetValues(A,1,&ip,1,&jp,&vp,INSERT_VALUES);
    MatSetValues(A,1,&jp,1,&ip,&vp,INSERT_VALUES);
  }
  for (kp=1; kp<np; kp++) {
    for (ip=kp*nps+np-1; ip<(kp+1)*nps-1; ip++) {
      vp    = 1; jp=ip-np;
      MatSetValues(A,1,&ip,1,&jp,&vp,INSERT_VALUES);
      MatSetValues(A,1,&jp,1,&ip,&vp,INSERT_VALUES);
    }
  }
  for (ip=1; ip<np-1; ip++) {
    vp    = 1; jp=ip-1;
    MatSetValues(A,1,&ip,1,&jp,&vp,INSERT_VALUES);
    MatSetValues(A,1,&jp,1,&ip,&vp,INSERT_VALUES);
  }
  for (kp=1; kp<nps; kp++) {
    for (ip=kp*np; ip<(kp+1)*np-1; ip++) {
      vp    = 1; jp=ip-1;
      MatSetValues(A,1,&ip,1,&jp,&vp,INSERT_VALUES);
      MatSetValues(A,1,&jp,1,&ip,&vp,INSERT_VALUES);
    }
  }
  for (ip=0; ip<np-2; ip++) {
    vp    = -4;
    MatSetValues(A,1,&ip,1,&ip,&vp,INSERT_VALUES);
  }
  ip=np-2; vp    = -3;
  MatSetValues(A,1,&ip,1,&ip,&vp,INSERT_VALUES);
  for (jp=2; jp<np+1; jp++) {
    for (ip=1; ip<np+1; ip++) {
      sp=np*(jp-1)+ip-2;
      if (jp==np&&ip==1) {vp    = -3;}
      else if (jp==np&&ip==np) {vp    = -3;}
      else if (jp==np&&ip>1&&ip<np) {vp    = -4;}
      else if (jp>1&&jp<np&&ip==1) {vp    = -4;}
      else if (jp>1&&jp<np&&ip==np) {vp    = -4;}
      else {vp    = -5;}
      MatSetValues(A,1,&sp,1,&sp,&vp,INSERT_VALUES);
    }
  }
  for (kp=2; kp<np+1; kp++) {
    for (jp=1; jp<np+1; jp++) {
      for (ip=1; ip<np+1; ip++) {
        sp=nps*(kp-1)+np*(jp-1)+ip-2;
        if (kp==np&&jp==1&&ip==1) {vp    = -3;}
        else if (kp==np&&jp==1&&ip==np) {vp    = -3;}
        else if (kp==np&&jp==1&&ip>1&&ip<np) {vp    = -4;}
        else if (kp==np&&jp==np&&ip==1) {vp    = -3;}
        else if (kp==np&&jp==np&&ip==np) {vp    = -3;}
        else if (kp==np&&jp==np&&ip>1&&ip<np) {vp    = -4;}
        else if (kp==np&&jp>1&&jp<np&&ip==1) {vp    = -4;}
        else if (kp==np&&jp>1&&jp<np&&ip==np) {vp    = -4;}
        else if (kp==np&&jp>1&&jp<np&&ip>1&&ip<np) {vp    = -5;}
        else if (kp>1&&kp<np&&jp==1&&ip==1) {vp    = -4;}
        else if (kp>1&&kp<np&&jp==1&&ip==np) {vp    = -4;}
        else if (kp>1&&kp<np&&jp==1&&ip>1&&ip<np) {vp    = -5;}
        else if (kp>1&&kp<np&&jp==np&&ip==1) {vp    = -4;}
        else if (kp>1&&kp<np&&jp==np&&ip==np) {vp    = -4;}
        else if (kp>1&&kp<np&&jp==np&&ip>1&&ip<np) {vp    = -5;}
        else if (kp>1&&kp<np&&jp>1&&jp<np&&ip==1) {vp    = -5;}
        else if (kp>1&&kp<np&&jp>1&&jp<np&&ip==np) {vp    = -5;}
        else {vp    = -6;}
        MatSetValues(A,1,&sp,1,&sp,&vp,INSERT_VALUES);
      }
    }
  }
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

  return 0;
}

int KSP_GENER(KSP ksp,Mat A)
{
  PC                 pc;

  KSPSetOperators(ksp,A,A);
  KSPSetType(ksp,KSPCG);
  KSPGetPC(ksp,&pc);
  PCSetType(pc,PCGAMG);
  PCGAMGSetType(pc,PCGAMGAGG);
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);
  KSPSetTolerances(ksp,PETSC_DEFAULT,1.e-5,PETSC_DEFAULT,2000);

  return 0;
}

int LI_SOLV(PetscInt npc,KSP ksp,Vec b,Vec x,long double *ig)
{
  PetscScalar        vp;
  PetscInt           ip,its;
  int i;

  KSPSolve(ksp,b,x);
  KSPGetIterationNumber(ksp,&its);
  for (ip=0; ip<npc-1; ip++) {
    VecGetValues(x,1,&ip,&vp);
    i=ip; ig[i]=vp;
  }
  PetscPrintf(PETSC_COMM_WORLD,"Iterations %D\n",its);
  VecView(x,PETSC_VIEWER_STDOUT_WORLD);

  return 0;
}

int DESTROY(KSP ksp,Mat A,Vec b,Vec x,Vec e)
{
  VecDestroy(&x); VecDestroy(&b); VecDestroy(&e);
  MatDestroy(&A); KSPDestroy(&ksp);
  PetscFinalize();

  return 0;
}
