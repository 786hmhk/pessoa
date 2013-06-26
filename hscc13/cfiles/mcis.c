
#include <stdio.h>
#include <math.h>

#include "mcis.h"

#include "../cudd-2.5.0/include/cudd.h"
#include "../cudd-2.5.0/include/cuddInt.h"
#include "../cudd-2.5.0/include/dddmp.h"



#include "mex.h"



/* global variables */
/* bdd manager */
DdManager* ddman=NULL;

/* root nodes of bdd represnting the set of boxes */
/* inner approximation */
DdNode* iSet=NULL;       
/* outer approximation */
DdNode* oSet=NULL;       

/* infinity nodes */
DdNode* MINF=NULL;
DdNode* PINF=NULL;

/* system data */
sys_t* sys;
/* polytopes */
polytope_t* P=NULL;



/*
 * Initialize cudd package 
 */
void mcis_init(sys_t *init_sys) {

  /* intitialize struct with system data */
  sys=init_sys;

  /* init cudd manager */
  ddman = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS,0);

  MINF=Cudd_ReadMinusInfinity(ddman);
  Cudd_Ref(MINF);

  PINF=Cudd_ReadPlusInfinity(ddman);
  Cudd_Ref(PINF);

  /* set rounding precision  of leaf nodes in the add*/
  if (sys->precision>0)
    Cudd_SetRoundPre(ddman,sys->precision);
  else
    sys->precision=1.0e-12;
  
  /* rounding type */
  Cudd_SetRoundType(ddman,DD_FLOOR);

  /* Cudd_DisableGarbageCollection(ddman);
   */

}



/* init polytopes */
void mcis_addPolytope(double *H,double *h,int num) {

  int i,j;

  /* dimension */
  int dim=sys->dim;

  /* allocate memory */
  polytope_t *Pnew=malloc(sizeof*Pnew);

  Pnew->H=malloc(num*dim*sizeof*Pnew->H);
  Pnew->h=malloc(num*sizeof*Pnew->h);
  Pnew->num=num;

  Pnew->fmin=NULL;

  Pnew->next=P;

  P=Pnew;

  for(i=0;i<num;i++) {
    for(j=0;j<dim;j++) 
      P->H[i+j*num]=H[i+j*num];
    P->h[i]=h[i];
  }
}


/* compute maximal controlled invariant set */
/* use Nbits for each dimension */
int mcis_compute(int Nbits) {
 
  /* initialization done ? */
  if(!ddman)
    return 0;

  /* no polytope given */
  if(!P)
    return 0;

  /* dimension */
  int dim=sys->dim;
  /* number of binary variables in each dimension */
  sys->Nbits=Nbits;

  /* add that store min values for intervals */ 
  DdNode *fmin; /* fmin stores  min{ c'x | x in Box(b1b2...) } */


  /* inner and outer polytope representation on bdd */
  DdNode* iPol=NULL;
  DdNode* oPol=NULL;

  DdNode** pol=malloc(2*sizeof*pol);

  /* helper variables */
  int i,j,k,l;

  DdNode* one;
  DdNode* oTmp;

  polytope_t *Q;

  DdNode *dd1,*dd2;

  /* one */
  one = Cudd_ReadOne(ddman);
  Cudd_Ref(one);

  /* init temp bdd for outer approximation */
  if (oSet!=NULL) {
    oTmp=Cudd_bddAnd(ddman,oSet,one);
    Cudd_Ref(oTmp);
  } else {
    oTmp = Cudd_ReadOne(ddman);
    Cudd_Ref(oTmp);
  }

  /*
   *
   * code starts 
   *
   */ 

  Q=P;

  /* 
   * init  add (fmin) for each polytope
   */
  
  /* loop over polytopes */
  while(Q) {  

    /* allocate memory to store fmin and offset for each 
     * inequality in the polytope */
    Q->fmin=malloc(Q->num*sizeof*Q->fmin);
    Q->offset=malloc(Q->num*sizeof*Q->offset);

    /* loop over number of inequalities */
    for (k=0;k<Q->num;k++){ 

      fmin=Cudd_addConst(ddman,0);
      Cudd_Ref(fmin);

      Q->offset[k]=0;

      for(j=0;j<dim;j++) {

        if (Q->H[j*Q->num+k]>=0)
          Q->offset[k]-=Q->H[j*Q->num+k];
        else
          Q->offset[k]+=Q->H[j*Q->num+k];
         
        dd1=Cudd_addConst(ddman,Q->H[j*Q->num+k]);
        Cudd_Ref(dd1);

        /* depending on the sign */
        if (Q->H[j*Q->num+k]>=0) 
          dd2=Cudd_addApply(ddman,Cudd_addMinus,fmin,dd1);
        else 
          dd2=Cudd_addApply(ddman,Cudd_addPlus,fmin,dd1);

        Cudd_Ref(dd2);
        Cudd_RecursiveDeref(ddman,fmin);
        fmin=dd2;

        Cudd_RecursiveDeref(ddman,dd1);
      }

      Q->fmin[k]=(void*)fmin;

    } /* inequalities */ 

    /* go over next polytope */
    Q=Q->next;

  }/* done init add's */


  /* clear previous results (iSet) */
  mcis_clearSet(FMIN);

  /* MAIN LOOP */
  /* loop over number bits */
  for(i=0;i<Nbits;i++){

    /* clear previous results (oSet) */
    mcis_clearSet(FMAX);

    if(DEBUG) 
      mexPrintf("Resolution %d.",i);
 
    /* loop over number of polytopes */
    l=0;
    Q=P;
    while(Q) {

      if(DEBUG) 
        mexPrintf("\nPolytope %d: Inequality: ",l++);

      /* init polytope bdds:
       * exclude all points outside of outer 
       * approximation (oSet) */
      if (iPol!=NULL)
        Cudd_RecursiveDeref(ddman,iPol);
      if (oPol!=NULL)
        Cudd_RecursiveDeref(ddman,oPol);

      iPol = Cudd_bddAnd(ddman,oTmp,one);
      Cudd_Ref(iPol);

      oPol = Cudd_bddAnd(ddman,oTmp,one);
      Cudd_Ref(oPol);
      /* end init */


      /* iterate over number of inequalities of the polytope */
      for(k=0;k<Q->num;k++) {

        if(DEBUG) {
          mexPrintf("%d ",k);
          mexEvalString("drawnow;");
        }


        fmin=(DdNode*)Q->fmin[k];

        /* update add representation with one more bit*/
        for(j=0;j<dim;j++) 
          /* fmin = min { c'x | x in box } */
          fmin=(DdNode*) mcis_refine((void*)fmin,Q->H[j*Q->num+k],j,i);

        Q->fmin[k]=(void*)fmin;
        /* Cudd_PrintDebug(ddman,fmin,0,2); */
        
        /* truncate fmin to -inf and +inf 
         * and update the bdds that represent
         * the inner and outer approximation 
         * of the polytope (iPol and oPol)*/
        pol[0]=iPol;
        pol[1]=oPol;

        mcis_truncate((void**)pol,Q,k,i);

        iPol=pol[0];
        oPol=pol[1];

      } /* end number of inequalities */

      /* update iSet and oSet */
      dd1=Cudd_bddOr(ddman,iPol,iSet);
      Cudd_Ref(dd1);
      Cudd_RecursiveDeref(ddman,iSet);
      iSet=dd1;

      dd1=Cudd_bddOr(ddman,oPol,oSet);
      Cudd_Ref(dd1);
      Cudd_RecursiveDeref(ddman,oSet);
      oSet=dd1;
      /* end update */

   
      Q=Q->next;

    } /* end loop over polytopes */

    if(DEBUG) 
      mexPrintf("\n");

    /* compute outer approximation of mcis */
#if 0
    mcis_iter(FMAX);


    Cudd_RecursiveDeref(ddman,oTmp);
    oTmp=Cudd_bddAnd(ddman,oSet,one);
    Cudd_Ref(oTmp);
#endif


  } /* end loop over Nbits */


  /* compute mcis */
#if 1
  mcis_iter(FMAX);
#endif
  mcis_iter(FMIN);


  

  /* clean up add inequalities */
  Cudd_RecursiveDeref(ddman,iPol);
  Cudd_RecursiveDeref(ddman,oPol);

  Q=P;
  while(Q) {

    for(i=0;i<Q->num;i++) 
      Cudd_RecursiveDeref(ddman,Q->fmin[i]);

    free(Q->fmin);
    free(Q->offset);

    Q=Q->next;
  } /* done cleaning */

  Cudd_RecursiveDeref(ddman,one);
  Cudd_RecursiveDeref(ddman,oTmp);

  free(pol);

  return 1;

}


/*
 * Computes the mcis until it converges
 */
void mcis_iter(int mm) {

  if (!ddman)
    return;
  
  /* helper */
  DdNode *dd1, *dd2;
  DdNode *cube, *var;
  DdNode *zero;
  DdNode *NewbddSet;
  DdNode *bddSet;

  if (mm==FMIN)
    bddSet=iSet;
  else
    bddSet=oSet;

  int dim=sys->dim;
  int Nbits=sys->Nbits;

  int permute[dim*Nbits];

  int i,j;

  int idx=0;

  zero=Cudd_Not(Cudd_ReadOne(ddman));
  Cudd_Ref(zero);

  /* create cube to convert first Nbits to don't cares */
  /* if several inputs, add also the corresponding Nbits at positions
   * of the controllability intervals */
  cube=Cudd_ReadOne(ddman);
  Cudd_Ref(cube);


  /* controllability intervals (mu) */
  for (j=0;j<sys->Nu;j++) {
    if (j>0)
      idx+=sys->mu[j-1];

    for (i=0;i<Nbits;i++) {

      var=Cudd_bddIthVar(ddman,idx+dim*i);
      dd1=Cudd_bddAnd(ddman,cube,var);
      Cudd_Ref(dd1);
 
      Cudd_RecursiveDeref(ddman,cube);
      cube=dd1;
    }
  }



  /* perpare the permutation array */

  for (i=0;i<dim*Nbits;i++) 
    permute[i]=i+1;
  for (i=0;i<Nbits;i++) 
    permute[(i+1)*dim-1]=i*dim;


  i=0;

  NewbddSet=Cudd_ReadOne(ddman);
  Cudd_Ref(NewbddSet);
  /* iterate until convergence */
  /* while new != old */
  while(!Cudd_EquivDC(ddman,NewbddSet,bddSet,zero)) {
    i++;

    Cudd_RecursiveDeref(ddman,NewbddSet);

    /* begin: pre computation */
    /* 1. shift */
    /*dd1=Cudd_bddExistAbstract(ddman,bddSet,cube);*/
    dd1=Cudd_bddPermute(ddman,bddSet,permute);
    Cudd_Ref(dd1);

    /* 2. set don't cares to 1st (and subsequent) dim */
    dd2=Cudd_bddExistAbstract(ddman,dd1,cube);
    /*dd2=Cudd_bddPermute(ddman,dd1,permute);*/
    Cudd_Ref(dd2);
    /* end: pre computation */


   /* begin: compute the intersection */
    NewbddSet=Cudd_bddAnd(ddman,bddSet,dd2);
    Cudd_Ref(NewbddSet);
    Cudd_RecursiveDeref(ddman,dd1);
    Cudd_RecursiveDeref(ddman,dd2);
    /* end: compute the intersection */

    /*
    dd1=Cudd_Not(NewbddSet);
    Cudd_Ref(dd1);
    
    dd2=Cudd_bddAnd(ddman,bddSet,dd1);
    Cudd_Ref(dd2);

   Cudd_PrintDebug(ddman,dd2,0,1); 
   Cudd_RecursiveDeref(ddman,dd1);
    Cudd_RecursiveDeref(ddman,dd2);
    */

    /* point to old iteration */
    dd1=bddSet;

    /* restart iteration */
    bddSet=NewbddSet;

    NewbddSet=dd1;

  }


  if(DEBUG) {
    mexPrintf("MCIS computed after %d iterations.\n",i);
    mexEvalString("drawnow;");
  }

  if (mm==FMIN) 
    iSet=bddSet;
  else
    oSet=bddSet;

  Cudd_RecursiveDeref(ddman,NewbddSet);
  Cudd_RecursiveDeref(ddman,cube);
  Cudd_RecursiveDeref(ddman,zero);
}


int mcis_query(int *query) {
  
  /* bdd form query cube array */
  DdNode *f;
  
  int i,j,k;

  int *cube=malloc(sys->Nbits*sys->dim*sizeof*cube);

  /* integer to binary */
  for(j=0;j<sys->dim;j++) {
    k=0;
    for(i=pow(2.0f,(double)(sys->Nbits-1)); i>0 ; i=i>>1)  {

      if (query[j] & i )
        cube[j+k*sys->dim] = 1;
      else
        cube[j+k*sys->dim] = 0;

      k++;

    }
  }

  /* for(j=0;j<sys->dim*sys->Nbits;j++)
   * printf("%d",cube[j]);
   */


  /* Cudd_PrintDebug(ddman,iSet,0,2); 
   */
  

  f=Cudd_CubeArrayToBdd(ddman,cube);
  Cudd_Ref(f);

  if  (Cudd_bddLeq(ddman,f,iSet)) {
    Cudd_RecursiveDeref(ddman,f);
    free(cube);
    return 1;
  }
  else {
    Cudd_RecursiveDeref(ddman,f);
    free(cube);
    return 0;
  }
}
  


int mcis_close(void) {

  /* is dd manager initialized */
  if (!ddman)
    return 0;

  /* FILE * fp; 
   * fp = fopen("bddInfo","w");
   * Cudd_DumpDot(ddman,1,&iSet,NULL,NULL,fp);
   * Cudd_PrintInfo(ddman,fp);
   * 
   * fclose(fp);
   */

  polytope_t *Pn;

  if (P) {

    free(P->H);
    free(P->h);

    Pn=P->next;

    free(P);

    P=Pn;

  }

  /*counts the number of references after dereferncing */
  int refcount;

  Cudd_RecursiveDeref(ddman,MINF);
  Cudd_RecursiveDeref(ddman,PINF);



  if (iSet)
    Cudd_RecursiveDeref(ddman,iSet);
  if (oSet) 
    Cudd_RecursiveDeref(ddman,oSet);

  refcount=Cudd_CheckZeroRef(ddman); 

  Cudd_Quit(ddman);

  ddman=NULL;

  return refcount;
 
}


/* get numbers of nodes inside the set */
long mcis_getNumCubes(int approx) {

  if (!ddman)
    return 0;

  if (!approx) {
    /* printf("Dag size  %d\n", Cudd_DagSize(iSet));
     * mexEvalString("drawnow;");
     */
    return Cudd_CountPathsToNonZero(iSet);
  }
  else {
    return Cudd_CountPathsToNonZero(oSet);
  }


}

/* function provides an array of char where each row 
 * corresponds to a node that is inside the 
 * inner or outer approximation
 */

char* mcis_getNodes(char* bddCubes, int approx, int m) {

  /* bdd not initialized*/
  if (ddman==NULL)
    return NULL;

  /* system specific numbers */
  int dim=sys->dim;
  int Nbits=sys->Nbits;

  /* get cubes (string of 0 and 1 for 
   * which the bdd evaluates to 1) */
  DdGen* gen;
  int **cube;
  double* value;

  /* inner or outer approximation */
  DdNode* set;

  if (!approx)
    set=iSet;
  else
    set=oSet;

  /* helper */
  int i,j,k;
  int suc;

  /* allocate memory */
  value=malloc(sizeof*value);
  cube=malloc(dim*Nbits*sizeof*cube);

  /* get first cube and prepare for iteration 
   * over all cubes */
  gen=Cudd_FirstCube(ddman,set,cube,value);
  if (gen)
    suc=1;
  else
    suc=0;


  j=0;
  k=0;

  while (suc) {

    if (fmod(k,m)==0) {

      /* copy cube to matlab workspace */
      for(i=0; i<dim*Nbits; i++) 
        bddCubes[j*dim*Nbits+i]=cube[0][i]+'0';

      j++;
    }

    suc=Cudd_NextCube(gen,cube,value);
    k++;

  }


  free(cube);
  free(value);
  
  return bddCubes;
}


/* mcis_truncate: nodes (cubes) that are ensured 
 * to be outside and inside the inequality
 * are determined and fmin as well as inner (iPol)
 * and outer (oPol) approximation of the polytope
 * are updated accordingly */
void mcis_truncate(void** pol, polytope_t* P, int k, int i) {

  /* helper */
  DdNode *dd1, *dd2, *dd3;

  DdNode *minf, *pinf;

  /* iPol=Pol[0] and oPol=Pol[1] */
  DdNode** Pol=(DdNode **)pol;

  /* rounding correction term */
  double rc=sys->Nbits*sys->dim*sys->precision;

  /* nodes that are contained in the set */
  dd1=Cudd_addNegate(ddman,P->fmin[k]);
  Cudd_Ref(dd1);
  dd2=Cudd_addBddThreshold(ddman,dd1,-P->h[k]-P->offset[k]*pow(2,-i)+rc);
  Cudd_Ref(dd2);
  Cudd_RecursiveDeref(ddman,dd1);

  /* update iPol with current inequality */
  dd3=Cudd_bddAnd(ddman,Pol[0],dd2);
  Cudd_Ref(dd3);
  Cudd_RecursiveDeref(ddman,Pol[0]);
  Pol[0]=dd3;
  /* done with iPol update */

  /* update fmin with -infty for nodes that are inside */
  dd1=Cudd_BddToAdd(ddman,dd2);
  Cudd_Ref(dd1);
  Cudd_RecursiveDeref(ddman,dd2);

  minf=Cudd_addApply(ddman,Cudd_addTimes,dd1,MINF);
  Cudd_Ref(minf);
  Cudd_RecursiveDeref(ddman,dd1);

  dd2=Cudd_addApply(ddman,Cudd_addPlus,P->fmin[k],minf);
  Cudd_Ref(dd2);
  Cudd_RecursiveDeref(ddman,minf);
  Cudd_RecursiveDeref(ddman,P->fmin[k]);
  P->fmin[k]=dd2;
  /* done with function update */



  /* nodes (cubes) that can be excluded the set */
  dd1=Cudd_addBddStrictThreshold(ddman,P->fmin[k],P->h[k]);
  Cudd_Ref(dd1);

  /* update oPol with current inequality */
  dd2=Cudd_Not(dd1);
  Cudd_Ref(dd2);

  dd3=Cudd_bddAnd(ddman,Pol[1],dd2);
  Cudd_Ref(dd3);
  Cudd_RecursiveDeref(ddman,dd2);
  Cudd_RecursiveDeref(ddman,Pol[1]);
  Pol[1]=dd3;

  /* done with oPol update */

  /* fmin +infty for nodes that are outside */
  dd2=Cudd_BddToAdd(ddman,dd1);
  Cudd_Ref(dd2);
  Cudd_RecursiveDeref(ddman,dd1);

  pinf=Cudd_addApply(ddman,Cudd_addTimes,dd2,PINF);
  Cudd_Ref(pinf);
  Cudd_RecursiveDeref(ddman,dd2);

  dd2=Cudd_addApply(ddman,Cudd_addPlus,P->fmin[k],pinf);
  Cudd_Ref(dd2);
  Cudd_RecursiveDeref(ddman,pinf);
  Cudd_RecursiveDeref(ddman,P->fmin[k]);
  P->fmin[k]=dd2;
  /* done with function update */

}



/* 
 * refine add by one binary variable in dimension d at given level
 * 
 * take care if function corresponds to min or max of a cell
 *
 */
void* mcis_refine(void* fv, double c, int d, int level) {

  /* helper */
  DdNode *dd1,*dd2;
  DdNode *con; 
  DdNode *var;
  
  DdNode *fun=(DdNode *)fv;

  /* 
   * f(b1...bi-1bi)=f(b1...bi-1)+1/2^(i+1)*c*bi
   */

  /* get pointer to negation of add variable */

  if(c<0) {

    /* var=Cudd_bddIthVar(ddman,d*sys->Nbits+level); */
    var=Cudd_bddIthVar(ddman,d+level*sys->dim);
    Cudd_Ref(var);

    dd2=Cudd_Not(var);
    Cudd_Ref(dd2);

    Cudd_RecursiveDeref(ddman,var);

    var=Cudd_BddToAdd(ddman,dd2);
    Cudd_Ref(var);

    Cudd_RecursiveDeref(ddman,dd2);

    c=-c;

  } else {

    /* var=Cudd_bddIthVar(ddman,d*sys->Nbits+level); */
    var=Cudd_addIthVar(ddman,d+level*sys->dim);
    Cudd_Ref(var);
  }

  con=Cudd_addConst(ddman,pow(2,-level)*c);
  Cudd_Ref(con);

  dd1=Cudd_addApply(ddman,Cudd_addTimes,var,con);
  Cudd_Ref(dd1);
 
  dd2=Cudd_addApply(ddman,Cudd_addPlus,dd1,fun);
  Cudd_Ref(dd2);

  Cudd_RecursiveDeref(ddman,fun);
  Cudd_RecursiveDeref(ddman,dd1);
  Cudd_RecursiveDeref(ddman,var);
  Cudd_RecursiveDeref(ddman,con);

  fun=dd2;


  return (void*)fun;

}

double mcis_getError(int bound) {

  /* is dd manager initialized */
  if (!ddman)
    return HUGE_VAL;

  int dim=sys->dim;
  int Nbits=sys->Nbits;



  /* number of boxes in inner and outer approx */
  double Ni=0;
  double No=0;
  double dN=0;


  /* get cubes (string of 0 and 1 for 
   * which the bdd evaluates to 1) */
  DdGen* gen;
  int **cube;
  double* value;

  /* helper */
  int i,j,k;
  int suc;
  DdNode *set;
  unsigned long n;


  /* allocate memory */
  value=malloc(sizeof*value);
  cube=malloc(dim*Nbits*sizeof*cube);




  
  for (j=0;j<2;j++) {
    
    /* init box counter */
    n=0;

    if (j==0)
      set=iSet; /* inner approximation */
    else
      set=oSet; /* outer approximation */

    /* get first cube and prepare for iteration over all cubes */
    gen=Cudd_FirstCube(ddman,set,cube,value);
    if (gen && Cudd_CountPathsToNonZero(set))
      suc=1;
    else
      suc=0;

    while (suc) {
      n++;
      /* account for don't cares */
      k=0;
      for(i=0; i<dim*Nbits; i++) {
        if (cube[0][i]==2) {
          n+=pow(2,k);
          k++;
        }
      }
      /* next cube */
      suc=Cudd_NextCube(gen,cube,value);
    }

    if (j==0)
      Ni=n; /* # cubes of inner approximation */
    else
      No=n; /* # cubes of outer approximation */

  }

  free(cube);
  free(value);


  dN=No-Ni;

  if (Ni==0)
    return HUGE_VAL;

  /* upper bound */
  if (bound)
    return dN/Ni;
  else /* lower bound */
    return dN/No;

}

void mcis_getInfo(int p) {

  /* is dd manager initialized */
  if (!ddman)
    return;

  FILE * fp; 

  /* print bdds? */
  if (p) {
    fp = fopen("innerApprox.dot","w");
    Cudd_DumpDot(ddman,1,&iSet,NULL,NULL,fp);
    Cudd_PrintInfo(ddman,fp);
    fclose(fp);
  
    fp = fopen("outerApprox.dot","w");
    Cudd_DumpDot(ddman,1,&oSet,NULL,NULL,fp);
    Cudd_PrintInfo(ddman,fp);
    fclose(fp);

  }

  mexPrintf("Number of unique constants: %d\n", Cudd_ReadKeys(ddman));
  mexPrintf("Number of unique constants: %d\n", ddman->constants.keys);

  mexPrintf("Dag size of inner approx: %d\n", Cudd_DagSize(iSet));
  mexPrintf("Dag size of outer approx: %d\n", Cudd_DagSize(oSet));

  mexPrintf("Num of cubes inner approx  %.0lf\n", Cudd_CountPathsToNonZero(iSet));
  mexPrintf("Num of cubes outer approx  %.0lf\n", Cudd_CountPathsToNonZero(oSet));


  mexPrintf("%ld KByte memory in use.\n", Cudd_ReadMemoryInUse(ddman)/1024);
  mexPrintf("Unique lookups %6.0lf\n",Cudd_ReadUniqueLookUps(ddman));
  mexPrintf("Unique links %6.0lf\n",Cudd_ReadUniqueLinks(ddman));
  mexPrintf("Num of garbage collections %d\n",Cudd_ReadGarbageCollections(ddman));
  mexPrintf("Time of garbage collections %ld\n",Cudd_ReadGarbageCollectionTime(ddman));
  mexPrintf("Cache hit ratio %.2lf\n",Cudd_ReadCacheHits(ddman)/Cudd_ReadCacheLookUps(ddman));
  mexPrintf("Rounding precision  %lf\n", Cudd_ReadRoundPre(ddman));
  mexEvalString("drawnow;");

}

void mcis_clearSet(int m) {

  if(!ddman)
    return;

  /* inner approximation */
  if (m==FMIN) {

    if (iSet)
      Cudd_RecursiveDeref(ddman,iSet);

    iSet = Cudd_Not(Cudd_ReadOne(ddman));
    Cudd_Ref(iSet);
  }

  /* outer approximation */
  if (m==FMAX) {

    if (oSet)
      Cudd_RecursiveDeref(ddman,oSet);

    oSet = Cudd_Not(Cudd_ReadOne(ddman));
    Cudd_Ref(oSet);
  }

}




