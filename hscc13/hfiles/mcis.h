#ifndef _MYMODEL_H
#define _MYMODEL_H

#ifdef __cplusplus
extern "C" {
#endif





/* macros */

#define FMIN 0
#define FMAX 1

#define DEBUG 1

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) > (b) ? (b) : (a))


typedef struct sys_ sys_t;
typedef struct polytope_ polytope_t;


/* polytope P={ x | Hx <= h } */
struct polytope_ {
  
  /* matrix and vector */
  double *H;
  double *h;

  /* number of inequalities */
  int num;

  /* pointer to add representing the half-space */
  void **fmin;

  /* to compute fmax from fmin */
  double* offset;

  /* pointer to next polytope */
  polytope_t *next;

};

/* sys */
struct sys_ { 
  /*state space dimension */ 
  int dim;   
  /* number of inputs */
  int Nu; 
  /* ctrb indices (shifted) */
  int *mu; 
  /*number of bits per dim */
  int Nbits;   
  /*rounding precision of leaf node values*/
  double precision;
};


void mcis_init(sys_t*);

int mcis_close(void);

void mcis_clearSet(int);

void mcis_addPolytope(double*,double*,int);

int mcis_compute(int);


int mcis_query(int*);

char* mcis_getNodes(char*,int,int);

long mcis_getNumCubes(int);

double mcis_getError(int);

void mcis_getInfo(int);


/* helper functions (internal) */
void mcis_iter(int);

void* mcis_refine(void*,double,int,int);

void mcis_truncate(void**,polytope_t*,int,int);



#ifdef __cplusplus
}
#endif
#endif


