/***********************************************************************

	PESSOA Version 1.4  
	------------------

Cyber-Physical Systems Laboratory
http://www.cyphylab.ee.ucla.edu
Author(s):	Pritam Roy pritam.roy@gmail.com

Dependencies:	CUDD library, MEX library
University of California, Los Angeles.
December 2010.

************************************************************************/

#include "mex.h"
#include "pessoa.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "cudd.h"
#include "dddmp.h"

using namespace std;

static DdManager* ddman;     // global ddmanager pointer variable


DdNode* createADDTree(DdManager* dd,int start, int bits){
	
	DdNode* v = Cudd_ReadZero(dd);
	Cudd_Ref(v);
	
	for(int i=start;i<(start+bits);i++){
		DdNode* tmp = Cudd_addIthVar(dd,i);
		Cudd_Ref(tmp);
		
		DdNode* cons = Cudd_addConst(dd,pow(2,bits-i+start-1));
		Cudd_Ref(cons);
		
		DdNode* tmp1 = Cudd_addApply(dd,Cudd_addTimes,tmp,cons);
		Cudd_Ref(tmp1);
		Cudd_RecursiveDeref(dd,tmp);
		Cudd_RecursiveDeref(dd,cons);
		DdNode* tmp2 = Cudd_addApply(dd,Cudd_addPlus,tmp1,v);
		Cudd_Ref(tmp2);
		Cudd_RecursiveDeref(dd,tmp1);
		Cudd_RecursiveDeref(dd,v);
		v = tmp2;
	}
	
	DdNode* cons1 = Cudd_addConst(dd,1);
	Cudd_Ref(cons1);
	
	DdNode* v1 = Cudd_addApply(dd,Cudd_addPlus,v,cons1);
	Cudd_Ref(v1);
	Cudd_RecursiveDeref(dd,v);
	Cudd_RecursiveDeref(dd,cons1);
	return v1;
}

DdNode* func_nondeter(DdManager* dd,int n,int m,double ad[],double bd[],double tau,double mu,double eta,
		 double vmax[],double minp[],double xset[],double uset[],int xoffset[],int uoffset[],int xoob[],int xbits[],int ubitsx[]){
	
	DdNode* one;		
	
	//Cudd_AutodynEnable(dd,CUDD_REORDER_SIFT);
	
	
	int totalxbits = 0;
	for(int i=0;i<n;i++)
		totalxbits += xbits[i];
	int ubits = 0 ; 
	for(int i=0;i<m;i++)
		ubits += ubitsx[i];
	
	
	int bits = 0;
	DdNode* bxtree = Cudd_ReadOne(dd);
	Cudd_Ref(bxtree);
	DdNode* xtreei[n];
	for(int i=0;i<n;i++){
		xtreei[i] = createADDTree(dd,bits,xbits[i]);
		Cudd_Ref(xtreei[i]);
		DdNode* tmp = Cudd_addBddInterval(dd,xtreei[i], (xset[2*i]/eta)+1-xoffset[i],(xset[2*i+1]/eta)+1-xoffset[i]);
		Cudd_Ref(tmp);
		DdNode* tmp1 = Cudd_bddAnd(dd,tmp,bxtree);
		Cudd_Ref(tmp1);
		Cudd_RecursiveDeref(dd,bxtree);
		Cudd_RecursiveDeref(dd,tmp);
		bxtree = tmp1;
		bits += xbits[i];
	}
	DdNode* butree = Cudd_ReadOne(dd);
	Cudd_Ref(butree);
	DdNode* utreei[m];
	for(int i=0;i<m;i++){
		utreei[i] = createADDTree(dd,bits,ubitsx[i]);
		Cudd_Ref(utreei[i]);
		DdNode* tmp = Cudd_addBddInterval(dd,utreei[i],(uset[2*i]/mu)+1-uoffset[i],(uset[2*i+1]/mu)+1-uoffset[i]); 
		Cudd_Ref(tmp);
		DdNode* tmp1 = Cudd_bddAnd(dd,tmp,butree);
		Cudd_Ref(tmp1);
		Cudd_RecursiveDeref(dd,butree);
		Cudd_RecursiveDeref(dd,tmp);
		butree = tmp1;
		bits+= ubitsx[i];
	}
	
	
	DdNode* bxu = Cudd_bddAnd(dd,bxtree,butree);
	Cudd_Ref(bxu);
	Cudd_RecursiveDeref(dd,bxtree);
	Cudd_RecursiveDeref(dd,butree);
	DdNode* xu = Cudd_BddToAdd(dd,bxu);
	Cudd_Ref(xu);
	
	DdNode* bnxtree = Cudd_ReadOne(dd);
	Cudd_Ref(bnxtree);
	DdNode* xoobBdd = Cudd_ReadOne(dd);
	Cudd_Ref(xoobBdd);
	DdNode* nxtreei[n];
	DdNode* nxset[n];
	for(int i=0;i<n;i++){
		nxtreei[i] = createADDTree(dd,bits,xbits[i]);
		Cudd_Ref(nxtreei[i]);
		DdNode* tmp = Cudd_addBddInterval(dd,nxtreei[i], (xset[2*i]/eta)+1-xoffset[i],(xset[2*i+1]/eta)+1-xoffset[i]);  
		Cudd_Ref(tmp);
		DdNode* nx  = Cudd_BddToAdd(dd,tmp);
		Cudd_Ref(nx);
		nxset[i] = Cudd_addApply(dd,Cudd_addTimes,nx,nxtreei[i]);
		
		Cudd_Ref(nxset[i]);
		
		DdNode* tmp1 = Cudd_bddAnd(dd,tmp,bnxtree);
		Cudd_Ref(tmp1);
		bnxtree = tmp1;
		DdNode* tmp2 = Cudd_addBddInterval(dd,nxtreei[i],xoob[i]+1,xoob[i]+1); 
		Cudd_Ref(tmp2);
		DdNode* tmp3 = Cudd_bddAnd(dd,tmp2,xoobBdd);
		Cudd_Ref(tmp3);
		Cudd_RecursiveDeref(dd,xoobBdd);
		xoobBdd = tmp3;
		bits += xbits[i];
	}
	
	DdNode* trout = Cudd_ReadLogicZero(dd);
	Cudd_Ref(trout);
	DdNode* lmaxset[n];
	DdNode* lminset[n];
	for(int j=0;j<n;j++){
		//create Ax
		DdNode* ax = Cudd_ReadZero(dd);
		Cudd_Ref(ax);
		for(int i=0;i<n;i++){
			DdNode* cons = Cudd_addConst(dd,ad[i+n*j]);
			Cudd_Ref(cons);
			DdNode* tmp = Cudd_addApply(dd,Cudd_addTimes,cons,xtreei[i]);
			Cudd_Ref(tmp);
			//Cudd_RecursiveDeref(dd,xtreei[i]);
			Cudd_RecursiveDeref(dd,cons);
			DdNode* tmp1 = Cudd_addApply(dd,Cudd_addPlus,tmp,ax);
			Cudd_Ref(tmp1);
			Cudd_RecursiveDeref(dd,ax);
			Cudd_RecursiveDeref(dd,tmp);
			ax = tmp1;
		}
		//create Bu
		DdNode* bu = Cudd_ReadZero(dd);
		Cudd_Ref(bu);
		for(int i=0;i<m;i++){
			DdNode* tmp = Cudd_addApply(dd,Cudd_addTimes,Cudd_addConst(dd,bd[i+m*j]),utreei[i]);
			Cudd_Ref(tmp);
			DdNode* tmp1 = Cudd_addApply(dd,Cudd_addPlus,tmp,bu);
			Cudd_Ref(tmp1);
			Cudd_RecursiveDeref(dd,bu);
			Cudd_RecursiveDeref(dd,tmp);
			bu = tmp1;
		}
		DdNode* axbu = Cudd_addApply(dd,Cudd_addPlus,ax,bu); 
		Cudd_Ref(axbu);
		Cudd_RecursiveDeref(dd,ax);
		Cudd_RecursiveDeref(dd,bu);
		DdNode* axbuset = Cudd_addApply(dd,Cudd_addTimes,xu,axbu);
		Cudd_Ref(axbuset);
		Cudd_RecursiveDeref(dd,axbu);
		double adsum = 0;
		double bdsum = 0;
		for(int i=0;i<n;i++)
			adsum += ad[i+n*j];
		for(int i=0;i<m;i++)
			bdsum += bd[i+m*j];
		// create minp + vmax
		double minpplus = minp[j] + vmax[j] - adsum -bdsum ;
		double minpminus = minp[j] - vmax[j] - adsum -bdsum ;
		//create maxrow 
		DdNode* minplusval = Cudd_addConst(dd,minpplus);
		Cudd_Ref(minplusval);
		DdNode* minminusval = Cudd_addConst(dd,minpminus);
		Cudd_Ref(minminusval);		
		DdNode* maxrow = Cudd_addApply(dd,Cudd_addPlus,axbuset,minplusval); 
		Cudd_Ref(maxrow);
		
		DdNode* minrow = Cudd_addApply(dd,Cudd_addPlus,axbuset,minminusval); 
		Cudd_Ref(minrow);
		Cudd_RecursiveDeref(dd,axbuset);
		DdNode* half = Cudd_addConst(dd,0.5);
		Cudd_Ref(half);
		
		DdNode* lmaxo = Cudd_addApply(dd,Cudd_addPlus,maxrow,half); 
		Cudd_Ref(lmaxo);
		
		DdNode* lmino = Cudd_addApply(dd,Cudd_addPlus,minrow,half); 
		Cudd_Ref(lmino);
		
		DdNode* lmaxrow = Cudd_addRoundOff(dd,lmaxo,0);
		Cudd_Ref(lmaxrow);
		DdNode* lminrow = Cudd_addRoundOff(dd,lmino,0);
		Cudd_Ref(lminrow);
		Cudd_RecursiveDeref(dd,maxrow);
		Cudd_RecursiveDeref(dd,minrow);
		DdNode* offsetB = Cudd_addConst(dd,xoffset[j]);
		Cudd_Ref(offsetB);
		DdNode* lmax = Cudd_addApply(dd,Cudd_addMinus,lmaxrow,offsetB);
		Cudd_Ref(lmax);
		DdNode* lmin = Cudd_addApply(dd,Cudd_addMinus,lminrow,offsetB);
		Cudd_Ref(lmin);
		Cudd_RecursiveDeref(dd,lmaxrow);
		Cudd_RecursiveDeref(dd,lminrow);
		lmaxset[j] = Cudd_addApply(dd,Cudd_addTimes,lmax,xu);
		Cudd_Ref(lmaxset[j]);
		lminset[j] = Cudd_addApply(dd,Cudd_addTimes,lmin,xu);
		Cudd_Ref(lminset[j]);
		Cudd_RecursiveDeref(dd,lmax);
		Cudd_RecursiveDeref(dd,lmin);
		DdNode* lminout = Cudd_addBddInterval(dd,lminset[j],-10000000,0); 
		Cudd_Ref(lminout);
		
		DdNode* lmaxout = Cudd_addBddInterval(dd,lmaxset[j],(2+(int)((xset[2*j+1]-xset[2*j])/eta)),10000000); 
		Cudd_Ref(lmaxout);
		
		DdNode* out =Cudd_bddOr(dd,lmaxout,lminout);
		Cudd_Ref(out);
		Cudd_RecursiveDeref(dd,lmaxout);
		Cudd_RecursiveDeref(dd,lminout);
		
		DdNode* outx = Cudd_bddAnd(dd,out,bxu);
		Cudd_Ref(outx);
		Cudd_RecursiveDeref(dd,out);
		
		DdNode* tmpout = Cudd_bddOr(dd,outx,trout);
		Cudd_Ref(tmpout);
		Cudd_RecursiveDeref(dd,trout);
		Cudd_RecursiveDeref(dd,outx);
		trout = tmpout;
	}
	
	DdNode* outTrans = Cudd_bddAnd(dd,trout,xoobBdd);
	Cudd_Ref(outTrans);
	Cudd_RecursiveDeref(dd,xoobBdd);
	
	DdNode* rest = Cudd_bddAnd(dd,bxu,Cudd_Not(trout));
	Cudd_Ref(rest);
	DdNode* inside = Cudd_ReadLogicZero(dd);
	Cudd_Ref(inside);
	if(rest != Cudd_ReadLogicZero(dd)){
		inside = Cudd_ReadOne(dd);
		Cudd_Ref(inside);
		DdNode* restadd = Cudd_BddToAdd(dd,rest);
		Cudd_Ref(restadd);
		
		for(int j=0;j<n;j++){
			DdNode* lmaxin =  Cudd_addApply(dd,Cudd_addTimes,lmaxset[j],restadd);
			Cudd_Ref(lmaxin); 
			Cudd_RecursiveDeref(dd,lmaxset[j]);
			
			DdNode* lminin =  Cudd_addApply(dd,Cudd_addTimes,lminset[j],restadd);
			Cudd_Ref(lminin);
			Cudd_RecursiveDeref(dd,lminset[j]);
			
			DdNode* nxlmaxin = Cudd_addApply(dd,Cudd_addMinus,lmaxin,nxset[j]);
			Cudd_Ref(nxlmaxin);
			Cudd_RecursiveDeref(dd,lmaxin);
			DdNode* nxlminin = Cudd_addApply(dd,Cudd_addMinus,nxset[j],lminin);
			Cudd_Ref(nxlminin);
			Cudd_RecursiveDeref(dd,nxset[j]);
			Cudd_RecursiveDeref(dd,lminin);
			DdNode* nxmaxb = Cudd_addBddInterval(dd,nxlmaxin,0,10000000); 
			Cudd_Ref(nxmaxb);
			Cudd_RecursiveDeref(dd,nxlmaxin);
			DdNode* nxminb = Cudd_addBddInterval(dd,nxlminin,0,10000000); 
			Cudd_Ref(nxminb);
			Cudd_RecursiveDeref(dd,nxlminin);
			DdNode* nxb =  Cudd_bddAnd(dd,nxmaxb,nxminb);
			Cudd_Ref(nxb);
			Cudd_RecursiveDeref(dd,nxmaxb);
			Cudd_RecursiveDeref(dd,nxminb);
			DdNode* nxb1 =  Cudd_bddAnd(dd,nxb,rest);
			Cudd_Ref(nxb1);
			Cudd_RecursiveDeref(dd,nxb);
			DdNode* tmpin = Cudd_bddAnd(dd,nxb1,inside);
			Cudd_Ref(tmpin);
			Cudd_RecursiveDeref(dd,inside);
			Cudd_RecursiveDeref(dd,nxb1);
			inside = tmpin;
		}
	}
	DdNode* trans = Cudd_bddOr(dd,inside,outTrans);
	Cudd_Ref(trans);
	Cudd_RecursiveDeref(dd,inside);
	Cudd_RecursiveDeref(dd,outTrans);
	return trans;		
}

DdNode* func_deter(DdManager* dd, int n,int m,double ad[],double bd[],double tau,double mu,double eta,
		  double vmax[],double minp[],double xset[],double uset[],int xoffset[],int uoffset[],int xoob[],int xbits[],int ubitsx[]){
	
	DdNode* one;		
	
	//Cudd_AutodynEnable(dd,CUDD_REORDER_SIFT);
	
	
	int totalxbits = 0;
	for(int i=0;i<n;i++)
		totalxbits += xbits[i];
	int ubits = 0 ; 
	for(int i=0;i<m;i++)
		ubits += ubitsx[i];
	
	
	int bits = 0;
	//create xtree on xvars  
	DdNode* bxtree = Cudd_ReadOne(dd);
	Cudd_Ref(bxtree);
	
	DdNode* xtreei[n];
	for(int i=0;i<n;i++){
		xtreei[i] = createADDTree(dd,bits,xbits[i]);
		Cudd_Ref(xtreei[i]);
		DdNode* tmp = Cudd_addBddInterval(dd,xtreei[i], (xset[2*i]/eta)+1-xoffset[i],(xset[2*i+1]/eta)+1-xoffset[i]); 
		Cudd_Ref(tmp);
		DdNode* tmp1 = Cudd_bddAnd(dd,tmp,bxtree);
		Cudd_Ref(tmp1);
		Cudd_RecursiveDeref(dd,bxtree);
		Cudd_RecursiveDeref(dd,tmp);
		bxtree = tmp1;
		bits += xbits[i];
	}
	//create utree on uvars
	DdNode* butree = Cudd_ReadOne(dd);
	Cudd_Ref(butree);
	DdNode* utreei[m];
	for(int i=0;i<m;i++){
		utreei[i] = createADDTree(dd,bits,ubitsx[i]);
		Cudd_Ref(utreei[i]);
		DdNode* tmp = Cudd_addBddInterval(dd,utreei[i],(uset[2*i]/mu)+1-uoffset[i],(uset[2*i+1]/mu)+1-uoffset[i]); 
		Cudd_Ref(tmp);
		DdNode* tmp1 = Cudd_bddAnd(dd,tmp,butree);
		Cudd_Ref(tmp1);
		Cudd_RecursiveDeref(dd,butree);
		Cudd_RecursiveDeref(dd,tmp);
		butree = tmp1;
		bits+= ubitsx[i];
	}
	
	
	DdNode* bxu = Cudd_bddAnd(dd,bxtree,butree);
	Cudd_Ref(bxu);
	Cudd_RecursiveDeref(dd,bxtree);
	Cudd_RecursiveDeref(dd,butree);
	DdNode* xu = Cudd_BddToAdd(dd,bxu);
	Cudd_Ref(xu);
	//create nxtree on nxvars
	
	DdNode* bnxtree = Cudd_ReadOne(dd);
	Cudd_Ref(bnxtree);
	DdNode* xoobBdd = Cudd_ReadOne(dd);
	Cudd_Ref(xoobBdd);
	DdNode* nxtreei[n];
	DdNode* nxset[n];
	for(int i=0;i<n;i++){
		nxtreei[i] = createADDTree(dd,bits,xbits[i]);
		Cudd_Ref(nxtreei[i]);
		DdNode* tmp = Cudd_addBddInterval(dd,nxtreei[i], (xset[2*i]/eta)+1-xoffset[i],(xset[2*i+1]/eta)+1-xoffset[i]); 
		Cudd_Ref(tmp);
		DdNode* nx  = Cudd_BddToAdd(dd,tmp);
		Cudd_Ref(nx);
		nxset[i] = Cudd_addApply(dd,Cudd_addTimes,nx,nxtreei[i]);
		
		Cudd_Ref(nxset[i]);
		
		DdNode* tmp1 = Cudd_bddAnd(dd,tmp,bnxtree);
		Cudd_Ref(tmp1);
		bnxtree = tmp1;
		DdNode* tmp2 = Cudd_addBddInterval(dd,nxtreei[i],xoob[i]+1,xoob[i]+1); 
		Cudd_Ref(tmp2);
		DdNode* tmp3 = Cudd_bddAnd(dd,tmp2,xoobBdd);
		Cudd_Ref(tmp3);
		Cudd_RecursiveDeref(dd,xoobBdd);
		xoobBdd = tmp3;
		bits += xbits[i];
	}
	
	DdNode* trout = Cudd_ReadLogicZero(dd);
	Cudd_Ref(trout);
	DdNode* lmaxset[n];
	DdNode* lminset[n];
	
	for(int j=0;j<n;j++){
		//create Ax
		
		DdNode* ax = Cudd_ReadZero(dd);
		Cudd_Ref(ax);
		for(int i=0;i<n;i++){
			DdNode* cons = Cudd_addConst(dd,ad[i+n*j]);
			Cudd_Ref(cons);
			DdNode* tmp = Cudd_addApply(dd,Cudd_addTimes,cons,xtreei[i]);
			Cudd_Ref(tmp);
			//Cudd_RecursiveDeref(dd,xtreei[i]);
			Cudd_RecursiveDeref(dd,cons);
			DdNode* tmp1 = Cudd_addApply(dd,Cudd_addPlus,tmp,ax);
			Cudd_Ref(tmp1);
			Cudd_RecursiveDeref(dd,ax);
			Cudd_RecursiveDeref(dd,tmp);
			ax = tmp1;
		}
		//create Bu
		DdNode* bu = Cudd_ReadZero(dd);
		Cudd_Ref(bu);
		for(int i=0;i<m;i++){
			DdNode* tmp = Cudd_addApply(dd,Cudd_addTimes,Cudd_addConst(dd,bd[i+m*j]),utreei[i]);
			Cudd_Ref(tmp);
			DdNode* tmp1 = Cudd_addApply(dd,Cudd_addPlus,tmp,bu);
			Cudd_Ref(tmp1);
			Cudd_RecursiveDeref(dd,bu);
			Cudd_RecursiveDeref(dd,tmp);
			bu = tmp1;
		}
		//create Ax+Bu
		DdNode* axbu = Cudd_addApply(dd,Cudd_addPlus,ax,bu); 
		Cudd_Ref(axbu);
		Cudd_RecursiveDeref(dd,ax);
		Cudd_RecursiveDeref(dd,bu);
		
		DdNode* axbuset = Cudd_addApply(dd,Cudd_addTimes,xu,axbu);
		Cudd_Ref(axbuset);
		Cudd_RecursiveDeref(dd,axbu);
		
		double adsum = 0;
		double bdsum = 0;
		for(int i=0;i<n;i++)
			adsum += ad[i+n*j];
		for(int i=0;i<m;i++)
			bdsum += bd[i+m*j];
		// create minp + vmax
		double minpplus = minp[j] - adsum -bdsum;// + 1 -0.5;
		
		//create maxrow 
		DdNode* maxrow = Cudd_addApply(dd,Cudd_addPlus,axbuset,Cudd_addConst(dd,minpplus)); 
		Cudd_Ref(maxrow);
		
		Cudd_RecursiveDeref(dd,axbuset);
		DdNode* half = Cudd_addConst(dd,0.5);
		Cudd_Ref(half);
		
		DdNode* lmaxo = Cudd_addApply(dd,Cudd_addPlus,maxrow,half); 
		Cudd_Ref(lmaxo);
		Cudd_RecursiveDeref(dd,maxrow);
		
		DdNode* lmaxrow = Cudd_addRoundOff(dd,lmaxo,0);
		Cudd_Ref(lmaxrow);
		
		DdNode* offsetB = Cudd_addConst(dd,xoffset[j]);
		Cudd_Ref(offsetB);
		DdNode* lmax = Cudd_addApply(dd,Cudd_addMinus,lmaxrow,offsetB);
		Cudd_Ref(lmax);
		
		Cudd_RecursiveDeref(dd,lmaxrow);
		
		lmaxset[j] = Cudd_addApply(dd,Cudd_addTimes,lmax,xu);
		Cudd_Ref(lmaxset[j]);
		
		Cudd_RecursiveDeref(dd,lmax);
		DdNode* lminout = Cudd_addBddInterval(dd,lmaxset[j],-100000000,0); 
		Cudd_Ref(lminout);
		
		DdNode* lmaxout = Cudd_addBddInterval(dd,lmaxset[j],2+(int)((xset[2*j+1]-xset[2*j])/eta),100000000); 
		Cudd_Ref(lmaxout);
		
		DdNode* out =Cudd_bddOr(dd,lmaxout,lminout);
		Cudd_Ref(out);
		Cudd_RecursiveDeref(dd,lmaxout);
		Cudd_RecursiveDeref(dd,lminout);
		
		DdNode* outx = Cudd_bddAnd(dd,out,bxu);
		Cudd_Ref(outx);
		Cudd_RecursiveDeref(dd,out);
		
		DdNode* tmpout = Cudd_bddOr(dd,outx,trout);
		Cudd_Ref(tmpout);
		Cudd_RecursiveDeref(dd,trout);
		Cudd_RecursiveDeref(dd,outx);
		trout = tmpout;
	}
	
	DdNode* outTrans = Cudd_bddAnd(dd,trout,xoobBdd);
	Cudd_Ref(outTrans);
	Cudd_RecursiveDeref(dd,xoobBdd);
	
	DdNode* rest = Cudd_bddAnd(dd,bxu,Cudd_Not(trout));
	Cudd_Ref(rest);
	DdNode* inside = Cudd_ReadLogicZero(dd);
	Cudd_Ref(inside);
	if(rest != Cudd_ReadLogicZero(dd)){
		inside = Cudd_ReadOne(dd);
		Cudd_Ref(inside);
		DdNode* restadd = Cudd_BddToAdd(dd,rest);
		Cudd_Ref(restadd);
		
		for(int j=0;j<n;j++){
			
			DdNode* lmaxin =  Cudd_addApply(dd,Cudd_addTimes,lmaxset[j],restadd);
			Cudd_Ref(lmaxin); 
			Cudd_RecursiveDeref(dd,lmaxset[j]);
			
						
			DdNode* nxlmaxin = Cudd_addApply(dd,Cudd_addMinus,lmaxin,nxset[j]);
			Cudd_Ref(nxlmaxin);
			Cudd_RecursiveDeref(dd,lmaxin);
			
			Cudd_RecursiveDeref(dd,nxset[j]);
			DdNode* nxmaxb = Cudd_addBddInterval(dd,nxlmaxin,0,0); 
			Cudd_Ref(nxmaxb);
			Cudd_RecursiveDeref(dd,nxlmaxin);
			
			DdNode* nxb1 =  Cudd_bddAnd(dd,nxmaxb,rest);
			Cudd_Ref(nxb1);
			Cudd_RecursiveDeref(dd,nxmaxb);
			DdNode* tmpin = Cudd_bddAnd(dd,nxb1,inside);
			Cudd_Ref(tmpin);
			Cudd_RecursiveDeref(dd,inside);
			Cudd_RecursiveDeref(dd,nxb1);
			inside = tmpin;
		}
	}
	DdNode* trans = Cudd_bddOr(dd,inside,outTrans);
	Cudd_Ref(trans);
	Cudd_RecursiveDeref(dd,inside);
	Cudd_RecursiveDeref(dd,outTrans);
	return trans;		
}
/************************************

 MEX-MAIN function "mexFunction"
 Initializes data (params_symb structure and pointer to params_symb)

*************************************/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{

	s_vector params_symb;
	int ret,ok;
	int verbose;
        mxArray *psv;
	double total_time;
	double ninputs,nstates;
	long nbatch,i,j,r,s;
	double totloops;
	long memuse1, memuse2, mem_res;
	int plot_ok;

	char *bddT;
	char suffix[5];
	char *bddT_name;

	//BDD Variables	
	short numVars;  // num of vars; if unknown set to 0
	short numVarsZ; // num of vars for ZBDDs; if unknown set to 0
	int numSlots; // default for CUDD package
	int cacheSize; // default for CUDD package
	//int maxCacheSize;   // default for CUDD package

	mexEvalString("t_start = tic;");
	//Copy data from Matlab workspace 
	psv=mexGetVariable("caller","params_symb");
	//Copy data to variables
	nbatch=(long)mxGetScalar(mexGetVariable("caller","nbatch"));
	totloops=mxGetScalar(mexGetVariable("caller","totloops"));
	params_symb.n=(int)mxGetScalar(mxGetField(psv,0,"n"));
	params_symb.m=(int)mxGetScalar(mxGetField(psv,0,"m"));
	params_symb.nume=(double *)mxGetPr(mxGetField(psv,0,"nume"));   
	params_symb.totbits=(int)mxGetScalar(mxGetField(psv,0,"totbits"));
	params_symb.nbitsloop=(int)mxGetScalar(mxGetField(psv,0,"nbitsloop"));
	params_symb.nbits=(double *)mxGetPr(mxGetField(psv,0,"nbits"));
	params_symb.deter=(int)mxGetScalar(mxGetField(psv,0,"deter"));
	params_symb.nbitsx=(int)mxGetScalar(mxGetField(psv,0,"nbitsx"));	
	
	//BDD INITIALIZATIONS		
	numVars=params_symb.totbits;  // num of vars; if unknown set to 0
	numVarsZ=0; // num of vars for ZBDDs; if unknown set to 0
	numSlots=CUDD_UNIQUE_SLOTS; // default for CUDD package
	cacheSize=CUDD_CACHE_SLOTS; // default for CUDD package
	long maxCacheSize=0; //10485760*2;   // default for CUDD package

	ddman = Cudd_Init(numVars, numVarsZ, numSlots, cacheSize, maxCacheSize); //maxCacheSize);
    
	//Each DdNode* needs to be properly init. Watch it for multi BDDs.
    nstates=1;
	for (i=0;i<params_symb.n;i++)
		nstates*=(params_symb.nume[i]+1);
	
	ninputs=1;
	for (i=params_symb.n;i<params_symb.n+params_symb.m;i++)
		ninputs*=(params_symb.nume[i]+1);
	
	// Print some extra-info
	mexPrintf("\nSymbolic model size: ");
	mexPrintf("%.0f", nstates);
	mexPrintf(" states; \n");
	mexPrintf("                     %.0f", ninputs);
	mexPrintf(" inputs. \n");
	
	int n = params_symb.n;
	int m = params_symb.m;
	double ad[n*n];
	double bd[n*m];
	double tau;
	double mu;
	double eta;
	double vmax[n];
	double minp[n];
	double xset[2*n];
	double uset[2*m];
	int xoffset[n];
	int uoffset[m];
	int xoob[n];
	int xbits[n];
	int ubits[m];
	
	const mxArray* Ad = prhs[2];
	double* adv = mxGetPr(Ad);
	int k=0;
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			ad[k]= adv[i+j*n];
			k++;
		}
	}

	const mxArray* Bd = prhs[3];
	double* bdv = mxGetPr(Bd);
	k=0;
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			bd[k]= bdv[i+j*m];
			k++;
		}
	}	
	const mxArray* minpv = prhs[4];
	double* min = mxGetPr(minpv);
	k=0;
	for(int i=0;i<n;i++){
		minp[k]= min[k];
		k++;
	}	
	const mxArray* vmaxv = prhs[5];
	double* max = mxGetPr(vmaxv);
	k=0;
	for(int i=0;i<n;i++){
		vmax[k]= max[k];
		k++;
	}	
	
	const mxArray* xo = prhs[6];
	double* xoff = mxGetPr(xo);
	k=0;
	for(int i=0;i<n;i++){
		xoffset[k]= (int)xoff[k];
		k++;
	}	
	
	const mxArray* uo = prhs[7];
	double* uoff = mxGetPr(uo);
	k=0;
	for(int i=0;i<m;i++){
		uoffset[k]= (int)uoff[k];
		k++;
	}	
	
	const mxArray* xoobv = prhs[8];
	double* xoobval = mxGetPr(xoobv);
	k=0;
	for(int i=0;i<n;i++){
		xoob[k]= (int)xoobval[k];
		k++;
	}	
	
	const mxArray* xsets = prhs[9];
	double* xsetval = mxGetPr(xsets);
	k=0;
	for(int i=0;i<n;i++){
		xset[k]= xsetval[i];
		k++;
		xset[k]= xsetval[i+n];
		k++;
	}	
	
	const mxArray* usets = prhs[10];
	double* usetval = mxGetPr(usets);
	k=0;
	for(int i=0;i<m;i++){
		uset[k]= usetval[i];
		k++;
		uset[k]= usetval[i+m];
		k++;
	}	
	
	for(int i=0;i<n;i++){
		xbits[i] = (int)params_symb.nbits[i];
	}
	for(int i=0;i<m;i++){
		ubits[i] = (int)params_symb.nbits[i+n];
	}	
	// Build model abstraction
	DdNode* trans;
	//= Cudd_ReadLogicZero(ddman);
	//Cudd_Ref(trans);
	
	tau = mxGetScalar(prhs[11]);
	mu = mxGetScalar(prhs[12]);
	eta = mxGetScalar(prhs[13]);
	
	if(params_symb.deter == 1) 
		trans = func_deter(ddman,n,m,ad,bd,tau,mu,eta,vmax,minp,xset,uset,xoffset,uoffset,xoob,xbits,ubits);
	else 
		trans = func_nondeter(ddman,n,m,ad,bd,tau,mu,eta,vmax,minp,xset,uset,xoffset,uoffset,xoob,xbits,ubits);
	
	Cudd_Ref(trans); // referenced twice in a row?
	// Save the abstract model
	bddT=mxArrayToString(prhs[0]);

	strcpy(suffix,".bdd");

	bddT_name=(char*)mxMalloc(strlen(bddT)+5);
	strcpy(bddT_name, bddT);
	strcat(bddT_name,suffix);
	mxFree(bddT);

	ok = Dddmp_cuddBddStore(ddman, NULL, trans, NULL,
					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
					bddT_name, NULL);
	if(ok)
		mexPrintf("Symbolic model successfully saved to '.bdd' file.\n");
	else
		mexPrintf("Symbolic model failed to be saved.\n");
			

	Cudd_Quit(ddman); //Maybe better to deref trans first?

	mexEvalString("t_end = toc(t_start);");
	total_time = mxGetScalar(mexGetVariable("caller","t_end"));
	if(total_time < 0.001)
		mexPrintf("Elapsed time: 0.001 seconds. \n");
	else		
	{	
		mexPrintf("Elapsed time: %.3f", total_time);
		mexPrintf(" seconds. \n", total_time);
	}

	mexPrintf("\n---------------------- Pessoa: Abstraction Terminated ------------------ \n");
}	
