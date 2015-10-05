#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>

#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP))


//this eliminates the apply() for difference of deviances. 
//we pass RX, which contains prepended column of 1's if appropriate
//and also the last column is the alternative model, and is ALSO
//DOUBLED. THIS IS VERY IMPORTANT. 
//SERIOUSLY REMEMBER THE ABOVE LINE
//RY consists of the genotyping data, which we will process row by row
//returns difference of deviance for the two rows
SEXP assoc(SEXP RX, SEXP RY, SEXP Rmi, SEXP Rtol){
    int *dimX, *dimY, maxiter = (int)(*REAL(Rmi));
    double *X, *Y, tol = (double)(*REAL(Rtol));
    int i, j, k, n, ind1, ind2, ind;
    int null;
    double devnull, devalt;
    int flag;

    dimX = getDims(RX);
    dimY = getDims(RY);
    PROTECT(RY = coerceVector(RY, REALSXP));
    PROTECT(RX = coerceVector(RX, REALSXP));
    Y = REAL(RY);
    X = REAL(RX);

    SEXP Rret; //returns difference of deviances
    double *ret;
    PROTECT(Rret = allocVector(REALSXP, dimY[0]));
    ret = REAL(Rret);

    //performs some checks...might need others?
    if(dimX[0] != (dimY[1]*2))
        Rprintf("dimension mismatch? crap\n");
    if(dimX[1] < 2)
        Rprintf("expecting constant term + another covariate?\n");

    //here we allocate temp memory
    int numblock = dimX[1]*64; //approximate? i assume ilnaev will return 64?
    double *b  = (double*) malloc(sizeof(double)*dimX[1]); //beta
    double *bl = (double*) malloc(sizeof(double)*dimX[1]); //beta last
    double *f  = (double*) malloc(sizeof(double)*dimX[1]); //tmp
    double *p  = (double*) malloc(sizeof(double)*dimX[0]); //mle
    double *y  = (double*) malloc(sizeof(double)*dimX[0]); //the doubled y

    double *w  = (double*) malloc(sizeof(double)*dimX[1]*dimX[1]);
    int* ipiv  = (int*)    malloc(sizeof(int)   *dimX[1]); //for inverting
    double *wo = (double*) malloc(sizeof(double)*numblock);
    double max; // check convergence

    int iter;
    double alpha = -1.0, zero = 0.0, one = 1.0;
    int ione = 1;
    int info=0;
    char tr = 'n';
    double tmp;
    

//scratch: Y: SNPs by samples. X: samples by sv(+trait)
//note: it's prudent evidentally to reset b for every set of null/alt
//fits.
    null = dimX[1] -1;
    //loop over all rows
    for(n = 0; n < dimY[0]; n++){
        //n marks beginning for Y
        //make the extended y
        ind = n;
        for(i = 0; i < dimY[1]; i++){
            if(Y[ind] == 0.0){
                y[i] = 0.0;
                y[i + dimY[1]] = 0.0;
            }
            else if(Y[ind] == 1.0){
                y[i] = 1.0;
                y[i + dimY[1]] = 0.0;
            }
            else if(Y[ind] == 2.0){
                y[i] = 1.0;
                y[i + dimY[1]] = 1.0;
            }
            ind += dimY[0];
        }

        //logistic regression for null model
        //IRLS
        iter = 1;
        for(i = 0; i < dimX[1]; i++) {
            b[i]  = 0;
            bl[i] = 0;
        }
        
        flag = 0;
        while(iter <= maxiter){
            ///////////////////////////////////////////////////////////////////////
            //p <- as.vector(1/(1 + exp(-X %*% b)))
            F77_CALL(dgemv)(&tr,dimX,&null,&alpha,X,dimX,b,&ione,&zero,p,&ione);
            for(i = 0; i < dimX[0]; i++) 
                p[i] = 1/(1+exp(p[i]));
            
            ///////////////////////////////////////////////////////////////////////
            //var.b <- solve(crossprod(X, p * (1 - p) * X))
            //
            //here, solve is inverting the matrix. 
            //p*(1-p) is applied to cols of X.
            //at the moment I am manually computing the crossprod
            //which is guaranteed to be symmetric
            for(i = 0; i < null; i++){ //rows
                for(j = i; j < null; j++){ //columns
                    ind1 = i*dimX[0]; //i-th col of X
                    ind2 = j*dimX[0]; //j-th col of X
                    ind  = null*i + j; //position on w
                    w[ind] = 0;
                    for(k = 0; k < dimX[0]; k++){ //loop over X'p(1-p)X
                        w[ind]+=X[ind1]*X[ind2]*p[k]*(1-p[k]);
                        ind1++;
                        ind2++;
                    }
                    if(i != j) //reflect it
                        w[null*j+i] = w[ind];
                }
            }


            //actually inverting here. remember to pay attention to includes
            F77_CALL(dgetrf)(&null,&null,w,&null,ipiv,&info);
            if(info != 0) {
                Rprintf("warning: dgetrf error (null model), NA used\n");
                Rprintf("n:%i info:%i iter:%i\n", n, info, iter);
                //error("dgetrf error (null model)\n");
                flag = 1;
                break;
            }
            F77_CALL(dgetri)(&null,w,&null,ipiv,wo,&numblock,&info);
            if(info != 0) {
                Rprintf("warning: dgetri error (null model), NA used\n");
                Rprintf("n:%i info:%i iter:%i\n", n, info, iter);
                //error("dgetri error (null model)\n");
                flag = 1;
                break;
            }    
    
            ///////////////////////////////////////////////////////////////////////
            //b <- b + var.b %*% crossprod(X, y - p)
            //use f to calculate crossprod(X,y-p) first.
            //then use dgemv
            ind  = 0; //since we are iterating over X in order
            for(i = 0; i < null; i++){ //cols of X, values of f
                f[i] = 0;
                for(j = 0; j < dimX[0]; j++){ //rows of X, values of y-p
                    f[i] += X[ind] * (y[j] - p[j]);
                    ind++;
                }
            }
    
            F77_CALL(dgemv)(&tr,&null,&null,&one,w,&null,f,&ione,&one,b,&ione);
            
            
            ///////////////////////////////////////////////////////////////////////
            //if (max(abs(b - b.last)/(abs(b.last) + 0.01*tol)) < tol) break
            //check to see if we need to break
            max = 0.0;
            for(i = 0; i < null; i++) {
                tmp = fabs(b[i] - bl[i])/(fabs(bl[i]) + 0.01*tol);
                if(tmp > max) max = tmp;
            }
    
            if(max < tol)
                break;
    
            
            ///////////////////////////////////////////////////////////////////////
            //b.last <- b
            //it <- it + 1
            for(i = 0; i < dimX[1]; i++) bl[i] = b[i];
            iter++; 
        }
        devnull = 0.0;
        if(flag == 0){
            for(i = 0; i < dimX[0]; i++) 
                devnull += (y[i]*log(p[i])) + ((1-y[i])*log(1-p[i]));
        }


        //logistic regression for alternative model
        //IRLS
        iter = 1;
        while(iter <= maxiter){
            ///////////////////////////////////////////////////////////////////////
            //p <- as.vector(1/(1 + exp(-X %*% b)))
            F77_CALL(dgemv)(&tr,dimX,dimX+1,&alpha,X,dimX,b,&ione,&zero,p,&ione);
            for(i = 0; i < dimX[0]; i++) 
                p[i] = 1/(1+exp(p[i]));
            
            ///////////////////////////////////////////////////////////////////////
            //var.b <- solve(crossprod(X, p * (1 - p) * X))
            //
            //here, solve is inverting the matrix. 
            //p*(1-p) is applied to cols of X.
            //at the moment I am manually computing the crossprod
            //which is guaranteed to be symmetric
            for(i = 0; i < dimX[1]; i++){ //rows
                for(j = i; j < dimX[1]; j++){ //columns
                    ind1 = i*dimX[0]; //i-th col of X
                    ind2 = j*dimX[0]; //j-th col of X
                    ind  = dimX[1]*i + j; //position on w
                    w[ind] = 0;
                    for(k = 0; k < dimX[0]; k++){ //loop over X'p(1-p)X
                        w[ind]+=X[ind1]*X[ind2]*p[k]*(1-p[k]);
                        ind1++;
                        ind2++;
                    }
                    if(i != j) //reflect it
                        w[dimX[1]*j+i] = w[ind];
                }
            }
            
            //actually inverting here. remember to pay attention to includes
            F77_CALL(dgetrf)(dimX+1,dimX+1,w,dimX+1,ipiv,&info);
            if(info != 0) {
                Rprintf("warning: dgetrf error (alt model), NA used\n");
                Rprintf("n:%i info:%i iter:%i\n", n, info, iter);
                //error("dgetrf error (alt model)\n");
                flag = 1;
                break;
            }
            F77_CALL(dgetri)(dimX+1,w,dimX+1,ipiv,wo,&numblock,&info);
            if(info != 0) {
                Rprintf("warning: dgetri error (alt model), NA used\n");
                Rprintf("n:%i info:%i iter:%i\n", n, info, iter);
                //error("dgetri error (alt model)\n");
                flag = 1;
                break;
            }
    
    
            ///////////////////////////////////////////////////////////////////////
            //b <- b + var.b %*% crossprod(X, y - p)
            //use f to calculate crossprod(X,y-p) first.
            //then use dgemv
            ind  = 0; //since we are iterating over X in order
            for(i = 0; i < dimX[1]; i++){ //cols of X, values of f
                f[i] = 0;
                for(j = 0; j < dimX[0]; j++){ //rows of X, values of y-p
                    f[i] += X[ind] * (y[j] - p[j]); //traversing Y
                    ind++; 
                }
            }
    
            F77_CALL(dgemv)(&tr,dimX+1,dimX+1,&one,w,dimX+1,f,&ione,&one,b,&ione);
            
            
            ///////////////////////////////////////////////////////////////////////
            //if (max(abs(b - b.last)/(abs(b.last) + 0.01*tol)) < tol) break
            //check to see if we need to break
            max = 0.0;
            for(i = 0; i < dimX[1]; i++) {
                tmp = fabs(b[i] - bl[i])/(fabs(bl[i]) + 0.01*tol);
                if(tmp > max) max = tmp;
            }
    
            if(max < tol)
                break;

            ///////////////////////////////////////////////////////////////////////
            //b.last <- b
            //it <- it + 1
            for(i = 0; i < dimX[1]; i++) bl[i] = b[i];
            iter++;
        }
        if(flag == 0){
            devalt = 0.0; //missing the -2
            for(i = 0; i < dimX[0]; i++) 
                devalt += (y[i]*log(p[i])) + ((1-y[i])*log(1-p[i]));
            ret[n] = -2*(devnull - devalt);
        }
        else{
            ret[n] =  R_NaReal;
        }
    }

    UNPROTECT(3);
    free(b);
    free(bl);
    free(f);
    free(p);
    free(y);
    free(w);
    free(ipiv);
    free(wo);
    return Rret;
}

