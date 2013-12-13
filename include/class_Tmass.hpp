////////////////////////////////////////////////////
///////////////                      ///////////////
///////////////     Tmassmatrix      ///////////////
///////////////                      ///////////////
////////////////////////////////////////////////////

template <typename MASS_MATRIX_TYPE, 
            typename MASK_MATRIX_TYPE, 
            typename STATE_MATRIX_TYPE, 
            typename STATE_MATRIX_CV_TYPE, 
	    typename LEAFCELLPOINTER_TYPE, 
	    typename CHILDREN_LEAFCELLPOINTER_TYPE, 
	    int SHAPEDIM, int STATEDIM, int CHILDDIM>
void Tmassmatrix<MASS_MATRIX_TYPE, MASK_MATRIX_TYPE, 
                 STATE_MATRIX_TYPE, STATE_MATRIX_CV_TYPE, 
                 LEAFCELLPOINTER_TYPE, CHILDREN_LEAFCELLPOINTER_TYPE, 
		 SHAPEDIM, STATEDIM, CHILDDIM>
::set_massmatrix (const Tmassmatrix & m) 
{
   int i,j,k;
   for (k=0;k<SHAPEDIM;k++) 
     for (i=0;i<SHAPEDIM;i++) { 
       A[k][i]    = m.A[k][i];
       for (j=0; j<CHILDDIM; j++) B[j][k][i]    = m.B[j][k][i];
       }
};

//////////////////////////////////////////////////////

template <typename MASS_MATRIX_TYPE, 
            typename MASK_MATRIX_TYPE, 
            typename STATE_MATRIX_TYPE,
	    typename STATE_MATRIX_CV_TYPE, 
	    typename LEAFCELLPOINTER_TYPE, 
	    typename CHILDREN_LEAFCELLPOINTER_TYPE, 
	    int SHAPEDIM, int STATEDIM, int CHILDDIM>
void Tmassmatrix<MASS_MATRIX_TYPE, MASK_MATRIX_TYPE, 
                 STATE_MATRIX_TYPE, STATE_MATRIX_CV_TYPE,        
                 LEAFCELLPOINTER_TYPE, CHILDREN_LEAFCELLPOINTER_TYPE, 
		 SHAPEDIM, STATEDIM, CHILDDIM>
::Cholesky_decomp () 
{
   int i,j,k;
   for (k=0;k<SHAPEDIM;k++) {
     for (j=0;j<k;j++) A[k][k] -= A[k][j]*A[k][j]*A[j][j];
     for (i=k+1;i<SHAPEDIM;i++) {
       for (j=0;j<k;j++) A[i][k] -= A[i][j]*A[j][j]*A[k][j];
       A[i][k]    /= A[k][k];
       }
     }
};

//////////////////////////////////////////////////////

template <typename MASS_MATRIX_TYPE, 
            typename MASK_MATRIX_TYPE, 
            typename STATE_MATRIX_TYPE, 
	    typename STATE_MATRIX_CV_TYPE, 
	    typename LEAFCELLPOINTER_TYPE, 
	    typename CHILDREN_LEAFCELLPOINTER_TYPE, 
	    int SHAPEDIM, int STATEDIM, int CHILDDIM>
void Tmassmatrix<MASS_MATRIX_TYPE, MASK_MATRIX_TYPE, 
                 STATE_MATRIX_TYPE, STATE_MATRIX_CV_TYPE, 
                 LEAFCELLPOINTER_TYPE, CHILDREN_LEAFCELLPOINTER_TYPE, 
		 SHAPEDIM, STATEDIM, CHILDDIM>
::Cholesky_solve (STATE_MATRIX_TYPE & w) 
{
   int i,j,k;
   for (i=0; i<STATEDIM; i++) {
     for (k=1;k<SHAPEDIM;k++) {
       for (j=0;j<k;j++) w[k][i] -= A[k][j]*w[j][i];
       }
     for (k=0;k<SHAPEDIM;k++) w[k][i] /= A[k][k];
     for (k=SHAPEDIM-2;k>-1;k--) {
       for (j=k+1;j<SHAPEDIM;j++) w[k][i] -= A[j][k]*w[j][i];
     }
   }
};


//////////////////////////////////////////////////////

template <typename MASS_MATRIX_TYPE, 
            typename MASK_MATRIX_TYPE, 
            typename STATE_MATRIX_TYPE, 
	    typename STATE_MATRIX_CV_TYPE, 
	    typename LEAFCELLPOINTER_TYPE, 
	    typename CHILDREN_LEAFCELLPOINTER_TYPE, 
	    int SHAPEDIM, int STATEDIM, int CHILDDIM>
void Tmassmatrix<MASS_MATRIX_TYPE, MASK_MATRIX_TYPE, 
                 STATE_MATRIX_TYPE, STATE_MATRIX_CV_TYPE, 
                 LEAFCELLPOINTER_TYPE, CHILDREN_LEAFCELLPOINTER_TYPE, 
		 SHAPEDIM, STATEDIM, CHILDDIM>
::Cholesky_subsolve (const int & subshapedim, STATE_MATRIX_TYPE & w) 
{
   int i,j,k;
   for (i=0; i<STATEDIM; i++) {
     for (k=1;k<subshapedim;k++) {
       for (j=0;j<k;j++) w[k][i] -= A[k][j]*w[j][i];
       }
     for (k=0;k<subshapedim;k++) w[k][i] /= A[k][k];
     for (k=subshapedim-2;k>-1;k--) {
       for (j=k+1;j<subshapedim;j++) w[k][i] -= A[j][k]*w[j][i];
     }
     
     for (k=subshapedim;k<SHAPEDIM;k++) w[k][i] = 0.0;
   }
};

//////////////////////////////////////////////////////

template <typename MASS_MATRIX_TYPE, 
            typename MASK_MATRIX_TYPE, 
            typename STATE_MATRIX_TYPE, 
	    typename STATE_MATRIX_CV_TYPE, 
	    typename LEAFCELLPOINTER_TYPE, 
	    typename CHILDREN_LEAFCELLPOINTER_TYPE, 
	    int SHAPEDIM, int STATEDIM, int CHILDDIM>
void Tmassmatrix<MASS_MATRIX_TYPE, MASK_MATRIX_TYPE, 
                 STATE_MATRIX_TYPE, STATE_MATRIX_CV_TYPE, 
                 LEAFCELLPOINTER_TYPE, CHILDREN_LEAFCELLPOINTER_TYPE, 
		 SHAPEDIM, STATEDIM, CHILDDIM>
::Cholesky_multiply (STATE_MATRIX_TYPE & w) 
{
   int i,j,k;
   for (i=0; i<STATEDIM; i++) {
     for (k=0;k<SHAPEDIM;k++) {
       for (j=k+1;j<SHAPEDIM;j++) w[k][i] += A[j][k]*w[j][i];
       w[k][i] *= A[k][k];
       }
     for (k=SHAPEDIM-1;k>-1;k--) {
       for (j=0;j<k;j++) w[k][i] += A[k][j]*w[j][i];
       }
     }
};

//////////////////////////////////////////////////////

template <typename MASS_MATRIX_TYPE, 
            typename MASK_MATRIX_TYPE, 
            typename STATE_MATRIX_TYPE, 
	    typename STATE_MATRIX_CV_TYPE, 
	    typename LEAFCELLPOINTER_TYPE, 
	    typename CHILDREN_LEAFCELLPOINTER_TYPE, 
	    int SHAPEDIM, int STATEDIM, int CHILDDIM>
void Tmassmatrix<MASS_MATRIX_TYPE, MASK_MATRIX_TYPE, 
                 STATE_MATRIX_TYPE, STATE_MATRIX_CV_TYPE, 
                 LEAFCELLPOINTER_TYPE, CHILDREN_LEAFCELLPOINTER_TYPE, 
		 SHAPEDIM, STATEDIM, CHILDDIM>
::Cholesky_submultiply (const int & subshapedim, STATE_MATRIX_TYPE & w) 
{
   int i,j,k;
   for (i=0; i<STATEDIM; i++) {
     for (k=0;k<subshapedim;k++) {
       for (j=k+1;j<subshapedim;j++) w[k][i] += A[j][k]*w[j][i];
       w[k][i] *= A[k][k];
       }
     for (k=subshapedim-1;k>-1;k--) {
       for (j=0;j<k;j++) w[k][i] += A[k][j]*w[j][i];
       }
     }
};

//////////////////////////////////////////////////////

template <typename MASS_MATRIX_TYPE, 
            typename MASK_MATRIX_TYPE, 
            typename STATE_MATRIX_TYPE, 
	    typename STATE_MATRIX_CV_TYPE, 
	    typename LEAFCELLPOINTER_TYPE, 
	    typename CHILDREN_LEAFCELLPOINTER_TYPE, 
	    int SHAPEDIM, int STATEDIM, int CHILDDIM>
void Tmassmatrix<MASS_MATRIX_TYPE, MASK_MATRIX_TYPE, 
                 STATE_MATRIX_TYPE, STATE_MATRIX_CV_TYPE, 
                 LEAFCELLPOINTER_TYPE, CHILDREN_LEAFCELLPOINTER_TYPE, 
		 SHAPEDIM, STATEDIM, CHILDDIM>
::Cholesky_multiply (const STATE_MATRIX_TYPE & v, STATE_MATRIX_TYPE & w) 
{
   int i,j,k;
   for (i=0; i<STATEDIM; i++) {
     for (k=0;k<SHAPEDIM;k++) w[k][i] = v[k][i];
     for (k=0;k<SHAPEDIM;k++) {
       for (j=k+1;j<SHAPEDIM;j++) w[k][i] += A[j][k]*w[j][i];
       w[k][i] *= A[k][k];
       }
     for (k=SHAPEDIM-1;k>-1;k--) {
       for (j=0;j<k;j++) w[k][i] += A[k][j]*w[j][i];
       }
     }
};

//////////////////////////////////////////////////////

template <typename MASS_MATRIX_TYPE, 
            typename MASK_MATRIX_TYPE, 
            typename STATE_MATRIX_TYPE, 
	    typename STATE_MATRIX_CV_TYPE, 
	    typename LEAFCELLPOINTER_TYPE, 
	    typename CHILDREN_LEAFCELLPOINTER_TYPE, 
	    int SHAPEDIM, int STATEDIM, int CHILDDIM>
void Tmassmatrix<MASS_MATRIX_TYPE, MASK_MATRIX_TYPE, 
                 STATE_MATRIX_TYPE, STATE_MATRIX_CV_TYPE, 
                 LEAFCELLPOINTER_TYPE, CHILDREN_LEAFCELLPOINTER_TYPE, 
		 SHAPEDIM, STATEDIM, CHILDDIM>
::Cholesky_submultiply (const int & subshapedim, const STATE_MATRIX_TYPE & v, 
                       STATE_MATRIX_TYPE & w) 
{
   int i,j,k;
   for (i=0; i<STATEDIM; i++) {
     for (k=0;k<subshapedim;k++) w[k][i] = v[k][i];
     for (k=0;k<subshapedim;k++) {
       for (j=k+1;j<subshapedim;j++) w[k][i] += A[j][k]*w[j][i];
       w[k][i] *= A[k][k];
       }
     for (k=subshapedim-1;k>-1;k--) {
       for (j=0;j<k;j++) w[k][i] += A[k][j]*w[j][i];
       }
     }
};

//////////////////////////////////////////////////////

template <typename MASS_MATRIX_TYPE, 
            typename MASK_MATRIX_TYPE, 
            typename STATE_MATRIX_TYPE, 
	    typename STATE_MATRIX_CV_TYPE, 
	    typename LEAFCELLPOINTER_TYPE, 
	    typename CHILDREN_LEAFCELLPOINTER_TYPE, 
	    int SHAPEDIM, int STATEDIM, int CHILDDIM>
void Tmassmatrix<MASS_MATRIX_TYPE, MASK_MATRIX_TYPE, 
                 STATE_MATRIX_TYPE, STATE_MATRIX_CV_TYPE, 
                 LEAFCELLPOINTER_TYPE, CHILDREN_LEAFCELLPOINTER_TYPE, 
		 SHAPEDIM, STATEDIM, CHILDDIM>
::Matrix_multiply (const STATE_MATRIX_TYPE & v, STATE_MATRIX_TYPE & w) 
{
   int i,j,k;
   for (i=0; i<STATEDIM; i++) {
     for (k=0;k<SHAPEDIM;k++) w[k][i] = 0.0;
     for (k=0;k<SHAPEDIM;k++) 
       for (j=0;j<SHAPEDIM;j++) w[k][i] += A[k][j]*v[j][i];
     }
};

//////////////////////////////////////////////////////

template <typename MASS_MATRIX_TYPE, 
            typename MASK_MATRIX_TYPE, 
            typename STATE_MATRIX_TYPE, 
	    typename STATE_MATRIX_CV_TYPE, 
	    typename LEAFCELLPOINTER_TYPE, 
	    typename CHILDREN_LEAFCELLPOINTER_TYPE, 
	    int SHAPEDIM, int STATEDIM, int CHILDDIM>
void Tmassmatrix<MASS_MATRIX_TYPE, MASK_MATRIX_TYPE, 
                 STATE_MATRIX_TYPE, STATE_MATRIX_CV_TYPE, 
                 LEAFCELLPOINTER_TYPE, CHILDREN_LEAFCELLPOINTER_TYPE, 
		 SHAPEDIM, STATEDIM, CHILDDIM>
::write () 
{
   int j,k;
   for (k=0;k<SHAPEDIM;k++) {
     for (j=0;j<SHAPEDIM;j++) cerr <<  A[k][j] << "   ";
     cerr << endl;
     }
};

//////////////////////////////////////////////////////

template <typename MASS_MATRIX_TYPE, 
            typename MASK_MATRIX_TYPE, 
            typename STATE_MATRIX_TYPE, 
	    typename STATE_MATRIX_CV_TYPE, 
	    typename LEAFCELLPOINTER_TYPE, 
	    typename CHILDREN_LEAFCELLPOINTER_TYPE, 
	    int SHAPEDIM, int STATEDIM, int CHILDDIM>
void Tmassmatrix<MASS_MATRIX_TYPE, MASK_MATRIX_TYPE, 
                 STATE_MATRIX_TYPE, STATE_MATRIX_CV_TYPE, 
                 LEAFCELLPOINTER_TYPE, CHILDREN_LEAFCELLPOINTER_TYPE, 
		 SHAPEDIM, STATEDIM, CHILDDIM>
::coarse_to_fine (const STATE_MATRIX_TYPE & v, 
                 CHILDREN_LEAFCELLPOINTER_TYPE & pvLC)
{
   int i,k,l,j;
   for (j=0; j<CHILDDIM; j++) 
     for (i=0; i<STATEDIM; i++) 
       for (l=0; l<SHAPEDIM; l++) {
         pvLC[j]->u[l][i] = 0.0; 
         for (k=0;k<SHAPEDIM;k++) 
           pvLC[j]->u[l][i] += B[j][k][l]*v[k][i];
       }
       
   for (j=0; j<CHILDDIM; j++) {
     Cholesky_solve (pvLC[j]->u);
     for (i=0; i<STATEDIM; i++) 
       for (l=0; l<SHAPEDIM; l++) pvLC[j]->u[l][i] *= 4.0; // Generalize the 4.0 !!!
     }  
};	

//////////////////////////////////////////////////////

template <typename MASS_MATRIX_TYPE, 
            typename MASK_MATRIX_TYPE, 
            typename STATE_MATRIX_TYPE, 
	    typename STATE_MATRIX_CV_TYPE, 
	    typename LEAFCELLPOINTER_TYPE, 
	    typename CHILDREN_LEAFCELLPOINTER_TYPE, 
	    int SHAPEDIM, int STATEDIM, int CHILDDIM>
void Tmassmatrix<MASS_MATRIX_TYPE, MASK_MATRIX_TYPE, 
                 STATE_MATRIX_TYPE, STATE_MATRIX_CV_TYPE, 
                 LEAFCELLPOINTER_TYPE, CHILDREN_LEAFCELLPOINTER_TYPE, 
		 SHAPEDIM, STATEDIM, CHILDDIM>
::coarse_to_fine_cv (const STATE_MATRIX_TYPE & v, 
                     STATE_MATRIX_CV_TYPE & w)
{
   int i,k,l,j;
   for (j=0; j<CHILDDIM; j++) 
     for (i=0; i<STATEDIM; i++) 
       for (l=0; l<SHAPEDIM; l++) {
         w[j][l][i] = 0.0; 
         for (k=0;k<SHAPEDIM;k++) 
           w[j][l][i] += B[j][k][l]*v[k][i];
       }
       
   for (j=0; j<CHILDDIM; j++) {
     Cholesky_solve (w[j]);
     for (i=0; i<STATEDIM; i++) 
       for (l=0; l<SHAPEDIM; l++) w[j][l][i] *= 4.0; // Generalize the 4.0 !!!
     }  
};	

//////////////////////////////////////////////////////

template <typename MASS_MATRIX_TYPE, 
            typename MASK_MATRIX_TYPE, 
            typename STATE_MATRIX_TYPE, 
	    typename STATE_MATRIX_CV_TYPE, 
	    typename LEAFCELLPOINTER_TYPE, 
	    typename CHILDREN_LEAFCELLPOINTER_TYPE, 
	    int SHAPEDIM, int STATEDIM, int CHILDDIM>
void Tmassmatrix<MASS_MATRIX_TYPE, MASK_MATRIX_TYPE, 
                 STATE_MATRIX_TYPE, STATE_MATRIX_CV_TYPE, 
                 LEAFCELLPOINTER_TYPE, CHILDREN_LEAFCELLPOINTER_TYPE, 
		 SHAPEDIM, STATEDIM, CHILDDIM>
::fine_to_coarse (const CHILDREN_LEAFCELLPOINTER_TYPE & pvLC, 
                 STATE_MATRIX_TYPE & v)
{
   int i,k,l,j;
   for (i=0; i<STATEDIM; i++) 
     for (k=0;k<SHAPEDIM;k++) {
       v[k][i] = 0.0;
       for (j=0; j<CHILDDIM; j++) 
         for (l=0; l<SHAPEDIM; l++) 
           v[k][i] += B[j][k][l]*pvLC[j]->u[l][i];
       }
   
   Cholesky_solve (v);
};	

//////////////////////////////////////////////////////

template <typename MASS_MATRIX_TYPE, 
            typename MASK_MATRIX_TYPE, 
            typename STATE_MATRIX_TYPE, 
	    typename STATE_MATRIX_CV_TYPE, 
	    typename LEAFCELLPOINTER_TYPE, 
	    typename CHILDREN_LEAFCELLPOINTER_TYPE, 
	    int SHAPEDIM, int STATEDIM, int CHILDDIM>
void Tmassmatrix<MASS_MATRIX_TYPE, MASK_MATRIX_TYPE, 
                 STATE_MATRIX_TYPE, STATE_MATRIX_CV_TYPE, 
                 LEAFCELLPOINTER_TYPE, CHILDREN_LEAFCELLPOINTER_TYPE, 
		 SHAPEDIM, STATEDIM, CHILDDIM>
::fine_to_coarse_cv (const STATE_MATRIX_CV_TYPE & w, 
                    STATE_MATRIX_TYPE & v)
{
   int i,k,l,j;
   for (i=0; i<STATEDIM; i++) 
     for (k=0;k<SHAPEDIM;k++) {
       v[k][i] = 0.0;
       for (j=0; j<CHILDDIM; j++) 
         for (l=0; l<SHAPEDIM; l++) 
           v[k][i] += B[j][k][l]*w[j][l][i];
       }
   
   Cholesky_solve (v);
};	

//////////////////////////////////////////////////
///////////////                    ///////////////
///////////////     Tqrmatrix      ///////////////
///////////////                    ///////////////
//////////////////////////////////////////////////

template <typename A_MATRIX_TYPE, 
            typename Q_MATRIX_TYPE,  
            typename R_MATRIX_TYPE,  
            typename x_VECTOR_TYPE,  
            typename b_VECTOR_TYPE,  
	    int mQR, int nQR>
void Tqrmatrix<A_MATRIX_TYPE, Q_MATRIX_TYPE, R_MATRIX_TYPE,  
               x_VECTOR_TYPE, b_VECTOR_TYPE, mQR, nQR>
::set_qrmatrix (const Tqrmatrix & qr) 
{
   int i,j;
   for (i=0; i<mQR; i++) 
     for (j=0; j<mQR; j++) Q[i][j] = qr.Q[i][j];

   for (i=0; i<nQR; i++) 
     for (j=0; j<nQR; j++) R[i][j] = qr.R[i][j];
};

//////////////////////////////////////////////////////

template <typename A_MATRIX_TYPE, 
            typename Q_MATRIX_TYPE,  
            typename R_MATRIX_TYPE,  
            typename x_VECTOR_TYPE,  
            typename b_VECTOR_TYPE,  
	    int mQR, int nQR>
void Tqrmatrix<A_MATRIX_TYPE, Q_MATRIX_TYPE, R_MATRIX_TYPE,  
               x_VECTOR_TYPE, b_VECTOR_TYPE, mQR, nQR>
::write () 
{
   int j,k;
   cerr << "Q:" << endl;
   for (k=0;k<mQR;k++) {
     for (j=0;j<mQR;j++) cerr <<  Q[k][j] << "   ";
     cerr << endl;
     }
   cerr << "R:" << endl;
   for (k=0;k<nQR;k++) {
     for (j=0;j<nQR;j++) cerr <<  R[k][j] << "   ";
     cerr << endl;
     }
};

//////////////////////////////////////////////////////

template <typename A_MATRIX_TYPE, 
            typename Q_MATRIX_TYPE,  
            typename R_MATRIX_TYPE,  
            typename x_VECTOR_TYPE,  
            typename b_VECTOR_TYPE,  
	    int mQR, int nQR>
void Tqrmatrix<A_MATRIX_TYPE, Q_MATRIX_TYPE, R_MATRIX_TYPE,  
               x_VECTOR_TYPE, b_VECTOR_TYPE, mQR, nQR>
::qr_decomposition (A_MATRIX_TYPE & A) 
{
   int i,j,k;
   double alpha,beta,gamma;
   b_VECTOR_TYPE v;
   
   for (i=0; i<mQR; i++) {
     for (j=0; j<mQR; j++) Q[i][j] = 0.0;
     Q[i][i] = 1.0;
     }
     
   for (i=0; i<nQR; i++) {
     alpha = 0;
     for (j=i; j<mQR; j++) alpha += sqr(A[j][i]);
     alpha = sqrt(alpha);
     if (A[i][i] < 0) alpha = -alpha;
      
     for (j=i; j<mQR; j++) v[j] = A[j][i];
     v[i] += alpha;
     
     beta = 0.0;
     for (j=i; j<mQR; j++) beta += sqr(v[j]);
     beta = 2.0/beta;
          
     for (k=i; k<nQR; k++) {
       gamma = 0.0;
       for (j=i; j<mQR; j++) gamma += v[j]*A[j][k];
       gamma *= beta;
       for (j=i; j<mQR; j++) A[j][k] -= gamma*v[j];
       }
     for (k=0; k<mQR; k++) {
       gamma = 0.0;
       for (j=i; j<mQR; j++) gamma += v[j]*Q[j][k];
       gamma *= beta;
       for (j=i; j<mQR; j++) Q[j][k] -= gamma*v[j];
       }
     }
   
   for (i=0; i<nQR; i++) {
     for (j=0; j<i; j++) R[i][j] = 0.0; 
     for (j=i; j<nQR; j++) R[i][j] = A[i][j];
     }
};

//////////////////////////////////////////////////////

template <typename A_MATRIX_TYPE, 
            typename Q_MATRIX_TYPE,  
            typename R_MATRIX_TYPE,  
            typename x_VECTOR_TYPE,  
            typename b_VECTOR_TYPE,  
	    int mQR, int nQR>
void Tqrmatrix<A_MATRIX_TYPE, Q_MATRIX_TYPE, R_MATRIX_TYPE,  
               x_VECTOR_TYPE, b_VECTOR_TYPE, mQR, nQR>
::leastsquare_solve (const b_VECTOR_TYPE & b, x_VECTOR_TYPE & x) 
{
   int i,j;
   double sum;
   for (i=0; i<nQR; i++) {
     x[i] = 0.0;
     for (j=0; j<mQR; j++) x[i] += Q[i][j]*b[j];
     }

   x[nQR-1] /= R[nQR-1][nQR-1];
   for (i=nQR-2; i>-1; i--) {
     sum = x[i];
     for (j=i+1; j<nQR; j++) sum -= R[i][j]*x[j];
     x[i] = sum/R[i][i];
     }
};

//////////////////////////////////////////////////////
