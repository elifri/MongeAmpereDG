
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
