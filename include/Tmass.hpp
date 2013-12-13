#ifndef TMASS_HPP
#define TMASS_HPP

template <typename CONFIG_TYPE, typename Tshape>
class Tmass { // call this class Tmass, if M and G are incorporated
public:  

   enum {     
     statedim       = Tshape::statedim, 
     shapedim       = Tshape::shapedim, 
     childdim       = CONFIG_TYPE::childdim  
     };
     
   typedef typename CONFIG_TYPE::value_type              value_type;
   typedef typename CONFIG_TYPE::space_type              space_type;
   typedef typename CONFIG_TYPE::leafcell_type           leafcell_type;
   typedef typename CONFIG_TYPE::leafcellptrvector_type  leafcellptrCV_type;

   typedef typename Tshape::baryc_type                   baryc_type;
   typedef typename Tshape::state_type                   state_type;
   typedef typename Tshape::Estate_type                  Estate_type;
   typedef typename Tshape::Eshape_type                  Eshape_type;
   typedef typename Tshape::EstateCV_type                EstateCV_type; // CV=children vector
   typedef typename Tshape::Emass_type                   Emass_type;// needed in other class?
   typedef typename Tshape::Emask_type                   Emask_type;// needed in other class?

   Emass_type A, A_full; // call it M_no_details
   Emask_type B; // call it G_no_details
   
   // typedef double test_type[9]; 
   
   /* incorporate in this class, or call this class Tmass and inherit Tmass from Tmass
   typedef value_type MSAmatrix_type[childdim*shapedim][childdim*shapedim];
  
   MSAmatrix_typ M;
   MSAmatrix_typ G;
   */
   
   void set_massmatrix (const Tmass & m);
   void Cholesky_decomp ();
   void Cholesky_solve (Estate_type & w);
   void Cholesky_solve (Eshape_type & w);
   void Cholesky_solve (Estate_type & w, const unsigned int & istate);
   void Cholesky_subsolve (const int & subshapedim, Estate_type & w);
   void Cholesky_multiply (Estate_type & w);
   void Cholesky_submultiply (const int & subshapedim, Estate_type & w);
   void Cholesky_multiply (const Estate_type & v, 
                           Estate_type & w);
   void Cholesky_submultiply (const int & subshapedim, 
                              const Estate_type & v, 
			      Estate_type & w);
   void Cholesky_multiply (const Estate_type & v, const unsigned int & istate, 
                           Eshape_type & w);
   void Matrix_multiply (const Estate_type & v, Estate_type & w);
   void write (); 
   void coarse_to_fine (const Estate_type & v, 
                        leafcellptrCV_type & pvLC);
   void coarse_to_fine_cv (const Estate_type & v, 
                           EstateCV_type & w);
   void fine_to_coarse (const leafcellptrCV_type & pvLC, 
                        Estate_type & v);
   void fine_to_coarse_cv (const EstateCV_type & w, 
                           Estate_type & v);
   void fine_to_coarse_max_cv (const EstateCV_type & w, 
                               Estate_type & v);
};

//////////////////////////////////////////////////////

template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::set_massmatrix (const Tmass & m) 
{
   int i,j,k;
   for (k=0;k<shapedim;k++) 
     for (i=0;i<shapedim;i++) { 
       A[k][i]    = m.A[k][i];
       for (j=0; j<childdim; j++) B[j][k][i]    = m.B[j][k][i];
       }
};

//////////////////////////////////////////////////////

template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::Cholesky_decomp () 
{
   int i,j,k;
   for (k=0;k<shapedim;k++) {
     for (j=0;j<k;j++) A[k][k] -= A[k][j]*A[k][j]*A[j][j];
     for (i=k+1;i<shapedim;i++) {
       for (j=0;j<k;j++) A[i][k] -= A[i][j]*A[j][j]*A[k][j];
       A[i][k]    /= A[k][k];
       }
     }
};

//////////////////////////////////////////////////////

template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::Cholesky_solve (Estate_type & w) 
{
   int i,j,k;
   for (i=0; i<statedim; i++) {
     for (k=1;k<shapedim;k++) {
       for (j=0;j<k;j++) w[k][i] -= A[k][j]*w[j][i];
       }
     for (k=0;k<shapedim;k++) w[k][i] /= A[k][k];
     for (k=shapedim-2;k>-1;k--) {
       for (j=k+1;j<shapedim;j++) w[k][i] -= A[j][k]*w[j][i];
     }
   }
};

//////////////////////////////////////////////////////

template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::Cholesky_solve (Eshape_type & w) 
{
   int j,k;
   for (k=1;k<shapedim;k++) {
     for (j=0;j<k;j++) w[k] -= A[k][j]*w[j];
     }
   for (k=0;k<shapedim;k++) w[k] /= A[k][k];
   for (k=shapedim-2;k>-1;k--) {
     for (j=k+1;j<shapedim;j++) w[k] -= A[j][k]*w[j];
     }
};

//////////////////////////////////////////////////////

template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::Cholesky_solve (Estate_type & w, const unsigned int & istate) 
{
   int j,k;
   for (k=1;k<shapedim;k++) {
     for (j=0;j<k;j++) w[k][istate] -= A[k][j]*w[j][istate];
     }
   for (k=0;k<shapedim;k++) w[k][istate] /= A[k][k];
   for (k=shapedim-2;k>-1;k--) {
     for (j=k+1;j<shapedim;j++) w[k][istate] -= A[j][k]*w[j][istate];
     }
};

//////////////////////////////////////////////////////

template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::Cholesky_subsolve (const int & subshapedim, Estate_type & w) 
{
   int i,j,k;
   for (i=0; i<statedim; i++) {
     for (k=1;k<subshapedim;k++) {
       for (j=0;j<k;j++) w[k][i] -= A[k][j]*w[j][i];
       }
     for (k=0;k<subshapedim;k++) w[k][i] /= A[k][k];
     for (k=subshapedim-2;k>-1;k--) {
       for (j=k+1;j<subshapedim;j++) w[k][i] -= A[j][k]*w[j][i];
     }
     
     for (k=subshapedim;k<shapedim;k++) w[k][i] = 0.0;
   }
};

//////////////////////////////////////////////////////

template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::Cholesky_multiply (Estate_type & w) 
{
   int i,j,k;
   for (i=0; i<statedim; i++) {
     for (k=0;k<shapedim;k++) {
       for (j=k+1;j<shapedim;j++) w[k][i] += A[j][k]*w[j][i];
       w[k][i] *= A[k][k];
       }
     for (k=shapedim-1;k>-1;k--) {
       for (j=0;j<k;j++) w[k][i] += A[k][j]*w[j][i];
       }
     }
};

//////////////////////////////////////////////////////

template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::Cholesky_submultiply (const int & subshapedim, Estate_type & w) 
{
   int i,j,k;
   for (i=0; i<statedim; i++) {
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

template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::Cholesky_multiply (const Estate_type & v, Estate_type & w) 
{
   int i,j,k;
   for (i=0; i<statedim; i++) {
     for (k=0;k<shapedim;k++) w[k][i] = v[k][i];
     for (k=0;k<shapedim;k++) {
       for (j=k+1;j<shapedim;j++) w[k][i] += A[j][k]*w[j][i];
       w[k][i] *= A[k][k];
       }
     for (k=shapedim-1;k>-1;k--) {
       for (j=0;j<k;j++) w[k][i] += A[k][j]*w[j][i];
       }
     }
};

//////////////////////////////////////////////////////

template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::Cholesky_multiply (const Estate_type & v, const unsigned int & istate, Eshape_type & w) 
{
   int j,k;
   for (k=0;k<shapedim;k++) w[k] = v[k][istate];
   for (k=0;k<shapedim;k++) {
     for (j=k+1;j<shapedim;j++) w[k] += A[j][k]*w[j];
     w[k] *= A[k][k];
     }
   for (k=shapedim-1;k>-1;k--) {
     for (j=0;j<k;j++) w[k] += A[k][j]*w[j];
     }
};

//////////////////////////////////////////////////////

template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::Cholesky_submultiply (const int & subshapedim, const Estate_type & v, 
                       Estate_type & w) 
{
   int i,j,k;
   for (i=0; i<statedim; i++) {
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

template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::Matrix_multiply (const Estate_type & v, Estate_type & w) 
{
   int i,j,k;
   for (i=0; i<statedim; i++) {
     for (k=0;k<shapedim;k++) w[k][i] = 0.0;
     for (k=0;k<shapedim;k++) 
       for (j=0;j<shapedim;j++) w[k][i] += A[k][j]*v[j][i];
     }
};

//////////////////////////////////////////////////////

template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::write () 
{
   int j,k;
   
   cerr << "A_full: " << endl;
   for (k=0;k<shapedim;k++) {
     for (j=0;j<shapedim;j++) cerr <<  A_full[k][j] << "   ";
     cerr << endl;
     }     
   cerr << "A: " << endl;
   for (k=0;k<shapedim;k++) {
     for (j=0;j<shapedim;j++) cerr <<  A[k][j] << "   ";
     cerr << endl;
     }
};

//////////////////////////////////////////////////////

template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::coarse_to_fine (const Estate_type & v, 
                 leafcellptrCV_type & pvLC)
{
   int i,k,l,j;
   for (j=0; j<childdim; j++) 
     for (i=0; i<statedim; i++) 
       for (l=0; l<shapedim; l++) {
         pvLC[j]->u[l][i] = 0.0; 
         for (k=0;k<shapedim;k++) 
           pvLC[j]->u[l][i] += B[j][k][l]*v[k][i];
       }
       
   for (j=0; j<childdim; j++) {
     Cholesky_solve (pvLC[j]->u);
     for (i=0; i<statedim; i++) 
       for (l=0; l<shapedim; l++) pvLC[j]->u[l][i] *= childdim; // Generalize the 4.0 !!!
     }  
};	

//////////////////////////////////////////////////////

template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::coarse_to_fine_cv (const Estate_type & v, 
                     EstateCV_type & w)
{
   int i,k,l,j;
   for (j=0; j<childdim; j++) 
     for (i=0; i<statedim; i++) 
       for (l=0; l<shapedim; l++) {
         w[j][l][i] = 0.0; 
         for (k=0;k<shapedim;k++) 
           w[j][l][i] += B[j][k][l]*v[k][i];
       }
       
   for (j=0; j<childdim; j++) {
     Cholesky_solve (w[j]);
     for (i=0; i<statedim; i++) 
       for (l=0; l<shapedim; l++) w[j][l][i] *= childdim; // Generalize the 4.0 !!!
     }  
};	

//////////////////////////////////////////////////////

template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::fine_to_coarse (const leafcellptrCV_type & pvLC, 
                 Estate_type & v)
{
   int i,k,l,j;
   for (i=0; i<statedim; i++) 
     for (k=0;k<shapedim;k++) {
       v[k][i] = 0.0;
       for (j=0; j<childdim; j++) 
         for (l=0; l<shapedim; l++) 
           v[k][i] += B[j][k][l]*pvLC[j]->u[l][i];
       }
   
   Cholesky_solve (v);
};	

//////////////////////////////////////////////////////

template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::fine_to_coarse_cv (const EstateCV_type & w, 
                    Estate_type & v)
{
   int i,k,l,j;
   for (i=0; i<statedim; i++) 
     for (k=0;k<shapedim;k++) {
       v[k][i] = 0.0;
       for (j=0; j<childdim; j++) 
         for (l=0; l<shapedim; l++) 
           v[k][i] += B[j][k][l]*w[j][l][i];
       }
   
   Cholesky_solve (v);
};	

//////////////////////////////////////////////////////

// Only for linear Lagrange-shapes at the moment, 
// this fine to coarse is meant to maintain maxima !!!

template <typename CONFIG_TYPE, typename Tshape>
void Tmass<CONFIG_TYPE, Tshape>
::fine_to_coarse_max_cv (const EstateCV_type & w, 
                    Estate_type & v)
{
   v[0][0] = w[1][0][0];
   v[1][0] = w[2][1][0];
   v[2][0] = w[3][2][0];
};	

#endif
