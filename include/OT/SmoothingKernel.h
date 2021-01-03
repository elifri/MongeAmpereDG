/*
 * SmoothingKernel.h
 *
 *  Created on: Feb 10, 2017
 *      Author: friebel
 */

#ifndef INCLUDE_OT_SMOOTHINGKERNEL_H_
#define INCLUDE_OT_SMOOTHINGKERNEL_H_

struct SmoothingKernel
{


  SmoothingKernel()
  {
    std::cout << " init smoothing kernel " << std::endl;

    smoothingKernelValues.resize(2*n_+1, 2*n_+1);
//////---------Gaussian Filter ---------

    /*
    for (int i = -n_; i <= n_; i++)
      for (int j = -n_; j <= n_; j++)
      {
        smoothingKernelValues.coeffRef(i+n_,j+n_) = 1./2./M_PI/sigmaSqr_*std::exp(-(i*i+j*j)/2./sigmaSqr_);
//        std::cout << " kernel " << i << ' ' << j << " " << smoothingKernelValues[i+n][j+n] << std::endl;
      }

*/

//-----------twice a moving average----------


    static_assert(n_== 2, "wrong filter dimension");
    for (int i = -n_; i <= n_; i++)
      for (int j = -n_; j <= n_; j++)
      {
        const int distance = std::abs(i) + std::abs(j);
        switch(distance)
        {
        case 0: smoothingKernelValues.coeffRef(i+n_,j+n_) = 9; break;
        case 1: smoothingKernelValues.coeffRef(i+n_,j+n_) = 6; break;
        case 2:
          if (i == 1 || i == -1)
            smoothingKernelValues.coeffRef(i+n_,j+n_) = 4;
          else
            smoothingKernelValues.coeffRef(i+n_,j+n_) = 3;
          break;
        case 3: smoothingKernelValues.coeffRef(i+n_,j+n_) = 2; break;
        case 4: smoothingKernelValues.coeffRef(i+n_,j+n_) = 1; break;
        }
      }


    //-----------simple moving average----------


/*
    static_assert(n_== 1, " wrong filter dimension");

        for (int i = -n_; i <= n_; i++)
          for (int j = -n_; j <= n_; j++)
          {
            smoothingKernelValues.coeffRef(i+n_,j+n_) = 1;
    //        std::cout << " kernel " << i << ' ' << j << " " << smoothingKernelValues[i+n][j+n] << std::endl;
          }
*/

    //-----------twice a moving average----------

/*

    static_assert(n_== 3, "wrong filter dimension");
    for (int i = -n_; i <= n_; i++)
      for (int j = -n_; j <= n_; j++)
      {
        const int distance = std::abs(i) + std::abs(j);
        switch(distance)
        {
        case 0: smoothingKernelValues.coeffRef(i+n_,j+n_) = 49; break;
        case 1: smoothingKernelValues.coeffRef(i+n_,j+n_) = 42; break;
        case 2:
          if (i == 1 || i == -1)
            smoothingKernelValues.coeffRef(i+n_,j+n_) = 36;
          else
            smoothingKernelValues.coeffRef(i+n_,j+n_) = 21;
          break;
        case 3:
          if (i == 0 || j == 0)
            smoothingKernelValues.coeffRef(i+n_,j+n_) = 7;
          else
            smoothingKernelValues.coeffRef(i+n_,j+n_) = 18;
          break;
        case 4:
          if (i == 2 || i == -2)
            smoothingKernelValues.coeffRef(i+n_,j+n_) = 9;
          else
            smoothingKernelValues.coeffRef(i+n_,j+n_) = 6;
          break;
        case 5: smoothingKernelValues.coeffRef(i+n_,j+n_) = 3; break;
        case 6: smoothingKernelValues.coeffRef(i+n_,j+n_) = 1; break;
        }
      }
*/


    //--------no smoothing------------------
/*
    static_assert(n_== 0, " wrong filter dimension");
    smoothingKernelValues.coeffRef(0,0) = 1.;
*/

    smoothingKernelValues /= smoothingKernelValues.sum();
    std::cout << "kernel " << smoothingKernelValues << " sum " << smoothingKernelValues.sum() << std::endl;
  }

//  double operator[](int i){ return smoothingKernelValues[i];}
  double operator()(int i, int j) const{ return smoothingKernelValues(i,j);}

  static const int n_ = 2;
  Config::DenseMatrixType smoothingKernelValues;
  static constexpr double sigmaSqr_=2;
};




#endif /* INCLUDE_OT_SMOOTHINGKERNEL_H_ */
