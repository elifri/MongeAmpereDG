/*
 * localfunction_utils.hh
 *
 *  Created on: Jul 11, 2017
 *      Author: friebel
 */

#ifndef INCLUDE_LOCALFUNCTIONS_LOCALFUNCTION_UTILS_HH_
#define INCLUDE_LOCALFUNCTIONS_LOCALFUNCTION_UTILS_HH_

namespace Dune
{

template<typename DomainLocal, typename BarycCoordType>
inline
void calcBarycCoords(const DomainLocal& x, BarycCoordType& baryc)
   {
     assert(x[0] < 1+1e-12 && x[1] <= 1+1e-12 && x[0] > -1e-12 && x[1] > -1e-12 && " These are not local coordinates");
/*
     const auto& b0 = geo_.corner(0);
     const auto& b1 = geo_.corner(1);
     const auto& b2 = geo_.corner(2);

     baryc[0] = (b1[0]*b2[1]-b1[0]*x[1]-b1[1]*b2[0]+b1[1]*x[0]+b2[0]*x[1]-b2[1]*x[0])/determinantBarycTrafo_;
     baryc[1] = -(b0[0]*b2[1]-b0[0]*x[1]-b0[1]*b2[0]+b0[1]*x[0]+b2[0]*x[1]-b2[1]*x[0])/determinantBarycTrafo_;
     baryc[2] = (b0[0]*b1[1]-b0[0]*x[1]-b0[1]*b1[0]+b0[1]*x[0]+b1[0]*x[1]-b1[1]*x[0])/determinantBarycTrafo_;
*/
     baryc[1] = x[0];
     baryc[2] = x[1];
     baryc[0] = 1 - x[0] -x[1];

     for (int i = 0; i < 3; i++)
       if (baryc[i] < 0)
       {
         baryc[i] = 0;
         baryc[(i+1)%3] = 1 - baryc[i%3] -baryc[(i+2)%3];
       }

//     std::cerr << std::setprecision(16);
//     std::cerr << " baryc sum " << baryc[0]+baryc[1]+baryc[2] << " > 0? " << !(baryc[0]+baryc[1]+baryc[2] > 1) << std::endl;

     assert(baryc[0] >= 0);
     assert(baryc[1] >= 0);
     assert(baryc[2] >= 0);
//     assert(!(baryc[0]+baryc[1]+baryc[2] > 1));

   }

}

#endif /* INCLUDE_LOCALFUNCTIONS_LOCALFUNCTION_UTILS_HH_ */
