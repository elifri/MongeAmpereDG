// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_SSplineLOCALINTERPOLATION_HH
#define DUNE_SSplineLOCALINTERPOLATION_HH

#include <vector>

namespace Dune
{
  template<class LB>
  class PS12SSplineInterpolation
  {
    /** \brief The number of degrees of freedom */
    enum {N = LB::N};

    /** \brief Export the element order */
    enum {k = LB::O};
    /** \brief Export the element order */
    enum {Dim = 2};

    const typename LB::Traits::RangeFieldType h = 1e-5;

  public:
    PS12SSplineInterpolation(const typename LB::Geometry& geo): geo_(geo) {}


    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      assert(false);
      out.resize(N);

      int k = 0;
      for (int i = 0; i < geo_.corners(); i++)
      {
        typename LB::Traits::RangeType tempvalue;
        typename LB::Traits::RangeFieldType value;
        f.evaluate(geo_.local(geo_.corner(i)), tempvalue);
        value = tempvalue;

        //set dofs associated with values at vertices
//        assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);
        out[k++] = value;

        //test if this was the right basis function
//        {
//          std::vector<FieldVector<double, 1> > functionValues(lFE.size());
//          lFE.localBasis().evaluateFunction(geo_.local(geo_.corner(i)), functionValues);
//          assert(std::abs(functionValues[k-1][0]-1) < 1e-10);
//        }

        typename LB::Traits::RangeType tempfPlus;
        typename LB::Traits::RangeFieldType fPlus;

        //set dofs associated with gradient values at vertices
        auto xValuePlus = geo_.local(geo_.corner(i));
        xValuePlus[0] += i % 2 == 0 ? h : - h;

//        assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);
        f.evaluate(xValuePlus,tempfPlus);
        fPlus = tempfPlus;
        out[k++] = i % 2 == 0 ? (fPlus-value) / h : -(fPlus-value) / h;

        xValuePlus = geo_.corner(i);
        xValuePlus[1] += i < 2 ? h : - h;

        f.evaluate(xValuePlus,tempfPlus);
        fPlus = tempfPlus;
        out[k++] = i < 2 ? (fPlus-value) / h : -(fPlus-value) / h;

//        //test if this were the right basis function
//        {
//          std::vector<FieldMatrix<double, 1, 2> > jacobianValues(lFE.size());
//          lFE.localBasis().evaluateJacobian(geo_.local(geo_.corner(i)), jacobianValues);
//          assert(std::abs(jacobianValues[k-2][0][0]-1) < 1e-10);
//          assert(std::abs(jacobianValues[k-1][0][1]-1) < 1e-10);
//        }

        k++;
      }

      assert(k == 12);

      const auto b0 = geo_.corner(0), b1 = geo_.corner(1), b2 = geo_.corner(2);

      std::vector<FieldVector<double, 2>> normals(geo_.corners());
      normals[0] = {b0[1]-b1[1] , - b0[0]+b1[0]}; if (normals[0]*(b2-b0) > 0) normals[0]*=-1;
      normals[1] = {b2[1]-b1[1] , - b2[0]+b1[0]}; if (normals[1]*(b0-b2) > 0) normals[1]*=-1;
      normals[2] = {b0[1]-b2[1] , - b0[0]+b2[0]}; if (normals[2]*(b1-b0) > 0) normals[2]*=-1;

      std::vector<FieldVector<double, 2>> edgemids(geo_.corners());
      edgemids[0] = geo_.corner(0)+geo_.corner(1); edgemids[0]/=2.;
      edgemids[1] = geo_.corner(1)+geo_.corner(2); edgemids[1]/=2.;
      edgemids[2] = geo_.corner(0)+geo_.corner(2); edgemids[2]/=2.;

      for (int i; i < geo_.corners(); i++) //loop over edges
      {
        // normal of center in face's reference element
        const FieldVector<double, Dim>& normal = normals[i];

        bool unit_pointUpwards;
        if (std::abs(normal[0]+normal[1])< 1e-12)
          unit_pointUpwards = (normal[1] > 0);
        else
          unit_pointUpwards = (normal[0]+normal[1] > 0);

        const auto& face_center = edgemids[i];

        FieldVector<double, 2> approxGradientF;

        typename LB::Traits::RangeType tempvalue, tempfPlus;
        typename LB::Traits::RangeFieldType value, fPlus;
        f.evaluate(face_center, tempvalue); value = tempvalue;

        //calculate finite difference in x0-direction
        auto xValuePlus = face_center;
        xValuePlus[0] += !(normal[0] > 0) ? h : - h;
        f.evaluate(xValuePlus,tempfPlus); fPlus = tempfPlus;
        approxGradientF[0] = !(normal[0] > 0)? (fPlus-value) / h : -(fPlus-value) / h;

        //calculate finite difference in x1-direction
        xValuePlus = face_center;
        xValuePlus[1] += !(normal[1] > 0) ? h : - h;
        f.evaluate(xValuePlus,tempfPlus); fPlus = tempfPlus;
        approxGradientF[1] = !(normal[1] > 0) ? (fPlus-value) / h : -(fPlus-value) / h;

        if (i == 0)
          k = 3;
        else
          if (i == 1)
            k = 11;
          else
            k = 7;

//        assert(lFE.localCoefficients().localKey(k).subEntity() == (unsigned int) i);
        out[k++] = unit_pointUpwards ? (approxGradientF*normal) : -(approxGradientF*normal);
      }
    }

    const typename LB::Geometry& geo_;
  };
}

#endif
