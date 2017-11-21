/*
 * imageOT.hpp
 *
 *  Created on: Apr 15, 2016
 *      Author: friebel
 */

#ifndef SRC_OT_IMAGEOT_HPP_
#define SRC_OT_IMAGEOT_HPP_


#include <CImg.h>

#include "config.h"
#include "ImageFunction.hpp"

template<typename Function, typename DerivativeFunction>
void print_image_OT(const Function& T, const DerivativeFunction& Jac_T,
    const ImageFunction& inputImage, const ImageFunction& targetImage,
    const std::string& outputFile, int pixel_width, int pixel_height, int scale = 6)
{
  const int originalHeight = inputImage.getOriginalImage().height();
  const int originalWidth = inputImage.getOriginalImage().width();
  cimg_library::CImg<double> transported(pixel_width,pixel_height);

  transported.fill(0);

  for (int i = 0; i < pixel_height*scale; i++)
  {
    for (int j = 0; j < pixel_height*scale; j++)
    {
      const double i_original = (double)i/pixel_width/scale*originalWidth;
      const double j_original = (double)j/pixel_height/scale*originalHeight;
      const auto x = inputImage.getCoordinate(i_original,j_original);
      auto transportedPixel = T(x);

      double transportedPixelvalue;
      inputImage.evaluate(x,transportedPixelvalue);


      typename DerivativeFunction::Hessian hessU;
      Jac_T.evaluateHess(x, hessU);
//      transportedPixelvalue/=(hessU[0][0]*hessU[1][1]-hessU[1][0]*hessU[0][1]);


      int i_transported, j_transported;
      targetImage.getPixel(transportedPixel, i_transported, j_transported);
      transported(i_transported, j_transported) += transportedPixelvalue;
      if (i % 100 == 0 && j % 10 == 0)
      {
        std::cerr << " i,j = (" << i<< "," <<  j << " ) -> x " << x << ", is transported to " << transportedPixel << std::endl;

        std::cerr << " add to (" << i_transported << "," <<  j_transported << " ) the value " << transportedPixelvalue << ", dethess was " << hessU[0][0]*hessU[1][1]-hessU[1][0]*hessU[0][1];
        std::cerr << " -> " << transported(i_transported, j_transported) << std::endl;
      }

    }
  }

  transported /= targetImage.factor_;
  transported /= 255.0 / (255.0 + targetImage.minValue());
  transported -= targetImage.minValue();

  transported.save(outputFile.c_str());
}



#endif /* SRC_OT_IMAGEOT_HPP_ */
