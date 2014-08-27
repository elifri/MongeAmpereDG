/*
 * utility.cpp
 *
 *  Created on: 07.01.2014
 *      Author: elisa
 */

#include "utility.hpp"

/** Calculate factorial by loop.
 *
 * @param n given number
 * @return factorial of n
 */
inline unsigned long int factorial(const unsigned long int &n)
{
  if (n != 0)
  {

    unsigned long int f = n;
    for (unsigned long int j = n - 1; j > 1; --j)
      f *= j;
    return f;

  }
  else
  {

    return 1;

  }
}

/** Calculate binomial coefficient "n choose k".
 *
 * @param n upper parameter
 * @param k lower parameter
 * @return binomial coefficient "n choose k"
 */
int binomial(const unsigned long int &n, const unsigned long int &k)
{

  if (k != 0)
  {

    unsigned long int b = n;
    for (unsigned long int j = 1; j < k; ++j)
      b *= n - j;

    return b / factorial(k);
  }
  else
  {
    return 1;
  }

}



/** Output separator line to standard output
 */
void outputseparator(){
  std::cout << "////////////////////////////////////////////////////" << std::endl;
}


////////////////////////////////////////////////////
///////////////                      ///////////////
///////////////   fixpoint format    ///////////////
///////////////                      ///////////////
////////////////////////////////////////////////////

// author James Kanze
// usage: cout << FFmt( 6 , 4 ) << f1 << FFmt( 8 , 5 ) << f2 << '\n' ;
// see also http://groups.google.de/group/de.comp.lang.iso-c++/browse_thread/thread/c3f95adf67bad7c1/9ef68ceb7959edbf?lnk=st&q=myOwner-%3Eflags&rnum=1&hl=de#9ef68ceb7959edbf

FFmt::FFmt(int width, int precision)
    :myWidth(width)
    , myPrecision(precision)
    , myOwner(NULL)
{}

FFmt::~FFmt()
{
  if (myOwner != NULL)
  {
    myOwner->flags(myOriginalFlags);
    myOwner->precision(myOriginalPrecision);
    myOwner->fill(myOriginalFill);
    myOwner->pword(myPointerIndex) = NULL;
  }
}

std::ostream & operator<<(std::ostream & dst, FFmt const &format)
{
  static int i = std::ios::xalloc();
  if (dst.pword(i) == NULL)
  {
    FFmt & f = const_cast < FFmt & >(format);
    dst.pword(i) = &dst;
    f.myOwner = &dst;
    f.myOriginalFlags = dst.flags();
    f.myOriginalFill = dst.fill();
    f.myOriginalPrecision = dst.precision();
    f.myPointerIndex = i;
  }
  dst.setf(std::ios::fixed, std::ios::floatfield);
  dst.precision(format.myPrecision);
  dst.width(format.myWidth);
  return dst;
}

#ifdef USE_PETSC

void check_petsc_error(PetscErrorCode ierr)
{
  if (ierr != 0)
  {
    cerr << "PETSc Error! " << ierr << endl;
    exit(ierr);
  }
}

#endif

template <class Cont>
Cont toupper(Cont cont)
{
  std::transform(cont.begin(), cont.end(), cont.begin(),
                 static_cast<int(*)(int)>(std::toupper));
  return cont;
}

bool string_contains_true(const std::string s) {
  std::string s2=toupper(s);
  return ( s2=="Y" || s2=="YES" || s2=="TRUE");
}


double fRand(const double fMin, const double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}
