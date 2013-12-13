#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <string>
#include <sstream>

#ifdef USE_PETSC
#include "petscksp.h"
#endif

/**
 * Macro for debugging output 
 * Advantage: If DEBUG is not set, no code is generated by the compiler.
 *            If level is to low, the if statement is ignored by the
 *            optimization techniques of the compiler whenever the 
 *            condition in the if statement can be evaluated at compile time.
 * 
 * @param level minimum level to display output
 * @param output text to output, use like ostream
 * 
 */
#ifdef DEBUG
  #define debugoutput(level,output) {if(DEBUG >= level) { cout << output; }}
  #define debugoutputcond(level,output,condition) {if(DEBUG >= level && (condition)) { cout << output; }}
#else
  #define debugoutput(level,output) {}
  #define debugoutputcond(level,output,condition) {}
#endif

const double Pi=4.0 * std::atan2(1.0, 1.0);

template < class T >
/** Convert number to string.
 * 
 * @param number number to convert
 * @return string containing number
 */
std::string NumberToString(const T number)
{
  std::stringstream str;

  str << number;
  if(str.fail())
  {
    throw("Conversion from number to string failed.");
  }
  std::string s(str.str());

  return s;
}

/** Calculate square of a double.
 * 
 * @param x given double
 * @return square of x
 */
inline double sqr(const double &x)
{
  return x * x;
}

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
inline int binomial(const unsigned long int &n, const unsigned long int &k)
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
  cout << "////////////////////////////////////////////////////" << endl;
}


////////////////////////////////////////////////////
///////////////                      ///////////////
///////////////   fixpoint format    ///////////////
///////////////                      ///////////////
////////////////////////////////////////////////////

// author James Kanze
// usage: cout << FFmt( 6 , 4 ) << f1 << FFmt( 8 , 5 ) << f2 << '\n' ;
// see also http://groups.google.de/group/de.comp.lang.iso-c++/browse_thread/thread/c3f95adf67bad7c1/9ef68ceb7959edbf?lnk=st&q=myOwner-%3Eflags&rnum=1&hl=de#9ef68ceb7959edbf

class FFmt
{
public:
  FFmt(int width, int precision);
  ~FFmt();
  friend ostream & operator<<(ostream &, FFmt const &);
private:
  int myWidth;
  int myPrecision;
  ios *myOwner;
  ios::fmtflags myOriginalFlags;
  int myOriginalPrecision;
  char myOriginalFill;
  int myPointerIndex;
};

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

ostream & operator<<(ostream & dst, FFmt const &format)
{
  static int i = ios::xalloc();
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
  dst.setf(ios::fixed, ios::floatfield);
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

 
#endif
