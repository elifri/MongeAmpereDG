/**
 *  \file timer.hpp
 *  \author Yasemin Hafizogullari and Andreas Platen
 *  \date 11.2011
 *  \brief
 */


#ifndef MAGICMIRRORTOOLS_TIMER_HPP
#define MAGICMIRRORTOOLS_TIMER_HPP



namespace mirror_problem
{


/*! \class Timer
 *  \brief Class for time measurement
 */
class Timer
{
  private:
    double tstart;

  public:
    /*! \brief Start time measurement */
    void   tic() { tstart = clock(); }
    /*! \brief Return measured time (does not stop the measurement) */
    double toc() { double time=clock()-tstart;
                   return time/CLOCKS_PER_SEC; }
};


} // end of namespace mirror_problem

#endif // MAGICMIRRORTOOLS_TIMER_HPP
