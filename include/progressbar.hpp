// copied from http://andreas.familie-steinel.de/trackbacks?article_id=475


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifdef __APPLE__
#define CLOCK_REALTIME 1
#include <mach/mach_time.h>
#endif


#ifndef PROGRESSBAR_HPP
#define PROGRESSBAR_HPP


class ProgressBar
{
    public:

        ProgressBar()
        {
            timer = NULL;
            memset(buffer, ' ', 80);
            size = 50;
        };

        void start()
        {
            if ( timer == NULL )
            {
                timer = (struct timespec*) malloc (sizeof(struct timespec));
            }
            my_clock_gettime(CLOCK_REALTIME, timer);
        }

        // returns computing time
        double stop()
        {
            set_buffer(100.0);
            int sec = (int) time_elapsed();
            int min = sec/60;
            printf("%5.1f%% %s took %4i min %2i sec\n", 100.0, buffer, min, int(sec%60));
            fflush(NULL);
            return time_elapsed();
        }

        double status(unsigned int current, unsigned int max)
        {
            float percent = (float) current / (float) max * 100.0;
            set_buffer(percent);
            float time = time_elapsed();
            int sec = (int) ((100.0 / percent * time) - time);
            int min = sec/60;
            if (abs(min) >= 1e4)
                printf("%5.1f%% %s >10000 min left\r", percent, buffer);
            else
                printf("%5.1f%% %s %4i min %2i sec left\r", percent, buffer, min, int(sec%60));
            fflush(NULL);

            return time;
        }

    private:

        void set_buffer(float percent)
        {
            buffer[0] = '[';
            unsigned int steps = floor( percent / 100 * size);
            unsigned int pos = 1;
            for (; pos <= steps; pos++)
            {
                buffer[pos] = '-';
            }
            pos++;
            for (; pos < size; pos++)
            {
                buffer[pos] = ' ';
            }

            buffer[size] = ']';
            buffer[size +1] = '\0';


        }

        float time_elapsed()
        {
            struct timespec *timer2 = (struct timespec*) malloc (sizeof(struct timespec));

            // save current time
            my_clock_gettime(CLOCK_REALTIME, timer2);

            // compute right time difference
            int sec = timer2->tv_sec - timer->tv_sec;
            int nsec = (timer2->tv_nsec - timer->tv_nsec);

            if (nsec < 0)
            {
                sec--;
                nsec = 1000000000 +nsec;
            }

            float ret = sec + (((float) nsec) / (1000*1000*1000));
            free(timer2);
            return ret;
        }

        int my_clock_gettime(int b, struct timespec *ts)
        {
#ifdef __APPLE__
                float conv = 0.0;
                uint64_t time = mach_absolute_time();
                mach_timebase_info_data_t info;
                kern_return_t err = mach_timebase_info( &info );
                if ( err == 0 )
                    conv = 1e-9 * (float) info.numer / (float) info.denom;

                uint64_t seconds = (uint64_t) (conv * (float) time);
                time -= (uint64_t) ((float) seconds / conv);

                ts->tv_sec = (time_t) seconds;
                ts->tv_nsec = (long) time;

                return 0;
#else
                return clock_gettime(b,ts);
#endif
        }

        struct timespec *timer;          // timer (global)
        char buffer[80];
        unsigned int size;

};

#endif // PROGRESSBAR_HPP
