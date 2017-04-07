
/**************************************************************************
*
* Cinderella -- COMPASS ONLINE FILTER
*
* Copyright (C) 2002-2004
*   Technische Universititaet Muenchen
*   Physik-Department
*   James Frank Strasse
*   D-85748 Garching
*   Germany
*
* Author(s) of this file: TN
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*
**************************************************************endofheader*/

/** \file
    Silicon Time Reconstruction.  Instructions for use outside of Cinderella
    may be found in silicon_timing_coral.h */

#include "silicon_timing.h"

using namespace std;

/** 
 * Reconstruct time information for all activated strips of a detector plane.
 * The 'noises' array contains the noise levels for all strips of the plane
 * (starting with strip zero at position zero).  Details of time
 * reconstruction are controlled through the 'options' structure.  The 'debug'
 * structure is intended for use inside Cinderella.  If debugging is not
 * desired, a NULL pointer may be given.  Input and output is handled through
 * the 'strip_hit' array of size 'strip_count'.  The digits 'a0f', 'a1f' and
 * 'a2f' are required as input, as output times 'time', 'time0', 'time1' and
 * errors 'el', 'eh', 'el0', 'eh0', 'el1', 'eh1' are provided.  Additionally,
 * consistency information is calculated for every strip and stored in the
 * field 'consistency'.
 **/
void
silicon_calculate_timings(strip_t * strip_hits, const unsigned int strip_count,
                          const double *const noises,
                          const time_reconstruction_option_t *const options,
                          clusterization_debug_t *debug) {

  unsigned int        i;

  for (i = 0; i < strip_count; i++) {
    strip_t             *strip = &strip_hits[i];
    const double         noise = noises[strip->wire];
    const double         a0 = strip->a0f;
    const double         a1 = strip->a1f;
    const double         a2 = strip->a2f;
    const double         r0 = strip->a0f / strip->a2f;
    const double         r1 = strip->a1f / strip->a2f;

    MSG(_DEBUG3, silicon, "\ttcsphase %.2f, tcscycle %.2f, delta %.2f, "
        "decay_time %.2f, min_early_error %.2f",
        options->tcsphase, options->tcscycle, options->delta,
        options->decay_time, options->min_early_error);

    /* fill a2 histogram */
    if (debug) {
      histogram_fill_raw(debug->a2, a2);

      /* some debug/devel histogramming on the ratio level */
      if (fpclassify(a2) != FP_ZERO) {

        histogram2_fill(debug->r0_phase, options->tcsphase, r0);
        histogram2_fill(debug->r1_phase, options->tcsphase, r1);

        /* fill banana plot histogram */
        histogram2_fill(debug->r0_r1, r1, r0);
      }
    }

    /* calculate single timing for r0 */
    single_timing(&strip->time0, &strip->el0, &strip->eh0, a0, a2, 0,
                  noise, options, 1);
    /* logging and histogramming */
    if (debug && isfinite(strip->time0)) {
      MSG(_DEBUG3, silicon, "\twire %i, t0 %.2f, el0 %.2f, eh0 %.2f",
          strip->wire, strip->time0, -strip->el0, strip->eh0);
      histogram_fill(debug->time0, strip->time0);
      histogram_fill(debug->clocktime0, strip->time0 - options->tcsphase);
      histogram_fill(debug->time0_error, -strip->el0);
      histogram_fill(debug->time0_error, strip->eh0);
      histogram_fill(debug->pull0, strip->time0
                     / error_towards(strip->time0, strip->el0, strip->eh0,
                                     0.));
      histogram2_fill(debug->err0_ratio, a0 / a2, -strip->el0);
      histogram2_fill(debug->err0_ratio, a0 / a2, strip->eh0);
      histogram2_fill(debug->err0_vs_time, strip->time0, -strip->el0);
      histogram2_fill(debug->err0_vs_time, strip->time0, strip->eh0);
      histogram2_fill(debug->r0_vs_clocktime,
                      strip->time0 - options->tcsphase, r0);
    }

    /* calculate single timing for r1 */
    single_timing(&strip->time1, &strip->el1, &strip->eh1, a1, a2, 1,
                  noise, options, 1);
    /* logging and histogramming */
    if (debug && isfinite(strip->time1)) {
      MSG(_DEBUG3, silicon, "\twire %i, t1 %.2f, el1 %.2f, eh1 %.2f",
          strip->wire, strip->time1, -strip->el1, strip->eh1);
      histogram_fill(debug->time1, strip->time1);
      histogram_fill(debug->clocktime1, strip->time1 - options->tcsphase);
      histogram_fill(debug->time1_error, -strip->el1);
      histogram_fill(debug->time1_error, strip->eh1);
      histogram_fill(debug->pull1, strip->time1
                     / error_towards(strip->time1, strip->el1, strip->eh1,
                                     0.));
      histogram2_fill(debug->err1_ratio, a1 / a2, -strip->el1);
      histogram2_fill(debug->err1_ratio, a1 / a2, strip->eh1);
      histogram2_fill(debug->err1_vs_time, strip->time1, -strip->el1);
      histogram2_fill(debug->err1_vs_time, strip->time1, strip->eh1);
      histogram2_fill(debug->r1_vs_clocktime,
                      strip->time1 - options->tcsphase, r1);
    }

    /* combine r0 and r1 timings to one joint timing */
    join_timing(&strip->time, &strip->el, &strip->eh,
                strip->time0, strip->el0, strip->eh0,
                strip->time1, strip->el1, strip->eh1);

    /* fill timing histogram */
    if (debug && isfinite(strip->time)) {
      debug->valid_strips++;
      MSG(_DEBUG3, silicon, "combined strip time: wire %u: %.2f -%.2f/+%.2f",
          strip->wire, strip->time, strip->el, strip->eh);
      histogram_fill(debug->time, strip->time);
      histogram_fill(debug->clocktime, strip->time - options->tcsphase);
      histogram_fill(debug->time_error, -strip->el);
      histogram_fill(debug->time_error, strip->eh);
      histogram_fill(debug->pull, strip->time
                     / error_towards(strip->time, strip->el, strip->eh, 0.));
      histogram2_fill(debug->err_vs_time, strip->time, -strip->el);
      histogram2_fill(debug->err_vs_time, strip->time, strip->eh);
    } else if (debug)
      MSG(_DEBUG3, silicon,
          "combined strip time: wire %u: %.2f -%.2f/+%.2f (BAD)",
          strip->wire, strip->time, strip->el, strip->eh);

    /* calculate intrinsic consistency of strip */
    strip->consistency =
        time_consistency(strip->time0, strip->el0, strip->eh0,
                         strip->time1, strip->el1, strip->eh1);

    /* fill intrinsic consistency histogram */
    if (debug && isfinite(strip->consistency)) {
      MSG(_DEBUG3, silicon, "intra-consistency %.2f", strip->consistency);
      histogram_fill(debug->ratio_consistency, strip->consistency);
      histogram2_fill(debug->time_vs_consistency, strip->time,
                      strip->consistency);
    }
  }
}

/* helper structure for tear_region */
typedef struct strip_amps_s {
  double a0sum;
  double a1sum;
  double a2sum;
  double a2wsum; /* a2sum weighted with wire number */
  double noisesum;
} strip_amps_t;

/* helper function for tear_region */
/*__A_CONST__*/ static inline strip_amps_t
clear_strip_amps(void)
{
  strip_amps_t a;
  a.a0sum = 0;
  a.a1sum = 0;
  a.a2sum = 0;
  a.a2wsum = 0;
  a.noisesum = 0;
  return a;
}

/* helper function for tear_region */
/*__A_PURE__*/ static inline strip_amps_t
add_strip_amp(strip_amps_t a, const strip_t * s, const double * noises, double scale_b)
{
  a.a0sum += scale_b * s->a0f;
  a.a1sum += scale_b * s->a1f;
  a.a2sum += scale_b * s->a2f;
  a.a2wsum += scale_b * s->a2f * s->wire;
  a.noisesum += scale_b * noises[s->wire];
  return a;
}

/* helper function for tear_region */
/*__A_CONST__*/ static inline strip_amps_t
add_strip_amps(strip_amps_t a, strip_amps_t b, double scale_b)
{
  a.a0sum += scale_b * b.a0sum;
  a.a1sum += scale_b * b.a1sum;
  a.a2sum += scale_b * b.a2sum;
  a.a2wsum += scale_b * b.a2wsum;
  a.noisesum += scale_b * b.noisesum;
  return a;
}

/** create the real clusters by tearing apart (if necessary) the neighbours */
void
tear_region(const strip_t *const region_strip,
            const unsigned int region_size,
            const double *const consistency,
            const double *const noises,
            const time_reconstruction_option_t *const time_recon,
            const clusterization_option_t *const options,
            cluster_writer_function_t output_function,
            void *reserved,
            clusterization_debug_t *debug) {
				

  unsigned int        start = 0;
  strip_amps_t        amp_sum = clear_strip_amps();
  // switch cut determines if noise hits from 3 and 4 strip clusters are allowed or cut
  bool cut = true;
  
  for (unsigned int i = 0; i < region_size; i++) {
    const strip_t      *const that = &region_strip[i];



    MSG(_DEBUG3, silicon, "tear_point_sigmas %.2f, min_cluster_a2sum %i, "
        "max_error %.2f, enable_position %i",
        options->tear_point_sigmas, options->min_cluster_a2sum,
        options->max_error, options->enable_position);

    if (i < region_size - 1)
      MSG(_DEBUG3, silicon, "Consistency %.2f at pos %i", consistency[i], i);

    /* increase digit sums */
    amp_sum = add_strip_amp(amp_sum, that, noises, 1.0);

	int clustersize=-1;		


	/* shall we tear at this point? */
    if (i == region_size - 1
        || consistency[i] > options->tear_point_sigmas
        || options->tear_point_sigmas < 0.
        ) {
      /* tear: calculate cluster time */
      double               time, el = NAN, eh = NAN, position, ep=0;
      /* histogramming of a2sum of cluster */
      if (debug) {
        histogram_fill_raw(debug->cluster_a2sum, amp_sum.a2sum);
        histogram2_fill(debug->cluster_a2sum_vs_phase, debug->tcsphase,
                        amp_sum.a2sum);
      }


      /* cut on sum of a2 */
      if (amp_sum.a2sum < options->min_cluster_a2sum)
        goto end_of_tear;

      if (region_size == 1 || i==start) {
        /* no need to calculate time for size=1 clusters */
        time = that->time;
        el = that->el;
        eh = that->eh;
        position = that->wire;
        ep = 1 / sqrt(12.);
      } else {
        /* calculate time from digit sums for region_size>1 clusters */
        double               t0, el0 = NAN, eh0 = NAN, t1, el1 = NAN, eh1 = NAN;

        single_timing(&t0, &el0, &eh0, amp_sum.a0sum, amp_sum.a2sum, 0, amp_sum.noisesum, time_recon, 1);
        single_timing(&t1, &el1, &eh1, amp_sum.a1sum, amp_sum.a2sum, 1, amp_sum.noisesum, time_recon, 1);
        join_timing(&time, &el, &eh, t0, el0, eh0, t1, el1, eh1);

        position = amp_sum.a2wsum / amp_sum.a2sum;
        ep = 1 / sqrt(12.);
                
      }

      MSG(_DEBUG3, silicon, "Joined %i neighbours for time %.2f -%.2f/+%.2f",
          i + 1 - start, time, el, eh);

      /* check whether time reconstruction was successful */
      if (!isfinite(time))
        goto end_of_tear;

      /* histogramming of cluster time */
      if (debug) {
        debug->valid_clusters++;
        histogram_fill_raw(debug->cluster_size, i - start + 1);
        histogram_fill_raw(debug->cluster_position, cind_iround(position));
        histogram_fill(debug->cluster_time, time);
        histogram_fill(debug->cluster_time_error, -el);
        histogram_fill(debug->cluster_time_error, eh);
        histogram_fill(debug->cluster_pull, time / MEAN(el, eh));
        histogram2_fill(debug->cluster_err_vs_time, time, -el);
        histogram2_fill(debug->cluster_err_vs_time, time, eh);
      }

      /* error cut (for noise suppression) */
      if (options->max_error > 0
          && error_towards(time, el, eh, 0.) > options->max_error) {
        MSG(_DEBUG3, silicon, "Rejected at error cut");
        goto end_of_tear;
      }

      /* histogramming of cluster time and position after error cut */
      if (debug) {
        histogram_fill(debug->cluster_time_ec, time);
        histogram_fill(debug->cluster_time_error_ec, -el);
        histogram_fill(debug->cluster_time_error_ec, eh);
        histogram_fill_raw(debug->cluster_position_ec, cind_iround(position));
        histogram_fill_raw(debug->cluster_a2sum_ec, amp_sum.a2sum);

        /* do triggercut just for histogramming */
        if (debug->cluster_a2sum_ec_tc) {
          const double         sigma = 5.;

          if (-eh * sigma < time && time < el * sigma)
            histogram_fill_raw(debug->cluster_a2sum_ec_tc, amp_sum.a2sum);
        }
      }

      /* if cluster size is >= 3, split it up if there is a minimum between two maxima
       * (i.e. if there are several peaks) */
      if ( i - start + 1 >= 3 ) {
        unsigned int  j;
        unsigned int  last_max = start,  last_min = 0, clsize;
        double        slope = 0;
        double        t0, el0 = NAN, eh0 = NAN, t1, el1 = NAN, eh1 = NAN;
        const strip_t *strip_min = &region_strip[start];
        const strip_t *strip_max = &region_strip[start];

        strip_amps_t  sum_a = clear_strip_amps();
        strip_amps_t  sum_b = clear_strip_amps();

        sum_a = add_strip_amp(sum_a, &region_strip[start], noises, 1.0);

        for (j = start + 1; j <= i; j++) {
          const strip_t *const last = &region_strip[j-1];
          const strip_t *const cur  = &region_strip[j];
          const strip_t       *next = NULL;

          if (j!=i)
            next = &region_strip[j+1];

          if ((int)last->a2f != (int)cur->a2f) slope = cur->a2f - last->a2f;

          sum_a = add_strip_amp(sum_a, cur, noises, 1.0);

          /* if local minimum: remember it was there */
          if (next && slope<0 && cur->a2f < next->a2f &&
              !(last_min > last_max && strip_min->a2f < cur->a2f) &&
              (strip_max->a2f - cur->a2f) > (noises[cur->wire] + noises[strip_max->wire])){
            last_min = j;
            strip_min = cur;

            sum_b = sum_a;
            sum_a = clear_strip_amps();
          }

          /* if local maximum or last strip in big cluster */
          if (j==i || (slope>0 && cur->a2f > next->a2f &&
                       (last_min > last_max || last_min == 0))) {

            /* if there was a peak and following minimum that we did not write yet
             * (i.e. always, except for last strip in region and falling edge of a signal): */
            if (last_min > start &&
                (cur->a2f - strip_min->a2f) > (noises[cur->wire] + noises[strip_min->wire])) {
              int start_shift = 0;

              /* 3 cases right now:
               * if the minimum was right next to one of the peaks (and there were at
               * least 2 strips between the peaks)
               * give it exclusively to that one (makes two cases for left&right)
               * otherwise, just give half of the minimum to each peak */

              clsize = last_min - start + 1;

              if (j - last_max > 2 && j - last_min == 1) {
                /* need to subtract min wire from peak. Also cluster size is one less
                 * as the min is not included */
                sum_b = add_strip_amp(sum_b, &region_strip[last_min], noises, -1.0);
                sum_a = add_strip_amp(sum_a, &region_strip[last_min], noises, +1.0);

                clsize = last_min - start;
              }
              else if (!(j - last_max > 2 && last_min - last_max == 1)) {
                /* here, halve the min strip */
                /* noise is not halved: we are only "guessing" how to split anyways */
                sum_b = add_strip_amp(sum_b, &region_strip[last_min], noises, -0.5);
                sum_b.noisesum += 0.5 * noises[region_strip[last_min].wire];
                sum_a = add_strip_amp(sum_a, &region_strip[last_min], noises, +0.5);
                sum_a.noisesum += 0.5 * noises[region_strip[last_min].wire];
              }
              else
                start_shift = 1;

              /* third case: easiest, split is already where we wanted to have it:
               * just need to shift "start" */

              /* in any of the cases, make time, position and output_function for first of the two peaks
               * using the variables for time&error from the "big cluster" should not do any harm here... */
              single_timing(&t0, &el0, &eh0, sum_b.a0sum, sum_b.a2sum, 0, sum_b.noisesum, time_recon, 2);
              single_timing(&t1, &el1, &eh1, sum_b.a1sum, sum_b.a2sum, 1, sum_b.noisesum, time_recon, 2);
              join_timing(&time, &el, &eh, t0, el0, eh0, t1, el1, eh1);

              position = sum_b.a2wsum / sum_b.a2sum;

              MSG(_DEBUG2, silicon, "Tearing multi-peak: j=%i, start=%i, last_min=%i, "
                  "cluster size %i, start_shift+last_min=%i",
                  j, start, last_min, clsize,last_min + start_shift);

              output_function(reserved, time, el, eh, position, ep,
                              clsize, region_strip[start].wire, sum_b.a2sum);

              sum_b = clear_strip_amps();
              start = last_min + start_shift;
              last_max = j;
              strip_max = cur;

            }
            else if (last_min > start &&
                     (cur->a2f - strip_min->a2f) <= (noises[cur->wire] + noises[strip_min->wire])
                     && (strip_max->a2f - strip_min->a2f) < (noises[strip_max->wire] + noises[strip_min->wire])) {
              sum_a = add_strip_amps(sum_a, sum_b, 1);
              sum_b = clear_strip_amps();
            }
            else if (!((cur->a2f - strip_min->a2f) <= (noises[cur->wire] + noises[strip_min->wire]))) {
              last_max = j;
              strip_max = cur;
            }



          }
        }

        amp_sum = add_strip_amps(sum_a, sum_b, 1);

        single_timing(&t0, &el0, &eh0, amp_sum.a0sum, amp_sum.a2sum, 0, amp_sum.noisesum, time_recon, 2);
        single_timing(&t1, &el1, &eh1, amp_sum.a1sum, amp_sum.a2sum, 1, amp_sum.noisesum, time_recon, 2);
        join_timing(&time, &el, &eh, t0, el0, eh0, t1, el1, eh1);

        position = amp_sum.a2wsum / amp_sum.a2sum;
      }
	  // mark certain clusters for studies, to be removed before checkin. //mleeb
	  
	  
	  
	 
	  clustersize = i + 1 - start;
	  

	  
      // left side small  s|B|B
	  if( (i -start +1) ==3 && region_strip[start].a2f < 15 && region_strip[start+1].a2f>14 &&  region_strip[start].a2f < region_strip[start+2].a2f)
	  {
		if(cut) {
	  	double               t0, el0 = NAN, eh0 = NAN, t1, el1 = NAN, eh1 = NAN;
		amp_sum = clear_strip_amps();
		amp_sum = add_strip_amp(amp_sum, &region_strip[start+1] , noises, 1.0);
		amp_sum = add_strip_amp(amp_sum, &region_strip[start+2] , noises, 1.0);
	        single_timing(&t0, &el0, &eh0, amp_sum.a0sum, amp_sum.a2sum, 0, amp_sum.noisesum, time_recon, 1);
	        single_timing(&t1, &el1, &eh1, amp_sum.a1sum, amp_sum.a2sum, 1, amp_sum.noisesum, time_recon, 1);
        
	        join_timing(&time, &el, &eh, t0, el0, eh0, t1, el1, eh1);

	        position = amp_sum.a2wsum / amp_sum.a2sum;
		}
        
        ep = 666;
        
     
//        cout << " LEFTSIDE CUT " << endl;
//        cout << " old one " << endl;
//        for (int k=start; k<i+1; k++) {
//			if(k<i) 	cout << region_strip[k].a2f << "  |  ";
//			else 		cout << region_strip[k].a2f << endl;
//		}	
        
        if(cut) {
	
        start=i-1;
        clustersize= i+1 -start;
		}
		
//  		cout << " new one " << endl;
//		for (int k=start; k<start+clustersize; k++) {
//			if(k<start+clustersize) 	cout << region_strip[k].a2f << "  |  ";
//			else 		cout << region_strip[k].a2f << endl;
//		}
                	  
	  }
	
	  // right side	small    B|B|s
	  else  if( (i -start +1) ==3 && region_strip[start+2].a2f < 15 && region_strip[start+1].a2f>14 &&  region_strip[start+2].a2f < region_strip[start].a2f)
	  {
		if(cut) {  
		double               t0, el0 = NAN, eh0 = NAN, t1, el1 = NAN, eh1 = NAN;
		amp_sum = clear_strip_amps();
		amp_sum = add_strip_amp(amp_sum, &region_strip[start] , noises, 1.0);
		amp_sum = add_strip_amp(amp_sum, &region_strip[start+1] , noises, 1.0);
 	        single_timing(&t0, &el0, &eh0, amp_sum.a0sum, amp_sum.a2sum, 0, amp_sum.noisesum, time_recon, 1);
        	single_timing(&t1, &el1, &eh1, amp_sum.a1sum, amp_sum.a2sum, 1, amp_sum.noisesum, time_recon, 1);
        
  	        join_timing(&time, &el, &eh, t0, el0, eh0, t1, el1, eh1);
	
        	position = amp_sum.a2wsum / amp_sum.a2sum;
		}
		
        ep = 666;
        
//        cout << " RIGHTSIDE CUT " << endl;
//        cout << " old one " << endl;
//        for (int k=start; k<i+1; k++) {
//			if(k<i) 	cout << region_strip[k].a2f << "  |  ";
//			else 		cout << region_strip[k].a2f << endl;
//		}
        
        if(cut) {
        start=i-2;
        clustersize= i+1 -start -1;
		}
		
//		cout << " new one " << endl;
//		for (int k=start; k<start+clustersize; k++) {
//			if(k<start+clustersize) 	cout << region_strip[k].a2f << "  |  ";
//			else 		cout << region_strip[k].a2f << endl;
//		}
        
	  }	
	
		  
	  // peak   s|B|B|s
	  else if( (i -start +1) ==4 && region_strip[start].a2f < 15 && region_strip[start+3].a2f < 15 && ( region_strip[start+1].a2f > 14 ||  region_strip[start+2].a2f>14) )
          {
		if(cut) {  
		double               t0, el0 = NAN, eh0 = NAN, t1, el1 = NAN, eh1 = NAN;
		amp_sum = clear_strip_amps();
		amp_sum = add_strip_amp(amp_sum, &region_strip[start+1] , noises, 1.0);
		amp_sum = add_strip_amp(amp_sum, &region_strip[start+2] , noises, 1.0);
	        single_timing(&t0, &el0, &eh0, amp_sum.a0sum, amp_sum.a2sum, 0, amp_sum.noisesum, time_recon, 1);
	        single_timing(&t1, &el1, &eh1, amp_sum.a1sum, amp_sum.a2sum, 1, amp_sum.noisesum, time_recon, 1);
        
	        join_timing(&time, &el, &eh, t0, el0, eh0, t1, el1, eh1);

        	position = amp_sum.a2wsum / amp_sum.a2sum;
		}
		
		
        ep = 666;  
        
        
//        cout << " PEAK CUT " << endl;
//        cout << " old one " << endl;
//        for (int k=start; k<i+1; k++) {
//			if(k<i) 	cout << region_strip[k].a2f << "  |  ";
//			else 		cout << region_strip[k].a2f << endl;
//		} 
         
        if(cut) {        
        start = i-2;
		clustersize= i+1 -start -1;  
		}
		
	
//		cout << " new one " << endl;
//		for (int k=start; k<start+clustersize; k++) {
//			if(k<start+clustersize) 	cout << region_strip[k].a2f << "  |  ";
//			else 		cout << region_strip[k].a2f << endl;
//		}

	  }
	
	  //double left	    s|s|B|B
	  else if( (i -start +1) ==4 && region_strip[start].a2f < 15 && region_strip[start+1].a2f < 15 && ( region_strip[start+2].a2f > 14 ||  region_strip[start+3].a2f>14) )
      	  {
		if(cut) {  
	  	double               t0, el0 = NAN, eh0 = NAN, t1, el1 = NAN, eh1 = NAN;
		amp_sum = clear_strip_amps();
		amp_sum = add_strip_amp(amp_sum, &region_strip[start+2] , noises, 1.0);
		amp_sum = add_strip_amp(amp_sum, &region_strip[start+3] , noises, 1.0);
	        single_timing(&t0, &el0, &eh0, amp_sum.a0sum, amp_sum.a2sum, 0, amp_sum.noisesum, time_recon, 1);
	        single_timing(&t1, &el1, &eh1, amp_sum.a1sum, amp_sum.a2sum, 1, amp_sum.noisesum, time_recon, 1);
        
	        join_timing(&time, &el, &eh, t0, el0, eh0, t1, el1, eh1);

	        position = amp_sum.a2wsum / amp_sum.a2sum;
		}
		
        ep = 666;
        
                
//        cout << "DOUBLE LEFTSIDE CUT " << endl;
//        cout << " old one " << endl;
//        for (int k=start; k<i+1; k++) {
//			if(k<i) 	cout << region_strip[k].a2f << "  |  ";
//			else 		cout << region_strip[k].a2f << endl;
//		}
        
        if(cut) {
        start =i-1;	  
		clustersize= i+1 -start;  
		}
		
//		cout << " new one " << endl;
//		for (int k=start; k<start+clustersize; k++) {
//			if(k<start+clustersize) 	cout << region_strip[k].a2f << "  |  ";
//			else 		cout << region_strip[k].a2f << endl;
//		}
		
	  }
	  
	  // double right	B|B|s|s
	  else if( (i -start +1) ==4 && region_strip[start+2].a2f < 15 && region_strip[start+3].a2f < 15 && ( region_strip[start].a2f > 14 ||  region_strip[start+1].a2f>14) )
      	  {
		  
		if(cut) {  
	  	double               t0, el0 = NAN, eh0 = NAN, t1, el1 = NAN, eh1 = NAN;
		amp_sum = clear_strip_amps();
		amp_sum = add_strip_amp(amp_sum, &region_strip[start] , noises, 1.0);
		amp_sum = add_strip_amp(amp_sum, &region_strip[start+1] , noises, 1.0);
	        single_timing(&t0, &el0, &eh0, amp_sum.a0sum, amp_sum.a2sum, 0, amp_sum.noisesum, time_recon, 1);
	        single_timing(&t1, &el1, &eh1, amp_sum.a1sum, amp_sum.a2sum, 1, amp_sum.noisesum, time_recon, 1);
        
	        join_timing(&time, &el, &eh, t0, el0, eh0, t1, el1, eh1);

	        position = amp_sum.a2wsum / amp_sum.a2sum;
		}
        
        ep = 666;	
        
        
        
//        cout << "DOUBLE RIGHTSIDE CUT " << endl;
//	  cout << " old one " << endl;
//        for (int k=start; k<i+1; k++) {
//			if(k<i) 	cout << region_strip[k].a2f << "  |  ";
//			else 		cout << region_strip[k].a2f << endl;
//		}
        
        if(cut) {
        start=i-3; 
		clustersize= i+1 -start -2;  
		}
		
//		cout << " new one " << endl;
//		for (int k=start; k<start+clustersize; k++) {
//			if(k<start+clustersize) 	cout << region_strip[k].a2f << "  |  ";
//			else 		cout << region_strip[k].a2f << endl;
//		}
		
	  }
	  
	  // shift time depending on clustersize
      if (i+1-start == 1 )      { time-=time_recon->region_size_shift[0]; }
      else if (i+1-start ==2 )  { time-=time_recon->region_size_shift[1]; }
	  
      output_function(reserved, time, el, eh, position, ep,
                     clustersize, region_strip[start].wire, amp_sum.a2sum);

      /* setup/clear structures for next cluster */
	  end_of_tear:
      start = i + 1;
      amp_sum = clear_strip_amps();
	//  carryover=false;
    }                                     /* end if tear condition */
  }                                       /* end for i */
  
}


/** 
 * Clusterization of the hits of a detector plane.  The 'strip_hit' array
 * (length 'strip_count') is assumed to have been processed with
 * silicon_calculate_timings().  Further on it is required to be sorted by
 * wire number (ascending).  The output of the clusters, is handled by calling
 * a user-defined cluster_writer function.  The 'reserved' parameter is
 * specified to be passed through to the cluster_writer function and may be
 * NULL.  For 'debug' a NULL pointer may be passed to disable debugging.
 **/
void
silicon_clusterize_plane(strip_t * strip_hits,
                         const unsigned int strip_count,
                         const double *const noises,
                         const time_reconstruction_option_t * const time_recon,
                         const clusterization_option_t *const options,
                         cluster_writer_function_t output_function,
                         void *reserved,
                         clusterization_debug_t *debug) {

  /* assume this intra-consistency if we cannot calculate it
     because we just have one (valid) ratio */
  const double         one_ratio_penalty = options->tear_point_sigmas * 0.25;

  unsigned int        i;

  /* size of neighbourhood region */
  unsigned int        region_size = 0;

  /* value of "i" that points to first strip of neighbourhood */
  unsigned int                 region_start = 0;

  /* total amount of clusters */
  int                 cluster_count = 1;

  /* overall digit sums for whole regions, calculated only for debug/devel
     purpose */
  double               region_a0sum = 0.;
  double               region_a1sum = 0.;
  double               region_a2sum = 0.;
  double               region_noisesum = 0.;

  /* array to hold consistency of neighbouring strips inside a region */
  double              *consistency;

  if (strip_count == 0)
    return;

  F_CALLOC(consistency, double, strip_count - 1);

  /* iterate over strips and identify regions of neighbouring strips */
  for (i = 0; i < strip_count; i++) {
	  
	 

	  
    unsigned int        j;
    strip_t            *that;
    strip_t            *next;

    that = &strip_hits[i];
    /* skip strips with invalid time values */
    if (!isfinite(that->time)) {
      if (region_start == i)
        region_start++;
      continue;
    }
    

    /* only set 'next' if there is a next strip with valid time */
    next = NULL;
    for (j = i + 1; j < strip_count; j++) {
      strip_t            *next_candidate = &strip_hits[j];

      /* sanity check for duplicate wire numbers */
      if (that->wire == next_candidate->wire) {
        MSG_T(_ERROR, silicon, "Duplicate digit: wire %u, a0 %.0f : %.0f, "
              "a1 %.0f : %.0f, a2 %.0f : %.0f",
              that->wire,
              that->a0f, next_candidate->a0f,
              that->a1f, next_candidate->a1f,
              that->a2f, next_candidate->a2f);
        if (debug)
          debug->error_count++;
        /* skip the second wire */
        continue;
      }

      if (isfinite(next_candidate->time)) {
        next = next_candidate;
        break;
      }
      
      
      
    }

    /* increase region size counter */
    region_size++;
    /* increase digit sums for debug plots */
    if (debug) {
      region_a0sum += that->a0f;
      region_a1sum += that->a1f;
      region_a2sum += that->a2f;
      region_noisesum += noises[that->wire];
    }

    if (next && that->wire == next->wire - 1) {
      /* we are INSIDE a neighbourhood: calculate consistency of neighbours */
      /* inter-consistency between this and next strip */
      const double inter_consist =
          time_consistency(that->time, that->el, that->eh,
                           next->time, next->el, next->eh);

      /* intra-consistency of this strip */
      const double that_intra_consist =
        isfinite(that->consistency) ? fabsf(that->consistency) : one_ratio_penalty;
      /* intra-consistency of next strip */
      const double next_intra_consist =
        isfinite(next->consistency) ? fabsf(next->consistency) : one_ratio_penalty;

      MSG(_DEBUG3, silicon,
          "Wire %i-%i, this_consist %.2f, next_consist %.2f, inter-consist %.2f",
          that->wire, next->wire, that_intra_consist, next_intra_consist,
          inter_consist);

      consistency[region_size - 1] =
          fabsf(inter_consist) * 2 + that_intra_consist + next_intra_consist;

      /* for debugging purposes: calculate neighbour joint timing */
      if (debug) {
        double               joint_time, joint_el = NAN, joint_eh = NAN;

        join_timing(&joint_time, &joint_el, &joint_eh,
                    that->time, that->el, that->eh,
                    next->time, next->el, next->eh);
                    
        histogram_fill(debug->neighbour_time, joint_time);
        histogram_fill(debug->neighbour_time_error, -joint_el);
        histogram_fill(debug->neighbour_time_error, joint_eh);

        histogram_fill(debug->neighbour_consistency, inter_consist);
        histogram2_fill(debug->neighbour_t0_t1, that->time, next->time);
        histogram2_fill(debug->neighbour_time_vs_consistency, joint_time,
                        inter_consist);
      }

    } else {
      /* we are at NEIGHBOURHOOD REGION BOUNDARY */
      /* histogramming of size of neighbourhood regions */
      if (debug)
        histogram_fill_raw(debug->neighbour_size, region_size);

      MSG(_DEBUG3, silicon, "Neighbourhood boundary: computing tear points "
          "for %i neighbours starting at wire %i", region_size, region_start);

      tear_region(&strip_hits[region_start],
                  region_size, consistency,
                  noises,
                  time_recon, options,
                  output_function, reserved,
                  debug);


      /* starting the next cluster */
      region_a0sum = 0.;
      region_a1sum = 0.;
      region_a2sum = 0.;
      region_noisesum = 0.;
      region_size = 0;
      region_start = i + 1;
      cluster_count++;
    }
  }
  free(consistency);
  
}


#if 0
{
  /* add assumed trigger width of 0.5ns */
  /** \todo move to config file and/or find a cleaner solution */
  el = sqrt(SQR(el) + SQR(0.5));
  eh = sqrt(SQR(eh) + SQR(0.5));


}
#endif
