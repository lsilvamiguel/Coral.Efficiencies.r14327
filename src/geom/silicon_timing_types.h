#ifndef _SILICON_TIMING_TYPES_H
#define _SILICON_TIMING_TYPES_H
#include <assert.h>

struct hist_s;
struct hist2_s;

/** structure for clustering */
struct strip_s {
  /** wire number */
  unsigned int        wire;
  /** digits (as read from APV) */
  double               a0f, a1f, a2f;
  /** joint reconstructed time (combined from time0 and time1) */
  double               time;
  /** asymmetric error of joint reconstructed time (downward) */
  double               el;
  /** asymmetric error of joint reconstructed time (upward) */
  double               eh;
  /** time reconstructed from ratio 0 */
  double               time0;
  /** time reconstructed from ratio 1 */
  double               time1;
  /** asymmetric error of time0 (downward) */
  double               el0;
  /** asymmetric error of time0 (upward) */
  double               eh0;
  /** asymmetric error of time1 (downward) */
  double               el1;
  /** asymmetric error of time1 (upward) */
  double               eh1;
  /** consistency between time0 and time1 expressed by their errors */
  double               consistency;
};
typedef struct strip_s strip_t;

/** parameter for Jan's function */
struct param_s {
  /** ignored */
  double               range_min, range_max;
  /** r0 */
  double               r0;
  /** t0 */
  double               t0;
  double               a;
  double               c;
  double               b;
  double               d;
  /** obsolete: array of precalculated function values */
  double              *precal;
  /** obsolete: width of one precal bin */
  double               width;
};
typedef struct param_s param_t;


struct eta_s{
	  double eta0;	// parameters for eta correction of cluster position. see phd thesis of dinkelbach
      double eta1;	// parameters for eta correction of cluster position. see phd thesis of dinkelbach
      double eta2;	// parameters for eta correction of cluster position. see phd thesis of dinkelbach
      double eta3;	// parameters for eta correction of cluster position. see phd thesis of dinkelbach
      double eta4;	// parameters for eta correction of cluster position. see phd thesis of dinkelbach
      double eta5;	// parameters for eta correction of cluster position. see phd thesis of dinkelbach
		
	};
typedef struct eta_s eta_t;



/** parameter for Micheal Wiesmann's amplitude parametrization */
struct amplitude_param_s {
  double               A0;
  double               t0;
  double               trise;
  double               tfall;
};
typedef struct amplitude_param_s amplitude_param_t;

struct calib_s {
  param_t             p0;
  param_t             p1;
};
typedef struct calib_s calib_t;

/** options/constants for timing reconstruction alogrithm */
struct time_reconstruction_option_s {
  /** TCS phase */
  double               tcsphase;

  /** Sampling interval of ADC in nanoseconds, should be set to 1000/38.8 ==
      25.77 */
  double               tcscycle;

  /** Minimal ratio cut: Ratios are required to be inside the interval
      [delta;r0-delta], where r0 denotes the maximum of the parametrisation
      function and is identical to the function parameter of the same name.
      A setting of 0.01 has been tested successfully. */
  double               delta;

  /** Assumed maximal decay time of the silicon signal (in ns).  300 seem
      to work well. */
  double               decay_time;

  /** Minimal downwards error of early particle ("late timing").  30 seem
      to work well. */
  double               min_early_error;

  /** Shift values for cluster sizes */
  double             region_size_shift[2];

  /** Parameters to Jan's approximation function for ratio 0 and ratio 1. */
  param_t            *param[2];
};
typedef struct time_reconstruction_option_s time_reconstruction_option_t;


/** options/constants for clusterization function */
struct clusterization_option_s {
  /** tear point in units of 4*sigma, set to -1 to tear at every point (useful for testing) */
  double               tear_point_sigmas;
  /** CUT on minimal a2sum of a cluster, set to zero to disable */
  unsigned int        min_cluster_a2sum;
  /** CUT on maximal error towards trigger (ns), set to -1 to disable */
  double               max_error;
  /** when not zero, cluster position is returned, too */
  unsigned int        enable_position;
};
typedef struct clusterization_option_s clusterization_option_t;


/** Debug structure with histogramming for clusterization. This is intended for
    use with Cinderella. Other applications (CORAL, etc.) probably want to
    supply a NULL pointer instead. */
struct clusterization_debug_s {
  /** counter for valid strips (only to be increased, never reset) */
  unsigned int        valid_strips;

  /** counter for valid clusters (only to be increased, never reset) */
  unsigned int        valid_clusters;

  /** error counter (only to be increased, never reset) */
  unsigned int        error_count;

  /** just for histogramming purposes */
  double               tcsphase;

  /* plethora of histograms */
  /* checked */
  struct hist_s      *time;
  struct hist_s      *clocktime;
  struct hist_s      *time_error;
  struct hist_s      *time0;
  struct hist_s      *clocktime0;
  struct hist_s      *time0_error;
  struct hist_s      *time1;
  struct hist_s      *clocktime1;
  struct hist_s      *time1_error;

  struct hist_s      *pull, *cluster_pull, *pull0, *pull1;

  struct hist2_s     *err_vs_time;
  struct hist2_s     *cluster_err_vs_time;
  struct hist2_s     *err0_vs_time;
  struct hist2_s     *err1_vs_time;

  struct hist2_s     *r0_r1;
  struct hist2_s     *r0_phase;
  struct hist2_s     *r1_phase;
  struct hist2_s     *r0_vs_clocktime;
  struct hist2_s     *r1_vs_clocktime;

  struct hist_s      *neighbour_size;

  struct hist_s      *cluster_size;
  struct hist_s      *cluster_time;
  struct hist_s      *cluster_time_error;

  struct hist_s      *cluster_time_ec;
  struct hist_s      *cluster_time_error_ec;

  struct hist_s      *a2;
  struct hist_s      *cluster_a2sum;
  struct hist_s      *cluster_a2sum_ec;
  struct hist_s      *cluster_a2sum_ec_tc;

  struct hist2_s     *cluster_a2sum_vs_phase;

  struct hist_s      *cluster_position;
  struct hist_s      *cluster_position_ec;

  /* unchecked */
  struct hist_s      *neighbour_time;
  struct hist_s      *neighbour_time_error;
  struct hist_s      *neighbour_consistency;
  struct hist_s      *ratio_consistency;

  struct hist2_s     *err0_ratio;
  struct hist2_s     *err1_ratio;
  struct hist2_s     *neighbour_t0_t1;
  struct hist2_s     *neighbour_time_vs_consistency;
  struct hist2_s     *time_vs_consistency;
};
typedef struct clusterization_debug_s clusterization_debug_t;



/** of value 'x' with asymmetric errors 'xel' and 'xeh' it is returned that
    error which is in direction of value 'y' */
#define error_towards(x, xel, xeh, y) ((x)<(y)?(xeh):(xel))


/** calculate absolute error of amplitude ratio ax/a2 using Gauss formula and assuming
    'noise' absolute error in the individual digits */
__A_CONST__ static inline double
ratio_error(double ax, double a2, double noise)
{
  const double         sqr_ax = ax * ax;
  const double         sqr_a2 = a2 * a2;

  /* The digitization error is 1/sqrt(12). */
  return (noise + 1 / sqrt(12.)) * sqrt(sqr_ax + sqr_a2) / sqr_a2;
}

/** deviation of two time points in units of their own errors (sigmas) */
__A_CONST__ static inline double
time_consistency(const double t0, const double el0, const double eh0,
                 const double t1, const double el1, const double eh1)
{

  /* return NaN if consistency cannot be computed */
  if (!std::isfinite(t0) || !std::isfinite(t1))
    return NAN;

  {
    const double         delta = t1 - t0;
    const double         e0 = error_towards(t0, el0, eh0, t1);
    const double         e1 = error_towards(t1, el1, eh1, t0);
    const double         consistency = delta / sqrt(SQR(e0) + SQR(e1));

    return consistency;
  }
}


/** randomize data to let histograms look better */
static inline double
digit_randomize(double a)
{
  double               R = rand();

  R /= RAND_MAX;
  R -= 0.5;

  return a + R;
}


/* Jan's functions */
/**
 p[0]: r0
 p[1]: t0
 p[2]: a
 p[3]: c
 p[4]: b
 p[5]: d

double eS(double x, double *p) {
  double xx=x-p[1];
  return 0.5*( (p[2]+p[3]) *       xx  + 
	       (p[2]-p[3]) * (sqrt(xx*xx+p[4]*p[4])-fabs(p[4])) ) + p[5];
}
double eSi(double s, double *p) {
  double apc=p[2]+p[3];
  double amc=p[2]-p[3];
  double f=s-p[5]+0.5*p[4]*amc;
  double ac=p[2]*p[3];
  double x=0.5/ac*(apc*f-amc*sqrt(f*f+p[4]*p[4]*ac));
  return x+p[1];
}
double APV_r(double x, double *p) {
  return p[0]*exp(-exp(-eS(x, p)));
}
double APV_t(double x, double *p) {
  return eSi(-log(-log(x/p[0])), p);
}


delta=25ns
N=amount of iterations
ratio=ratio
pp=parameters corresponding to ratio

double APV_SI_ampl(double t, int ratio,
		       double *pp,
		       double delta, int N) {
  t+=delta;
  double ampl = 1;
  double tau2 = delta/log(pp[0]);
  for (int i=0; i<N; i++) {
    ampl *= (ratio==0) ? APV_r(t+double(i)*delta, pp)
                       : APV_r(t+2*double(i)*delta, pp);
  }
  return ampl * exp(-(t+double(N)*delta)/tau2);
}
*/

__A_PURE__ static inline double
eS(const double t, const param_t *p)
{
  const double         apc = p->a + p->c;
  const double         amc = p->a - p->c;
  const double         bb = p->b * p->b;

  return 0.5 * (apc * t + amc * (sqrt(t * t + bb) - fabs(p->b))) + p->d;
}

__A_PURE__ static inline double
APV_r(const double t, const param_t *p)
{
  return p->r0 * exp(-exp(-eS(-t - p->t0, p)));
}

__A_PURE__ static inline double
eSi(const double s, const param_t *p)
{
  const double         apc = p->a + p->c;
  const double         amc = p->a - p->c;
  const double         f = s - p->d + 0.5 * p->b * amc;
  const double         ac = p->a * p->c;
  const double         x = 0.5 / ac
      * (apc * f - amc * sqrt(f * f + p->b * p->b * ac));
  return x + p->t0;
}

__A_PURE__ static inline double
APV_t(const double r, const param_t *p)
{
  return -eSi(-log(-log(r / p->r0)), p);
}


__A_PURE__ static inline double
APV_SI_ampl(double time, const unsigned int ratio,
            const param_t *p, double delta_, unsigned int N)
{
  unsigned int        i;
  const double         delta = delta_ * (ratio + 1);
  const double         tau2 = delta / log(p->r0);
  double               ampl = 1;

  for (i = 0; i < N; i++)
    ampl *= APV_r(-(time + i * delta), p);

  return ampl * exp(-(time + N * delta) / tau2);
}

__A_PURE__ static inline double
amplitude(const double time, const amplitude_param_t *p)
{
  const double         tbar = -(time - p->t0);
  const double         amp =
      p->A0 * exp(-tbar / p->tfall) * (1 - exp(-tbar / p->trise));
  return (amp > 0.) ? amp : 0.;
}

/** Either runs the function directly or takes a value from the lookup table, 
 depending on configuration. */
#define APV_TIME(x,p) \
  (p->precal \
  ? p->precal[cind_iround(((x)-p->range_min)/p->width)] \
  : APV_t(x, p) )


/** 
 * Calculates time and asymmetric error of a single ratio in the variables
 * 'time', 'err_up' and 'err_down'.  Input is ratio 'r', its error 'r_err',
 * the parameter array 'p' and the 'tcsphase'.  On error 'time' is set NaN.
 * When 'time' is not NaN, it is guaranteed that 'el' and 'eh' also are finite
 * and not zero.
**/
static inline void
single_timing(double *time, double *err_down, double *err_up,
              const double ax, const double a2, const unsigned int ratio,
              const double noise,
              const time_reconstruction_option_t *const options,
              const int region_size)
{
  // size dependent timing shift only defined for size 1 and 2
  assert(region_size <= 2);
  assert(region_size >= 1);
  if (!std::isnormal(a2)) {
    *time = NAN;
    return;
  }

  {
    const param_t      *p = options->param[ratio];
    const double         r = ax / a2;
    const double         r_err = ratio_error(ax, a2, noise);

    MSG(_DEBUG3, silicon, "single_timing: \tinput: "
        "a%i=%.2f, a2=%.2f, r%i=%.2f +/- %.3f (noise=%.2f)",
        ratio, ax, a2, ratio, r, r_err, noise);
    MSG(_DEBUG3, silicon, "single_timing: \tparams: "
        "r0=%.2f, t0=%.2f, a=%.2f, b=%.2f, c=%.2f, d=%.2f, %s precal.",
        p->r0, p->t0, p->a, p->b, p->c, p->d, (p->precal ? "with" : "no"));

    /* calculate time */
    if (options->delta < r && r < p->r0 - options->delta)
      *time = APV_TIME(r, p);
    else {
      *time = NAN;
      return;
    }

    /* check sanity of calculated time */
    if (!std::isfinite(*time)) {
      MSG_T(_ERROR, silicon, "Time finity assertion failed: "
	    "r = %.f ==> time = %.2fns", r, *time);
      return;
    }

    /* calculate error */
    {
      const double         rl = r - r_err;
      const double         rh = r + r_err;

      /* since ratio function is decreasing, we need to evaluate it at (r -
         r_err) to get upwards error */
      if (rl < options->delta)
        /* if 'rl' is out of range, the upwards error is determined by the
           length of the TCS cycle: if the particle had been later, the ratio 
           would have been larger and we would have noticed it */
        /* But if the failure to calculate the upwards error from the ratio
           function stems from small amplitudes, the error cannot be smaller
           than the difference to the time of infinitely small ratio:
           APV_TIME(delta) */
        *err_up = MAX(options->tcscycle * (ratio + 1) / sqrt(12.),
                      APV_TIME(options->delta, p) - *time);
      else {
        MSG(_DEBUG3, silicon, "rl: %.3f, t(rl): %.5f", rl, APV_TIME(rl, p));
        *err_up = APV_TIME(rl, p) - *time;
      }

      /* at the opposite we need to evaluate the ratio function at (r +
         r_err) to get downwards error */
      if (rh > p->r0 - options->delta)
        /* if 'rh' is out of range, the downwards error is determined by the
           decay time of the signal */
        *err_down =
            MAX(options->decay_time + *time,
                options->min_early_error * sqrt(12.)) / sqrt(12.);
      else {
        MSG(_DEBUG3, silicon, "rh: %.3f, t(rh): %.5f", rh, APV_TIME(rh, p));
        *err_down = *time - APV_TIME(rh, p);
      }
    }

    /* add cluster size dependent time shift */
    *time += options->region_size_shift[region_size - 1];

    /* add tcsphase */
    *time += options->tcsphase;

    MSG(_DEBUG3, silicon,
        "single_timing: r%i=%.2f, t=%.2fns -%.2fns/+%.2fns, tcsphase %.2fns",
        ratio, r, *time, *err_down, *err_up, options->tcsphase);

    /* check for infinity or zeroness of error (should not be necessary) */
    if (!std::isnormal(*err_up) || !std::isnormal(*err_down) || *err_down < 0.
        || *err_up < 0.) {
      MSG_T(_ERROR, silicon,
            "Time error assertion failed: noise = %.1f, "
            "r = %.f / %.f = %.2f +/- %.2f ==> time = %.2fns -%.2fns/+%.2fns",
            noise, ax, a2, r, r_err, *time, *err_down, *err_up);
      *time = NAN;
      return;
    }
  }
  return;
}


/**
 * Joins two time values including asymmetric error. Returns NaN for 'time'
 * when both input timings are NaN or errors are zero. If 'time' is finite, it
 * is guaranteed for 'el' and 'eh' also to be finite and not zero.
 **/
static inline void
join_timing(double *time, double *el, double *eh,
            const double t0, const double el0, const double eh0,
            const double t1, const double el1, const double eh1)
{

  if (std::isfinite(t0) && std::isfinite(t1)) {
    /* calculate joint time and error using Gauss formula */
    const double         e0 = error_towards(t0, el0, eh0, t1);
    const double         e1 = error_towards(t1, el1, eh1, t0);
    const double         e0_2 = 1 / SQR(e0);
    const double         e1_2 = 1 / SQR(e1);
    const double         el0_2 = 1 / SQR(el0);
    const double         eh0_2 = 1 / SQR(eh0);
    const double         el1_2 = 1 / SQR(el1);
    const double         eh1_2 = 1 / SQR(eh1);

    *time = (t0 * e0_2 + t1 * e1_2) / (e0_2 + e1_2);
    *el = 1 / sqrt(el0_2 + el1_2);
    *eh = 1 / sqrt(eh0_2 + eh1_2);
  } else if (std::isfinite(t0)) {
    *time = t0;
    *el = el0;
    *eh = eh0;
  } else if (std::isfinite(t1)) {
    *time = t1;
    *el = el1;
    *eh = eh1;
  } else {
    *time = NAN;
    return;
  }

  if (!std::isnormal(*el) || !std::isnormal(*eh))
    *time = NAN;
}

#ifdef __cplusplus
  extern "C" {
#endif

/** 
 * Callback function to store one silicon cluster.  'reserved' is custom
 * identifier passed on to the callback function by
 * silicon_clusterize_plane().  The time of the cluster is returned via 'time'
 * with asymmetric errors 'el' and 'eh'.  The cluster position is returned via
 * 'postion' and error 'ep' in units of wire numbers.  Note that 'ep' is not
 * implemented reasonably, yet.
 **/
typedef void        (*cluster_writer_function_t) (void *reserved, double time,
                                                  double el, double eh,
                                                  double position, double ep,
                                                  unsigned int cluster_size,
                                                  unsigned int first_strip,
						  double a2sum);

void
                    silicon_calculate_timings(strip_t * strip_hits, const unsigned int strip_count,
                                              const double *const noises, const
                                              time_reconstruction_option_t
                                              *const time_recon,
                                              clusterization_debug_t *debug);

void
                    silicon_clusterize_plane(strip_t * strip_hits, const unsigned int strip_count,
                                             const double *const noises, const
                                             time_reconstruction_option_t
                                             *const time_recon, const clusterization_option_t
                                             *const options, cluster_writer_function_t
                                             output_function, void *reserved,
                                             clusterization_debug_t *debug);

#ifdef __cplusplus
}
#endif

#endif
