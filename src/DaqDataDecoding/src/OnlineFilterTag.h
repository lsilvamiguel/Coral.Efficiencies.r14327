#ifndef _FILTER_EVENT_H_
#define _FILTER_EVENT_H_

/**
 *
 * The type field of the filter information blocks describes which
 * kind of information is contained in each outblock.
 *
 **/
typedef enum {
  ot_decisions = 1,
  ot_calib_t0 = 2,
  ot_filter_version = 3,
  ot_filter_config = 4,
  ot_missing_events = 5,
  ot_eob_filter = 6
} outblock_type_t;

/**
 *
 * Each filter is identified externally by this ID, e.g. in connection
 * with the filter decision output block.
 *
 **/
typedef enum {
  gentoc = 1,
  gentriggertime = 2,
  bms = 3,
  scifi = 4,
  silicon = 5
} filter_id_t;

/*************************************************************
 * first define the structures which don't contain bitfields *
 *************************************************************/

struct output_missing_events_s {
  int32               missing;
  int32               discarded;
  int32               not_written;
  int32               non_physics;
  int32               too_big;
  int32               rejected;
};
typedef struct output_missing_events_s output_missing_events_t;

struct output_eob_filter_s {
  uint32              accepted;
  uint32              rejected;
  uint32              direct_writeout;
  uint32              max_filter;
  uint32              max_trigger;
};
typedef struct output_eob_filter_s output_eob_filter_t;

/** \def MAX_FILTERTIME
 *
 * There currently is a 14bit wide field to store the time each filter
 * has used in the output_decision_t. To avoid rollover problems, the
 * stored value is limited by this define.
 *
 **/
#define MAX_FILTERTIME 0x3fff

/***********************************************************************
 *                                                                     *
 * !!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT NOTICE !!!!!!!!!!!!!!!!!!!!!!!! *
 *                                                                     *
 * Unfortunately the order in which bit fields are filled is           *
 * implementation dependent, so that we have to provide two versions   *
 * of each declaration. The positive side is that the swapped          *
 * (non-x86) order is easier to read as the highest bit positions have *
 * to be put first.                                                    *
 *                                                                     *
 * !!!!!!!!! IF YOU CHANGE SOMETHING, DO SO IN BOTH VARIANTS !!!!!!!!! *
 *                                                                     *
 ***********************************************************************/
#if CONFIG_SWAP_BITFIELDS

/***********************************************************************/

/* first come the declarations for non-x86-like platforms */

struct outblock_s {
  uint32              monitoring:1;
  uint32              accepted:1;
  uint32              reserved:6;
  outblock_type_t     type:8;
  uint32              size:16;
};
typedef struct outblock_s outblock_t;

/**
 *
 * Holds the decision of one filter. filter_id is placed at the
 * beginning and followed by the reserved field so that it can easily
 * be expanded in the future. This means that the reserved bits MUST
 * be set to zero.
 *
 **/
struct output_decision_s {
  uint32              accepted:1;
  uint32              decided:1;
  uint32              time:14;
  uint32              reserved:8;
  filter_id_t         filter_id:8;
};
typedef struct output_decision_s output_decision_t;

/**********************************************************************/
#else /* CONFIG_SWAP_BITFIELDS */

/**********************************************************************/

struct outblock_s {
  uint32              size:16;
  outblock_type_t     type:8;
  uint32              reserved:6;
  uint32              accepted:1;
  uint32              monitoring:1;
};
typedef struct outblock_s outblock_t;

/**
 *
 * Holds the decision of one filter. filter_id is placed at the
 * beginning and followed by the reserved field so that it can easily
 * be expanded in the future. This means that the reserved bits MUST
 * be set to zero.
 *
 **/
struct output_decision_s {
  filter_id_t         filter_id:8;
  uint32              reserved:8;
  uint32              time:14;
  uint32              decided:1;
  uint32              accepted:1;
};
typedef struct output_decision_s output_decision_t;

/**********************************************************************/
#endif /* CONFIG_SWAP_BITFIELDS */

/**********************************************************************/


#endif
