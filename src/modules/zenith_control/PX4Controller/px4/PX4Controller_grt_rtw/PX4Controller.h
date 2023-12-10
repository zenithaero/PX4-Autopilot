/*
 * PX4Controller.h
 *
 * Code generation for model "PX4Controller".
 *
 * Model version              : 1.278
 * Simulink Coder version : 9.7 (R2022a) 13-Nov-2021
 * C++ source code generated on : Sun Dec 10 03:06:45 2023
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_PX4Controller_h_
#define RTW_HEADER_PX4Controller_h_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "PX4Controller_types.h"

extern "C" {

#include "rtGetInf.h"

}
#include <cstring>

extern "C" {

#include "rt_nonfinite.h"

}
#include "zero_crossing_types.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
#define rtmGetContStateDisabled(rtm)   ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
#define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
#define rtmGetContStates(rtm)          ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
#define rtmSetContStates(rtm, val)     ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
#define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
#define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetIntgData
#define rtmGetIntgData(rtm)            ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
#define rtmSetIntgData(rtm, val)       ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
#define rtmGetOdeF(rtm)                ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
#define rtmSetOdeF(rtm, val)           ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
#define rtmGetOdeY(rtm)                ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
#define rtmSetOdeY(rtm, val)           ((rtm)->odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
#define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
#define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
#define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
#define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
#define rtmGetZCCacheNeedsReset(rtm)   ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
#define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
#define rtmGetdX(rtm)                  ((rtm)->derivs)
#endif

#ifndef rtmSetdX
#define rtmSetdX(rtm, val)             ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

/* Block signals (default storage) */
struct B_PX4Controller_T {
  real_T Gain7;                        /* '<S75>/Gain7' */
  real_T steErr;                       /* '<S77>/Sum' */
  real_T VectorConcatenate[6];         /* '<S13>/Vector Concatenate' */
  real_T Gain[6];                      /* '<S15>/Gain' */
  real_T Switch1[4];                   /* '<S11>/Switch1' */
  real_T Switch2[6];                   /* '<S11>/Switch2' */
  real_T Switch3[2];                   /* '<S11>/Switch3' */
  real_T Switch4[2];                   /* '<S11>/Switch4' */
  real_T Product[2];                   /* '<S34>/Product' */
  real_T Product_g;                    /* '<S79>/Product' */
  real_T steDotErr;                    /* '<S77>/Sum1' */
  boolean_T OR;                        /* '<S73>/OR' */
  boolean_T OR_k;                      /* '<S26>/OR' */
};

/* Block states (default storage) for system '<Root>' */
struct DW_PX4Controller_T {
  real_T Memory_1_PreviousInput[4];    /* '<S1>/Memory' */
  real_T Memory_2_PreviousInput[6];    /* '<S1>/Memory' */
  real_T Memory_3_PreviousInput[2];    /* '<S1>/Memory' */
  real_T Memory_4_PreviousInput[2];    /* '<S1>/Memory' */
  uint32_T Prelookup4_DWORK1;          /* '<S8>/Prelookup4' */
  uint32_T Prelookup6_DWORK1;          /* '<S8>/Prelookup6' */
  uint32_T Prelookup1_DWORK1;          /* '<S8>/Prelookup1' */
  uint32_T Prelookup2_DWORK1;          /* '<S8>/Prelookup2' */
  uint32_T Prelookup3_DWORK1;          /* '<S8>/Prelookup3' */
};

/* Continuous states (default storage) */
struct X_PX4Controller_T {
  real_T Integrator1_CSTATE;           /* '<S76>/Integrator1' */
  real_T Integrator_CSTATE[2];         /* '<S3>/Integrator' */
  real_T Lowpass_CSTATE;               /* '<S77>/Low pass' */
  real_T Integrator_CSTATE_f;          /* '<S77>/Integrator' */
};

/* State derivatives (default storage) */
struct XDot_PX4Controller_T {
  real_T Integrator1_CSTATE;           /* '<S76>/Integrator1' */
  real_T Integrator_CSTATE[2];         /* '<S3>/Integrator' */
  real_T Lowpass_CSTATE;               /* '<S77>/Low pass' */
  real_T Integrator_CSTATE_f;          /* '<S77>/Integrator' */
};

/* State disabled  */
struct XDis_PX4Controller_T {
  boolean_T Integrator1_CSTATE;        /* '<S76>/Integrator1' */
  boolean_T Integrator_CSTATE[2];      /* '<S3>/Integrator' */
  boolean_T Lowpass_CSTATE;            /* '<S77>/Low pass' */
  boolean_T Integrator_CSTATE_f;       /* '<S77>/Integrator' */
};

/* Zero-crossing (trigger) state */
struct PrevZCX_PX4Controller_T {
  ZCSigState Integrator1_Reset_ZCE;    /* '<S76>/Integrator1' */
  ZCSigState Integrator_Reset_ZCE;     /* '<S3>/Integrator' */
  ZCSigState Integrator_Reset_ZCE_l;   /* '<S77>/Integrator' */
};

#ifndef ODE4_INTG
#define ODE4_INTG

/* ODE4 Integration Data */
struct ODE4_IntgData {
  real_T *y;                           /* output */
  real_T *f[4];                        /* derivatives */
};

#endif

/* External inputs (root inport signals with default storage) */
struct ExtU_PX4Controller_T {
  CmdBus CmdBusIn;                     /* '<Root>/CmdBusIn' */
  StateBus StateBus_m;                 /* '<Root>/StateBus' */
};

/* External outputs (root outports fed by signals with default storage) */
struct ExtY_PX4Controller_T {
  CmdBus CmdBusOut;                    /* '<Root>/CmdBusOut' */
  ActBus ActBus_e;                     /* '<Root>/ActBus' */
};

/* Parameters (default storage) */
struct P_PX4Controller_T_ {
  real_T SubsystemReference_fmMask[6];
                                    /* Mask Parameter: SubsystemReference_fmMask
                                     * Referenced by:
                                     *   '<S15>/Gain'
                                     *   '<S15>/Gain1'
                                     */
  real_T SubsystemReference1_fmMask[6];
                                   /* Mask Parameter: SubsystemReference1_fmMask
                                    * Referenced by:
                                    *   '<S16>/Gain'
                                    *   '<S16>/Gain1'
                                    */
  real_T SubsystemReference2_fmMask[6];
                                   /* Mask Parameter: SubsystemReference2_fmMask
                                    * Referenced by:
                                    *   '<S17>/Gain'
                                    *   '<S17>/Gain1'
                                    */
  real_T SubsystemReference3_table[36];
                                    /* Mask Parameter: SubsystemReference3_table
                                     * Referenced by: '<S31>/Interpolation Using Prelookup'
                                     */
  real_T SubsystemReference3_table_m[36];
                                  /* Mask Parameter: SubsystemReference3_table_m
                                   * Referenced by: '<S71>/Interpolation Using Prelookup'
                                   */
  real_T SubsystemReference3_table_p[48];
                                  /* Mask Parameter: SubsystemReference3_table_p
                                   * Referenced by: '<S19>/Interpolation Using Prelookup'
                                   */
  real_T SubsystemReference4_table[72];
                                    /* Mask Parameter: SubsystemReference4_table
                                     * Referenced by: '<S20>/Interpolation Using Prelookup'
                                     */
  real_T SubsystemReference5_table[24];
                                    /* Mask Parameter: SubsystemReference5_table
                                     * Referenced by: '<S21>/Interpolation Using Prelookup'
                                     */
  real_T SubsystemReference6_table[24];
                                    /* Mask Parameter: SubsystemReference6_table
                                     * Referenced by: '<S22>/Interpolation Using Prelookup'
                                     */
  real_T Gain_Gain;                    /* Expression: 10
                                        * Referenced by: '<S12>/Gain'
                                        */
  real_T Gain_Gain_d;               /* Expression: config.gains.allocation.maxFx
                                     * Referenced by: '<S70>/Gain'
                                     */
  real_T Gain_Gain_g;               /* Expression: config.gains.allocation.maxFx
                                     * Referenced by: '<S27>/Gain'
                                     */
  real_T Saturation_UpperSat;          /* Expression: inf
                                        * Referenced by: '<S77>/Saturation'
                                        */
  real_T Saturation_LowerSat;          /* Expression: 0.1
                                        * Referenced by: '<S77>/Saturation'
                                        */
  real_T InterpolationUsingPrelookup_Tab[4];/* Expression: squeeze(table)
                                             * Referenced by: '<S81>/Interpolation Using Prelookup'
                                             */
  real_T InterpolationUsingPrelookup_T_n[4];/* Expression: squeeze(table)
                                             * Referenced by: '<S83>/Interpolation Using Prelookup'
                                             */
  real_T InterpolationUsingPrelookup_T_c[4];/* Expression: squeeze(table)
                                             * Referenced by: '<S82>/Interpolation Using Prelookup'
                                             */
  real_T M_Gain;                       /* Expression: config.aircraft.M
                                        * Referenced by: '<S77>/M'
                                        */
  real_T Switch1_Threshold;            /* Expression: 0
                                        * Referenced by: '<S27>/Switch1'
                                        */
  real_T Switch1_Threshold_g;          /* Expression: 0
                                        * Referenced by: '<S70>/Switch1'
                                        */
  real_T Gain1_Gain;                   /* Expression: 10
                                        * Referenced by: '<S12>/Gain1'
                                        */
  real_T Gain2_Gain;                   /* Expression: -10
                                        * Referenced by: '<S12>/Gain2'
                                        */
  real_T Gain3_Gain;                   /* Expression: 5
                                        * Referenced by: '<S12>/Gain3'
                                        */
  real_T Gain1_Gain_f;                 /* Computed Parameter: Gain1_Gain_f
                                        * Referenced by: '<S7>/Gain1'
                                        */
  real_T Gain2_Gain_k;                 /* Expression: -1
                                        * Referenced by: '<S7>/Gain2'
                                        */
  real_T Gain3_Gain_c;                 /* Expression: 0
                                        * Referenced by: '<S7>/Gain3'
                                        */
  real_T InterpolationUsingPrelookup__cs[12];/* Expression: squeeze(table)
                                              * Referenced by: '<S72>/Interpolation Using Prelookup'
                                              */
  real_T Gain_Gain_a;                  /* Computed Parameter: Gain_Gain_a
                                        * Referenced by: '<S7>/Gain'
                                        */
  real_T Saturation_UpperSat_c[3];     /* Expression: config.gains.rate.max
                                        * Referenced by: '<S7>/Saturation'
                                        */
  real_T Saturation_LowerSat_d[3];     /* Expression: -config.gains.rate.max
                                        * Referenced by: '<S7>/Saturation'
                                        */
  real_T _Gain;                        /* Computed Parameter: _Gain
                                        * Referenced by: '<S7>/ '
                                        */
  real_T _Gain_m;                      /* Computed Parameter: _Gain_m
                                        * Referenced by: '<S7>/   '
                                        */
  real_T Saturation1_UpperSat;       /* Expression: config.gains.attitude.max(2)
                                      * Referenced by: '<S74>/Saturation1'
                                      */
  real_T Saturation1_LowerSat;      /* Expression: -config.gains.attitude.max(2)
                                     * Referenced by: '<S74>/Saturation1'
                                     */
  real_T Saturation1_UpperSat_g;     /* Expression: config.gains.attitude.max(1)
                                      * Referenced by: '<S68>/Saturation1'
                                      */
  real_T Saturation1_LowerSat_b;    /* Expression: -config.gains.attitude.max(1)
                                     * Referenced by: '<S68>/Saturation1'
                                     */
  real_T Saturation_UpperSat_p[3];     /* Expression: config.gains.rate.max
                                        * Referenced by: '<S28>/Saturation'
                                        */
  real_T Saturation_LowerSat_c[3];     /* Expression: -config.gains.rate.max
                                        * Referenced by: '<S28>/Saturation'
                                        */
  real_T Gain1_Gain_p;               /* Expression: config.gains.attitude.max(1)
                                      * Referenced by: '<S27>/Gain1'
                                      */
  real_T Gain2_Gain_p;              /* Expression: -config.gains.attitude.max(2)
                                     * Referenced by: '<S27>/Gain2'
                                     */
  real_T Gain3_Gain_b;               /* Expression: config.gains.attitude.max(3)
                                      * Referenced by: '<S27>/Gain3'
                                      */
  real_T Gain_Gain_a0;                 /* Expression: sqrt(2)
                                        * Referenced by: '<S66>/Gain'
                                        */
  real_T Gain_Gain_k;                  /* Expression: 0.5 * 1.225
                                        * Referenced by: '<S8>/Gain'
                                        */
  real_T Prelookup4_BreakpointsData[2];
                       /* Expression: 1/2 * 1.224 * config.gains.points.tas .^ 2
                        * Referenced by: '<S8>/Prelookup4'
                        */
  real_T Prelookup6_BreakpointsData[2];/* Expression: config.trim.points.h
                                        * Referenced by: '<S8>/Prelookup6'
                                        */
  real_T InterpolationUsingPrelookup_T_h[4];/* Expression: squeeze(table)
                                             * Referenced by: '<S69>/Interpolation Using Prelookup'
                                             */
  real_T Saturation_UpperSat_m;        /* Expression: inf
                                        * Referenced by: '<S66>/Saturation'
                                        */
  real_T Saturation_LowerSat_k;        /* Expression: config.gains.L1.minL1
                                        * Referenced by: '<S66>/Saturation'
                                        */
  real_T Gain_Gain_do;                 /* Expression: -1
                                        * Referenced by: '<S67>/Gain'
                                        */
  real_T Saturation_UpperSat_a;        /* Expression: inf
                                        * Referenced by: '<S67>/Saturation'
                                        */
  real_T Saturation_LowerSat_e;        /* Expression: 0
                                        * Referenced by: '<S67>/Saturation'
                                        */
  real_T Gain_Gain_dn;                 /* Expression: 2
                                        * Referenced by: '<S65>/Gain'
                                        */
  real_T uG_Gain;                      /* Expression: 1/config.gains.tecs.G
                                        * Referenced by: '<S65>/1//G'
                                        */
  real_T Gain2_Gain_h;                 /* Computed Parameter: Gain2_Gain_h
                                        * Referenced by: '<S6>/Gain2'
                                        */
  real_T Integrator1_IC;               /* Expression: 0
                                        * Referenced by: '<S76>/Integrator1'
                                        */
  real_T Constant_Value;               /* Expression: config.gains.tecs.balSpeW
                                        * Referenced by: '<S76>/Constant'
                                        */
  real_T Gain_Gain_f;                  /* Computed Parameter: Gain_Gain_f
                                        * Referenced by: '<S75>/Gain'
                                        */
  real_T Multiply_Gain;                /* Expression: config.gains.tecs.G
                                        * Referenced by: '<S75>/Multiply'
                                        */
  real_T Constant1_Value;              /* Expression: config.gains.tecs.balSkeW
                                        * Referenced by: '<S76>/Constant1'
                                        */
  real_T Gain2_Gain_o;                 /* Computed Parameter: Gain2_Gain_o
                                        * Referenced by: '<S75>/Gain2'
                                        */
  real_T Multiply5_Gain;               /* Expression: 0.5
                                        * Referenced by: '<S75>/Multiply5'
                                        */
  real_T Gain2_Gain_c;                 /* Expression: -1
                                        * Referenced by: '<S37>/Gain2'
                                        */
  real_T Gain1_Gain_h;                 /* Computed Parameter: Gain1_Gain_h
                                        * Referenced by: '<S75>/Gain1'
                                        */
  real_T Multiply1_Gain;               /* Expression: config.gains.tecs.G
                                        * Referenced by: '<S75>/Multiply1'
                                        */
  real_T Gain3_Gain_j;                 /* Computed Parameter: Gain3_Gain_j
                                        * Referenced by: '<S75>/Gain3'
                                        */
  real_T Multiply3_Gain;               /* Expression: 0.5
                                        * Referenced by: '<S75>/Multiply3'
                                        */
  real_T InterpolationUsingPrelookup_T_j[4];/* Expression: squeeze(table)
                                             * Referenced by: '<S80>/Interpolation Using Prelookup'
                                             */
  real_T Gain6_Gain;                   /* Computed Parameter: Gain6_Gain
                                        * Referenced by: '<S75>/Gain6'
                                        */
  real_T Multiply4_Gain;               /* Expression: config.gains.tecs.G
                                        * Referenced by: '<S75>/Multiply4'
                                        */
  real_T Gain5_Gain;                   /* Computed Parameter: Gain5_Gain
                                        * Referenced by: '<S75>/Gain5'
                                        */
  real_T Gain3_Gain_p;                 /* Expression: -1
                                        * Referenced by: '<S37>/Gain3'
                                        */
  real_T Gain4_Gain;                   /* Computed Parameter: Gain4_Gain
                                        * Referenced by: '<S75>/Gain4'
                                        */
  real_T Multiply2_Gain;               /* Expression: config.gains.tecs.G
                                        * Referenced by: '<S75>/Multiply2'
                                        */
  real_T Constant_Value_b;             /* Expression: 0
                                        * Referenced by: '<S37>/Constant'
                                        */
  real_T Gain7_Gain;                   /* Computed Parameter: Gain7_Gain
                                        * Referenced by: '<S75>/Gain7'
                                        */
  real_T InterpolationUsingPrelookup_T_b[4];/* Expression: squeeze(table)
                                             * Referenced by: '<S78>/Interpolation Using Prelookup'
                                             */
  real_T Gain1_Gain_hh;                /* Expression: config.gains.tecs.G
                                        * Referenced by: '<S76>/Gain1'
                                        */
  real_T Saturation_UpperSat_l;        /* Expression: inf
                                        * Referenced by: '<S76>/Saturation'
                                        */
  real_T Saturation_LowerSat_b;        /* Expression: 0.1
                                        * Referenced by: '<S76>/Saturation'
                                        */
  real_T Constant2_Value;              /* Expression: 0
                                        * Referenced by: '<S76>/Constant2'
                                        */
  real_T Switch7_Threshold;            /* Expression: 0
                                        * Referenced by: '<S27>/Switch7'
                                        */
  real_T Gain1_Gain_l;                 /* Expression: config.gains.rate.max(1)
                                        * Referenced by: '<S70>/Gain1'
                                        */
  real_T Gain2_Gain_m;                 /* Expression: -config.gains.rate.max(2)
                                        * Referenced by: '<S70>/Gain2'
                                        */
  real_T Gain3_Gain_k;                 /* Expression: config.gains.rate.max(3)
                                        * Referenced by: '<S70>/Gain3'
                                        */
  real_T Integrator_IC;                /* Expression: 0
                                        * Referenced by: '<S3>/Integrator'
                                        */
  real_T Prelookup1_BreakpointsData[2];
                        /* Expression: 1/2 * 1.224 * config.trim.points.tas .^ 2
                         * Referenced by: '<S8>/Prelookup1'
                         */
  real_T Prelookup2_BreakpointsData[3];/* Expression: config.trim.points.phi
                                        * Referenced by: '<S8>/Prelookup2'
                                        */
  real_T Prelookup3_BreakpointsData[2];/* Expression: config.trim.points.h
                                        * Referenced by: '<S8>/Prelookup3'
                                        */
  real_T Constant_Value_h[3];          /* Expression: 0:size(table, 4) - 1
                                        * Referenced by: '<S31>/Constant'
                                        */
  real_T _Gain_n;                      /* Computed Parameter: _Gain_n
                                        * Referenced by: '<S3>/   '
                                        */
  real_T Switch2_Threshold;            /* Expression: 0
                                        * Referenced by: '<S27>/Switch2'
                                        */
  real_T Saturation1_UpperSat_e[2];/* Expression: config.gains.attitude.max(1:2)
                                    * Referenced by: '<S3>/Saturation1'
                                    */
  real_T Saturation1_LowerSat_k[2];
                                  /* Expression: -config.gains.attitude.max(1:2)
                                   * Referenced by: '<S3>/Saturation1'
                                   */
  real_T Gain1_Gain_g;                 /* Computed Parameter: Gain1_Gain_g
                                        * Referenced by: '<S3>/Gain1'
                                        */
  real_T Gain2_Gain_c2;                /* Computed Parameter: Gain2_Gain_c2
                                        * Referenced by: '<S3>/Gain2'
                                        */
  real_T Constant_Value_bd[2];         /* Expression: 0:size(table, 4) - 1
                                        * Referenced by: '<S35>/Constant'
                                        */
  real_T InterpolationUsingPrelookup_T_i[8];/* Expression: squeeze(table)
                                             * Referenced by: '<S35>/Interpolation Using Prelookup'
                                             */
  real_T Gain4_Gain_g;                 /* Computed Parameter: Gain4_Gain_g
                                        * Referenced by: '<S3>/Gain4'
                                        */
  real_T Gain5_Gain_i;                 /* Expression: -1
                                        * Referenced by: '<S3>/Gain5'
                                        */
  real_T Constant_Value_k[2];          /* Expression: 0:size(table, 4) - 1
                                        * Referenced by: '<S33>/Constant'
                                        */
  real_T InterpolationUsingPrelookup_T_l[8];/* Expression: squeeze(table)
                                             * Referenced by: '<S33>/Interpolation Using Prelookup'
                                             */
  real_T Gain1_Gain_e;                 /* Expression: config.gains.tecs.G
                                        * Referenced by: '<S32>/Gain1'
                                        */
  real_T Saturation_UpperSat_ag;       /* Expression: inf
                                        * Referenced by: '<S32>/Saturation'
                                        */
  real_T Saturation_LowerSat_kd;       /* Expression: 0.1
                                        * Referenced by: '<S32>/Saturation'
                                        */
  real_T Gain3_Gain_g;                 /* Computed Parameter: Gain3_Gain_g
                                        * Referenced by: '<S3>/Gain3'
                                        */
  real_T Gain6_Gain_m;                 /* Expression: 0
                                        * Referenced by: '<S3>/Gain6'
                                        */
  real_T Switch7_Threshold_k;          /* Expression: 0
                                        * Referenced by: '<S70>/Switch7'
                                        */
  real_T Constant_Value_i[9];          /* Expression: config.aircraft.I
                                        * Referenced by: '<S7>/Constant'
                                        */
  real_T Constant_Value_a[3];          /* Expression: 0:size(table, 4) - 1
                                        * Referenced by: '<S71>/Constant'
                                        */
  real_T Switch2_Threshold_l;          /* Expression: 0
                                        * Referenced by: '<S70>/Switch2'
                                        */
  real_T Constant_Value_f[3];          /* Expression: 0:size(table, 4) - 1
                                        * Referenced by: '<S72>/Constant'
                                        */
  real_T Switch7_Threshold_e;          /* Expression: 0
                                        * Referenced by: '<S12>/Switch7'
                                        */
  real_T Lowpass_A;                    /* Computed Parameter: Lowpass_A
                                        * Referenced by: '<S77>/Low pass'
                                        */
  real_T Lowpass_C;                    /* Computed Parameter: Lowpass_C
                                        * Referenced by: '<S77>/Low pass'
                                        */
  real_T Integrator_IC_a;              /* Expression: 0
                                        * Referenced by: '<S77>/Integrator'
                                        */
  real_T Switch1_Threshold_l;          /* Expression: 0
                                        * Referenced by: '<S12>/Switch1'
                                        */
  real_T Constant_Value_hr[6];         /* Expression: zeros(6, 1)
                                        * Referenced by: '<S2>/Constant'
                                        */
  real_T Constant_Value_n[84];         /* Expression: 0:size(table, 4) - 1
                                        * Referenced by: '<S18>/Constant'
                                        */
  real_T InterpolationUsingPrelookup_T_p[336];/* Expression: squeeze(table)
                                               * Referenced by: '<S18>/Interpolation Using Prelookup'
                                               */
  real_T Constant_Value_bp[84];        /* Expression: zeros(14, 6)
                                        * Referenced by: '<S14>/Constant'
                                        */
  real_T Constant_Value_np[4];         /* Expression: 0:size(table, 4) - 1
                                        * Referenced by: '<S25>/Constant'
                                        */
  real_T InterpolationUsingPrelookup_T_m[16];/* Expression: squeeze(table)
                                              * Referenced by: '<S25>/Interpolation Using Prelookup'
                                              */
  real_T Constant1_Value_e[2];         /* Expression: [0 0]
                                        * Referenced by: '<S13>/Constant1'
                                        */
  real_T Constant_Value_e[4];          /* Expression: 0:size(table, 4) - 1
                                        * Referenced by: '<S19>/Constant'
                                        */
  real_T Constant_Value_m[6];          /* Expression: 0:size(table, 4) - 1
                                        * Referenced by: '<S20>/Constant'
                                        */
  real_T Constant_Value_e4[2];         /* Expression: 0:size(table, 4) - 1
                                        * Referenced by: '<S21>/Constant'
                                        */
  real_T Constant_Value_ns[2];         /* Expression: 0:size(table, 4) - 1
                                        * Referenced by: '<S22>/Constant'
                                        */
  real_T Saturation_UpperSat_pe[14];   /* Expression: config.act.limits(2, :)
                                        * Referenced by: '<S15>/Saturation'
                                        */
  real_T Saturation_LowerSat_bc[14];   /* Expression: config.act.limits(1, :)
                                        * Referenced by: '<S15>/Saturation'
                                        */
  real_T Constant_Value_ke[2];         /* Expression: 0:size(table, 4) - 1
                                        * Referenced by: '<S23>/Constant'
                                        */
  real_T InterpolationUsingPrelookup__pv[8];/* Expression: squeeze(table)
                                             * Referenced by: '<S23>/Interpolation Using Prelookup'
                                             */
  real_T Saturation_UpperSat_cb[14];   /* Expression: config.act.limits(2, :)
                                        * Referenced by: '<S16>/Saturation'
                                        */
  real_T Saturation_LowerSat_c1[14];   /* Expression: config.act.limits(1, :)
                                        * Referenced by: '<S16>/Saturation'
                                        */
  real_T Constant_Value_l[16];         /* Expression: 0:size(table, 4) - 1
                                        * Referenced by: '<S24>/Constant'
                                        */
  real_T InterpolationUsingPrelookup__m0[64];/* Expression: squeeze(table)
                                              * Referenced by: '<S24>/Interpolation Using Prelookup'
                                              */
  real_T Saturation_UpperSat_l4[14];   /* Expression: config.act.limits(2, :)
                                        * Referenced by: '<S17>/Saturation'
                                        */
  real_T Saturation_LowerSat_ey[14];   /* Expression: config.act.limits(1, :)
                                        * Referenced by: '<S17>/Saturation'
                                        */
  real_T Memory_1_InitialCondition;    /* Expression: 0
                                        * Referenced by: '<S1>/Memory'
                                        */
  real_T Memory_2_InitialCondition;    /* Expression: 0
                                        * Referenced by: '<S1>/Memory'
                                        */
  real_T Memory_3_InitialCondition;    /* Expression: 0
                                        * Referenced by: '<S1>/Memory'
                                        */
  real_T Memory_4_InitialCondition;    /* Expression: 0
                                        * Referenced by: '<S1>/Memory'
                                        */
  real_T Constant_Value_j[4];          /* Expression: [1 1 1 1]
                                        * Referenced by: '<S11>/Constant'
                                        */
  real_T Constant1_Value_h[6];         /* Expression: [1 1 1 -1 -1 -1]
                                        * Referenced by: '<S11>/Constant1'
                                        */
  real_T Constant2_Value_k[2];         /* Expression: [-1 -1]
                                        * Referenced by: '<S11>/Constant2'
                                        */
  real_T Constant3_Value[2];           /* Expression: [1 1]
                                        * Referenced by: '<S11>/Constant3'
                                        */
  real_T Switch1_Threshold_m;          /* Expression: 0
                                        * Referenced by: '<S11>/Switch1'
                                        */
  real_T Switch2_Threshold_o;          /* Expression: 0
                                        * Referenced by: '<S11>/Switch2'
                                        */
  real_T Switch3_Threshold;            /* Expression: 0
                                        * Referenced by: '<S11>/Switch3'
                                        */
  real_T Switch4_Threshold;            /* Expression: 0
                                        * Referenced by: '<S11>/Switch4'
                                        */
  real_T Constant_Value_c[2];          /* Expression: 0:size(table, 4) - 1
                                        * Referenced by: '<S34>/Constant'
                                        */
  real_T InterpolationUsingPrelookup__ho[8];/* Expression: squeeze(table)
                                             * Referenced by: '<S34>/Interpolation Using Prelookup'
                                             */
  real_T InterpolationUsingPrelookup__mw[4];/* Expression: squeeze(table)
                                             * Referenced by: '<S79>/Interpolation Using Prelookup'
                                             */
  uint32_T InterpolationUsingPrelookup_dim[3];
                          /* Computed Parameter: InterpolationUsingPrelookup_dim
                           * Referenced by: '<S31>/Interpolation Using Prelookup'
                           */
  uint32_T InterpolationUsingPrelookup_d_f[3];
                          /* Computed Parameter: InterpolationUsingPrelookup_d_f
                           * Referenced by: '<S71>/Interpolation Using Prelookup'
                           */
  uint32_T InterpolationUsingPrelookup_d_d[3];
                          /* Computed Parameter: InterpolationUsingPrelookup_d_d
                           * Referenced by: '<S19>/Interpolation Using Prelookup'
                           */
  uint32_T InterpolationUsingPrelookup_d_o[3];
                          /* Computed Parameter: InterpolationUsingPrelookup_d_o
                           * Referenced by: '<S20>/Interpolation Using Prelookup'
                           */
  uint32_T InterpolationUsingPrelookup_d_m[3];
                          /* Computed Parameter: InterpolationUsingPrelookup_d_m
                           * Referenced by: '<S21>/Interpolation Using Prelookup'
                           */
  uint32_T InterpolationUsingPrelookup__mi[3];
                          /* Computed Parameter: InterpolationUsingPrelookup__mi
                           * Referenced by: '<S22>/Interpolation Using Prelookup'
                           */
};

/* Real-time Model Data Structure */
struct tag_RTM_PX4Controller_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X_PX4Controller_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[5];
  real_T odeF[4][5];
  ODE4_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Block parameters (default storage) */
#ifdef __cplusplus

extern "C" {

#endif

  extern P_PX4Controller_T PX4Controller_P;

#ifdef __cplusplus

}
#endif

/* Block signals (default storage) */
#ifdef __cplusplus

extern "C" {

#endif

  extern struct B_PX4Controller_T PX4Controller_B;

#ifdef __cplusplus

}
#endif

/* Continuous states (default storage) */
extern X_PX4Controller_T PX4Controller_X;

/* Block states (default storage) */
extern struct DW_PX4Controller_T PX4Controller_DW;

/* Zero-crossing (trigger) state */
extern PrevZCX_PX4Controller_T PX4Controller_PrevZCX;

#ifdef __cplusplus

extern "C" {

#endif

  /* External inputs (root inport signals with default storage) */
  extern struct ExtU_PX4Controller_T PX4Controller_U;

  /* External outputs (root outports fed by signals with default storage) */
  extern struct ExtY_PX4Controller_T PX4Controller_Y;

#ifdef __cplusplus

}
#endif

#ifdef __cplusplus

extern "C" {

#endif

  /* Model entry point functions */
  extern void PX4Controller_initialize(void);
  extern void PX4Controller_step(void);
  extern void PX4Controller_terminate(void);

#ifdef __cplusplus

}
#endif

/* Real-time Model object */
#ifdef __cplusplus

extern "C" {

#endif

  extern RT_MODEL_PX4Controller_T *const PX4Controller_M;

#ifdef __cplusplus

}
#endif

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'PX4Controller'
 * '<S1>'   : 'PX4Controller/Controller'
 * '<S2>'   : 'PX4Controller/Controller/Allocation'
 * '<S3>'   : 'PX4Controller/Controller/Attitude'
 * '<S4>'   : 'PX4Controller/Controller/Baselines'
 * '<S5>'   : 'PX4Controller/Controller/Estimator'
 * '<S6>'   : 'PX4Controller/Controller/L1'
 * '<S7>'   : 'PX4Controller/Controller/Rate'
 * '<S8>'   : 'PX4Controller/Controller/SchedBus'
 * '<S9>'   : 'PX4Controller/Controller/TECS'
 * '<S10>'  : 'PX4Controller/Controller/Allocation/ Trim'
 * '<S11>'  : 'PX4Controller/Controller/Allocation/Actuation override'
 * '<S12>'  : 'PX4Controller/Controller/Allocation/FM Override'
 * '<S13>'  : 'PX4Controller/Controller/Allocation/FMCmd'
 * '<S14>'  : 'PX4Controller/Controller/Allocation/Subsystem'
 * '<S15>'  : 'PX4Controller/Controller/Allocation/Subsystem Reference'
 * '<S16>'  : 'PX4Controller/Controller/Allocation/Subsystem Reference1'
 * '<S17>'  : 'PX4Controller/Controller/Allocation/Subsystem Reference2'
 * '<S18>'  : 'PX4Controller/Controller/Allocation/Subsystem Reference3'
 * '<S19>'  : 'PX4Controller/Controller/Allocation/ Trim/Subsystem Reference3'
 * '<S20>'  : 'PX4Controller/Controller/Allocation/ Trim/Subsystem Reference4'
 * '<S21>'  : 'PX4Controller/Controller/Allocation/ Trim/Subsystem Reference5'
 * '<S22>'  : 'PX4Controller/Controller/Allocation/ Trim/Subsystem Reference6'
 * '<S23>'  : 'PX4Controller/Controller/Allocation/Subsystem/Subsystem Reference1'
 * '<S24>'  : 'PX4Controller/Controller/Allocation/Subsystem/Subsystem Reference2'
 * '<S25>'  : 'PX4Controller/Controller/Allocation/Subsystem/Subsystem Reference3'
 * '<S26>'  : 'PX4Controller/Controller/Attitude/Disable'
 * '<S27>'  : 'PX4Controller/Controller/Attitude/Manual override'
 * '<S28>'  : 'PX4Controller/Controller/Attitude/Saturation'
 * '<S29>'  : 'PX4Controller/Controller/Attitude/Signals'
 * '<S30>'  : 'PX4Controller/Controller/Attitude/Subsystem Reference'
 * '<S31>'  : 'PX4Controller/Controller/Attitude/Subsystem Reference3'
 * '<S32>'  : 'PX4Controller/Controller/Attitude/Turn coordinator'
 * '<S33>'  : 'PX4Controller/Controller/Attitude/attitude.Kd'
 * '<S34>'  : 'PX4Controller/Controller/Attitude/attitude.Ki'
 * '<S35>'  : 'PX4Controller/Controller/Attitude/attitude.Kp'
 * '<S36>'  : 'PX4Controller/Controller/Attitude/Signals/Subsystem Reference1'
 * '<S37>'  : 'PX4Controller/Controller/Estimator/Computed'
 * '<S38>'  : 'PX4Controller/Controller/Estimator/Subsystem Reference'
 * '<S39>'  : 'PX4Controller/Controller/Estimator/Subsystem Reference1'
 * '<S40>'  : 'PX4Controller/Controller/Estimator/Subsystem Reference2'
 * '<S41>'  : 'PX4Controller/Controller/Estimator/Subsystem Reference3'
 * '<S42>'  : 'PX4Controller/Controller/Estimator/Subsystem Reference4'
 * '<S43>'  : 'PX4Controller/Controller/Estimator/Subsystem1'
 * '<S44>'  : 'PX4Controller/Controller/Estimator/Subsystem10'
 * '<S45>'  : 'PX4Controller/Controller/Estimator/Subsystem5'
 * '<S46>'  : 'PX4Controller/Controller/Estimator/Subsystem6'
 * '<S47>'  : 'PX4Controller/Controller/Estimator/Subsystem7'
 * '<S48>'  : 'PX4Controller/Controller/Estimator/Subsystem8'
 * '<S49>'  : 'PX4Controller/Controller/Estimator/Subsystem9'
 * '<S50>'  : 'PX4Controller/Controller/Estimator/Subsystem Reference/Subsystem'
 * '<S51>'  : 'PX4Controller/Controller/Estimator/Subsystem Reference/Subsystem1'
 * '<S52>'  : 'PX4Controller/Controller/Estimator/Subsystem Reference/Subsystem2'
 * '<S53>'  : 'PX4Controller/Controller/Estimator/Subsystem Reference1/Subsystem'
 * '<S54>'  : 'PX4Controller/Controller/Estimator/Subsystem Reference1/Subsystem1'
 * '<S55>'  : 'PX4Controller/Controller/Estimator/Subsystem Reference1/Subsystem2'
 * '<S56>'  : 'PX4Controller/Controller/Estimator/Subsystem Reference2/Subsystem'
 * '<S57>'  : 'PX4Controller/Controller/Estimator/Subsystem Reference2/Subsystem1'
 * '<S58>'  : 'PX4Controller/Controller/Estimator/Subsystem Reference2/Subsystem2'
 * '<S59>'  : 'PX4Controller/Controller/Estimator/Subsystem Reference3/Subsystem'
 * '<S60>'  : 'PX4Controller/Controller/Estimator/Subsystem Reference3/Subsystem1'
 * '<S61>'  : 'PX4Controller/Controller/Estimator/Subsystem Reference3/Subsystem2'
 * '<S62>'  : 'PX4Controller/Controller/Estimator/Subsystem Reference4/Subsystem'
 * '<S63>'  : 'PX4Controller/Controller/Estimator/Subsystem Reference4/Subsystem1'
 * '<S64>'  : 'PX4Controller/Controller/Estimator/Subsystem Reference4/Subsystem2'
 * '<S65>'  : 'PX4Controller/Controller/L1/Controller'
 * '<S66>'  : 'PX4Controller/Controller/L1/L1'
 * '<S67>'  : 'PX4Controller/Controller/L1/Line'
 * '<S68>'  : 'PX4Controller/Controller/L1/Saturation'
 * '<S69>'  : 'PX4Controller/Controller/L1/L1/gains.L1.wnInv'
 * '<S70>'  : 'PX4Controller/Controller/Rate/Manual override'
 * '<S71>'  : 'PX4Controller/Controller/Rate/Subsystem Reference3'
 * '<S72>'  : 'PX4Controller/Controller/Rate/rate.Kp'
 * '<S73>'  : 'PX4Controller/Controller/TECS/Disable'
 * '<S74>'  : 'PX4Controller/Controller/TECS/Saturation'
 * '<S75>'  : 'PX4Controller/Controller/TECS/Specific energy '
 * '<S76>'  : 'PX4Controller/Controller/TECS/Specific energy balance'
 * '<S77>'  : 'PX4Controller/Controller/TECS/Total specific energy'
 * '<S78>'  : 'PX4Controller/Controller/TECS/Specific energy balance/tecs.balKd'
 * '<S79>'  : 'PX4Controller/Controller/TECS/Specific energy balance/tecs.balKi'
 * '<S80>'  : 'PX4Controller/Controller/TECS/Specific energy balance/tecs.balKp'
 * '<S81>'  : 'PX4Controller/Controller/TECS/Total specific energy/gains.tecs.totKi'
 * '<S82>'  : 'PX4Controller/Controller/TECS/Total specific energy/tecs.totKd'
 * '<S83>'  : 'PX4Controller/Controller/TECS/Total specific energy/tecs.totKp'
 */
#endif                                 /* RTW_HEADER_PX4Controller_h_ */
