/*
 * PX4Controller.cpp
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

#include "PX4Controller.h"
#include "rtwtypes.h"
#include <cmath>
#include "PX4Controller_private.h"
#include <cstring>
#include "rt_defines.h"
#include "zero_crossing_types.h"

extern "C" {

#include "rt_nonfinite.h"

}
/* Block signals (default storage) */
  B_PX4Controller_T PX4Controller_B;

/* Continuous states */
X_PX4Controller_T PX4Controller_X;

/* Block states (default storage) */
DW_PX4Controller_T PX4Controller_DW;

/* Previous zero-crossings (trigger) states */
PrevZCX_PX4Controller_T PX4Controller_PrevZCX;

/* External inputs (root inport signals with default storage) */
ExtU_PX4Controller_T PX4Controller_U;

/* External outputs (root outports fed by signals with default storage) */
ExtY_PX4Controller_T PX4Controller_Y;

/* Real-time model */
RT_MODEL_PX4Controller_T PX4Controller_M_{ };

RT_MODEL_PX4Controller_T *const PX4Controller_M{ &PX4Controller_M_ };

uint32_T plook_bincp(real_T u, const real_T bp[], uint32_T maxIndex, real_T
                     *fraction, uint32_T *prevIndex)
{
  uint32_T bpIndex;

  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Clip'
     Use previous index: 'on'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u <= bp[0U]) {
    bpIndex = 0U;
    *fraction = 0.0;
  } else if (u < bp[maxIndex]) {
    bpIndex = binsearch_u32d_prevIdx(u, bp, *prevIndex, maxIndex);
    *fraction = (u - bp[bpIndex]) / (bp[bpIndex + 1U] - bp[bpIndex]);
  } else {
    bpIndex = maxIndex - 1U;
    *fraction = 1.0;
  }

  *prevIndex = bpIndex;
  return bpIndex;
}

real_T intrp2d_l_pw(const uint32_T bpIndex[], const real_T frac[], const real_T
                    table[], const uint32_T stride)
{
  real_T yL_0d0;
  real_T yL_0d1;
  uint32_T offset_1d;

  /* Column-major Interpolation 2-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  offset_1d = bpIndex[1U] * stride + bpIndex[0U];
  yL_0d0 = table[offset_1d];
  yL_0d0 += (table[offset_1d + 1U] - yL_0d0) * frac[0U];
  offset_1d += stride;
  yL_0d1 = table[offset_1d];
  return (((table[offset_1d + 1U] - yL_0d1) * frac[0U] + yL_0d1) - yL_0d0) *
    frac[1U] + yL_0d0;
}

real_T intrp3d_l_pw(const uint32_T bpIndex[], const real_T frac[], const real_T
                    table[], const uint32_T stride[])
{
  real_T yL_0d0;
  real_T yL_1d;
  real_T yL_2d;
  uint32_T offset_0d;
  uint32_T offset_2d;

  /* Column-major Interpolation 3-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  offset_2d = (bpIndex[2U] * stride[2U] + bpIndex[1U] * stride[1U]) + bpIndex[0U];
  yL_0d0 = table[offset_2d];
  yL_1d = (table[offset_2d + 1U] - yL_0d0) * frac[0U] + yL_0d0;
  offset_0d = offset_2d + stride[1U];
  yL_0d0 = table[offset_0d];
  yL_2d = (((table[offset_0d + 1U] - yL_0d0) * frac[0U] + yL_0d0) - yL_1d) *
    frac[1U] + yL_1d;
  offset_2d += stride[2U];
  yL_0d0 = table[offset_2d];
  yL_1d = (table[offset_2d + 1U] - yL_0d0) * frac[0U] + yL_0d0;
  offset_0d = offset_2d + stride[1U];
  yL_0d0 = table[offset_0d];
  return (((((table[offset_0d + 1U] - yL_0d0) * frac[0U] + yL_0d0) - yL_1d) *
           frac[1U] + yL_1d) - yL_2d) * frac[2U] + yL_2d;
}

uint32_T binsearch_u32d_prevIdx(real_T u, const real_T bp[], uint32_T startIndex,
  uint32_T maxIndex)
{
  uint32_T bpIndex;
  uint32_T found;
  uint32_T iLeft;
  uint32_T iRght;

  /* Binary Search using Previous Index */
  bpIndex = startIndex;
  iLeft = 0U;
  iRght = maxIndex;
  found = 0U;
  while (found == 0U) {
    if (u < bp[bpIndex]) {
      iRght = bpIndex - 1U;
      bpIndex = ((bpIndex + iLeft) - 1U) >> 1U;
    } else if (u < bp[bpIndex + 1U]) {
      found = 1U;
    } else {
      iLeft = bpIndex + 1U;
      bpIndex = ((bpIndex + iRght) + 1U) >> 1U;
    }
  }

  return bpIndex;
}

/*
 * This function updates continuous states using the ODE4 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T t { rtsiGetT(si) };

  time_T tnew { rtsiGetSolverStopTime(si) };

  time_T h { rtsiGetStepSize(si) };

  real_T *x { rtsiGetContStates(si) };

  ODE4_IntgData *id { static_cast<ODE4_IntgData *>(rtsiGetSolverData(si)) };

  real_T *y { id->y };

  real_T *f0 { id->f[0] };

  real_T *f1 { id->f[1] };

  real_T *f2 { id->f[2] };

  real_T *f3 { id->f[3] };

  real_T temp;
  int_T i;
  int_T nXc { 5 };

  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) std::memcpy(y, x,
                     static_cast<uint_T>(nXc)*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  PX4Controller_derivatives();

  /* f1 = f(t + (h/2), y + (h/2)*f0) */
  temp = 0.5 * h;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f0[i]);
  }

  rtsiSetT(si, t + temp);
  rtsiSetdX(si, f1);
  PX4Controller_step();
  PX4Controller_derivatives();

  /* f2 = f(t + (h/2), y + (h/2)*f1) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f1[i]);
  }

  rtsiSetdX(si, f2);
  PX4Controller_step();
  PX4Controller_derivatives();

  /* f3 = f(t + h, y + h*f2) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (h*f2[i]);
  }

  rtsiSetT(si, tnew);
  rtsiSetdX(si, f3);
  PX4Controller_step();
  PX4Controller_derivatives();

  /* tnew = t + h
     ynew = y + (h/6)*(f0 + 2*f1 + 2*f2 + 2*f3) */
  temp = h / 6.0;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + temp*(f0[i] + 2.0*f1[i] + 2.0*f2[i] + f3[i]);
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

real_T rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  if (std::isnan(u0) || std::isnan(u1)) {
    y = (rtNaN);
  } else if (std::isinf(u0) && std::isinf(u1)) {
    int32_T u0_0;
    int32_T u1_0;
    if (u0 > 0.0) {
      u0_0 = 1;
    } else {
      u0_0 = -1;
    }

    if (u1 > 0.0) {
      u1_0 = 1;
    } else {
      u1_0 = -1;
    }

    y = std::atan2(static_cast<real_T>(u0_0), static_cast<real_T>(u1_0));
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = std::atan2(u0, u1);
  }

  return y;
}

/* Model step function */
void PX4Controller_step(void)
{
  real_T rtb_Assignment2[84];
  real_T rtb_InterpolationUsingPreloo_dw[84];
  real_T rtb_InterpolationUsingPrelook_c[16];
  real_T rtb_Saturation_az[14];
  real_T rtb_Saturation_j[14];
  real_T rtb_Subtract1[6];
  real_T rtb_Subtract1_tmp[6];
  real_T rtb_Subtract1_tmp_0[6];
  real_T tmp[6];
  real_T frac_2[3];
  real_T frac_5[3];
  real_T frac_c[3];
  real_T frac_d[3];
  real_T frac_e[3];
  real_T frac_f[3];
  real_T rtb_Switch7_e[3];
  real_T frac[2];
  real_T frac_0[2];
  real_T frac_1[2];
  real_T frac_3[2];
  real_T frac_4[2];
  real_T frac_6[2];
  real_T frac_7[2];
  real_T frac_8[2];
  real_T frac_9[2];
  real_T frac_a[2];
  real_T frac_b[2];
  real_T frac_g[2];
  real_T frac_h[2];
  real_T frac_i[2];
  real_T frac_j[2];
  real_T rtb_Gain1_a;
  real_T rtb_Gain1_g;
  real_T rtb_Gain2;
  real_T rtb_Gain2_m;
  real_T rtb_Gain6;
  real_T rtb_InterpolationUsingPrelook_0;
  real_T rtb_InterpolationUsingPrelook_4;
  real_T rtb_InterpolationUsingPrelook_5;
  real_T rtb_InterpolationUsingPrelook_6;
  real_T rtb_InterpolationUsingPrelook_k;
  real_T rtb_Switch1_o;
  real_T rtb_Switch7_idx_0;
  real_T rtb_Switch7_idx_1;
  real_T rtb_Tan1;
  real_T rtb_f1;
  real_T rtb_f1_m;
  real_T rtb_f2;
  real_T rtb_f2_p;
  real_T rtb_sebErr;
  real_T rtb_skeCmd;
  real_T rtb_skeDot;
  real_T rtb_skeDotCmd;
  real_T rtb_spe;
  real_T rtb_speCmd;
  real_T rtb_speDot;
  real_T rtb_speDotCmd;
  real_T rtb_tecsthetaCmd;
  uint32_T bpIndex_2[4];
  uint32_T bpIndex_5[4];
  uint32_T bpIndex_c[4];
  uint32_T bpIndex_d[4];
  uint32_T bpIndex_e[4];
  uint32_T bpIndex_f[4];
  uint32_T bpIndex_3[3];
  uint32_T bpIndex_4[3];
  uint32_T bpIndex_6[3];
  uint32_T bpIndex_8[3];
  uint32_T bpIndex_9[3];
  uint32_T bpIndex_g[3];
  uint32_T bpIndex_h[3];
  uint32_T bpIndex_i[3];
  uint32_T bpIndex[2];
  uint32_T bpIndex_0[2];
  uint32_T bpIndex_1[2];
  uint32_T bpIndex_7[2];
  uint32_T bpIndex_a[2];
  uint32_T bpIndex_b[2];
  uint32_T bpIndex_j[2];
  uint32_T rtb_k1;
  uint32_T rtb_k1_g;
  uint32_T rtb_k2;
  uint32_T rtb_k2_e;
  uint32_T rtb_k3;
  boolean_T didZcEventOccur;
  if (rtmIsMajorTimeStep(PX4Controller_M)) {
    /* set solver stop time */
    if (!(PX4Controller_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&PX4Controller_M->solverInfo,
                            ((PX4Controller_M->Timing.clockTickH0 + 1) *
        PX4Controller_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&PX4Controller_M->solverInfo,
                            ((PX4Controller_M->Timing.clockTick0 + 1) *
        PX4Controller_M->Timing.stepSize0 + PX4Controller_M->Timing.clockTickH0 *
        PX4Controller_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(PX4Controller_M)) {
    PX4Controller_M->Timing.t[0] = rtsiGetT(&PX4Controller_M->solverInfo);
  }

  /* Gain: '<S27>/Gain1' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_Gain1_g = PX4Controller_P.Gain1_Gain_p * PX4Controller_U.CmdBusIn.rc[1];

  /* Gain: '<S27>/Gain2' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_Tan1 = PX4Controller_P.Gain2_Gain_p * PX4Controller_U.CmdBusIn.rc[2];

  /* Gain: '<S27>/Gain3' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_Gain6 = PX4Controller_P.Gain3_Gain_b * PX4Controller_U.CmdBusIn.rc[3];

  /* DotProduct: '<S66>/Dot Product' incorporates:
   *  Inport: '<Root>/StateBus'
   */
  rtb_InterpolationUsingPrelook_0 = PX4Controller_U.StateBus_m.vNed[0] *
    PX4Controller_U.StateBus_m.vNed[0] + PX4Controller_U.StateBus_m.vNed[1] *
    PX4Controller_U.StateBus_m.vNed[1];

  /* Gain: '<S8>/Gain' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   *  Math: '<S8>/Square'
   */
  rtb_f1_m = PX4Controller_U.CmdBusIn.tasCmd * PX4Controller_U.CmdBusIn.tasCmd *
    PX4Controller_P.Gain_Gain_k;

  /* PreLookup: '<S8>/Prelookup4' */
  rtb_k1 = plook_bincp(rtb_f1_m, PX4Controller_P.Prelookup4_BreakpointsData, 1U,
                       &rtb_f1, &PX4Controller_DW.Prelookup4_DWORK1);

  /* PreLookup: '<S8>/Prelookup6' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_k2 = plook_bincp(PX4Controller_U.CmdBusIn.hCmd,
                       PX4Controller_P.Prelookup6_BreakpointsData, 1U, &rtb_f2,
                       &PX4Controller_DW.Prelookup6_DWORK1);

  /* Interpolation_n-D: '<S69>/Interpolation Using Prelookup' */
  frac[0] = rtb_f1;
  frac[1] = rtb_f2;
  bpIndex[0] = rtb_k1;
  bpIndex[1] = rtb_k2;

  /* Product: '<S69>/Product' incorporates:
   *  DotProduct: '<S66>/Dot Product'
   *  Gain: '<S66>/Gain'
   *  Interpolation_n-D: '<S69>/Interpolation Using Prelookup'
   *  Sqrt: '<S66>/Square Root'
   */
  rtb_skeDot = PX4Controller_P.Gain_Gain_a0 * std::sqrt
    (rtb_InterpolationUsingPrelook_0) * intrp2d_l_pw(bpIndex, frac,
    PX4Controller_P.InterpolationUsingPrelookup_T_h, 2U);

  /* Saturate: '<S66>/Saturation' */
  if (rtb_skeDot > PX4Controller_P.Saturation_UpperSat_m) {
    rtb_skeDot = PX4Controller_P.Saturation_UpperSat_m;
  } else if (rtb_skeDot < PX4Controller_P.Saturation_LowerSat_k) {
    rtb_skeDot = PX4Controller_P.Saturation_LowerSat_k;
  }

  /* End of Saturate: '<S66>/Saturation' */

  /* Trigonometry: '<S67>/Cos1' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_spe = std::sin(PX4Controller_U.CmdBusIn.trackLine[0]);

  /* Trigonometry: '<S67>/Cos' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_skeCmd = std::cos(PX4Controller_U.CmdBusIn.trackLine[0]);

  /* SignalConversion generated from: '<S67>/Dot Product1' incorporates:
   *  Gain: '<S67>/Gain'
   */
  frac[1] = PX4Controller_P.Gain_Gain_do * rtb_skeCmd;

  /* DotProduct: '<S67>/Dot Product1' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   *  Inport: '<Root>/StateBus'
   *  SignalConversion generated from: '<S67>/Dot Product1'
   *  Sum: '<S67>/Minus'
   */
  rtb_Gain2 = (PX4Controller_U.StateBus_m.xNed[0] -
               PX4Controller_U.CmdBusIn.trackLine[1]) * rtb_spe +
    (PX4Controller_U.StateBus_m.xNed[1] - PX4Controller_U.CmdBusIn.trackLine[2])
    * frac[1];

  /* Sum: '<S67>/Minus1' incorporates:
   *  DotProduct: '<S67>/Dot Product1'
   *  Math: '<S67>/Square'
   *  Math: '<S67>/Square2'
   */
  rtb_f2_p = rtb_skeDot * rtb_skeDot - rtb_Gain2 * rtb_Gain2;

  /* Saturate: '<S67>/Saturation' */
  if (rtb_f2_p > PX4Controller_P.Saturation_UpperSat_a) {
    rtb_f2_p = PX4Controller_P.Saturation_UpperSat_a;
  } else if (rtb_f2_p < PX4Controller_P.Saturation_LowerSat_e) {
    rtb_f2_p = PX4Controller_P.Saturation_LowerSat_e;
  }

  /* End of Saturate: '<S67>/Saturation' */

  /* Trigonometry: '<S65>/Cos' incorporates:
   *  Inport: '<Root>/StateBus'
   *  Trigonometry: '<S30>/sincos'
   *  Trigonometry: '<S36>/sincos'
   */
  rtb_InterpolationUsingPrelook_6 = std::cos
    (PX4Controller_U.StateBus_m.eulerBody[1]);

  /* Gain: '<S6>/Gain2' incorporates:
   *  DotProduct: '<S66>/Dot Product'
   *  DotProduct: '<S67>/Dot Product1'
   *  DotProduct: '<S67>/Dot Product2'
   *  DotProduct: '<S67>/Dot Product3'
   *  Gain: '<S65>/1//G'
   *  Gain: '<S65>/Gain'
   *  Inport: '<Root>/StateBus'
   *  Product: '<S65>/Divide'
   *  Product: '<S65>/Product'
   *  SignalConversion generated from: '<S67>/Dot Product1'
   *  Sqrt: '<S67>/Square Root'
   *  Sum: '<S67>/Minus2'
   *  Trigonometry: '<S65>/Atan'
   *  Trigonometry: '<S65>/Cos'
   *  Trigonometry: '<S65>/Sin'
   *  Trigonometry: '<S67>/Atan1'
   *  Trigonometry: '<S67>/Atan2'
   */
  rtb_Gain2 = std::atan(std::sin(rt_atan2d_snf(rtb_spe *
    PX4Controller_U.StateBus_m.vNed[0] + frac[1] *
    PX4Controller_U.StateBus_m.vNed[1], rtb_skeCmd *
    PX4Controller_U.StateBus_m.vNed[0] + rtb_spe *
    PX4Controller_U.StateBus_m.vNed[1]) + rt_atan2d_snf(rtb_Gain2, std::sqrt
    (rtb_f2_p))) * (1.0 / rtb_skeDot) * rtb_InterpolationUsingPrelook_0 *
                        PX4Controller_P.Gain_Gain_dn * PX4Controller_P.uG_Gain *
                        rtb_InterpolationUsingPrelook_6) *
    PX4Controller_P.Gain2_Gain_h;

  /* Logic: '<S73>/OR' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  PX4Controller_B.OR = ((PX4Controller_U.CmdBusIn.manualAttitude != 0.0) ||
                        (PX4Controller_U.CmdBusIn.manualRate != 0.0) ||
                        (PX4Controller_U.CmdBusIn.manualFM != 0.0) ||
                        (PX4Controller_U.CmdBusIn.manualActuation != 0.0));

  /* Integrator: '<S76>/Integrator1' */
  if (rtsiIsModeUpdateTimeStep(&PX4Controller_M->solverInfo)) {
    didZcEventOccur = (((PX4Controller_PrevZCX.Integrator1_Reset_ZCE ==
                         POS_ZCSIG) != PX4Controller_B.OR) &&
                       (PX4Controller_PrevZCX.Integrator1_Reset_ZCE !=
                        UNINITIALIZED_ZCSIG));
    PX4Controller_PrevZCX.Integrator1_Reset_ZCE = PX4Controller_B.OR;

    /* evaluate zero-crossings and the level of the reset signal */
    if (didZcEventOccur || PX4Controller_B.OR) {
      PX4Controller_X.Integrator1_CSTATE = PX4Controller_P.Integrator1_IC;
    }
  }

  /* Gain: '<S75>/Multiply' incorporates:
   *  Gain: '<S75>/Gain'
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_speCmd = PX4Controller_P.Gain_Gain_f * PX4Controller_U.CmdBusIn.hCmd *
    PX4Controller_P.Multiply_Gain;

  /* Gain: '<S75>/Gain2' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_skeCmd = PX4Controller_P.Gain2_Gain_o * PX4Controller_U.CmdBusIn.tasCmd;

  /* Gain: '<S75>/Multiply5' incorporates:
   *  Math: '<S75>/Square1'
   */
  rtb_skeCmd = rtb_skeCmd * rtb_skeCmd * PX4Controller_P.Multiply5_Gain;

  /* Gain: '<S75>/Multiply1' incorporates:
   *  Gain: '<S37>/Gain2'
   *  Gain: '<S75>/Gain1'
   *  Inport: '<Root>/StateBus'
   */
  rtb_spe = PX4Controller_P.Gain2_Gain_c * PX4Controller_U.StateBus_m.xNed[2] *
    PX4Controller_P.Gain1_Gain_h * PX4Controller_P.Multiply1_Gain;

  /* Gain: '<S75>/Gain3' incorporates:
   *  Inport: '<Root>/StateBus'
   */
  rtb_skeDot = PX4Controller_P.Gain3_Gain_j * PX4Controller_U.StateBus_m.tas;

  /* Gain: '<S75>/Multiply3' incorporates:
   *  Math: '<S75>/Square'
   */
  rtb_InterpolationUsingPrelook_0 = rtb_skeDot * rtb_skeDot *
    PX4Controller_P.Multiply3_Gain;

  /* Sum: '<S76>/Sum' incorporates:
   *  Constant: '<S76>/Constant'
   *  Constant: '<S76>/Constant1'
   *  Product: '<S76>/Product'
   *  Product: '<S76>/Product1'
   *  Product: '<S76>/Product2'
   *  Product: '<S76>/Product3'
   *  Sum: '<S76>/seb'
   *  Sum: '<S76>/sebCmd'
   */
  rtb_sebErr = (PX4Controller_P.Constant_Value * rtb_speCmd -
                PX4Controller_P.Constant1_Value * rtb_skeCmd) -
    (PX4Controller_P.Constant_Value * rtb_spe - PX4Controller_P.Constant1_Value *
     rtb_InterpolationUsingPrelook_0);

  /* Interpolation_n-D: '<S80>/Interpolation Using Prelookup' */
  frac_0[0] = rtb_f1;
  frac_0[1] = rtb_f2;
  bpIndex_0[0] = rtb_k1;
  bpIndex_0[1] = rtb_k2;

  /* Gain: '<S75>/Multiply4' incorporates:
   *  Gain: '<S75>/Gain6'
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_speDotCmd = PX4Controller_P.Gain6_Gain * PX4Controller_U.CmdBusIn.hDotCmd *
    PX4Controller_P.Multiply4_Gain;

  /* Product: '<S75>/Product1' incorporates:
   *  Gain: '<S75>/Gain5'
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_skeDotCmd = PX4Controller_P.Gain5_Gain *
    PX4Controller_U.CmdBusIn.tasDotCmd * rtb_skeDot;

  /* Gain: '<S75>/Multiply2' incorporates:
   *  Gain: '<S37>/Gain3'
   *  Gain: '<S75>/Gain4'
   *  Inport: '<Root>/StateBus'
   */
  rtb_speDot = PX4Controller_P.Gain3_Gain_p * PX4Controller_U.StateBus_m.vNed[2]
    * PX4Controller_P.Gain4_Gain * PX4Controller_P.Multiply2_Gain;
  if (rtmIsMajorTimeStep(PX4Controller_M)) {
    /* Gain: '<S75>/Gain7' incorporates:
     *  Constant: '<S37>/Constant'
     */
    PX4Controller_B.Gain7 = PX4Controller_P.Gain7_Gain *
      PX4Controller_P.Constant_Value_b;
  }

  /* Product: '<S75>/Product' */
  rtb_skeDot *= PX4Controller_B.Gain7;

  /* Interpolation_n-D: '<S78>/Interpolation Using Prelookup' */
  frac_1[0] = rtb_f1;
  frac_1[1] = rtb_f2;
  bpIndex_1[0] = rtb_k1;
  bpIndex_1[1] = rtb_k2;

  /* Gain: '<S76>/Gain1' incorporates:
   *  Inport: '<Root>/StateBus'
   */
  rtb_f2_p = PX4Controller_P.Gain1_Gain_hh * PX4Controller_U.StateBus_m.tas;

  /* Saturate: '<S76>/Saturation' */
  if (rtb_f2_p > PX4Controller_P.Saturation_UpperSat_l) {
    rtb_f2_p = PX4Controller_P.Saturation_UpperSat_l;
  } else if (rtb_f2_p < PX4Controller_P.Saturation_LowerSat_b) {
    rtb_f2_p = PX4Controller_P.Saturation_LowerSat_b;
  }

  /* End of Saturate: '<S76>/Saturation' */

  /* Sum: '<S76>/Add1' incorporates:
   *  Constant: '<S76>/Constant'
   *  Constant: '<S76>/Constant1'
   *  Constant: '<S76>/Constant2'
   *  Integrator: '<S76>/Integrator1'
   *  Interpolation_n-D: '<S78>/Interpolation Using Prelookup'
   *  Interpolation_n-D: '<S80>/Interpolation Using Prelookup'
   *  Product: '<S76>/Divide'
   *  Product: '<S76>/Product4'
   *  Product: '<S76>/Product5'
   *  Product: '<S76>/Product6'
   *  Product: '<S76>/Product7'
   *  Product: '<S78>/Product'
   *  Product: '<S80>/Product'
   *  Sum: '<S76>/Add'
   *  Sum: '<S76>/Sum1'
   *  Sum: '<S76>/sebDot'
   *  Sum: '<S76>/sebDotCmd '
   */
  rtb_tecsthetaCmd = (((PX4Controller_P.Constant_Value * rtb_speDotCmd -
                        PX4Controller_P.Constant1_Value * rtb_skeDotCmd) -
                       (PX4Controller_P.Constant_Value * rtb_speDot -
                        PX4Controller_P.Constant1_Value * rtb_skeDot)) *
                      intrp2d_l_pw(bpIndex_1, frac_1,
    PX4Controller_P.InterpolationUsingPrelookup_T_b, 2U) + (rtb_sebErr *
    intrp2d_l_pw(bpIndex_0, frac_0,
                 PX4Controller_P.InterpolationUsingPrelookup_T_j, 2U) +
    PX4Controller_X.Integrator1_CSTATE)) / rtb_f2_p +
    PX4Controller_P.Constant2_Value;

  /* Switch: '<S27>/Switch7' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  if (PX4Controller_U.CmdBusIn.manualAttitude >
      PX4Controller_P.Switch7_Threshold) {
    rtb_Switch7_idx_0 = rtb_Gain1_g;
    rtb_Switch7_idx_1 = rtb_Tan1;

    /* Outport: '<Root>/CmdBusOut' */
    PX4Controller_Y.CmdBusOut.eulerCmd[2] = rtb_Gain6;
  } else {
    /* Outport: '<Root>/CmdBusOut' incorporates:
     *  Assignment: '<S74>/Assignment'
     *  SignalConversion generated from: '<S68>/Assignment'
     */
    PX4Controller_Y.CmdBusOut.eulerCmd[2] = PX4Controller_U.CmdBusIn.eulerCmd[2];

    /* Assignment: '<S68>/Assignment' */
    rtb_Switch7_idx_0 = rtb_Gain2;

    /* Assignment: '<S74>/Assignment' */
    rtb_Switch7_idx_1 = rtb_tecsthetaCmd;
  }

  /* End of Switch: '<S27>/Switch7' */

  /* Gain: '<S70>/Gain1' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_Gain1_a = PX4Controller_P.Gain1_Gain_l * PX4Controller_U.CmdBusIn.rc[1];

  /* Gain: '<S70>/Gain2' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_Gain2_m = PX4Controller_P.Gain2_Gain_m * PX4Controller_U.CmdBusIn.rc[2];

  /* Gain: '<S70>/Gain3' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_Switch1_o = PX4Controller_P.Gain3_Gain_k * PX4Controller_U.CmdBusIn.rc[3];

  /* Logic: '<S26>/OR' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  PX4Controller_B.OR_k = ((PX4Controller_U.CmdBusIn.manualRate != 0.0) ||
    (PX4Controller_U.CmdBusIn.manualFM != 0.0) ||
    (PX4Controller_U.CmdBusIn.manualActuation != 0.0));

  /* Integrator: '<S3>/Integrator' */
  if (rtsiIsModeUpdateTimeStep(&PX4Controller_M->solverInfo)) {
    didZcEventOccur = (((PX4Controller_PrevZCX.Integrator_Reset_ZCE == POS_ZCSIG)
                        != PX4Controller_B.OR_k) &&
                       (PX4Controller_PrevZCX.Integrator_Reset_ZCE !=
                        UNINITIALIZED_ZCSIG));
    PX4Controller_PrevZCX.Integrator_Reset_ZCE = PX4Controller_B.OR_k;

    /* evaluate zero-crossings and the level of the reset signal */
    if (didZcEventOccur || PX4Controller_B.OR_k) {
      PX4Controller_X.Integrator_CSTATE[0] = PX4Controller_P.Integrator_IC;
      PX4Controller_X.Integrator_CSTATE[1] = PX4Controller_P.Integrator_IC;
    }
  }

  /* PreLookup: '<S8>/Prelookup1' */
  rtb_k1_g = plook_bincp(rtb_f1_m, PX4Controller_P.Prelookup1_BreakpointsData,
    1U, &rtb_f1_m, &PX4Controller_DW.Prelookup1_DWORK1);

  /* PreLookup: '<S8>/Prelookup2' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_k2_e = plook_bincp(PX4Controller_U.CmdBusIn.eulerCmd[0],
    PX4Controller_P.Prelookup2_BreakpointsData, 2U, &rtb_f2_p,
    &PX4Controller_DW.Prelookup2_DWORK1);

  /* PreLookup: '<S8>/Prelookup3' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_k3 = plook_bincp(PX4Controller_U.CmdBusIn.hCmd,
                       PX4Controller_P.Prelookup3_BreakpointsData, 1U,
                       &rtb_InterpolationUsingPrelook_k,
                       &PX4Controller_DW.Prelookup3_DWORK1);

  /* Interpolation_n-D: '<S31>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S31>/Constant'
   */
  frac_2[0] = rtb_f1_m;
  frac_2[1] = rtb_f2_p;
  frac_2[2] = rtb_InterpolationUsingPrelook_k;
  bpIndex_2[0] = rtb_k1_g;
  bpIndex_2[1] = rtb_k2_e;
  bpIndex_2[2] = rtb_k3;
  if (PX4Controller_P.Constant_Value_h[0] > 2.0) {
    bpIndex_2[3] = 2U;
  } else if (PX4Controller_P.Constant_Value_h[0] >= 0.0) {
    bpIndex_2[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_h[0]);
  } else {
    bpIndex_2[3] = 0U;
  }

  rtb_InterpolationUsingPrelook_5 = intrp3d_l_pw(bpIndex_2, frac_2,
    &PX4Controller_P.SubsystemReference3_table[12U * bpIndex_2[3]],
    PX4Controller_P.InterpolationUsingPrelookup_dim);

  /* Gain: '<S3>/   ' */
  rtb_InterpolationUsingPrelook_4 = PX4Controller_P._Gain_n *
    rtb_InterpolationUsingPrelook_5;

  /* Interpolation_n-D: '<S31>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S31>/Constant'
   */
  if (PX4Controller_P.Constant_Value_h[1] > 2.0) {
    bpIndex_2[3] = 2U;
  } else if (PX4Controller_P.Constant_Value_h[1] >= 0.0) {
    bpIndex_2[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_h[1]);
  } else {
    bpIndex_2[3] = 0U;
  }

  rtb_InterpolationUsingPrelook_5 = intrp3d_l_pw(bpIndex_2, frac_2,
    &PX4Controller_P.SubsystemReference3_table[12U * bpIndex_2[3]],
    PX4Controller_P.InterpolationUsingPrelookup_dim);

  /* Switch: '<S27>/Switch2' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  if (PX4Controller_U.CmdBusIn.manualAttitude >
      PX4Controller_P.Switch2_Threshold) {
    frac_2[0] = rtb_Gain1_g;
    frac_2[1] = rtb_Tan1;
    frac_2[2] = rtb_Gain6;
  } else {
    /* Assignment: '<S74>/Assignment1' incorporates:
     *  SignalConversion generated from: '<S68>/Assignment1'
     */
    frac_2[2] = PX4Controller_U.CmdBusIn.eulerSat[2];

    /* Saturate: '<S68>/Saturation1' */
    if (rtb_Gain2 > PX4Controller_P.Saturation1_UpperSat_g) {
      /* Assignment: '<S68>/Assignment1' incorporates:
       *  Assignment: '<S74>/Assignment1'
       */
      frac_2[0] = PX4Controller_P.Saturation1_UpperSat_g;
    } else if (rtb_Gain2 < PX4Controller_P.Saturation1_LowerSat_b) {
      /* Assignment: '<S68>/Assignment1' incorporates:
       *  Assignment: '<S74>/Assignment1'
       */
      frac_2[0] = PX4Controller_P.Saturation1_LowerSat_b;
    } else {
      /* Assignment: '<S68>/Assignment1' incorporates:
       *  Assignment: '<S74>/Assignment1'
       */
      frac_2[0] = rtb_Gain2;
    }

    /* End of Saturate: '<S68>/Saturation1' */

    /* Saturate: '<S74>/Saturation1' */
    if (rtb_tecsthetaCmd > PX4Controller_P.Saturation1_UpperSat) {
      /* Assignment: '<S74>/Assignment1' */
      frac_2[1] = PX4Controller_P.Saturation1_UpperSat;
    } else if (rtb_tecsthetaCmd < PX4Controller_P.Saturation1_LowerSat) {
      /* Assignment: '<S74>/Assignment1' */
      frac_2[1] = PX4Controller_P.Saturation1_LowerSat;
    } else {
      /* Assignment: '<S74>/Assignment1' */
      frac_2[1] = rtb_tecsthetaCmd;
    }

    /* End of Saturate: '<S74>/Saturation1' */
  }

  /* End of Switch: '<S27>/Switch2' */

  /* Interpolation_n-D: '<S35>/Interpolation Using Prelookup' */
  frac_3[0] = rtb_f1;
  frac_3[1] = rtb_f2;
  bpIndex_3[0] = rtb_k1;
  bpIndex_3[1] = rtb_k2;

  /* Saturate: '<S3>/Saturation1' */
  if (frac_2[0] > PX4Controller_P.Saturation1_UpperSat_e[0]) {
    rtb_Gain2 = PX4Controller_P.Saturation1_UpperSat_e[0];
  } else if (frac_2[0] < PX4Controller_P.Saturation1_LowerSat_k[0]) {
    rtb_Gain2 = PX4Controller_P.Saturation1_LowerSat_k[0];
  } else {
    rtb_Gain2 = frac_2[0];
  }

  /* Sum: '<S3>/Subtract' incorporates:
   *  Gain: '<S3>/Gain1'
   *  Gain: '<S3>/Gain2'
   *  Inport: '<Root>/StateBus'
   */
  frac[0] = (PX4Controller_P.Gain1_Gain_g * rtb_Gain2 +
             rtb_InterpolationUsingPrelook_4) - PX4Controller_P.Gain2_Gain_c2 *
    PX4Controller_U.StateBus_m.eulerBody[0];

  /* Interpolation_n-D: '<S35>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S35>/Constant'
   */
  if (PX4Controller_P.Constant_Value_bd[0] > 1.0) {
    bpIndex_3[2] = 1U;
  } else if (PX4Controller_P.Constant_Value_bd[0] >= 0.0) {
    bpIndex_3[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_bd[0]);
  } else {
    bpIndex_3[2] = 0U;
  }

  /* Gain: '<S3>/Gain2' incorporates:
   *  Interpolation_n-D: '<S35>/Interpolation Using Prelookup'
   */
  frac_0[0] = intrp2d_l_pw(bpIndex_3, frac_3,
    &PX4Controller_P.InterpolationUsingPrelookup_T_i[bpIndex_3[2] << 2], 2U);

  /* Saturate: '<S3>/Saturation1' */
  if (frac_2[1] > PX4Controller_P.Saturation1_UpperSat_e[1]) {
    rtb_Gain2 = PX4Controller_P.Saturation1_UpperSat_e[1];
  } else if (frac_2[1] < PX4Controller_P.Saturation1_LowerSat_k[1]) {
    rtb_Gain2 = PX4Controller_P.Saturation1_LowerSat_k[1];
  } else {
    rtb_Gain2 = frac_2[1];
  }

  /* Sum: '<S3>/Subtract' incorporates:
   *  Gain: '<S3>/   '
   *  Gain: '<S3>/Gain1'
   *  Gain: '<S3>/Gain2'
   *  Inport: '<Root>/StateBus'
   */
  frac[1] = (PX4Controller_P._Gain_n * rtb_InterpolationUsingPrelook_5 +
             PX4Controller_P.Gain1_Gain_g * rtb_Gain2) -
    PX4Controller_P.Gain2_Gain_c2 * PX4Controller_U.StateBus_m.eulerBody[1];

  /* Interpolation_n-D: '<S35>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S35>/Constant'
   */
  if (PX4Controller_P.Constant_Value_bd[1] > 1.0) {
    bpIndex_3[2] = 1U;
  } else if (PX4Controller_P.Constant_Value_bd[1] >= 0.0) {
    bpIndex_3[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_bd[1]);
  } else {
    bpIndex_3[2] = 0U;
  }

  rtb_Gain2 = intrp2d_l_pw(bpIndex_3, frac_3,
    &PX4Controller_P.InterpolationUsingPrelookup_T_i[bpIndex_3[2] << 2], 2U);

  /* Trigonometry: '<S36>/sincos' incorporates:
   *  Inport: '<Root>/StateBus'
   *  Trigonometry: '<S30>/sincos'
   */
  rtb_InterpolationUsingPrelook_5 = std::sin
    (PX4Controller_U.StateBus_m.eulerBody[0]);
  rtb_InterpolationUsingPrelook_4 = std::cos
    (PX4Controller_U.StateBus_m.eulerBody[0]);

  /* Interpolation_n-D: '<S33>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S33>/Constant'
   */
  frac_4[0] = rtb_f1;
  frac_4[1] = rtb_f2;
  bpIndex_4[0] = rtb_k1;
  bpIndex_4[1] = rtb_k2;
  if (PX4Controller_P.Constant_Value_k[0] > 1.0) {
    bpIndex_4[2] = 1U;
  } else if (PX4Controller_P.Constant_Value_k[0] >= 0.0) {
    bpIndex_4[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_k[0]);
  } else {
    bpIndex_4[2] = 0U;
  }

  rtb_Gain6 = intrp2d_l_pw(bpIndex_4, frac_4,
    &PX4Controller_P.InterpolationUsingPrelookup_T_l[bpIndex_4[2] << 2], 2U);

  /* Trigonometry: '<S36>/sincos' incorporates:
   *  Inport: '<Root>/StateBus'
   *  Trigonometry: '<S30>/sincos'
   */
  rtb_Gain1_g = std::sin(PX4Controller_U.StateBus_m.eulerBody[1]);

  /* Sum: '<S3>/Add' incorporates:
   *  Fcn: '<S36>/phidot'
   *  Gain: '<S3>/Gain4'
   *  Gain: '<S3>/Gain5'
   *  Inport: '<Root>/StateBus'
   *  Integrator: '<S3>/Integrator'
   *  Product: '<S33>/Product'
   *  Product: '<S35>/Product'
   *  Trigonometry: '<S36>/sincos'
   */
  rtb_Tan1 = ((rtb_InterpolationUsingPrelook_5 *
               PX4Controller_U.StateBus_m.wBody[1] +
               rtb_InterpolationUsingPrelook_4 *
               PX4Controller_U.StateBus_m.wBody[2]) * (rtb_Gain1_g /
    rtb_InterpolationUsingPrelook_6) + PX4Controller_U.StateBus_m.wBody[0]) *
    PX4Controller_P.Gain4_Gain_g * PX4Controller_P.Gain5_Gain_i * rtb_Gain6 +
    (frac[0] * frac_0[0] + PX4Controller_X.Integrator_CSTATE[0]);

  /* Interpolation_n-D: '<S33>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S33>/Constant'
   */
  if (PX4Controller_P.Constant_Value_k[1] > 1.0) {
    bpIndex_4[2] = 1U;
  } else if (PX4Controller_P.Constant_Value_k[1] >= 0.0) {
    bpIndex_4[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_k[1]);
  } else {
    bpIndex_4[2] = 0U;
  }

  rtb_Gain6 = intrp2d_l_pw(bpIndex_4, frac_4,
    &PX4Controller_P.InterpolationUsingPrelookup_T_l[bpIndex_4[2] << 2], 2U);

  /* Sum: '<S3>/Add' incorporates:
   *  Fcn: '<S36>/thetadot'
   *  Gain: '<S3>/Gain4'
   *  Gain: '<S3>/Gain5'
   *  Inport: '<Root>/StateBus'
   *  Integrator: '<S3>/Integrator'
   *  Product: '<S33>/Product'
   *  Product: '<S35>/Product'
   *  Trigonometry: '<S36>/sincos'
   */
  rtb_Gain2 = (rtb_InterpolationUsingPrelook_4 *
               PX4Controller_U.StateBus_m.wBody[1] -
               rtb_InterpolationUsingPrelook_5 *
               PX4Controller_U.StateBus_m.wBody[2]) *
    PX4Controller_P.Gain4_Gain_g * PX4Controller_P.Gain5_Gain_i * rtb_Gain6 +
    (frac[1] * rtb_Gain2 + PX4Controller_X.Integrator_CSTATE[1]);

  /* Saturate: '<S32>/Saturation' incorporates:
   *  Inport: '<Root>/StateBus'
   */
  if (PX4Controller_U.StateBus_m.tas > PX4Controller_P.Saturation_UpperSat_ag) {
    rtb_tecsthetaCmd = PX4Controller_P.Saturation_UpperSat_ag;
  } else if (PX4Controller_U.StateBus_m.tas <
             PX4Controller_P.Saturation_LowerSat_kd) {
    rtb_tecsthetaCmd = PX4Controller_P.Saturation_LowerSat_kd;
  } else {
    rtb_tecsthetaCmd = PX4Controller_U.StateBus_m.tas;
  }

  /* End of Saturate: '<S32>/Saturation' */

  /* Gain: '<S3>/Gain6' incorporates:
   *  Gain: '<S32>/Gain1'
   *  Gain: '<S3>/Gain3'
   *  Product: '<S32>/Divide1'
   *  Sum: '<S32>/Add'
   *  Trigonometry: '<S32>/Tan'
   *  Trigonometry: '<S32>/Tan1'
   */
  rtb_Gain6 = (std::tan(rtb_Switch7_idx_0) + std::cos(rtb_Switch7_idx_1)) *
    PX4Controller_P.Gain1_Gain_e / rtb_tecsthetaCmd *
    PX4Controller_P.Gain3_Gain_g * PX4Controller_P.Gain6_Gain_m;

  /* Fcn: '<S30>/phidot' */
  rtb_Gain1_g = rtb_Tan1 - rtb_Gain6 * rtb_Gain1_g;

  /* Fcn: '<S30>/thetadot' */
  rtb_Tan1 = rtb_Gain6 * rtb_InterpolationUsingPrelook_5 *
    rtb_InterpolationUsingPrelook_6 + rtb_InterpolationUsingPrelook_4 *
    rtb_Gain2;

  /* Fcn: '<S30>/psidot' */
  rtb_Gain6 = rtb_Gain6 * rtb_InterpolationUsingPrelook_4 *
    rtb_InterpolationUsingPrelook_6 - rtb_InterpolationUsingPrelook_5 *
    rtb_Gain2;

  /* Switch: '<S70>/Switch7' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  if (PX4Controller_U.CmdBusIn.manualRate > PX4Controller_P.Switch7_Threshold_k)
  {
    /* Outport: '<Root>/CmdBusOut' */
    PX4Controller_Y.CmdBusOut.wCmd[0] = rtb_Gain1_a;
    PX4Controller_Y.CmdBusOut.wCmd[1] = rtb_Gain2_m;
    PX4Controller_Y.CmdBusOut.wCmd[2] = rtb_Switch1_o;
  } else {
    /* Outport: '<Root>/CmdBusOut' */
    PX4Controller_Y.CmdBusOut.wCmd[0] = rtb_Gain1_g;
    PX4Controller_Y.CmdBusOut.wCmd[1] = rtb_Tan1;
    PX4Controller_Y.CmdBusOut.wCmd[2] = rtb_Gain6;
  }

  /* End of Switch: '<S70>/Switch7' */

  /* Interpolation_n-D: '<S71>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S71>/Constant'
   */
  frac_5[0] = rtb_f1_m;
  frac_5[1] = rtb_f2_p;
  frac_5[2] = rtb_InterpolationUsingPrelook_k;
  bpIndex_5[0] = rtb_k1_g;
  bpIndex_5[1] = rtb_k2_e;
  bpIndex_5[2] = rtb_k3;
  if (PX4Controller_P.Constant_Value_a[0] > 2.0) {
    bpIndex_5[3] = 2U;
  } else if (PX4Controller_P.Constant_Value_a[0] >= 0.0) {
    bpIndex_5[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_a[0]);
  } else {
    bpIndex_5[3] = 0U;
  }

  rtb_InterpolationUsingPrelook_4 = intrp3d_l_pw(bpIndex_5, frac_5,
    &PX4Controller_P.SubsystemReference3_table_m[12U * bpIndex_5[3]],
    PX4Controller_P.InterpolationUsingPrelookup_d_f);
  if (PX4Controller_P.Constant_Value_a[1] > 2.0) {
    bpIndex_5[3] = 2U;
  } else if (PX4Controller_P.Constant_Value_a[1] >= 0.0) {
    bpIndex_5[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_a[1]);
  } else {
    bpIndex_5[3] = 0U;
  }

  rtb_InterpolationUsingPrelook_5 = intrp3d_l_pw(bpIndex_5, frac_5,
    &PX4Controller_P.SubsystemReference3_table_m[12U * bpIndex_5[3]],
    PX4Controller_P.InterpolationUsingPrelookup_d_f);
  if (PX4Controller_P.Constant_Value_a[2] > 2.0) {
    bpIndex_5[3] = 2U;
  } else if (PX4Controller_P.Constant_Value_a[2] >= 0.0) {
    bpIndex_5[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_a[2]);
  } else {
    bpIndex_5[3] = 0U;
  }

  rtb_InterpolationUsingPrelook_6 = intrp3d_l_pw(bpIndex_5, frac_5,
    &PX4Controller_P.SubsystemReference3_table_m[12U * bpIndex_5[3]],
    PX4Controller_P.InterpolationUsingPrelookup_d_f);

  /* End of Interpolation_n-D: '<S71>/Interpolation Using Prelookup' */

  /* Switch: '<S70>/Switch2' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   *  Saturate: '<S28>/Saturation'
   */
  if (PX4Controller_U.CmdBusIn.manualRate > PX4Controller_P.Switch2_Threshold_l)
  {
    rtb_Gain1_g = rtb_Gain1_a;
    rtb_Tan1 = rtb_Gain2_m;
    rtb_Gain6 = rtb_Switch1_o;
  } else {
    if (rtb_Gain1_g > PX4Controller_P.Saturation_UpperSat_p[0]) {
      /* Saturate: '<S28>/Saturation' */
      rtb_Gain1_g = PX4Controller_P.Saturation_UpperSat_p[0];
    } else if (rtb_Gain1_g < PX4Controller_P.Saturation_LowerSat_c[0]) {
      /* Saturate: '<S28>/Saturation' */
      rtb_Gain1_g = PX4Controller_P.Saturation_LowerSat_c[0];
    }

    /* Saturate: '<S28>/Saturation' */
    if (rtb_Tan1 > PX4Controller_P.Saturation_UpperSat_p[1]) {
      rtb_Tan1 = PX4Controller_P.Saturation_UpperSat_p[1];
    } else if (rtb_Tan1 < PX4Controller_P.Saturation_LowerSat_c[1]) {
      rtb_Tan1 = PX4Controller_P.Saturation_LowerSat_c[1];
    }

    if (rtb_Gain6 > PX4Controller_P.Saturation_UpperSat_p[2]) {
      rtb_Gain6 = PX4Controller_P.Saturation_UpperSat_p[2];
    } else if (rtb_Gain6 < PX4Controller_P.Saturation_LowerSat_c[2]) {
      rtb_Gain6 = PX4Controller_P.Saturation_LowerSat_c[2];
    }
  }

  /* End of Switch: '<S70>/Switch2' */

  /* Switch: '<S12>/Switch7' incorporates:
   *  Constant: '<S72>/Constant'
   *  Constant: '<S7>/Constant'
   *  Gain: '<S12>/Gain1'
   *  Gain: '<S12>/Gain2'
   *  Gain: '<S12>/Gain3'
   *  Inport: '<Root>/CmdBusIn'
   *  Interpolation_n-D: '<S72>/Interpolation Using Prelookup'
   *  Product: '<S7>/Matrix Multiply'
   */
  if (PX4Controller_U.CmdBusIn.manualFM > PX4Controller_P.Switch7_Threshold_e) {
    rtb_Switch7_e[0] = PX4Controller_P.Gain1_Gain * PX4Controller_U.CmdBusIn.rc
      [1];
    rtb_Switch7_e[1] = PX4Controller_P.Gain2_Gain * PX4Controller_U.CmdBusIn.rc
      [2];
    rtb_Switch7_e[2] = PX4Controller_P.Gain3_Gain * PX4Controller_U.CmdBusIn.rc
      [3];
  } else {
    /* Interpolation_n-D: '<S72>/Interpolation Using Prelookup' */
    frac_6[0] = rtb_f1;
    frac_6[1] = rtb_f2;
    bpIndex_6[0] = rtb_k1;
    bpIndex_6[1] = rtb_k2;
    if (PX4Controller_P.Constant_Value_f[0] > 2.0) {
      bpIndex_6[2] = 2U;
    } else if (PX4Controller_P.Constant_Value_f[0] >= 0.0) {
      bpIndex_6[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_f[0]);
    } else {
      bpIndex_6[2] = 0U;
    }

    /* Interpolation_n-D: '<S72>/Interpolation Using Prelookup' incorporates:
     *  Constant: '<S72>/Constant'
     */
    rtb_Gain2 = intrp2d_l_pw(bpIndex_6, frac_6,
      &PX4Controller_P.InterpolationUsingPrelookup__cs[bpIndex_6[2] << 2], 2U);

    /* Saturate: '<S7>/Saturation' incorporates:
     *  Gain: '<S7>/   '
     */
    if (rtb_Gain1_g > PX4Controller_P.Saturation_UpperSat_c[0]) {
      rtb_Switch1_o = PX4Controller_P.Saturation_UpperSat_c[0];
    } else if (rtb_Gain1_g < PX4Controller_P.Saturation_LowerSat_d[0]) {
      rtb_Switch1_o = PX4Controller_P.Saturation_LowerSat_d[0];
    } else {
      rtb_Switch1_o = rtb_Gain1_g;
    }

    /* Sum: '<S7>/Add' incorporates:
     *  Gain: '<S7>/ '
     *  Gain: '<S7>/   '
     *  Gain: '<S7>/Gain'
     *  Gain: '<S7>/Gain1'
     *  Gain: '<S7>/Gain2'
     *  Gain: '<S7>/Gain3'
     *  Inport: '<Root>/StateBus'
     *  Product: '<S72>/Product'
     *  Sum: '<S7>/Subtract'
     */
    rtb_Switch1_o = ((PX4Controller_P._Gain_m * rtb_InterpolationUsingPrelook_4
                      + PX4Controller_P._Gain * rtb_Switch1_o) -
                     PX4Controller_P.Gain_Gain_a *
                     PX4Controller_U.StateBus_m.wBody[0]) * rtb_Gain2 +
      PX4Controller_P.Gain1_Gain_f * PX4Controller_U.StateBus_m.wDotBody[0] *
      PX4Controller_P.Gain2_Gain_k * PX4Controller_P.Gain3_Gain_c;
    if (PX4Controller_P.Constant_Value_f[1] > 2.0) {
      bpIndex_6[2] = 2U;
    } else if (PX4Controller_P.Constant_Value_f[1] >= 0.0) {
      bpIndex_6[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_f[1]);
    } else {
      bpIndex_6[2] = 0U;
    }

    /* Interpolation_n-D: '<S72>/Interpolation Using Prelookup' incorporates:
     *  Constant: '<S72>/Constant'
     */
    rtb_Gain2 = intrp2d_l_pw(bpIndex_6, frac_6,
      &PX4Controller_P.InterpolationUsingPrelookup__cs[bpIndex_6[2] << 2], 2U);

    /* Saturate: '<S7>/Saturation' incorporates:
     *  Gain: '<S7>/   '
     */
    if (rtb_Tan1 > PX4Controller_P.Saturation_UpperSat_c[1]) {
      rtb_Gain1_a = PX4Controller_P.Saturation_UpperSat_c[1];
    } else if (rtb_Tan1 < PX4Controller_P.Saturation_LowerSat_d[1]) {
      rtb_Gain1_a = PX4Controller_P.Saturation_LowerSat_d[1];
    } else {
      rtb_Gain1_a = rtb_Tan1;
    }

    /* Sum: '<S7>/Add' incorporates:
     *  Gain: '<S7>/ '
     *  Gain: '<S7>/   '
     *  Gain: '<S7>/Gain'
     *  Gain: '<S7>/Gain1'
     *  Gain: '<S7>/Gain2'
     *  Gain: '<S7>/Gain3'
     *  Inport: '<Root>/StateBus'
     *  Product: '<S72>/Product'
     *  Sum: '<S7>/Subtract'
     */
    rtb_Gain1_a = ((PX4Controller_P._Gain_m * rtb_InterpolationUsingPrelook_5 +
                    PX4Controller_P._Gain * rtb_Gain1_a) -
                   PX4Controller_P.Gain_Gain_a *
                   PX4Controller_U.StateBus_m.wBody[1]) * rtb_Gain2 +
      PX4Controller_P.Gain1_Gain_f * PX4Controller_U.StateBus_m.wDotBody[1] *
      PX4Controller_P.Gain2_Gain_k * PX4Controller_P.Gain3_Gain_c;
    if (PX4Controller_P.Constant_Value_f[2] > 2.0) {
      bpIndex_6[2] = 2U;
    } else if (PX4Controller_P.Constant_Value_f[2] >= 0.0) {
      bpIndex_6[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_f[2]);
    } else {
      bpIndex_6[2] = 0U;
    }

    /* Interpolation_n-D: '<S72>/Interpolation Using Prelookup' incorporates:
     *  Constant: '<S72>/Constant'
     */
    rtb_Gain2 = intrp2d_l_pw(bpIndex_6, frac_6,
      &PX4Controller_P.InterpolationUsingPrelookup__cs[bpIndex_6[2] << 2], 2U);

    /* Saturate: '<S7>/Saturation' incorporates:
     *  Gain: '<S7>/   '
     */
    if (rtb_Gain6 > PX4Controller_P.Saturation_UpperSat_c[2]) {
      rtb_Gain2_m = PX4Controller_P.Saturation_UpperSat_c[2];
    } else if (rtb_Gain6 < PX4Controller_P.Saturation_LowerSat_d[2]) {
      rtb_Gain2_m = PX4Controller_P.Saturation_LowerSat_d[2];
    } else {
      rtb_Gain2_m = rtb_Gain6;
    }

    /* Sum: '<S7>/Add' incorporates:
     *  Gain: '<S7>/ '
     *  Gain: '<S7>/   '
     *  Gain: '<S7>/Gain'
     *  Gain: '<S7>/Gain1'
     *  Gain: '<S7>/Gain2'
     *  Gain: '<S7>/Gain3'
     *  Inport: '<Root>/StateBus'
     *  Product: '<S72>/Product'
     *  Sum: '<S7>/Subtract'
     */
    rtb_Gain2 = ((PX4Controller_P._Gain_m * rtb_InterpolationUsingPrelook_6 +
                  PX4Controller_P._Gain * rtb_Gain2_m) -
                 PX4Controller_P.Gain_Gain_a * PX4Controller_U.StateBus_m.wBody
                 [2]) * rtb_Gain2 + PX4Controller_P.Gain1_Gain_f *
      PX4Controller_U.StateBus_m.wDotBody[2] * PX4Controller_P.Gain2_Gain_k *
      PX4Controller_P.Gain3_Gain_c;
    for (int32_T i{0}; i < 3; i++) {
      rtb_Switch7_e[i] = (PX4Controller_P.Constant_Value_i[i + 3] * rtb_Gain1_a
                          + PX4Controller_P.Constant_Value_i[i] * rtb_Switch1_o)
        + PX4Controller_P.Constant_Value_i[i + 6] * rtb_Gain2;
    }
  }

  /* End of Switch: '<S12>/Switch7' */

  /* Sum: '<S77>/Sum' incorporates:
   *  Sum: '<S77>/ste'
   *  Sum: '<S77>/steCmd'
   */
  PX4Controller_B.steErr = (rtb_speCmd + rtb_skeCmd) - (rtb_spe +
    rtb_InterpolationUsingPrelook_0);

  /* Integrator: '<S77>/Integrator' */
  if (rtsiIsModeUpdateTimeStep(&PX4Controller_M->solverInfo)) {
    didZcEventOccur = (((PX4Controller_PrevZCX.Integrator_Reset_ZCE_l ==
                         POS_ZCSIG) != PX4Controller_B.OR) &&
                       (PX4Controller_PrevZCX.Integrator_Reset_ZCE_l !=
                        UNINITIALIZED_ZCSIG));
    PX4Controller_PrevZCX.Integrator_Reset_ZCE_l = PX4Controller_B.OR;

    /* evaluate zero-crossings and the level of the reset signal */
    if (didZcEventOccur || PX4Controller_B.OR) {
      PX4Controller_X.Integrator_CSTATE_f = PX4Controller_P.Integrator_IC_a;
    }
  }

  /* Switch: '<S12>/Switch1' incorporates:
   *  Gain: '<S12>/Gain'
   *  Inport: '<Root>/CmdBusIn'
   *  Switch: '<S27>/Switch1'
   *  Switch: '<S70>/Switch1'
   */
  if (PX4Controller_U.CmdBusIn.manualFM > PX4Controller_P.Switch1_Threshold_l) {
    rtb_Switch1_o = PX4Controller_P.Gain_Gain * PX4Controller_U.CmdBusIn.rc[0];
  } else if (PX4Controller_U.CmdBusIn.manualRate >
             PX4Controller_P.Switch1_Threshold_g) {
    /* Switch: '<S70>/Switch1' incorporates:
     *  Gain: '<S70>/Gain'
     */
    rtb_Switch1_o = PX4Controller_P.Gain_Gain_d * PX4Controller_U.CmdBusIn.rc[0];
  } else if (PX4Controller_U.CmdBusIn.manualAttitude >
             PX4Controller_P.Switch1_Threshold) {
    /* Switch: '<S27>/Switch1' incorporates:
     *  Gain: '<S27>/Gain'
     *  Switch: '<S70>/Switch1'
     */
    rtb_Switch1_o = PX4Controller_P.Gain_Gain_g * PX4Controller_U.CmdBusIn.rc[0];
  } else {
    /* Interpolation_n-D: '<S81>/Interpolation Using Prelookup' incorporates:
     *  Switch: '<S27>/Switch1'
     *  Switch: '<S70>/Switch1'
     */
    frac_7[0] = rtb_f1;
    frac_7[1] = rtb_f2;
    bpIndex_7[0] = rtb_k1;
    bpIndex_7[1] = rtb_k2;

    /* Interpolation_n-D: '<S83>/Interpolation Using Prelookup' incorporates:
     *  Switch: '<S27>/Switch1'
     *  Switch: '<S70>/Switch1'
     */
    frac_a[0] = rtb_f1;
    frac_a[1] = rtb_f2;
    bpIndex_a[0] = rtb_k1;
    bpIndex_a[1] = rtb_k2;

    /* Interpolation_n-D: '<S82>/Interpolation Using Prelookup' incorporates:
     *  Switch: '<S27>/Switch1'
     *  Switch: '<S70>/Switch1'
     */
    frac_b[0] = rtb_f1;
    frac_b[1] = rtb_f2;
    bpIndex_b[0] = rtb_k1;
    bpIndex_b[1] = rtb_k2;

    /* Saturate: '<S77>/Saturation' incorporates:
     *  Inport: '<Root>/StateBus'
     *  Switch: '<S27>/Switch1'
     *  Switch: '<S70>/Switch1'
     */
    if (PX4Controller_U.StateBus_m.tas > PX4Controller_P.Saturation_UpperSat) {
      rtb_tecsthetaCmd = PX4Controller_P.Saturation_UpperSat;
    } else if (PX4Controller_U.StateBus_m.tas <
               PX4Controller_P.Saturation_LowerSat) {
      rtb_tecsthetaCmd = PX4Controller_P.Saturation_LowerSat;
    } else {
      rtb_tecsthetaCmd = PX4Controller_U.StateBus_m.tas;
    }

    /* End of Saturate: '<S77>/Saturation' */

    /* Switch: '<S27>/Switch1' incorporates:
     *  Gain: '<S77>/M'
     *  Integrator: '<S77>/Integrator'
     *  Interpolation_n-D: '<S81>/Interpolation Using Prelookup'
     *  Interpolation_n-D: '<S82>/Interpolation Using Prelookup'
     *  Interpolation_n-D: '<S83>/Interpolation Using Prelookup'
     *  Product: '<S77>/Divide'
     *  Product: '<S81>/Product'
     *  Product: '<S82>/Product'
     *  Product: '<S83>/Product'
     *  Sum: '<S77>/Sum2'
     *  Switch: '<S70>/Switch1'
     *  TransferFcn: '<S77>/Low pass'
     */
    rtb_Switch1_o = ((PX4Controller_P.Lowpass_C * PX4Controller_X.Lowpass_CSTATE
                      * intrp2d_l_pw(bpIndex_b, frac_b,
      PX4Controller_P.InterpolationUsingPrelookup_T_c, 2U) +
                      PX4Controller_B.steErr * intrp2d_l_pw(bpIndex_a, frac_a,
      PX4Controller_P.InterpolationUsingPrelookup_T_n, 2U)) +
                     PX4Controller_X.Integrator_CSTATE_f * intrp2d_l_pw
                     (bpIndex_7, frac_7,
                      PX4Controller_P.InterpolationUsingPrelookup_Tab, 2U)) *
      PX4Controller_P.M_Gain / rtb_tecsthetaCmd;
  }

  /* End of Switch: '<S12>/Switch1' */

  /* Interpolation_n-D: '<S18>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S18>/Constant'
   */
  frac_8[0] = rtb_f1;
  frac_8[1] = rtb_f2;
  bpIndex_8[0] = rtb_k1;
  bpIndex_8[1] = rtb_k2;
  for (int32_T i{0}; i < 84; i++) {
    if (PX4Controller_P.Constant_Value_n[i] > 83.0) {
      bpIndex_8[2] = 83U;
    } else if (PX4Controller_P.Constant_Value_n[i] >= 0.0) {
      bpIndex_8[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_n[i]);
    } else {
      bpIndex_8[2] = 0U;
    }

    rtb_InterpolationUsingPreloo_dw[i] = intrp2d_l_pw(bpIndex_8, frac_8,
      &PX4Controller_P.InterpolationUsingPrelookup_T_p[bpIndex_8[2] << 2], 2U);
  }

  /* End of Interpolation_n-D: '<S18>/Interpolation Using Prelookup' */

  /* Interpolation_n-D: '<S25>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S25>/Constant'
   */
  frac_9[0] = rtb_f1;
  frac_9[1] = rtb_f2;
  bpIndex_9[0] = rtb_k1;
  bpIndex_9[1] = rtb_k2;
  if (PX4Controller_P.Constant_Value_np[0] > 3.0) {
    bpIndex_9[2] = 3U;
  } else if (PX4Controller_P.Constant_Value_np[0] >= 0.0) {
    bpIndex_9[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_np[0]);
  } else {
    bpIndex_9[2] = 0U;
  }

  rtb_InterpolationUsingPrelook_0 = intrp2d_l_pw(bpIndex_9, frac_9,
    &PX4Controller_P.InterpolationUsingPrelookup_T_m[bpIndex_9[2] << 2], 2U);
  if (PX4Controller_P.Constant_Value_np[1] > 3.0) {
    bpIndex_9[2] = 3U;
  } else if (PX4Controller_P.Constant_Value_np[1] >= 0.0) {
    bpIndex_9[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_np[1]);
  } else {
    bpIndex_9[2] = 0U;
  }

  rtb_spe = intrp2d_l_pw(bpIndex_9, frac_9,
    &PX4Controller_P.InterpolationUsingPrelookup_T_m[bpIndex_9[2] << 2], 2U);
  if (PX4Controller_P.Constant_Value_np[2] > 3.0) {
    bpIndex_9[2] = 3U;
  } else if (PX4Controller_P.Constant_Value_np[2] >= 0.0) {
    bpIndex_9[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_np[2]);
  } else {
    bpIndex_9[2] = 0U;
  }

  rtb_skeCmd = intrp2d_l_pw(bpIndex_9, frac_9,
    &PX4Controller_P.InterpolationUsingPrelookup_T_m[bpIndex_9[2] << 2], 2U);
  if (PX4Controller_P.Constant_Value_np[3] > 3.0) {
    bpIndex_9[2] = 3U;
  } else if (PX4Controller_P.Constant_Value_np[3] >= 0.0) {
    bpIndex_9[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_np[3]);
  } else {
    bpIndex_9[2] = 0U;
  }

  rtb_speCmd = intrp2d_l_pw(bpIndex_9, frac_9,
    &PX4Controller_P.InterpolationUsingPrelookup_T_m[bpIndex_9[2] << 2], 2U);

  /* End of Interpolation_n-D: '<S25>/Interpolation Using Prelookup' */

  /* Assignment: '<S14>/Assignment' incorporates:
   *  Assignment: '<S14>/Assignment2'
   *  Constant: '<S14>/Constant'
   *  Interpolation_n-D: '<S19>/Interpolation Using Prelookup'
   */
  std::memcpy(&rtb_Assignment2[0], &PX4Controller_P.Constant_Value_bp[0], 84U *
              sizeof(real_T));
  rtb_Assignment2[0] = rtb_InterpolationUsingPrelook_0;
  rtb_Assignment2[1] = rtb_spe;
  rtb_Assignment2[2] = rtb_skeCmd;
  rtb_Assignment2[3] = rtb_speCmd;

  /* SignalConversion generated from: '<S13>/Vector Concatenate' incorporates:
   *  Concatenate: '<S13>/Vector Concatenate'
   */
  PX4Controller_B.VectorConcatenate[0] = rtb_Switch1_o;
  if (rtmIsMajorTimeStep(PX4Controller_M)) {
    /* Constant: '<S13>/Constant1' incorporates:
     *  Concatenate: '<S13>/Vector Concatenate'
     */
    PX4Controller_B.VectorConcatenate[1] = PX4Controller_P.Constant1_Value_e[0];
    PX4Controller_B.VectorConcatenate[2] = PX4Controller_P.Constant1_Value_e[1];
    for (int32_T i{0}; i < 6; i++) {
      /* Gain: '<S15>/Gain' incorporates:
       *  Constant: '<S2>/Constant'
       */
      PX4Controller_B.Gain[i] = PX4Controller_P.SubsystemReference_fmMask[i] *
        PX4Controller_P.Constant_Value_hr[i];
    }
  }

  /* SignalConversion generated from: '<S13>/Vector Concatenate' incorporates:
   *  Concatenate: '<S13>/Vector Concatenate'
   */
  PX4Controller_B.VectorConcatenate[3] = rtb_Switch7_e[0];
  PX4Controller_B.VectorConcatenate[4] = rtb_Switch7_e[1];
  PX4Controller_B.VectorConcatenate[5] = rtb_Switch7_e[2];

  /* Sum: '<S15>/Subtract' incorporates:
   *  Gain: '<S15>/Gain1'
   */
  for (int32_T i{0}; i < 6; i++) {
    rtb_Subtract1[i] = PX4Controller_P.SubsystemReference_fmMask[i] *
      PX4Controller_B.VectorConcatenate[i] - PX4Controller_B.Gain[i];
  }

  /* End of Sum: '<S15>/Subtract' */

  /* Product: '<S15>/Matrix Multiply2' incorporates:
   *  Assignment: '<S14>/Assignment2'
   */
  std::memset(&rtb_Saturation_az[0], 0, 14U * sizeof(real_T));
  for (int32_T i_0{0}; i_0 < 6; i_0++) {
    for (int32_T i{0}; i < 14; i++) {
      rtb_Saturation_az[i] += rtb_Assignment2[14 * i_0 + i] * rtb_Subtract1[i_0];
    }
  }

  /* End of Product: '<S15>/Matrix Multiply2' */

  /* Interpolation_n-D: '<S19>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S19>/Constant'
   */
  frac_c[0] = rtb_f1_m;
  frac_c[1] = rtb_f2_p;
  frac_c[2] = rtb_InterpolationUsingPrelook_k;
  bpIndex_c[0] = rtb_k1_g;
  bpIndex_c[1] = rtb_k2_e;
  bpIndex_c[2] = rtb_k3;
  if (PX4Controller_P.Constant_Value_e[0] > 3.0) {
    bpIndex_c[3] = 3U;
  } else if (PX4Controller_P.Constant_Value_e[0] >= 0.0) {
    bpIndex_c[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_e[0]);
  } else {
    bpIndex_c[3] = 0U;
  }

  /* SignalConversion generated from: '<S15>/Matrix Multiply' incorporates:
   *  Interpolation_n-D: '<S19>/Interpolation Using Prelookup'
   */
  rtb_Saturation_j[0] = intrp3d_l_pw(bpIndex_c, frac_c,
    &PX4Controller_P.SubsystemReference3_table_p[12U * bpIndex_c[3]],
    PX4Controller_P.InterpolationUsingPrelookup_d_d);

  /* Interpolation_n-D: '<S19>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S19>/Constant'
   */
  if (PX4Controller_P.Constant_Value_e[1] > 3.0) {
    bpIndex_c[3] = 3U;
  } else if (PX4Controller_P.Constant_Value_e[1] >= 0.0) {
    bpIndex_c[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_e[1]);
  } else {
    bpIndex_c[3] = 0U;
  }

  /* SignalConversion generated from: '<S15>/Matrix Multiply' incorporates:
   *  Interpolation_n-D: '<S19>/Interpolation Using Prelookup'
   */
  rtb_Saturation_j[1] = intrp3d_l_pw(bpIndex_c, frac_c,
    &PX4Controller_P.SubsystemReference3_table_p[12U * bpIndex_c[3]],
    PX4Controller_P.InterpolationUsingPrelookup_d_d);

  /* Interpolation_n-D: '<S19>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S19>/Constant'
   */
  if (PX4Controller_P.Constant_Value_e[2] > 3.0) {
    bpIndex_c[3] = 3U;
  } else if (PX4Controller_P.Constant_Value_e[2] >= 0.0) {
    bpIndex_c[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_e[2]);
  } else {
    bpIndex_c[3] = 0U;
  }

  /* SignalConversion generated from: '<S15>/Matrix Multiply' incorporates:
   *  Interpolation_n-D: '<S19>/Interpolation Using Prelookup'
   */
  rtb_Saturation_j[2] = intrp3d_l_pw(bpIndex_c, frac_c,
    &PX4Controller_P.SubsystemReference3_table_p[12U * bpIndex_c[3]],
    PX4Controller_P.InterpolationUsingPrelookup_d_d);

  /* Interpolation_n-D: '<S19>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S19>/Constant'
   */
  if (PX4Controller_P.Constant_Value_e[3] > 3.0) {
    bpIndex_c[3] = 3U;
  } else if (PX4Controller_P.Constant_Value_e[3] >= 0.0) {
    bpIndex_c[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_e[3]);
  } else {
    bpIndex_c[3] = 0U;
  }

  /* SignalConversion generated from: '<S15>/Matrix Multiply' incorporates:
   *  Interpolation_n-D: '<S19>/Interpolation Using Prelookup'
   */
  rtb_Saturation_j[3] = intrp3d_l_pw(bpIndex_c, frac_c,
    &PX4Controller_P.SubsystemReference3_table_p[12U * bpIndex_c[3]],
    PX4Controller_P.InterpolationUsingPrelookup_d_d);

  /* Interpolation_n-D: '<S20>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S20>/Constant'
   */
  frac_d[0] = rtb_f1_m;
  frac_d[1] = rtb_f2_p;
  frac_d[2] = rtb_InterpolationUsingPrelook_k;
  bpIndex_d[0] = rtb_k1_g;
  bpIndex_d[1] = rtb_k2_e;
  bpIndex_d[2] = rtb_k3;
  for (int32_T i{0}; i < 6; i++) {
    if (PX4Controller_P.Constant_Value_m[i] > 5.0) {
      bpIndex_d[3] = 5U;
    } else if (PX4Controller_P.Constant_Value_m[i] >= 0.0) {
      bpIndex_d[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_m[i]);
    } else {
      bpIndex_d[3] = 0U;
    }

    rtb_Subtract1[i] = intrp3d_l_pw(bpIndex_d, frac_d,
      &PX4Controller_P.SubsystemReference4_table[12U * bpIndex_d[3]],
      PX4Controller_P.InterpolationUsingPrelookup_d_o);
  }

  /* End of Interpolation_n-D: '<S20>/Interpolation Using Prelookup' */

  /* Interpolation_n-D: '<S21>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S21>/Constant'
   */
  frac_e[0] = rtb_f1_m;
  frac_e[1] = rtb_f2_p;
  frac_e[2] = rtb_InterpolationUsingPrelook_k;
  bpIndex_e[0] = rtb_k1_g;
  bpIndex_e[1] = rtb_k2_e;
  bpIndex_e[2] = rtb_k3;
  if (PX4Controller_P.Constant_Value_e4[0] > 1.0) {
    bpIndex_e[3] = 1U;
  } else if (PX4Controller_P.Constant_Value_e4[0] >= 0.0) {
    bpIndex_e[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_e4[0]);
  } else {
    bpIndex_e[3] = 0U;
  }

  frac_0[0] = intrp3d_l_pw(bpIndex_e, frac_e,
    &PX4Controller_P.SubsystemReference5_table[12U * bpIndex_e[3]],
    PX4Controller_P.InterpolationUsingPrelookup_d_m);
  if (PX4Controller_P.Constant_Value_e4[1] > 1.0) {
    bpIndex_e[3] = 1U;
  } else if (PX4Controller_P.Constant_Value_e4[1] >= 0.0) {
    bpIndex_e[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_e4[1]);
  } else {
    bpIndex_e[3] = 0U;
  }

  frac_0[1] = intrp3d_l_pw(bpIndex_e, frac_e,
    &PX4Controller_P.SubsystemReference5_table[12U * bpIndex_e[3]],
    PX4Controller_P.InterpolationUsingPrelookup_d_m);

  /* End of Interpolation_n-D: '<S21>/Interpolation Using Prelookup' */

  /* Interpolation_n-D: '<S22>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S22>/Constant'
   */
  frac_f[0] = rtb_f1_m;
  frac_f[1] = rtb_f2_p;
  frac_f[2] = rtb_InterpolationUsingPrelook_k;
  bpIndex_f[0] = rtb_k1_g;
  bpIndex_f[1] = rtb_k2_e;
  bpIndex_f[2] = rtb_k3;
  if (PX4Controller_P.Constant_Value_ns[0] > 1.0) {
    bpIndex_f[3] = 1U;
  } else if (PX4Controller_P.Constant_Value_ns[0] >= 0.0) {
    bpIndex_f[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_ns[0]);
  } else {
    bpIndex_f[3] = 0U;
  }

  frac_3[0] = intrp3d_l_pw(bpIndex_f, frac_f,
    &PX4Controller_P.SubsystemReference6_table[12U * bpIndex_f[3]],
    PX4Controller_P.InterpolationUsingPrelookup__mi);
  if (PX4Controller_P.Constant_Value_ns[1] > 1.0) {
    bpIndex_f[3] = 1U;
  } else if (PX4Controller_P.Constant_Value_ns[1] >= 0.0) {
    bpIndex_f[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_ns[1]);
  } else {
    bpIndex_f[3] = 0U;
  }

  frac_3[1] = intrp3d_l_pw(bpIndex_f, frac_f,
    &PX4Controller_P.SubsystemReference6_table[12U * bpIndex_f[3]],
    PX4Controller_P.InterpolationUsingPrelookup__mi);

  /* End of Interpolation_n-D: '<S22>/Interpolation Using Prelookup' */

  /* SignalConversion generated from: '<S15>/Matrix Multiply' */
  for (int32_T i{0}; i < 6; i++) {
    rtb_Saturation_j[i + 4] = rtb_Subtract1[i];
  }

  rtb_Saturation_j[10] = frac_0[0];
  rtb_Saturation_j[12] = frac_3[0];
  rtb_Saturation_j[11] = frac_0[1];
  rtb_Saturation_j[13] = frac_3[1];

  /* Saturate: '<S15>/Saturation' incorporates:
   *  Sum: '<S15>/Sum1'
   */
  for (int32_T i{0}; i < 14; i++) {
    /* Sum: '<S15>/Sum1' */
    rtb_f2_p = rtb_Saturation_az[i] + rtb_Saturation_j[i];

    /* Saturate: '<S15>/Saturation' incorporates:
     *  Sum: '<S15>/Sum1'
     */
    rtb_f1_m = PX4Controller_P.Saturation_LowerSat_bc[i];
    rtb_InterpolationUsingPrelook_k = PX4Controller_P.Saturation_UpperSat_pe[i];
    if (rtb_f2_p > rtb_InterpolationUsingPrelook_k) {
      rtb_f2_p = rtb_InterpolationUsingPrelook_k;
    } else if (rtb_f2_p < rtb_f1_m) {
      rtb_f2_p = rtb_f1_m;
    }

    rtb_Saturation_az[i] = rtb_f2_p;
  }

  /* End of Saturate: '<S15>/Saturation' */
  for (int32_T i{0}; i < 6; i++) {
    /* Product: '<S15>/Matrix Multiply3' incorporates:
     *  Product: '<S16>/Matrix Multiply'
     */
    rtb_Subtract1_tmp[i] = 0.0;
  }

  for (int32_T i_0{0}; i_0 < 14; i_0++) {
    for (int32_T i{0}; i < 6; i++) {
      /* Product: '<S15>/Matrix Multiply3' incorporates:
       *  Interpolation_n-D: '<S18>/Interpolation Using Prelookup'
       */
      rtb_Subtract1_tmp[i] += rtb_InterpolationUsingPreloo_dw[6 * i_0 + i] *
        rtb_Saturation_az[i_0];
    }
  }

  /* Sum: '<S15>/Subtract1' incorporates:
   *  Constant: '<S2>/Constant'
   *  Interpolation_n-D: '<S18>/Interpolation Using Prelookup'
   *  Product: '<S15>/Matrix Multiply'
   *  Product: '<S15>/Matrix Multiply3'
   */
  for (int32_T i{0}; i < 6; i++) {
    rtb_tecsthetaCmd = 0.0;
    for (int32_T i_0{0}; i_0 < 14; i_0++) {
      rtb_tecsthetaCmd += rtb_InterpolationUsingPreloo_dw[6 * i_0 + i] *
        rtb_Saturation_j[i_0];
    }

    rtb_Subtract1[i] = (PX4Controller_P.Constant_Value_hr[i] +
                        rtb_Subtract1_tmp[i]) - rtb_tecsthetaCmd;
  }

  /* End of Sum: '<S15>/Subtract1' */

  /* Interpolation_n-D: '<S23>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S23>/Constant'
   */
  frac_g[0] = rtb_f1;
  frac_g[1] = rtb_f2;
  bpIndex_g[0] = rtb_k1;
  bpIndex_g[1] = rtb_k2;
  if (PX4Controller_P.Constant_Value_ke[0] > 1.0) {
    bpIndex_g[2] = 1U;
  } else if (PX4Controller_P.Constant_Value_ke[0] >= 0.0) {
    bpIndex_g[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_ke[0]);
  } else {
    bpIndex_g[2] = 0U;
  }

  frac_0[0] = intrp2d_l_pw(bpIndex_g, frac_g,
    &PX4Controller_P.InterpolationUsingPrelookup__pv[bpIndex_g[2] << 2], 2U);
  if (PX4Controller_P.Constant_Value_ke[1] > 1.0) {
    bpIndex_g[2] = 1U;
  } else if (PX4Controller_P.Constant_Value_ke[1] >= 0.0) {
    bpIndex_g[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_ke[1]);
  } else {
    bpIndex_g[2] = 0U;
  }

  frac_0[1] = intrp2d_l_pw(bpIndex_g, frac_g,
    &PX4Controller_P.InterpolationUsingPrelookup__pv[bpIndex_g[2] << 2], 2U);

  /* End of Interpolation_n-D: '<S23>/Interpolation Using Prelookup' */

  /* Assignment: '<S14>/Assignment1' incorporates:
   *  Assignment: '<S14>/Assignment2'
   *  Constant: '<S14>/Constant'
   *  Interpolation_n-D: '<S34>/Interpolation Using Prelookup'
   */
  std::memcpy(&rtb_Assignment2[0], &PX4Controller_P.Constant_Value_bp[0], 84U *
              sizeof(real_T));
  rtb_Assignment2[82] = frac_0[0];
  rtb_Assignment2[83] = frac_0[1];

  /* Sum: '<S16>/Subtract' incorporates:
   *  Gain: '<S16>/Gain'
   *  Gain: '<S16>/Gain1'
   */
  for (int32_T i{0}; i < 6; i++) {
    tmp[i] = PX4Controller_P.SubsystemReference1_fmMask[i] *
      PX4Controller_B.VectorConcatenate[i] -
      PX4Controller_P.SubsystemReference1_fmMask[i] * rtb_Subtract1[i];
  }

  /* End of Sum: '<S16>/Subtract' */
  for (int32_T i{0}; i < 14; i++) {
    /* Sum: '<S16>/Sum1' incorporates:
     *  Assignment: '<S14>/Assignment2'
     *  Product: '<S16>/Matrix Multiply2'
     */
    rtb_tecsthetaCmd = 0.0;
    for (int32_T i_0{0}; i_0 < 6; i_0++) {
      rtb_tecsthetaCmd += rtb_Assignment2[14 * i_0 + i] * tmp[i_0];
    }

    rtb_f2_p = rtb_tecsthetaCmd + rtb_Saturation_az[i];

    /* End of Sum: '<S16>/Sum1' */

    /* Saturate: '<S16>/Saturation' */
    rtb_f1_m = PX4Controller_P.Saturation_LowerSat_c1[i];
    rtb_InterpolationUsingPrelook_k = PX4Controller_P.Saturation_UpperSat_cb[i];
    if (rtb_f2_p > rtb_InterpolationUsingPrelook_k) {
      rtb_Saturation_j[i] = rtb_InterpolationUsingPrelook_k;
    } else if (rtb_f2_p < rtb_f1_m) {
      rtb_Saturation_j[i] = rtb_f1_m;
    } else {
      rtb_Saturation_j[i] = rtb_f2_p;
    }

    /* End of Saturate: '<S16>/Saturation' */
  }

  for (int32_T i{0}; i < 6; i++) {
    /* Product: '<S16>/Matrix Multiply3' incorporates:
     *  Product: '<S17>/Matrix Multiply'
     */
    rtb_Subtract1_tmp_0[i] = 0.0;
  }

  for (int32_T i_0{0}; i_0 < 14; i_0++) {
    for (int32_T i{0}; i < 6; i++) {
      /* Product: '<S16>/Matrix Multiply3' incorporates:
       *  Interpolation_n-D: '<S18>/Interpolation Using Prelookup'
       */
      rtb_Subtract1_tmp_0[i] += rtb_InterpolationUsingPreloo_dw[6 * i_0 + i] *
        rtb_Saturation_j[i_0];
    }
  }

  /* Sum: '<S16>/Subtract1' incorporates:
   *  Product: '<S16>/Matrix Multiply'
   *  Product: '<S16>/Matrix Multiply3'
   */
  for (int32_T i{0}; i < 6; i++) {
    rtb_Subtract1[i] = (rtb_Subtract1[i] + rtb_Subtract1_tmp_0[i]) -
      rtb_Subtract1_tmp[i];
  }

  /* End of Sum: '<S16>/Subtract1' */

  /* Interpolation_n-D: '<S24>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S24>/Constant'
   */
  frac_h[0] = rtb_f1;
  frac_h[1] = rtb_f2;
  bpIndex_h[0] = rtb_k1;
  bpIndex_h[1] = rtb_k2;
  for (int32_T i{0}; i < 16; i++) {
    if (PX4Controller_P.Constant_Value_l[i] > 15.0) {
      bpIndex_h[2] = 15U;
    } else if (PX4Controller_P.Constant_Value_l[i] >= 0.0) {
      bpIndex_h[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_l[i]);
    } else {
      bpIndex_h[2] = 0U;
    }

    rtb_InterpolationUsingPrelook_c[i] = intrp2d_l_pw(bpIndex_h, frac_h,
      &PX4Controller_P.InterpolationUsingPrelookup__m0[bpIndex_h[2] << 2], 2U);
  }

  /* End of Interpolation_n-D: '<S24>/Interpolation Using Prelookup' */

  /* Assignment: '<S14>/Assignment2' incorporates:
   *  Constant: '<S14>/Constant'
   *  Interpolation_n-D: '<S24>/Interpolation Using Prelookup'
   */
  std::memcpy(&rtb_Assignment2[0], &PX4Controller_P.Constant_Value_bp[0], 84U *
              sizeof(real_T));
  for (int32_T i{0}; i < 2; i++) {
    std::memcpy(&rtb_Assignment2[i * 14 + 46],
                &rtb_InterpolationUsingPrelook_c[i << 3], sizeof(real_T) << 3U);
  }

  /* End of Assignment: '<S14>/Assignment2' */

  /* Sum: '<S17>/Subtract' incorporates:
   *  Gain: '<S17>/Gain'
   *  Gain: '<S17>/Gain1'
   */
  for (int32_T i{0}; i < 6; i++) {
    tmp[i] = PX4Controller_P.SubsystemReference2_fmMask[i] *
      PX4Controller_B.VectorConcatenate[i] -
      PX4Controller_P.SubsystemReference2_fmMask[i] * rtb_Subtract1[i];
  }

  /* End of Sum: '<S17>/Subtract' */
  for (int32_T i{0}; i < 14; i++) {
    /* Sum: '<S17>/Sum1' incorporates:
     *  Assignment: '<S14>/Assignment2'
     *  Product: '<S17>/Matrix Multiply2'
     */
    rtb_tecsthetaCmd = 0.0;
    for (int32_T i_0{0}; i_0 < 6; i_0++) {
      rtb_tecsthetaCmd += rtb_Assignment2[14 * i_0 + i] * tmp[i_0];
    }

    rtb_f2_p = rtb_tecsthetaCmd + rtb_Saturation_j[i];

    /* End of Sum: '<S17>/Sum1' */

    /* Saturate: '<S17>/Saturation' */
    rtb_f1_m = PX4Controller_P.Saturation_LowerSat_ey[i];
    rtb_InterpolationUsingPrelook_k = PX4Controller_P.Saturation_UpperSat_l4[i];
    if (rtb_f2_p > rtb_InterpolationUsingPrelook_k) {
      rtb_Saturation_az[i] = rtb_InterpolationUsingPrelook_k;
    } else if (rtb_f2_p < rtb_f1_m) {
      rtb_Saturation_az[i] = rtb_f1_m;
    } else {
      rtb_Saturation_az[i] = rtb_f2_p;
    }

    /* End of Saturate: '<S17>/Saturation' */
  }

  /* Sum: '<S17>/Subtract1' incorporates:
   *  Interpolation_n-D: '<S18>/Interpolation Using Prelookup'
   *  Product: '<S17>/Matrix Multiply'
   *  Product: '<S17>/Matrix Multiply3'
   */
  for (int32_T i{0}; i < 6; i++) {
    rtb_tecsthetaCmd = 0.0;
    for (int32_T i_0{0}; i_0 < 14; i_0++) {
      rtb_tecsthetaCmd += rtb_InterpolationUsingPreloo_dw[6 * i_0 + i] *
        rtb_Saturation_az[i_0];
    }

    rtb_Subtract1[i] = (rtb_Subtract1[i] + rtb_tecsthetaCmd) -
      rtb_Subtract1_tmp_0[i];
  }

  /* End of Sum: '<S17>/Subtract1' */

  /* Outport: '<Root>/CmdBusOut' incorporates:
   *  BusCreator generated from: '<Root>/CmdBusOut'
   *  Inport: '<Root>/CmdBusIn'
   */
  PX4Controller_Y.CmdBusOut.rc[0] = PX4Controller_U.CmdBusIn.rc[0];
  PX4Controller_Y.CmdBusOut.rc[1] = PX4Controller_U.CmdBusIn.rc[1];
  PX4Controller_Y.CmdBusOut.rc[2] = PX4Controller_U.CmdBusIn.rc[2];
  PX4Controller_Y.CmdBusOut.rc[3] = PX4Controller_U.CmdBusIn.rc[3];
  PX4Controller_Y.CmdBusOut.hCmd = PX4Controller_U.CmdBusIn.hCmd;
  PX4Controller_Y.CmdBusOut.hDotCmd = PX4Controller_U.CmdBusIn.hDotCmd;
  PX4Controller_Y.CmdBusOut.tasCmd = PX4Controller_U.CmdBusIn.tasCmd;
  PX4Controller_Y.CmdBusOut.tasDotCmd = PX4Controller_U.CmdBusIn.tasDotCmd;
  PX4Controller_Y.CmdBusOut.manualAttitude =
    PX4Controller_U.CmdBusIn.manualAttitude;
  PX4Controller_Y.CmdBusOut.manualRate = PX4Controller_U.CmdBusIn.manualRate;
  PX4Controller_Y.CmdBusOut.manualFM = PX4Controller_U.CmdBusIn.manualFM;
  PX4Controller_Y.CmdBusOut.FCmd = rtb_Switch1_o;
  PX4Controller_Y.CmdBusOut.manualActuation =
    PX4Controller_U.CmdBusIn.manualActuation;
  PX4Controller_Y.CmdBusOut.trackLine[0] = PX4Controller_U.CmdBusIn.trackLine[0];
  PX4Controller_Y.CmdBusOut.eulerCmd[0] = rtb_Switch7_idx_0;
  PX4Controller_Y.CmdBusOut.MCmd[0] = rtb_Switch7_e[0];
  PX4Controller_Y.CmdBusOut.eulerSat[0] = frac_2[0];
  PX4Controller_Y.CmdBusOut.wSat[0] = rtb_Gain1_g;
  PX4Controller_Y.CmdBusOut.MSat[0] = rtb_Subtract1[3];
  PX4Controller_Y.CmdBusOut.trackLine[1] = PX4Controller_U.CmdBusIn.trackLine[1];
  PX4Controller_Y.CmdBusOut.eulerCmd[1] = rtb_Switch7_idx_1;
  PX4Controller_Y.CmdBusOut.MCmd[1] = rtb_Switch7_e[1];
  PX4Controller_Y.CmdBusOut.eulerSat[1] = frac_2[1];
  PX4Controller_Y.CmdBusOut.wSat[1] = rtb_Tan1;
  PX4Controller_Y.CmdBusOut.MSat[1] = rtb_Subtract1[4];
  PX4Controller_Y.CmdBusOut.trackLine[2] = PX4Controller_U.CmdBusIn.trackLine[2];
  PX4Controller_Y.CmdBusOut.MCmd[2] = rtb_Switch7_e[2];
  PX4Controller_Y.CmdBusOut.eulerSat[2] = frac_2[2];
  PX4Controller_Y.CmdBusOut.wSat[2] = rtb_Gain6;
  PX4Controller_Y.CmdBusOut.MSat[2] = rtb_Subtract1[5];
  PX4Controller_Y.CmdBusOut.FSat = rtb_Subtract1[0];
  if (rtmIsMajorTimeStep(PX4Controller_M)) {
    /* Memory generated from: '<S1>/Memory' */
    PX4Controller_Y.ActBus_e.motors[0] =
      PX4Controller_DW.Memory_1_PreviousInput[0];
    PX4Controller_Y.ActBus_e.motors[1] =
      PX4Controller_DW.Memory_1_PreviousInput[1];
    PX4Controller_Y.ActBus_e.motors[2] =
      PX4Controller_DW.Memory_1_PreviousInput[2];
    PX4Controller_Y.ActBus_e.motors[3] =
      PX4Controller_DW.Memory_1_PreviousInput[3];

    /* Memory generated from: '<S1>/Memory' */
    for (int32_T i{0}; i < 6; i++) {
      PX4Controller_Y.ActBus_e.ailerons[i] =
        PX4Controller_DW.Memory_2_PreviousInput[i];
    }

    /* Memory generated from: '<S1>/Memory' */
    PX4Controller_Y.ActBus_e.elevators[0] =
      PX4Controller_DW.Memory_3_PreviousInput[0];

    /* Memory generated from: '<S1>/Memory' */
    PX4Controller_Y.ActBus_e.rudders[0] =
      PX4Controller_DW.Memory_4_PreviousInput[0];

    /* Memory generated from: '<S1>/Memory' */
    PX4Controller_Y.ActBus_e.elevators[1] =
      PX4Controller_DW.Memory_3_PreviousInput[1];

    /* Memory generated from: '<S1>/Memory' */
    PX4Controller_Y.ActBus_e.rudders[1] =
      PX4Controller_DW.Memory_4_PreviousInput[1];
  }

  /* Switch: '<S11>/Switch1' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  if (PX4Controller_U.CmdBusIn.manualActuation >
      PX4Controller_P.Switch1_Threshold_m) {
    /* Switch: '<S11>/Switch1' incorporates:
     *  Constant: '<S11>/Constant'
     *  Product: '<S11>/Product'
     */
    PX4Controller_B.Switch1[0] = PX4Controller_U.CmdBusIn.rc[0] *
      PX4Controller_P.Constant_Value_j[0];
    PX4Controller_B.Switch1[1] = PX4Controller_U.CmdBusIn.rc[0] *
      PX4Controller_P.Constant_Value_j[1];
    PX4Controller_B.Switch1[2] = PX4Controller_U.CmdBusIn.rc[0] *
      PX4Controller_P.Constant_Value_j[2];
    PX4Controller_B.Switch1[3] = PX4Controller_U.CmdBusIn.rc[0] *
      PX4Controller_P.Constant_Value_j[3];
  } else {
    /* Switch: '<S11>/Switch1' */
    PX4Controller_B.Switch1[0] = rtb_Saturation_az[0];
    PX4Controller_B.Switch1[1] = rtb_Saturation_az[1];
    PX4Controller_B.Switch1[2] = rtb_Saturation_az[2];
    PX4Controller_B.Switch1[3] = rtb_Saturation_az[3];
  }

  /* End of Switch: '<S11>/Switch1' */

  /* Switch: '<S11>/Switch2' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  if (PX4Controller_U.CmdBusIn.manualActuation >
      PX4Controller_P.Switch2_Threshold_o) {
    /* Switch: '<S11>/Switch2' incorporates:
     *  Constant: '<S11>/Constant1'
     *  Product: '<S11>/Product1'
     */
    for (int32_T i{0}; i < 6; i++) {
      PX4Controller_B.Switch2[i] = PX4Controller_U.CmdBusIn.rc[1] *
        PX4Controller_P.Constant1_Value_h[i];
    }
  } else {
    /* Switch: '<S11>/Switch2' */
    for (int32_T i{0}; i < 6; i++) {
      PX4Controller_B.Switch2[i] = rtb_Saturation_az[i + 4];
    }
  }

  /* End of Switch: '<S11>/Switch2' */

  /* Switch: '<S11>/Switch3' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  if (PX4Controller_U.CmdBusIn.manualActuation >
      PX4Controller_P.Switch3_Threshold) {
    /* Switch: '<S11>/Switch3' incorporates:
     *  Constant: '<S11>/Constant2'
     *  Product: '<S11>/Product2'
     */
    PX4Controller_B.Switch3[0] = PX4Controller_P.Constant2_Value_k[0] *
      PX4Controller_U.CmdBusIn.rc[2];
    PX4Controller_B.Switch3[1] = PX4Controller_P.Constant2_Value_k[1] *
      PX4Controller_U.CmdBusIn.rc[2];
  } else {
    /* Switch: '<S11>/Switch3' */
    PX4Controller_B.Switch3[0] = rtb_Saturation_az[10];
    PX4Controller_B.Switch3[1] = rtb_Saturation_az[11];
  }

  /* End of Switch: '<S11>/Switch3' */

  /* Switch: '<S11>/Switch4' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  if (PX4Controller_U.CmdBusIn.manualActuation >
      PX4Controller_P.Switch4_Threshold) {
    /* Switch: '<S11>/Switch4' incorporates:
     *  Constant: '<S11>/Constant3'
     *  Product: '<S11>/Product3'
     */
    PX4Controller_B.Switch4[0] = PX4Controller_P.Constant3_Value[0] *
      PX4Controller_U.CmdBusIn.rc[3];
    PX4Controller_B.Switch4[1] = PX4Controller_P.Constant3_Value[1] *
      PX4Controller_U.CmdBusIn.rc[3];
  } else {
    /* Switch: '<S11>/Switch4' */
    PX4Controller_B.Switch4[0] = rtb_Saturation_az[12];
    PX4Controller_B.Switch4[1] = rtb_Saturation_az[13];
  }

  /* End of Switch: '<S11>/Switch4' */

  /* Interpolation_n-D: '<S34>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S34>/Constant'
   */
  frac_i[0] = rtb_f1;
  frac_i[1] = rtb_f2;
  bpIndex_i[0] = rtb_k1;
  bpIndex_i[1] = rtb_k2;
  if (PX4Controller_P.Constant_Value_c[0] > 1.0) {
    bpIndex_i[2] = 1U;
  } else if (PX4Controller_P.Constant_Value_c[0] >= 0.0) {
    bpIndex_i[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_c[0]);
  } else {
    bpIndex_i[2] = 0U;
  }

  rtb_Gain2 = intrp2d_l_pw(bpIndex_i, frac_i,
    &PX4Controller_P.InterpolationUsingPrelookup__ho[bpIndex_i[2] << 2], 2U);

  /* Product: '<S34>/Product' */
  PX4Controller_B.Product[0] = frac[0] * rtb_Gain2;

  /* Interpolation_n-D: '<S34>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S34>/Constant'
   */
  if (PX4Controller_P.Constant_Value_c[1] > 1.0) {
    bpIndex_i[2] = 1U;
  } else if (PX4Controller_P.Constant_Value_c[1] >= 0.0) {
    bpIndex_i[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_c[1]);
  } else {
    bpIndex_i[2] = 0U;
  }

  rtb_Gain2 = intrp2d_l_pw(bpIndex_i, frac_i,
    &PX4Controller_P.InterpolationUsingPrelookup__ho[bpIndex_i[2] << 2], 2U);

  /* Product: '<S34>/Product' */
  PX4Controller_B.Product[1] = frac[1] * rtb_Gain2;

  /* Interpolation_n-D: '<S79>/Interpolation Using Prelookup' */
  frac_j[0] = rtb_f1;
  frac_j[1] = rtb_f2;
  bpIndex_j[0] = rtb_k1;
  bpIndex_j[1] = rtb_k2;

  /* Product: '<S79>/Product' incorporates:
   *  Interpolation_n-D: '<S79>/Interpolation Using Prelookup'
   */
  PX4Controller_B.Product_g = rtb_sebErr * intrp2d_l_pw(bpIndex_j, frac_j,
    PX4Controller_P.InterpolationUsingPrelookup__mw, 2U);

  /* Sum: '<S77>/Sum1' incorporates:
   *  Sum: '<S77>/steDot'
   *  Sum: '<S77>/steDotCmd '
   */
  PX4Controller_B.steDotErr = (rtb_speDotCmd + rtb_skeDotCmd) - (rtb_speDot +
    rtb_skeDot);
  if (rtmIsMajorTimeStep(PX4Controller_M)) {
    if (rtmIsMajorTimeStep(PX4Controller_M)) {
      /* Update for Memory generated from: '<S1>/Memory' */
      PX4Controller_DW.Memory_1_PreviousInput[0] = PX4Controller_B.Switch1[0];
      PX4Controller_DW.Memory_1_PreviousInput[1] = PX4Controller_B.Switch1[1];
      PX4Controller_DW.Memory_1_PreviousInput[2] = PX4Controller_B.Switch1[2];
      PX4Controller_DW.Memory_1_PreviousInput[3] = PX4Controller_B.Switch1[3];

      /* Update for Memory generated from: '<S1>/Memory' */
      for (int32_T i{0}; i < 6; i++) {
        PX4Controller_DW.Memory_2_PreviousInput[i] = PX4Controller_B.Switch2[i];
      }

      /* Update for Memory generated from: '<S1>/Memory' */
      PX4Controller_DW.Memory_3_PreviousInput[0] = PX4Controller_B.Switch3[0];

      /* Update for Memory generated from: '<S1>/Memory' */
      PX4Controller_DW.Memory_4_PreviousInput[0] = PX4Controller_B.Switch4[0];

      /* Update for Memory generated from: '<S1>/Memory' */
      PX4Controller_DW.Memory_3_PreviousInput[1] = PX4Controller_B.Switch3[1];

      /* Update for Memory generated from: '<S1>/Memory' */
      PX4Controller_DW.Memory_4_PreviousInput[1] = PX4Controller_B.Switch4[1];
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(PX4Controller_M)) {
    rt_ertODEUpdateContinuousStates(&PX4Controller_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++PX4Controller_M->Timing.clockTick0)) {
      ++PX4Controller_M->Timing.clockTickH0;
    }

    PX4Controller_M->Timing.t[0] = rtsiGetSolverStopTime
      (&PX4Controller_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.004s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.004, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       * Timer of this task consists of two 32 bit unsigned integers.
       * The two integers represent the low bits Timing.clockTick1 and the high bits
       * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
       */
      PX4Controller_M->Timing.clockTick1++;
      if (!PX4Controller_M->Timing.clockTick1) {
        PX4Controller_M->Timing.clockTickH1++;
      }
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void PX4Controller_derivatives(void)
{
  XDot_PX4Controller_T *_rtXdot;
  _rtXdot = ((XDot_PX4Controller_T *) PX4Controller_M->derivs);

  /* Derivatives for Integrator: '<S76>/Integrator1' */
  if (!PX4Controller_B.OR) {
    _rtXdot->Integrator1_CSTATE = PX4Controller_B.Product_g;
  } else {
    /* level reset is active */
    _rtXdot->Integrator1_CSTATE = 0.0;
  }

  /* End of Derivatives for Integrator: '<S76>/Integrator1' */

  /* Derivatives for Integrator: '<S3>/Integrator' */
  if (!PX4Controller_B.OR_k) {
    _rtXdot->Integrator_CSTATE[0] = PX4Controller_B.Product[0];
    _rtXdot->Integrator_CSTATE[1] = PX4Controller_B.Product[1];
  } else {
    /* level reset is active */
    _rtXdot->Integrator_CSTATE[0] = 0.0;
    _rtXdot->Integrator_CSTATE[1] = 0.0;
  }

  /* End of Derivatives for Integrator: '<S3>/Integrator' */

  /* Derivatives for TransferFcn: '<S77>/Low pass' */
  _rtXdot->Lowpass_CSTATE = PX4Controller_P.Lowpass_A *
    PX4Controller_X.Lowpass_CSTATE;
  _rtXdot->Lowpass_CSTATE += PX4Controller_B.steDotErr;

  /* Derivatives for Integrator: '<S77>/Integrator' */
  if (!PX4Controller_B.OR) {
    _rtXdot->Integrator_CSTATE_f = PX4Controller_B.steErr;
  } else {
    /* level reset is active */
    _rtXdot->Integrator_CSTATE_f = 0.0;
  }

  /* End of Derivatives for Integrator: '<S77>/Integrator' */
}

/* Model initialize function */
void PX4Controller_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  PX4Controller_P.Saturation_UpperSat = rtInf;
  PX4Controller_P.Saturation_UpperSat_m = rtInf;
  PX4Controller_P.Saturation_UpperSat_a = rtInf;
  PX4Controller_P.Saturation_UpperSat_l = rtInf;
  PX4Controller_P.Saturation_UpperSat_ag = rtInf;

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&PX4Controller_M->solverInfo,
                          &PX4Controller_M->Timing.simTimeStep);
    rtsiSetTPtr(&PX4Controller_M->solverInfo, &rtmGetTPtr(PX4Controller_M));
    rtsiSetStepSizePtr(&PX4Controller_M->solverInfo,
                       &PX4Controller_M->Timing.stepSize0);
    rtsiSetdXPtr(&PX4Controller_M->solverInfo, &PX4Controller_M->derivs);
    rtsiSetContStatesPtr(&PX4Controller_M->solverInfo, (real_T **)
                         &PX4Controller_M->contStates);
    rtsiSetNumContStatesPtr(&PX4Controller_M->solverInfo,
      &PX4Controller_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&PX4Controller_M->solverInfo,
      &PX4Controller_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&PX4Controller_M->solverInfo,
      &PX4Controller_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&PX4Controller_M->solverInfo,
      &PX4Controller_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&PX4Controller_M->solverInfo, (&rtmGetErrorStatus
      (PX4Controller_M)));
    rtsiSetRTModelPtr(&PX4Controller_M->solverInfo, PX4Controller_M);
  }

  rtsiSetSimTimeStep(&PX4Controller_M->solverInfo, MAJOR_TIME_STEP);
  PX4Controller_M->intgData.y = PX4Controller_M->odeY;
  PX4Controller_M->intgData.f[0] = PX4Controller_M->odeF[0];
  PX4Controller_M->intgData.f[1] = PX4Controller_M->odeF[1];
  PX4Controller_M->intgData.f[2] = PX4Controller_M->odeF[2];
  PX4Controller_M->intgData.f[3] = PX4Controller_M->odeF[3];
  PX4Controller_M->contStates = ((X_PX4Controller_T *) &PX4Controller_X);
  rtsiSetSolverData(&PX4Controller_M->solverInfo, static_cast<void *>
                    (&PX4Controller_M->intgData));
  rtsiSetIsMinorTimeStepWithModeChange(&PX4Controller_M->solverInfo, false);
  rtsiSetSolverName(&PX4Controller_M->solverInfo,"ode4");
  rtmSetTPtr(PX4Controller_M, &PX4Controller_M->Timing.tArray[0]);
  PX4Controller_M->Timing.stepSize0 = 0.004;

  /* block I/O */
  (void) std::memset((static_cast<void *>(&PX4Controller_B)), 0,
                     sizeof(B_PX4Controller_T));

  /* states (continuous) */
  {
    (void) std::memset(static_cast<void *>(&PX4Controller_X), 0,
                       sizeof(X_PX4Controller_T));
  }

  /* states (dwork) */
  (void) std::memset(static_cast<void *>(&PX4Controller_DW), 0,
                     sizeof(DW_PX4Controller_T));

  /* external inputs */
  (void)std::memset(&PX4Controller_U, 0, sizeof(ExtU_PX4Controller_T));

  /* external outputs */
  (void)std::memset(&PX4Controller_Y, 0, sizeof(ExtY_PX4Controller_T));
  PX4Controller_PrevZCX.Integrator1_Reset_ZCE = UNINITIALIZED_ZCSIG;
  PX4Controller_PrevZCX.Integrator_Reset_ZCE = UNINITIALIZED_ZCSIG;
  PX4Controller_PrevZCX.Integrator_Reset_ZCE_l = UNINITIALIZED_ZCSIG;

  /* InitializeConditions for Integrator: '<S76>/Integrator1' */
  PX4Controller_X.Integrator1_CSTATE = PX4Controller_P.Integrator1_IC;

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  PX4Controller_X.Integrator_CSTATE[0] = PX4Controller_P.Integrator_IC;
  PX4Controller_X.Integrator_CSTATE[1] = PX4Controller_P.Integrator_IC;

  /* InitializeConditions for TransferFcn: '<S77>/Low pass' */
  PX4Controller_X.Lowpass_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S77>/Integrator' */
  PX4Controller_X.Integrator_CSTATE_f = PX4Controller_P.Integrator_IC_a;

  /* InitializeConditions for Memory generated from: '<S1>/Memory' */
  PX4Controller_DW.Memory_1_PreviousInput[0] =
    PX4Controller_P.Memory_1_InitialCondition;
  PX4Controller_DW.Memory_1_PreviousInput[1] =
    PX4Controller_P.Memory_1_InitialCondition;
  PX4Controller_DW.Memory_1_PreviousInput[2] =
    PX4Controller_P.Memory_1_InitialCondition;
  PX4Controller_DW.Memory_1_PreviousInput[3] =
    PX4Controller_P.Memory_1_InitialCondition;

  /* InitializeConditions for Memory generated from: '<S1>/Memory' */
  for (int32_T i{0}; i < 6; i++) {
    PX4Controller_DW.Memory_2_PreviousInput[i] =
      PX4Controller_P.Memory_2_InitialCondition;
  }

  /* InitializeConditions for Memory generated from: '<S1>/Memory' */
  PX4Controller_DW.Memory_3_PreviousInput[0] =
    PX4Controller_P.Memory_3_InitialCondition;

  /* InitializeConditions for Memory generated from: '<S1>/Memory' */
  PX4Controller_DW.Memory_4_PreviousInput[0] =
    PX4Controller_P.Memory_4_InitialCondition;

  /* InitializeConditions for Memory generated from: '<S1>/Memory' */
  PX4Controller_DW.Memory_3_PreviousInput[1] =
    PX4Controller_P.Memory_3_InitialCondition;

  /* InitializeConditions for Memory generated from: '<S1>/Memory' */
  PX4Controller_DW.Memory_4_PreviousInput[1] =
    PX4Controller_P.Memory_4_InitialCondition;
}

/* Model terminate function */
void PX4Controller_terminate(void)
{
  /* (no terminate code required) */
}
