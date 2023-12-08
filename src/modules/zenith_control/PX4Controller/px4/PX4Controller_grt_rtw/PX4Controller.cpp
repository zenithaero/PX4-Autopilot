/*
 * PX4Controller.cpp
 *
 * Code generation for model "PX4Controller".
 *
 * Model version              : 1.277
 * Simulink Coder version : 9.7 (R2022a) 13-Nov-2021
 * C++ source code generated on : Fri Dec  8 03:03:21 2023
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
  real_T rtb_InterpolationUsingPrelook_g[84];
  real_T rtb_InterpolationUsingPrelook_l[84];
  real_T rtb_Divide2[14];
  real_T rtb_Sum2_d[14];
  real_T rtb_Sum2_f[14];
  real_T rtb_InterpolationUsingPrelook_d[6];
  real_T rtb_MatrixMultiply1_k[6];
  real_T rtb_Sum2_a[6];
  real_T frac_2[3];
  real_T frac_5[3];
  real_T frac_b[3];
  real_T frac_c[3];
  real_T frac_d[3];
  real_T frac_e[3];
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
  real_T frac_f[2];
  real_T frac_g[2];
  real_T frac_h[2];
  real_T frac_i[2];
  real_T frac_j[2];
  real_T frac_k[2];
  real_T frac_l[2];
  real_T rtb_Assignment1_idx_0;
  real_T rtb_Assignment1_idx_1;
  real_T rtb_Assignment1_idx_2;
  real_T rtb_Gain2;
  real_T rtb_InterpolationUsingPrelook_k;
  real_T rtb_Saturation1;
  real_T rtb_Saturation7;
  real_T rtb_Saturation7_g;
  real_T rtb_Switch1_b;
  real_T rtb_Switch7_d_idx_0;
  real_T rtb_Switch7_d_idx_0_0;
  real_T rtb_Switch7_d_idx_1;
  real_T rtb_Switch7_d_idx_2;
  real_T rtb_Switch7_idx_0;
  real_T rtb_Switch7_idx_1;
  real_T rtb_f1;
  real_T rtb_f1_m;
  real_T rtb_f2;
  real_T rtb_f2_p;
  real_T rtb_f3;
  real_T rtb_sebErr;
  real_T rtb_skeCmd;
  real_T rtb_skeDot;
  real_T rtb_skeDotCmd;
  real_T rtb_ske_0;
  real_T rtb_spe;
  real_T rtb_speDot;
  real_T rtb_speDotCmd;
  real_T rtb_tecsthetaCmd;
  uint32_T bpIndex_2[4];
  uint32_T bpIndex_5[4];
  uint32_T bpIndex_b[4];
  uint32_T bpIndex_c[4];
  uint32_T bpIndex_d[4];
  uint32_T bpIndex_e[4];
  uint32_T bpIndex_3[3];
  uint32_T bpIndex_4[3];
  uint32_T bpIndex_6[3];
  uint32_T bpIndex_9[3];
  uint32_T bpIndex_f[3];
  uint32_T bpIndex_g[3];
  uint32_T bpIndex_h[3];
  uint32_T bpIndex_i[3];
  uint32_T bpIndex_j[3];
  uint32_T bpIndex_k[3];
  uint32_T bpIndex[2];
  uint32_T bpIndex_0[2];
  uint32_T bpIndex_1[2];
  uint32_T bpIndex_7[2];
  uint32_T bpIndex_8[2];
  uint32_T bpIndex_a[2];
  uint32_T bpIndex_l[2];
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

  /* DotProduct: '<S67>/Dot Product' incorporates:
   *  Inport: '<Root>/StateBus'
   */
  rtb_ske_0 = PX4Controller_U.StateBus_m.vNed[0] *
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

  /* Interpolation_n-D: '<S70>/Interpolation Using Prelookup' */
  frac[0] = rtb_f1;
  frac[1] = rtb_f2;
  bpIndex[0] = rtb_k1;
  bpIndex[1] = rtb_k2;

  /* Product: '<S70>/Product' incorporates:
   *  DotProduct: '<S67>/Dot Product'
   *  Gain: '<S67>/Gain'
   *  Interpolation_n-D: '<S70>/Interpolation Using Prelookup'
   *  Sqrt: '<S67>/Square Root'
   */
  rtb_skeDot = PX4Controller_P.Gain_Gain_a * std::sqrt(rtb_ske_0) * intrp2d_l_pw
    (bpIndex, frac, PX4Controller_P.InterpolationUsingPrelookup_T_h, 2U);

  /* Saturate: '<S67>/Saturation' */
  if (rtb_skeDot > PX4Controller_P.Saturation_UpperSat_m) {
    rtb_skeDot = PX4Controller_P.Saturation_UpperSat_m;
  } else if (rtb_skeDot < PX4Controller_P.Saturation_LowerSat_k) {
    rtb_skeDot = PX4Controller_P.Saturation_LowerSat_k;
  }

  /* End of Saturate: '<S67>/Saturation' */

  /* Trigonometry: '<S68>/Cos1' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_spe = std::sin(PX4Controller_U.CmdBusIn.trackLine[0]);

  /* Trigonometry: '<S68>/Cos' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_skeCmd = std::cos(PX4Controller_U.CmdBusIn.trackLine[0]);

  /* SignalConversion generated from: '<S68>/Dot Product1' incorporates:
   *  Gain: '<S68>/Gain'
   */
  frac[1] = PX4Controller_P.Gain_Gain_do * rtb_skeCmd;

  /* DotProduct: '<S68>/Dot Product1' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   *  Inport: '<Root>/StateBus'
   *  SignalConversion generated from: '<S68>/Dot Product1'
   *  Sum: '<S68>/Minus'
   */
  rtb_Switch7_d_idx_0 = (PX4Controller_U.StateBus_m.xNed[0] -
    PX4Controller_U.CmdBusIn.trackLine[1]) * rtb_spe +
    (PX4Controller_U.StateBus_m.xNed[1] - PX4Controller_U.CmdBusIn.trackLine[2])
    * frac[1];

  /* Sum: '<S68>/Minus1' incorporates:
   *  DotProduct: '<S68>/Dot Product1'
   *  Math: '<S68>/Square'
   *  Math: '<S68>/Square2'
   */
  rtb_f3 = rtb_skeDot * rtb_skeDot - rtb_Switch7_d_idx_0 * rtb_Switch7_d_idx_0;

  /* Saturate: '<S68>/Saturation' */
  if (rtb_f3 > PX4Controller_P.Saturation_UpperSat_a) {
    rtb_f3 = PX4Controller_P.Saturation_UpperSat_a;
  } else if (rtb_f3 < PX4Controller_P.Saturation_LowerSat_e) {
    rtb_f3 = PX4Controller_P.Saturation_LowerSat_e;
  }

  /* End of Saturate: '<S68>/Saturation' */

  /* Trigonometry: '<S66>/Cos' incorporates:
   *  Inport: '<Root>/StateBus'
   *  Trigonometry: '<S31>/sincos'
   *  Trigonometry: '<S37>/sincos'
   */
  rtb_InterpolationUsingPrelook_k = std::cos
    (PX4Controller_U.StateBus_m.eulerBody[1]);

  /* Gain: '<S6>/Gain2' incorporates:
   *  DotProduct: '<S67>/Dot Product'
   *  DotProduct: '<S68>/Dot Product1'
   *  DotProduct: '<S68>/Dot Product2'
   *  DotProduct: '<S68>/Dot Product3'
   *  Gain: '<S66>/1//G'
   *  Gain: '<S66>/Gain'
   *  Inport: '<Root>/StateBus'
   *  Product: '<S66>/Divide'
   *  Product: '<S66>/Product'
   *  SignalConversion generated from: '<S68>/Dot Product1'
   *  Sqrt: '<S68>/Square Root'
   *  Sum: '<S68>/Minus2'
   *  Trigonometry: '<S66>/Atan'
   *  Trigonometry: '<S66>/Cos'
   *  Trigonometry: '<S66>/Sin'
   *  Trigonometry: '<S68>/Atan1'
   *  Trigonometry: '<S68>/Atan2'
   */
  rtb_Gain2 = std::atan(std::sin(rt_atan2d_snf(rtb_spe *
    PX4Controller_U.StateBus_m.vNed[0] + frac[1] *
    PX4Controller_U.StateBus_m.vNed[1], rtb_skeCmd *
    PX4Controller_U.StateBus_m.vNed[0] + rtb_spe *
    PX4Controller_U.StateBus_m.vNed[1]) + rt_atan2d_snf(rtb_Switch7_d_idx_0, std::
    sqrt(rtb_f3))) * (1.0 / rtb_skeDot) * rtb_ske_0 *
                        PX4Controller_P.Gain_Gain_dn * PX4Controller_P.uG_Gain *
                        rtb_InterpolationUsingPrelook_k) *
    PX4Controller_P.Gain2_Gain_h;

  /* Logic: '<S74>/OR' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  PX4Controller_B.OR = ((PX4Controller_U.CmdBusIn.manualAttitude != 0.0) ||
                        (PX4Controller_U.CmdBusIn.manualRate != 0.0) ||
                        (PX4Controller_U.CmdBusIn.manualFM != 0.0) ||
                        (PX4Controller_U.CmdBusIn.manualActuation != 0.0));

  /* Integrator: '<S77>/Integrator1' */
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

  rtb_f3 = PX4Controller_X.Integrator1_CSTATE;

  /* Gain: '<S76>/Multiply' incorporates:
   *  Gain: '<S76>/Gain'
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_Saturation1 = PX4Controller_P.Gain_Gain_f * PX4Controller_U.CmdBusIn.hCmd *
    PX4Controller_P.Multiply_Gain;

  /* Gain: '<S76>/Gain2' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_skeCmd = PX4Controller_P.Gain2_Gain_o * PX4Controller_U.CmdBusIn.tasCmd;

  /* Gain: '<S76>/Multiply5' incorporates:
   *  Math: '<S76>/Square1'
   */
  rtb_skeCmd = rtb_skeCmd * rtb_skeCmd * PX4Controller_P.Multiply5_Gain;

  /* Gain: '<S76>/Multiply1' incorporates:
   *  Gain: '<S38>/Gain2'
   *  Gain: '<S76>/Gain1'
   *  Inport: '<Root>/StateBus'
   */
  rtb_spe = PX4Controller_P.Gain2_Gain_c * PX4Controller_U.StateBus_m.xNed[2] *
    PX4Controller_P.Gain1_Gain_h * PX4Controller_P.Multiply1_Gain;

  /* Gain: '<S76>/Gain3' incorporates:
   *  Inport: '<Root>/StateBus'
   */
  rtb_skeDot = PX4Controller_P.Gain3_Gain_j * PX4Controller_U.StateBus_m.tas;

  /* Gain: '<S76>/Multiply3' incorporates:
   *  Math: '<S76>/Square'
   */
  rtb_ske_0 = rtb_skeDot * rtb_skeDot * PX4Controller_P.Multiply3_Gain;

  /* Sum: '<S77>/Sum' incorporates:
   *  Constant: '<S77>/Constant'
   *  Constant: '<S77>/Constant1'
   *  Product: '<S77>/Product'
   *  Product: '<S77>/Product1'
   *  Product: '<S77>/Product2'
   *  Product: '<S77>/Product3'
   *  Sum: '<S77>/seb'
   *  Sum: '<S77>/sebCmd'
   */
  rtb_sebErr = (PX4Controller_P.Constant_Value * rtb_Saturation1 -
                PX4Controller_P.Constant1_Value * rtb_skeCmd) -
    (PX4Controller_P.Constant_Value * rtb_spe - PX4Controller_P.Constant1_Value *
     rtb_ske_0);

  /* Interpolation_n-D: '<S81>/Interpolation Using Prelookup' */
  frac_0[0] = rtb_f1;
  frac_0[1] = rtb_f2;
  bpIndex_0[0] = rtb_k1;
  bpIndex_0[1] = rtb_k2;

  /* Gain: '<S76>/Multiply4' incorporates:
   *  Gain: '<S76>/Gain6'
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_speDotCmd = PX4Controller_P.Gain6_Gain * PX4Controller_U.CmdBusIn.hDotCmd *
    PX4Controller_P.Multiply4_Gain;

  /* Product: '<S76>/Product1' incorporates:
   *  Gain: '<S76>/Gain5'
   *  Inport: '<Root>/CmdBusIn'
   */
  rtb_skeDotCmd = PX4Controller_P.Gain5_Gain *
    PX4Controller_U.CmdBusIn.tasDotCmd * rtb_skeDot;

  /* Gain: '<S76>/Multiply2' incorporates:
   *  Gain: '<S38>/Gain3'
   *  Gain: '<S76>/Gain4'
   *  Inport: '<Root>/StateBus'
   */
  rtb_speDot = PX4Controller_P.Gain3_Gain_p * PX4Controller_U.StateBus_m.vNed[2]
    * PX4Controller_P.Gain4_Gain * PX4Controller_P.Multiply2_Gain;
  if (rtmIsMajorTimeStep(PX4Controller_M)) {
    /* Gain: '<S76>/Gain7' incorporates:
     *  Constant: '<S38>/Constant'
     */
    PX4Controller_B.Gain7 = PX4Controller_P.Gain7_Gain *
      PX4Controller_P.Constant_Value_b;
  }

  /* Product: '<S76>/Product' */
  rtb_skeDot *= PX4Controller_B.Gain7;

  /* Interpolation_n-D: '<S79>/Interpolation Using Prelookup' */
  frac_1[0] = rtb_f1;
  frac_1[1] = rtb_f2;
  bpIndex_1[0] = rtb_k1;
  bpIndex_1[1] = rtb_k2;

  /* Gain: '<S77>/Gain1' incorporates:
   *  Inport: '<Root>/StateBus'
   */
  rtb_f2_p = PX4Controller_P.Gain1_Gain_hh * PX4Controller_U.StateBus_m.tas;

  /* Saturate: '<S77>/Saturation' */
  if (rtb_f2_p > PX4Controller_P.Saturation_UpperSat_l) {
    rtb_f2_p = PX4Controller_P.Saturation_UpperSat_l;
  } else if (rtb_f2_p < PX4Controller_P.Saturation_LowerSat_b) {
    rtb_f2_p = PX4Controller_P.Saturation_LowerSat_b;
  }

  /* End of Saturate: '<S77>/Saturation' */

  /* Sum: '<S77>/Add1' incorporates:
   *  Constant: '<S77>/Constant'
   *  Constant: '<S77>/Constant1'
   *  Constant: '<S77>/Constant2'
   *  Integrator: '<S77>/Integrator1'
   *  Interpolation_n-D: '<S79>/Interpolation Using Prelookup'
   *  Interpolation_n-D: '<S81>/Interpolation Using Prelookup'
   *  Product: '<S77>/Divide'
   *  Product: '<S77>/Product4'
   *  Product: '<S77>/Product5'
   *  Product: '<S77>/Product6'
   *  Product: '<S77>/Product7'
   *  Product: '<S79>/Product'
   *  Product: '<S81>/Product'
   *  Sum: '<S77>/Add'
   *  Sum: '<S77>/Sum1'
   *  Sum: '<S77>/sebDot'
   *  Sum: '<S77>/sebDotCmd '
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

  /* Switch: '<S28>/Switch7' incorporates:
   *  Gain: '<S28>/Gain1'
   *  Gain: '<S28>/Gain2'
   *  Inport: '<Root>/CmdBusIn'
   */
  if (PX4Controller_U.CmdBusIn.manualAttitude >
      PX4Controller_P.Switch7_Threshold) {
    rtb_Switch7_idx_0 = PX4Controller_P.Gain1_Gain_p *
      PX4Controller_U.CmdBusIn.rc[1];
    rtb_Switch7_idx_1 = PX4Controller_P.Gain2_Gain_p *
      PX4Controller_U.CmdBusIn.rc[2];

    /* Outport: '<Root>/CmdBusOut' incorporates:
     *  Gain: '<S28>/Gain1'
     *  Gain: '<S28>/Gain2'
     *  Gain: '<S28>/Gain3'
     */
    PX4Controller_Y.CmdBusOut.eulerCmd[2] = PX4Controller_P.Gain3_Gain_b *
      PX4Controller_U.CmdBusIn.rc[3];
  } else {
    /* Outport: '<Root>/CmdBusOut' incorporates:
     *  Assignment: '<S75>/Assignment'
     *  SignalConversion generated from: '<S69>/Assignment'
     */
    PX4Controller_Y.CmdBusOut.eulerCmd[2] = PX4Controller_U.CmdBusIn.eulerCmd[2];

    /* Assignment: '<S69>/Assignment' */
    rtb_Switch7_idx_0 = rtb_Gain2;

    /* Assignment: '<S75>/Assignment' */
    rtb_Switch7_idx_1 = rtb_tecsthetaCmd;
  }

  /* End of Switch: '<S28>/Switch7' */

  /* Logic: '<S27>/OR' incorporates:
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
                       PX4Controller_P.Prelookup3_BreakpointsData, 1U, &rtb_f3,
                       &PX4Controller_DW.Prelookup3_DWORK1);

  /* Interpolation_n-D: '<S32>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S32>/Constant'
   */
  frac_2[0] = rtb_f1_m;
  frac_2[1] = rtb_f2_p;
  frac_2[2] = rtb_f3;
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

  rtb_Switch1_b = intrp3d_l_pw(bpIndex_2, frac_2,
    &PX4Controller_P.SubsystemReference3_table[12U * bpIndex_2[3]],
    PX4Controller_P.InterpolationUsingPrelookup_dim);

  /* Gain: '<S3>/   ' */
  rtb_Saturation7_g = PX4Controller_P._Gain * rtb_Switch1_b;

  /* Interpolation_n-D: '<S32>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S32>/Constant'
   */
  if (PX4Controller_P.Constant_Value_h[1] > 2.0) {
    bpIndex_2[3] = 2U;
  } else if (PX4Controller_P.Constant_Value_h[1] >= 0.0) {
    bpIndex_2[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_h[1]);
  } else {
    bpIndex_2[3] = 0U;
  }

  rtb_Switch1_b = intrp3d_l_pw(bpIndex_2, frac_2,
    &PX4Controller_P.SubsystemReference3_table[12U * bpIndex_2[3]],
    PX4Controller_P.InterpolationUsingPrelookup_dim);

  /* Interpolation_n-D: '<S36>/Interpolation Using Prelookup' */
  frac_3[0] = rtb_f1;
  frac_3[1] = rtb_f2;
  bpIndex_3[0] = rtb_k1;
  bpIndex_3[1] = rtb_k2;

  /* Saturate: '<S3>/Saturation1' */
  if (rtb_Switch7_idx_0 > PX4Controller_P.Saturation1_UpperSat[0]) {
    rtb_Saturation7 = PX4Controller_P.Saturation1_UpperSat[0];
  } else if (rtb_Switch7_idx_0 < PX4Controller_P.Saturation1_LowerSat[0]) {
    rtb_Saturation7 = PX4Controller_P.Saturation1_LowerSat[0];
  } else {
    rtb_Saturation7 = rtb_Switch7_idx_0;
  }

  /* Sum: '<S3>/Subtract' incorporates:
   *  Gain: '<S3>/Gain1'
   *  Gain: '<S3>/Gain2'
   *  Inport: '<Root>/StateBus'
   */
  frac[0] = (PX4Controller_P.Gain1_Gain_g * rtb_Saturation7 + rtb_Saturation7_g)
    - PX4Controller_P.Gain2_Gain_c2 * PX4Controller_U.StateBus_m.eulerBody[0];

  /* Interpolation_n-D: '<S36>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S36>/Constant'
   */
  if (PX4Controller_P.Constant_Value_bd[0] > 1.0) {
    bpIndex_3[2] = 1U;
  } else if (PX4Controller_P.Constant_Value_bd[0] >= 0.0) {
    bpIndex_3[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_bd[0]);
  } else {
    bpIndex_3[2] = 0U;
  }

  /* Gain: '<S3>/Gain2' incorporates:
   *  Interpolation_n-D: '<S36>/Interpolation Using Prelookup'
   */
  frac_0[0] = intrp2d_l_pw(bpIndex_3, frac_3,
    &PX4Controller_P.InterpolationUsingPrelookup_T_i[bpIndex_3[2] << 2], 2U);

  /* Saturate: '<S3>/Saturation1' */
  if (rtb_Switch7_idx_1 > PX4Controller_P.Saturation1_UpperSat[1]) {
    rtb_Saturation7_g = PX4Controller_P.Saturation1_UpperSat[1];
  } else if (rtb_Switch7_idx_1 < PX4Controller_P.Saturation1_LowerSat[1]) {
    rtb_Saturation7_g = PX4Controller_P.Saturation1_LowerSat[1];
  } else {
    rtb_Saturation7_g = rtb_Switch7_idx_1;
  }

  /* Sum: '<S3>/Subtract' incorporates:
   *  Gain: '<S3>/   '
   *  Gain: '<S3>/Gain1'
   *  Gain: '<S3>/Gain2'
   *  Inport: '<Root>/StateBus'
   */
  frac[1] = (PX4Controller_P._Gain * rtb_Switch1_b +
             PX4Controller_P.Gain1_Gain_g * rtb_Saturation7_g) -
    PX4Controller_P.Gain2_Gain_c2 * PX4Controller_U.StateBus_m.eulerBody[1];

  /* Interpolation_n-D: '<S36>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S36>/Constant'
   */
  if (PX4Controller_P.Constant_Value_bd[1] > 1.0) {
    bpIndex_3[2] = 1U;
  } else if (PX4Controller_P.Constant_Value_bd[1] >= 0.0) {
    bpIndex_3[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_bd[1]);
  } else {
    bpIndex_3[2] = 0U;
  }

  rtb_Switch7_d_idx_0 = intrp2d_l_pw(bpIndex_3, frac_3,
    &PX4Controller_P.InterpolationUsingPrelookup_T_i[bpIndex_3[2] << 2], 2U);

  /* Trigonometry: '<S37>/sincos' incorporates:
   *  Inport: '<Root>/StateBus'
   *  Trigonometry: '<S31>/sincos'
   */
  rtb_Switch7_d_idx_1 = std::sin(PX4Controller_U.StateBus_m.eulerBody[0]);
  rtb_Switch7_d_idx_2 = std::cos(PX4Controller_U.StateBus_m.eulerBody[0]);

  /* Interpolation_n-D: '<S34>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S34>/Constant'
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

  rtb_Switch1_b = intrp2d_l_pw(bpIndex_4, frac_4,
    &PX4Controller_P.InterpolationUsingPrelookup_T_l[bpIndex_4[2] << 2], 2U);

  /* Trigonometry: '<S37>/sincos' incorporates:
   *  Inport: '<Root>/StateBus'
   *  Trigonometry: '<S31>/sincos'
   */
  rtb_Saturation7_g = std::sin(PX4Controller_U.StateBus_m.eulerBody[1]);

  /* Sum: '<S3>/Add' incorporates:
   *  Fcn: '<S37>/phidot'
   *  Gain: '<S3>/Gain4'
   *  Gain: '<S3>/Gain5'
   *  Inport: '<Root>/StateBus'
   *  Integrator: '<S3>/Integrator'
   *  Product: '<S34>/Product'
   *  Product: '<S36>/Product'
   *  Trigonometry: '<S37>/sincos'
   */
  rtb_Saturation7 = ((rtb_Switch7_d_idx_1 * PX4Controller_U.StateBus_m.wBody[1]
                      + rtb_Switch7_d_idx_2 * PX4Controller_U.StateBus_m.wBody[2])
                     * (rtb_Saturation7_g / rtb_InterpolationUsingPrelook_k) +
                     PX4Controller_U.StateBus_m.wBody[0]) *
    PX4Controller_P.Gain4_Gain_g * PX4Controller_P.Gain5_Gain_i * rtb_Switch1_b
    + (frac[0] * frac_0[0] + PX4Controller_X.Integrator_CSTATE[0]);

  /* Interpolation_n-D: '<S34>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S34>/Constant'
   */
  if (PX4Controller_P.Constant_Value_k[1] > 1.0) {
    bpIndex_4[2] = 1U;
  } else if (PX4Controller_P.Constant_Value_k[1] >= 0.0) {
    bpIndex_4[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_k[1]);
  } else {
    bpIndex_4[2] = 0U;
  }

  rtb_Switch1_b = intrp2d_l_pw(bpIndex_4, frac_4,
    &PX4Controller_P.InterpolationUsingPrelookup_T_l[bpIndex_4[2] << 2], 2U);

  /* Sum: '<S3>/Add' incorporates:
   *  Fcn: '<S37>/thetadot'
   *  Gain: '<S3>/Gain4'
   *  Gain: '<S3>/Gain5'
   *  Inport: '<Root>/StateBus'
   *  Integrator: '<S3>/Integrator'
   *  Product: '<S34>/Product'
   *  Product: '<S36>/Product'
   *  Trigonometry: '<S37>/sincos'
   */
  rtb_Switch7_d_idx_0 = (rtb_Switch7_d_idx_2 * PX4Controller_U.StateBus_m.wBody
    [1] - rtb_Switch7_d_idx_1 * PX4Controller_U.StateBus_m.wBody[2]) *
    PX4Controller_P.Gain4_Gain_g * PX4Controller_P.Gain5_Gain_i * rtb_Switch1_b
    + (frac[1] * rtb_Switch7_d_idx_0 + PX4Controller_X.Integrator_CSTATE[1]);

  /* Saturate: '<S33>/Saturation' incorporates:
   *  Inport: '<Root>/StateBus'
   */
  if (PX4Controller_U.StateBus_m.tas > PX4Controller_P.Saturation_UpperSat_ag) {
    rtb_Assignment1_idx_2 = PX4Controller_P.Saturation_UpperSat_ag;
  } else if (PX4Controller_U.StateBus_m.tas <
             PX4Controller_P.Saturation_LowerSat_kd) {
    rtb_Assignment1_idx_2 = PX4Controller_P.Saturation_LowerSat_kd;
  } else {
    rtb_Assignment1_idx_2 = PX4Controller_U.StateBus_m.tas;
  }

  /* End of Saturate: '<S33>/Saturation' */

  /* Gain: '<S3>/Gain6' incorporates:
   *  Gain: '<S33>/Gain1'
   *  Gain: '<S3>/Gain3'
   *  Product: '<S33>/Divide1'
   *  Sum: '<S33>/Add'
   *  Trigonometry: '<S33>/Tan'
   *  Trigonometry: '<S33>/Tan1'
   */
  rtb_Switch1_b = (std::tan(rtb_Switch7_idx_0) + std::cos(rtb_Switch7_idx_1)) *
    PX4Controller_P.Gain1_Gain_e / rtb_Assignment1_idx_2 *
    PX4Controller_P.Gain3_Gain_g * PX4Controller_P.Gain6_Gain_m;

  /* Fcn: '<S31>/phidot' */
  rtb_Saturation7 -= rtb_Switch1_b * rtb_Saturation7_g;

  /* Fcn: '<S31>/thetadot' */
  rtb_Saturation7_g = rtb_Switch1_b * rtb_Switch7_d_idx_1 *
    rtb_InterpolationUsingPrelook_k + rtb_Switch7_d_idx_2 * rtb_Switch7_d_idx_0;

  /* Fcn: '<S31>/psidot' */
  rtb_InterpolationUsingPrelook_k = rtb_Switch1_b * rtb_Switch7_d_idx_2 *
    rtb_InterpolationUsingPrelook_k - rtb_Switch7_d_idx_1 * rtb_Switch7_d_idx_0;

  /* Switch: '<S71>/Switch7' incorporates:
   *  Gain: '<S71>/Gain1'
   *  Gain: '<S71>/Gain2'
   *  Gain: '<S71>/Gain3'
   *  Inport: '<Root>/CmdBusIn'
   */
  if (PX4Controller_U.CmdBusIn.manualRate > PX4Controller_P.Switch7_Threshold_k)
  {
    rtb_Switch7_d_idx_0 = PX4Controller_P.Gain1_Gain_l *
      PX4Controller_U.CmdBusIn.rc[1];
    rtb_Switch7_d_idx_1 = PX4Controller_P.Gain2_Gain_m *
      PX4Controller_U.CmdBusIn.rc[2];
    rtb_Switch7_d_idx_2 = PX4Controller_P.Gain3_Gain_k *
      PX4Controller_U.CmdBusIn.rc[3];
  } else {
    rtb_Switch7_d_idx_0 = rtb_Saturation7;
    rtb_Switch7_d_idx_1 = rtb_Saturation7_g;
    rtb_Switch7_d_idx_2 = rtb_InterpolationUsingPrelook_k;
  }

  /* End of Switch: '<S71>/Switch7' */

  /* Interpolation_n-D: '<S72>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S72>/Constant'
   */
  frac_5[0] = rtb_f1_m;
  frac_5[1] = rtb_f2_p;
  frac_5[2] = rtb_f3;
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

  rtb_Assignment1_idx_0 = intrp3d_l_pw(bpIndex_5, frac_5,
    &PX4Controller_P.SubsystemReference3_table_m[12U * bpIndex_5[3]],
    PX4Controller_P.InterpolationUsingPrelookup_d_f);
  if (PX4Controller_P.Constant_Value_a[1] > 2.0) {
    bpIndex_5[3] = 2U;
  } else if (PX4Controller_P.Constant_Value_a[1] >= 0.0) {
    bpIndex_5[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_a[1]);
  } else {
    bpIndex_5[3] = 0U;
  }

  rtb_Assignment1_idx_1 = intrp3d_l_pw(bpIndex_5, frac_5,
    &PX4Controller_P.SubsystemReference3_table_m[12U * bpIndex_5[3]],
    PX4Controller_P.InterpolationUsingPrelookup_d_f);
  if (PX4Controller_P.Constant_Value_a[2] > 2.0) {
    bpIndex_5[3] = 2U;
  } else if (PX4Controller_P.Constant_Value_a[2] >= 0.0) {
    bpIndex_5[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_a[2]);
  } else {
    bpIndex_5[3] = 0U;
  }

  rtb_Assignment1_idx_2 = intrp3d_l_pw(bpIndex_5, frac_5,
    &PX4Controller_P.SubsystemReference3_table_m[12U * bpIndex_5[3]],
    PX4Controller_P.InterpolationUsingPrelookup_d_f);

  /* End of Interpolation_n-D: '<S72>/Interpolation Using Prelookup' */

  /* Interpolation_n-D: '<S73>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S73>/Constant'
   */
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

  rtb_Switch1_b = intrp2d_l_pw(bpIndex_6, frac_6,
    &PX4Controller_P.InterpolationUsingPrelookup__cs[bpIndex_6[2] << 2], 2U);

  /* Saturate: '<S7>/Saturation' incorporates:
   *  Sum: '<S7>/Subtract'
   */
  if (rtb_Switch7_d_idx_0 > PX4Controller_P.Saturation_UpperSat_c[0]) {
    rtb_Switch7_d_idx_0_0 = PX4Controller_P.Saturation_UpperSat_c[0];
  } else if (rtb_Switch7_d_idx_0 < PX4Controller_P.Saturation_LowerSat_d[0]) {
    rtb_Switch7_d_idx_0_0 = PX4Controller_P.Saturation_LowerSat_d[0];
  } else {
    rtb_Switch7_d_idx_0_0 = rtb_Switch7_d_idx_0;
  }

  /* Sum: '<S7>/Add' incorporates:
   *  Gain: '<S7>/ '
   *  Gain: '<S7>/   '
   *  Gain: '<S7>/Gain'
   *  Gain: '<S7>/Gain1'
   *  Gain: '<S7>/Gain2'
   *  Gain: '<S7>/Gain3'
   *  Inport: '<Root>/StateBus'
   *  Product: '<S73>/Product'
   *  Sum: '<S7>/Subtract'
   */
  rtb_Assignment1_idx_0 = ((PX4Controller_P._Gain_m * rtb_Assignment1_idx_0 +
    PX4Controller_P._Gain_o * rtb_Switch7_d_idx_0_0) -
    PX4Controller_P.Gain_Gain_at * PX4Controller_U.StateBus_m.wBody[0]) *
    rtb_Switch1_b + PX4Controller_P.Gain1_Gain_f *
    PX4Controller_U.StateBus_m.wDotBody[0] * PX4Controller_P.Gain2_Gain_k *
    PX4Controller_P.Gain3_Gain_c;

  /* Interpolation_n-D: '<S73>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S73>/Constant'
   */
  if (PX4Controller_P.Constant_Value_f[1] > 2.0) {
    bpIndex_6[2] = 2U;
  } else if (PX4Controller_P.Constant_Value_f[1] >= 0.0) {
    bpIndex_6[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_f[1]);
  } else {
    bpIndex_6[2] = 0U;
  }

  rtb_Switch1_b = intrp2d_l_pw(bpIndex_6, frac_6,
    &PX4Controller_P.InterpolationUsingPrelookup__cs[bpIndex_6[2] << 2], 2U);

  /* Saturate: '<S7>/Saturation' incorporates:
   *  Sum: '<S7>/Subtract'
   */
  if (rtb_Switch7_d_idx_1 > PX4Controller_P.Saturation_UpperSat_c[1]) {
    rtb_Switch7_d_idx_0_0 = PX4Controller_P.Saturation_UpperSat_c[1];
  } else if (rtb_Switch7_d_idx_1 < PX4Controller_P.Saturation_LowerSat_d[1]) {
    rtb_Switch7_d_idx_0_0 = PX4Controller_P.Saturation_LowerSat_d[1];
  } else {
    rtb_Switch7_d_idx_0_0 = rtb_Switch7_d_idx_1;
  }

  /* Sum: '<S7>/Add' incorporates:
   *  Gain: '<S7>/ '
   *  Gain: '<S7>/   '
   *  Gain: '<S7>/Gain'
   *  Gain: '<S7>/Gain1'
   *  Gain: '<S7>/Gain2'
   *  Gain: '<S7>/Gain3'
   *  Inport: '<Root>/StateBus'
   *  Product: '<S73>/Product'
   *  Sum: '<S7>/Subtract'
   */
  rtb_Assignment1_idx_1 = ((PX4Controller_P._Gain_m * rtb_Assignment1_idx_1 +
    PX4Controller_P._Gain_o * rtb_Switch7_d_idx_0_0) -
    PX4Controller_P.Gain_Gain_at * PX4Controller_U.StateBus_m.wBody[1]) *
    rtb_Switch1_b + PX4Controller_P.Gain1_Gain_f *
    PX4Controller_U.StateBus_m.wDotBody[1] * PX4Controller_P.Gain2_Gain_k *
    PX4Controller_P.Gain3_Gain_c;

  /* Interpolation_n-D: '<S73>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S73>/Constant'
   */
  if (PX4Controller_P.Constant_Value_f[2] > 2.0) {
    bpIndex_6[2] = 2U;
  } else if (PX4Controller_P.Constant_Value_f[2] >= 0.0) {
    bpIndex_6[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_f[2]);
  } else {
    bpIndex_6[2] = 0U;
  }

  rtb_Switch1_b = intrp2d_l_pw(bpIndex_6, frac_6,
    &PX4Controller_P.InterpolationUsingPrelookup__cs[bpIndex_6[2] << 2], 2U);

  /* Saturate: '<S7>/Saturation' incorporates:
   *  Sum: '<S7>/Subtract'
   */
  if (rtb_Switch7_d_idx_2 > PX4Controller_P.Saturation_UpperSat_c[2]) {
    rtb_Switch7_d_idx_0_0 = PX4Controller_P.Saturation_UpperSat_c[2];
  } else if (rtb_Switch7_d_idx_2 < PX4Controller_P.Saturation_LowerSat_d[2]) {
    rtb_Switch7_d_idx_0_0 = PX4Controller_P.Saturation_LowerSat_d[2];
  } else {
    rtb_Switch7_d_idx_0_0 = rtb_Switch7_d_idx_2;
  }

  /* Sum: '<S7>/Add' incorporates:
   *  Gain: '<S7>/ '
   *  Gain: '<S7>/   '
   *  Gain: '<S7>/Gain'
   *  Gain: '<S7>/Gain1'
   *  Gain: '<S7>/Gain2'
   *  Gain: '<S7>/Gain3'
   *  Inport: '<Root>/StateBus'
   *  Product: '<S73>/Product'
   *  Sum: '<S7>/Subtract'
   */
  rtb_Switch1_b = ((PX4Controller_P._Gain_m * rtb_Assignment1_idx_2 +
                    PX4Controller_P._Gain_o * rtb_Switch7_d_idx_0_0) -
                   PX4Controller_P.Gain_Gain_at *
                   PX4Controller_U.StateBus_m.wBody[2]) * rtb_Switch1_b +
    PX4Controller_P.Gain1_Gain_f * PX4Controller_U.StateBus_m.wDotBody[2] *
    PX4Controller_P.Gain2_Gain_k * PX4Controller_P.Gain3_Gain_c;

  /* Product: '<S7>/Matrix Multiply' incorporates:
   *  Constant: '<S7>/Constant'
   */
  for (int32_T i{0}; i < 3; i++) {
    frac_2[i] = (PX4Controller_P.Constant_Value_i[i + 3] * rtb_Assignment1_idx_1
                 + PX4Controller_P.Constant_Value_i[i] * rtb_Assignment1_idx_0)
      + PX4Controller_P.Constant_Value_i[i + 6] * rtb_Switch1_b;
  }

  /* End of Product: '<S7>/Matrix Multiply' */

  /* Sum: '<S78>/Sum' incorporates:
   *  Sum: '<S78>/ste'
   *  Sum: '<S78>/steCmd'
   */
  PX4Controller_B.steErr = (rtb_Saturation1 + rtb_skeCmd) - (rtb_spe + rtb_ske_0);

  /* Integrator: '<S78>/Integrator' */
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

  /* Switch: '<S71>/Switch1' incorporates:
   *  Gain: '<S71>/Gain'
   *  Inport: '<Root>/CmdBusIn'
   *  Switch: '<S28>/Switch1'
   */
  if (PX4Controller_U.CmdBusIn.manualRate > PX4Controller_P.Switch1_Threshold_g)
  {
    rtb_Switch1_b = PX4Controller_P.Gain_Gain_d * PX4Controller_U.CmdBusIn.rc[0];
  } else if (PX4Controller_U.CmdBusIn.manualAttitude >
             PX4Controller_P.Switch1_Threshold) {
    /* Switch: '<S28>/Switch1' incorporates:
     *  Gain: '<S28>/Gain'
     */
    rtb_Switch1_b = PX4Controller_P.Gain_Gain_g * PX4Controller_U.CmdBusIn.rc[0];
  } else {
    /* Interpolation_n-D: '<S82>/Interpolation Using Prelookup' incorporates:
     *  Switch: '<S28>/Switch1'
     */
    frac_7[0] = rtb_f1;
    frac_7[1] = rtb_f2;
    bpIndex_7[0] = rtb_k1;
    bpIndex_7[1] = rtb_k2;

    /* Interpolation_n-D: '<S84>/Interpolation Using Prelookup' incorporates:
     *  Switch: '<S28>/Switch1'
     */
    frac_8[0] = rtb_f1;
    frac_8[1] = rtb_f2;
    bpIndex_8[0] = rtb_k1;
    bpIndex_8[1] = rtb_k2;

    /* Interpolation_n-D: '<S83>/Interpolation Using Prelookup' incorporates:
     *  Switch: '<S28>/Switch1'
     */
    frac_a[0] = rtb_f1;
    frac_a[1] = rtb_f2;
    bpIndex_a[0] = rtb_k1;
    bpIndex_a[1] = rtb_k2;

    /* Saturate: '<S78>/Saturation' incorporates:
     *  Inport: '<Root>/StateBus'
     *  Switch: '<S28>/Switch1'
     */
    if (PX4Controller_U.StateBus_m.tas > PX4Controller_P.Saturation_UpperSat) {
      rtb_Assignment1_idx_2 = PX4Controller_P.Saturation_UpperSat;
    } else if (PX4Controller_U.StateBus_m.tas <
               PX4Controller_P.Saturation_LowerSat) {
      rtb_Assignment1_idx_2 = PX4Controller_P.Saturation_LowerSat;
    } else {
      rtb_Assignment1_idx_2 = PX4Controller_U.StateBus_m.tas;
    }

    /* End of Saturate: '<S78>/Saturation' */

    /* Switch: '<S28>/Switch1' incorporates:
     *  Gain: '<S78>/M'
     *  Integrator: '<S78>/Integrator'
     *  Interpolation_n-D: '<S82>/Interpolation Using Prelookup'
     *  Interpolation_n-D: '<S83>/Interpolation Using Prelookup'
     *  Interpolation_n-D: '<S84>/Interpolation Using Prelookup'
     *  Product: '<S78>/Divide'
     *  Product: '<S82>/Product'
     *  Product: '<S83>/Product'
     *  Product: '<S84>/Product'
     *  Sum: '<S78>/Sum2'
     *  TransferFcn: '<S78>/Low pass'
     */
    rtb_Switch1_b = ((PX4Controller_P.Lowpass_C * PX4Controller_X.Lowpass_CSTATE
                      * intrp2d_l_pw(bpIndex_a, frac_a,
      PX4Controller_P.InterpolationUsingPrelookup_T_c, 2U) +
                      PX4Controller_B.steErr * intrp2d_l_pw(bpIndex_8, frac_8,
      PX4Controller_P.InterpolationUsingPrelookup_T_n, 2U)) +
                     PX4Controller_X.Integrator_CSTATE_f * intrp2d_l_pw
                     (bpIndex_7, frac_7,
                      PX4Controller_P.InterpolationUsingPrelookup_Tab, 2U)) *
      PX4Controller_P.M_Gain / rtb_Assignment1_idx_2;
  }

  /* End of Switch: '<S71>/Switch1' */

  /* Outport: '<Root>/CmdBusOut' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   *  SignalConversion generated from: '<S69>/Assignment1'
   */
  PX4Controller_Y.CmdBusOut.eulerSat[2] = PX4Controller_U.CmdBusIn.eulerSat[2];

  /* Saturate: '<S69>/Saturation1' */
  if (rtb_Gain2 > PX4Controller_P.Saturation1_UpperSat_g) {
    /* Assignment: '<S69>/Assignment1' incorporates:
     *  Outport: '<Root>/CmdBusOut'
     */
    PX4Controller_Y.CmdBusOut.eulerSat[0] =
      PX4Controller_P.Saturation1_UpperSat_g;
  } else if (rtb_Gain2 < PX4Controller_P.Saturation1_LowerSat_b) {
    /* Assignment: '<S69>/Assignment1' incorporates:
     *  Outport: '<Root>/CmdBusOut'
     */
    PX4Controller_Y.CmdBusOut.eulerSat[0] =
      PX4Controller_P.Saturation1_LowerSat_b;
  } else {
    /* Assignment: '<S69>/Assignment1' incorporates:
     *  Outport: '<Root>/CmdBusOut'
     */
    PX4Controller_Y.CmdBusOut.eulerSat[0] = rtb_Gain2;
  }

  /* End of Saturate: '<S69>/Saturation1' */

  /* Saturate: '<S75>/Saturation1' */
  if (rtb_tecsthetaCmd > PX4Controller_P.Saturation1_UpperSat_k) {
    /* Assignment: '<S75>/Assignment1' incorporates:
     *  Outport: '<Root>/CmdBusOut'
     */
    PX4Controller_Y.CmdBusOut.eulerSat[1] =
      PX4Controller_P.Saturation1_UpperSat_k;
  } else if (rtb_tecsthetaCmd < PX4Controller_P.Saturation1_LowerSat_n) {
    /* Assignment: '<S75>/Assignment1' incorporates:
     *  Outport: '<Root>/CmdBusOut'
     */
    PX4Controller_Y.CmdBusOut.eulerSat[1] =
      PX4Controller_P.Saturation1_LowerSat_n;
  } else {
    /* Assignment: '<S75>/Assignment1' incorporates:
     *  Outport: '<Root>/CmdBusOut'
     */
    PX4Controller_Y.CmdBusOut.eulerSat[1] = rtb_tecsthetaCmd;
  }

  /* End of Saturate: '<S75>/Saturation1' */

  /* Saturate: '<S29>/Saturation' */
  if (rtb_Saturation7 > PX4Controller_P.Saturation_UpperSat_p[0]) {
    /* Outport: '<Root>/CmdBusOut' */
    PX4Controller_Y.CmdBusOut.wSat[0] = PX4Controller_P.Saturation_UpperSat_p[0];
  } else if (rtb_Saturation7 < PX4Controller_P.Saturation_LowerSat_c[0]) {
    /* Outport: '<Root>/CmdBusOut' */
    PX4Controller_Y.CmdBusOut.wSat[0] = PX4Controller_P.Saturation_LowerSat_c[0];
  } else {
    /* Outport: '<Root>/CmdBusOut' */
    PX4Controller_Y.CmdBusOut.wSat[0] = rtb_Saturation7;
  }

  if (rtb_Saturation7_g > PX4Controller_P.Saturation_UpperSat_p[1]) {
    /* Outport: '<Root>/CmdBusOut' */
    PX4Controller_Y.CmdBusOut.wSat[1] = PX4Controller_P.Saturation_UpperSat_p[1];
  } else if (rtb_Saturation7_g < PX4Controller_P.Saturation_LowerSat_c[1]) {
    /* Outport: '<Root>/CmdBusOut' */
    PX4Controller_Y.CmdBusOut.wSat[1] = PX4Controller_P.Saturation_LowerSat_c[1];
  } else {
    /* Outport: '<Root>/CmdBusOut' */
    PX4Controller_Y.CmdBusOut.wSat[1] = rtb_Saturation7_g;
  }

  if (rtb_InterpolationUsingPrelook_k > PX4Controller_P.Saturation_UpperSat_p[2])
  {
    /* Outport: '<Root>/CmdBusOut' */
    PX4Controller_Y.CmdBusOut.wSat[2] = PX4Controller_P.Saturation_UpperSat_p[2];
  } else if (rtb_InterpolationUsingPrelook_k <
             PX4Controller_P.Saturation_LowerSat_c[2]) {
    /* Outport: '<Root>/CmdBusOut' */
    PX4Controller_Y.CmdBusOut.wSat[2] = PX4Controller_P.Saturation_LowerSat_c[2];
  } else {
    /* Outport: '<Root>/CmdBusOut' */
    PX4Controller_Y.CmdBusOut.wSat[2] = rtb_InterpolationUsingPrelook_k;
  }

  /* End of Saturate: '<S29>/Saturation' */

  /* Interpolation_n-D: '<S22>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S22>/Constant'
   */
  frac_9[0] = rtb_f1;
  frac_9[1] = rtb_f2;
  bpIndex_9[0] = rtb_k1;
  bpIndex_9[1] = rtb_k2;
  for (int32_T i{0}; i < 84; i++) {
    if (PX4Controller_P.Constant_Value_n[i] > 83.0) {
      bpIndex_9[2] = 83U;
    } else if (PX4Controller_P.Constant_Value_n[i] >= 0.0) {
      bpIndex_9[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_n[i]);
    } else {
      bpIndex_9[2] = 0U;
    }

    rtb_InterpolationUsingPrelook_g[i] = intrp2d_l_pw(bpIndex_9, frac_9,
      &PX4Controller_P.InterpolationUsingPrelookup__cj[bpIndex_9[2] << 2], 2U);
  }

  /* End of Interpolation_n-D: '<S22>/Interpolation Using Prelookup' */
  if (rtmIsMajorTimeStep(PX4Controller_M)) {
    for (int32_T i{0}; i < 14; i++) {
      /* Selector: '<S14>/Selector5' incorporates:
       *  Constant: '<S14>/allLimits'
       */
      PX4Controller_B.Selector5[i] = PX4Controller_P.allLimits_Value[(i << 1) +
        1];
    }
  }

  /* Interpolation_n-D: '<S17>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S17>/Constant'
   */
  frac_b[0] = rtb_f1_m;
  frac_b[1] = rtb_f2_p;
  frac_b[2] = rtb_f3;
  bpIndex_b[0] = rtb_k1_g;
  bpIndex_b[1] = rtb_k2_e;
  bpIndex_b[2] = rtb_k3;
  if (PX4Controller_P.Constant_Value_e[0] > 3.0) {
    bpIndex_b[3] = 3U;
  } else if (PX4Controller_P.Constant_Value_e[0] >= 0.0) {
    bpIndex_b[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_e[0]);
  } else {
    bpIndex_b[3] = 0U;
  }

  rtb_spe = intrp3d_l_pw(bpIndex_b, frac_b,
    &PX4Controller_P.SubsystemReference3_table_p[12U * bpIndex_b[3]],
    PX4Controller_P.InterpolationUsingPrelookup_d_d);
  if (PX4Controller_P.Constant_Value_e[1] > 3.0) {
    bpIndex_b[3] = 3U;
  } else if (PX4Controller_P.Constant_Value_e[1] >= 0.0) {
    bpIndex_b[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_e[1]);
  } else {
    bpIndex_b[3] = 0U;
  }

  rtb_skeCmd = intrp3d_l_pw(bpIndex_b, frac_b,
    &PX4Controller_P.SubsystemReference3_table_p[12U * bpIndex_b[3]],
    PX4Controller_P.InterpolationUsingPrelookup_d_d);
  if (PX4Controller_P.Constant_Value_e[2] > 3.0) {
    bpIndex_b[3] = 3U;
  } else if (PX4Controller_P.Constant_Value_e[2] >= 0.0) {
    bpIndex_b[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_e[2]);
  } else {
    bpIndex_b[3] = 0U;
  }

  rtb_Gain2 = intrp3d_l_pw(bpIndex_b, frac_b,
    &PX4Controller_P.SubsystemReference3_table_p[12U * bpIndex_b[3]],
    PX4Controller_P.InterpolationUsingPrelookup_d_d);
  if (PX4Controller_P.Constant_Value_e[3] > 3.0) {
    bpIndex_b[3] = 3U;
  } else if (PX4Controller_P.Constant_Value_e[3] >= 0.0) {
    bpIndex_b[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_e[3]);
  } else {
    bpIndex_b[3] = 0U;
  }

  rtb_Saturation1 = intrp3d_l_pw(bpIndex_b, frac_b,
    &PX4Controller_P.SubsystemReference3_table_p[12U * bpIndex_b[3]],
    PX4Controller_P.InterpolationUsingPrelookup_d_d);

  /* End of Interpolation_n-D: '<S17>/Interpolation Using Prelookup' */

  /* Interpolation_n-D: '<S18>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S18>/Constant'
   */
  frac_c[0] = rtb_f1_m;
  frac_c[1] = rtb_f2_p;
  frac_c[2] = rtb_f3;
  bpIndex_c[0] = rtb_k1_g;
  bpIndex_c[1] = rtb_k2_e;
  bpIndex_c[2] = rtb_k3;
  for (int32_T i{0}; i < 6; i++) {
    if (PX4Controller_P.Constant_Value_m[i] > 5.0) {
      bpIndex_c[3] = 5U;
    } else if (PX4Controller_P.Constant_Value_m[i] >= 0.0) {
      bpIndex_c[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_m[i]);
    } else {
      bpIndex_c[3] = 0U;
    }

    rtb_MatrixMultiply1_k[i] = intrp3d_l_pw(bpIndex_c, frac_c,
      &PX4Controller_P.SubsystemReference4_table[12U * bpIndex_c[3]],
      PX4Controller_P.InterpolationUsingPrelookup_d_o);
  }

  /* End of Interpolation_n-D: '<S18>/Interpolation Using Prelookup' */

  /* Interpolation_n-D: '<S19>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S19>/Constant'
   */
  frac_d[0] = rtb_f1_m;
  frac_d[1] = rtb_f2_p;
  frac_d[2] = rtb_f3;
  bpIndex_d[0] = rtb_k1_g;
  bpIndex_d[1] = rtb_k2_e;
  bpIndex_d[2] = rtb_k3;
  if (PX4Controller_P.Constant_Value_e4[0] > 1.0) {
    bpIndex_d[3] = 1U;
  } else if (PX4Controller_P.Constant_Value_e4[0] >= 0.0) {
    bpIndex_d[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_e4[0]);
  } else {
    bpIndex_d[3] = 0U;
  }

  frac_0[0] = intrp3d_l_pw(bpIndex_d, frac_d,
    &PX4Controller_P.SubsystemReference5_table[12U * bpIndex_d[3]],
    PX4Controller_P.InterpolationUsingPrelookup_d_m);
  if (PX4Controller_P.Constant_Value_e4[1] > 1.0) {
    bpIndex_d[3] = 1U;
  } else if (PX4Controller_P.Constant_Value_e4[1] >= 0.0) {
    bpIndex_d[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_e4[1]);
  } else {
    bpIndex_d[3] = 0U;
  }

  frac_0[1] = intrp3d_l_pw(bpIndex_d, frac_d,
    &PX4Controller_P.SubsystemReference5_table[12U * bpIndex_d[3]],
    PX4Controller_P.InterpolationUsingPrelookup_d_m);

  /* End of Interpolation_n-D: '<S19>/Interpolation Using Prelookup' */

  /* Interpolation_n-D: '<S20>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S20>/Constant'
   */
  frac_e[0] = rtb_f1_m;
  frac_e[1] = rtb_f2_p;
  frac_e[2] = rtb_f3;
  bpIndex_e[0] = rtb_k1_g;
  bpIndex_e[1] = rtb_k2_e;
  bpIndex_e[2] = rtb_k3;
  if (PX4Controller_P.Constant_Value_ns[0] > 1.0) {
    bpIndex_e[3] = 1U;
  } else if (PX4Controller_P.Constant_Value_ns[0] >= 0.0) {
    bpIndex_e[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_ns[0]);
  } else {
    bpIndex_e[3] = 0U;
  }

  frac_3[0] = intrp3d_l_pw(bpIndex_e, frac_e,
    &PX4Controller_P.SubsystemReference6_table[12U * bpIndex_e[3]],
    PX4Controller_P.InterpolationUsingPrelookup__mi);
  if (PX4Controller_P.Constant_Value_ns[1] > 1.0) {
    bpIndex_e[3] = 1U;
  } else if (PX4Controller_P.Constant_Value_ns[1] >= 0.0) {
    bpIndex_e[3] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_ns[1]);
  } else {
    bpIndex_e[3] = 0U;
  }

  frac_3[1] = intrp3d_l_pw(bpIndex_e, frac_e,
    &PX4Controller_P.SubsystemReference6_table[12U * bpIndex_e[3]],
    PX4Controller_P.InterpolationUsingPrelookup__mi);

  /* End of Interpolation_n-D: '<S20>/Interpolation Using Prelookup' */

  /* Interpolation_n-D: '<S21>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S21>/Constant'
   */
  frac_f[0] = rtb_f1;
  frac_f[1] = rtb_f2;
  bpIndex_f[0] = rtb_k1;
  bpIndex_f[1] = rtb_k2;
  for (int32_T i{0}; i < 84; i++) {
    if (PX4Controller_P.Constant_Value_fu[i] > 83.0) {
      bpIndex_f[2] = 83U;
    } else if (PX4Controller_P.Constant_Value_fu[i] >= 0.0) {
      bpIndex_f[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_fu[i]);
    } else {
      bpIndex_f[2] = 0U;
    }

    rtb_InterpolationUsingPrelook_l[i] = intrp2d_l_pw(bpIndex_f, frac_f,
      &PX4Controller_P.InterpolationUsingPrelookup__nf[bpIndex_f[2] << 2], 2U);
  }

  /* End of Interpolation_n-D: '<S21>/Interpolation Using Prelookup' */

  /* Switch: '<S12>/Switch1' incorporates:
   *  Gain: '<S12>/Gain'
   *  Inport: '<Root>/CmdBusIn'
   */
  if (PX4Controller_U.CmdBusIn.manualFM > PX4Controller_P.Switch1_Threshold_l) {
    rtb_InterpolationUsingPrelook_k = PX4Controller_P.Gain_Gain *
      PX4Controller_U.CmdBusIn.rc[0];
  } else {
    rtb_InterpolationUsingPrelook_k = rtb_Switch1_b;
  }

  /* End of Switch: '<S12>/Switch1' */

  /* Switch: '<S12>/Switch7' incorporates:
   *  Gain: '<S12>/Gain1'
   *  Gain: '<S12>/Gain2'
   *  Gain: '<S12>/Gain3'
   *  Inport: '<Root>/CmdBusIn'
   */
  if (PX4Controller_U.CmdBusIn.manualFM > PX4Controller_P.Switch7_Threshold_e) {
    rtb_Assignment1_idx_0 = PX4Controller_P.Gain1_Gain *
      PX4Controller_U.CmdBusIn.rc[1];
    rtb_Assignment1_idx_1 = PX4Controller_P.Gain2_Gain *
      PX4Controller_U.CmdBusIn.rc[2];
    rtb_Assignment1_idx_2 = PX4Controller_P.Gain3_Gain *
      PX4Controller_U.CmdBusIn.rc[3];
  } else {
    rtb_Assignment1_idx_0 = frac_2[0];
    rtb_Assignment1_idx_1 = frac_2[1];
    rtb_Assignment1_idx_2 = frac_2[2];
  }

  /* End of Switch: '<S12>/Switch7' */

  /* Product: '<S13>/Product' incorporates:
   *  Constant: '<S13>/Constant5'
   *  Constant: '<S13>/Constant6'
   */
  rtb_InterpolationUsingPrelook_d[0] = rtb_InterpolationUsingPrelook_k *
    PX4Controller_P.Constant6_Value[0];
  rtb_InterpolationUsingPrelook_d[1] = PX4Controller_P.Constant5_Value[0] *
    PX4Controller_P.Constant6_Value[1];
  rtb_InterpolationUsingPrelook_d[2] = PX4Controller_P.Constant5_Value[1] *
    PX4Controller_P.Constant6_Value[2];
  rtb_InterpolationUsingPrelook_d[3] = rtb_Assignment1_idx_0 *
    PX4Controller_P.Constant6_Value[3];
  rtb_InterpolationUsingPrelook_d[4] = rtb_Assignment1_idx_1 *
    PX4Controller_P.Constant6_Value[4];
  rtb_InterpolationUsingPrelook_d[5] = rtb_Assignment1_idx_2 *
    PX4Controller_P.Constant6_Value[5];

  /* Product: '<S14>/Matrix Multiply' incorporates:
   *  Interpolation_n-D: '<S26>/Interpolation Using Prelookup'
   */
  std::memset(&rtb_Sum2_d[0], 0, 14U * sizeof(real_T));
  for (int32_T i_0{0}; i_0 < 6; i_0++) {
    for (int32_T i{0}; i < 14; i++) {
      rtb_Sum2_d[i] += rtb_InterpolationUsingPrelook_l[14 * i_0 + i] *
        rtb_InterpolationUsingPrelook_d[i_0];
    }
  }

  /* End of Product: '<S14>/Matrix Multiply' */

  /* Saturate: '<S14>/Saturation5' */
  if (rtb_Sum2_d[0] > PX4Controller_P.Saturation5_UpperSat) {
    rtb_ske_0 = PX4Controller_P.Saturation5_UpperSat;
  } else if (rtb_Sum2_d[0] < PX4Controller_P.Saturation5_LowerSat) {
    rtb_ske_0 = PX4Controller_P.Saturation5_LowerSat;
  } else {
    rtb_ske_0 = rtb_Sum2_d[0];
  }

  /* Product: '<S14>/Divide' incorporates:
   *  Product: '<S16>/Divide2'
   *  Selector: '<S14>/Selector5'
   *  Sum: '<S14>/Sum'
   */
  rtb_Divide2[0] = (PX4Controller_B.Selector5[0] - rtb_spe) / rtb_ske_0;

  /* Saturate: '<S14>/Saturation5' */
  if (rtb_Sum2_d[1] > PX4Controller_P.Saturation5_UpperSat) {
    rtb_ske_0 = PX4Controller_P.Saturation5_UpperSat;
  } else if (rtb_Sum2_d[1] < PX4Controller_P.Saturation5_LowerSat) {
    rtb_ske_0 = PX4Controller_P.Saturation5_LowerSat;
  } else {
    rtb_ske_0 = rtb_Sum2_d[1];
  }

  /* Product: '<S14>/Divide' incorporates:
   *  Product: '<S16>/Divide2'
   *  Selector: '<S14>/Selector5'
   *  Sum: '<S14>/Sum'
   */
  rtb_Divide2[1] = (PX4Controller_B.Selector5[1] - rtb_skeCmd) / rtb_ske_0;

  /* Saturate: '<S14>/Saturation5' */
  if (rtb_Sum2_d[2] > PX4Controller_P.Saturation5_UpperSat) {
    rtb_ske_0 = PX4Controller_P.Saturation5_UpperSat;
  } else if (rtb_Sum2_d[2] < PX4Controller_P.Saturation5_LowerSat) {
    rtb_ske_0 = PX4Controller_P.Saturation5_LowerSat;
  } else {
    rtb_ske_0 = rtb_Sum2_d[2];
  }

  /* Product: '<S14>/Divide' incorporates:
   *  Product: '<S16>/Divide2'
   *  Selector: '<S14>/Selector5'
   *  Sum: '<S14>/Sum'
   */
  rtb_Divide2[2] = (PX4Controller_B.Selector5[2] - rtb_Gain2) / rtb_ske_0;

  /* Saturate: '<S14>/Saturation5' */
  if (rtb_Sum2_d[3] > PX4Controller_P.Saturation5_UpperSat) {
    rtb_ske_0 = PX4Controller_P.Saturation5_UpperSat;
  } else if (rtb_Sum2_d[3] < PX4Controller_P.Saturation5_LowerSat) {
    rtb_ske_0 = PX4Controller_P.Saturation5_LowerSat;
  } else {
    rtb_ske_0 = rtb_Sum2_d[3];
  }

  /* Product: '<S14>/Divide' incorporates:
   *  Product: '<S16>/Divide2'
   *  Selector: '<S14>/Selector5'
   *  Sum: '<S14>/Sum'
   */
  rtb_Divide2[3] = (PX4Controller_B.Selector5[3] - rtb_Saturation1) / rtb_ske_0;
  for (int32_T i{0}; i < 6; i++) {
    /* Saturate: '<S14>/Saturation5' */
    rtb_f3 = rtb_Sum2_d[i + 4];
    if (rtb_f3 > PX4Controller_P.Saturation5_UpperSat) {
      rtb_f3 = PX4Controller_P.Saturation5_UpperSat;
    } else if (rtb_f3 < PX4Controller_P.Saturation5_LowerSat) {
      rtb_f3 = PX4Controller_P.Saturation5_LowerSat;
    }

    rtb_Divide2[i + 4] = (PX4Controller_B.Selector5[i + 4] -
                          rtb_MatrixMultiply1_k[i]) / rtb_f3;
  }

  /* Saturate: '<S14>/Saturation5' */
  if (rtb_Sum2_d[10] > PX4Controller_P.Saturation5_UpperSat) {
    rtb_ske_0 = PX4Controller_P.Saturation5_UpperSat;
  } else if (rtb_Sum2_d[10] < PX4Controller_P.Saturation5_LowerSat) {
    rtb_ske_0 = PX4Controller_P.Saturation5_LowerSat;
  } else {
    rtb_ske_0 = rtb_Sum2_d[10];
  }

  /* Product: '<S14>/Divide' incorporates:
   *  Product: '<S16>/Divide2'
   *  Selector: '<S14>/Selector5'
   *  Sum: '<S14>/Sum'
   */
  rtb_Divide2[10] = (PX4Controller_B.Selector5[10] - frac_0[0]) / rtb_ske_0;

  /* Saturate: '<S14>/Saturation5' */
  if (rtb_Sum2_d[12] > PX4Controller_P.Saturation5_UpperSat) {
    rtb_ske_0 = PX4Controller_P.Saturation5_UpperSat;
  } else if (rtb_Sum2_d[12] < PX4Controller_P.Saturation5_LowerSat) {
    rtb_ske_0 = PX4Controller_P.Saturation5_LowerSat;
  } else {
    rtb_ske_0 = rtb_Sum2_d[12];
  }

  /* Product: '<S14>/Divide' incorporates:
   *  Product: '<S16>/Divide2'
   *  Selector: '<S14>/Selector5'
   *  Sum: '<S14>/Sum'
   */
  rtb_Divide2[12] = (PX4Controller_B.Selector5[12] - frac_3[0]) / rtb_ske_0;

  /* Saturate: '<S14>/Saturation5' */
  if (rtb_Sum2_d[11] > PX4Controller_P.Saturation5_UpperSat) {
    rtb_ske_0 = PX4Controller_P.Saturation5_UpperSat;
  } else if (rtb_Sum2_d[11] < PX4Controller_P.Saturation5_LowerSat) {
    rtb_ske_0 = PX4Controller_P.Saturation5_LowerSat;
  } else {
    rtb_ske_0 = rtb_Sum2_d[11];
  }

  /* Product: '<S14>/Divide' incorporates:
   *  Product: '<S16>/Divide2'
   *  Selector: '<S14>/Selector5'
   *  Sum: '<S14>/Sum'
   */
  rtb_Divide2[11] = (PX4Controller_B.Selector5[11] - frac_0[1]) / rtb_ske_0;

  /* Saturate: '<S14>/Saturation5' */
  if (rtb_Sum2_d[13] > PX4Controller_P.Saturation5_UpperSat) {
    rtb_ske_0 = PX4Controller_P.Saturation5_UpperSat;
  } else if (rtb_Sum2_d[13] < PX4Controller_P.Saturation5_LowerSat) {
    rtb_ske_0 = PX4Controller_P.Saturation5_LowerSat;
  } else {
    rtb_ske_0 = rtb_Sum2_d[13];
  }

  /* Product: '<S14>/Divide' incorporates:
   *  Product: '<S16>/Divide2'
   *  Selector: '<S14>/Selector5'
   *  Sum: '<S14>/Sum'
   */
  rtb_Divide2[13] = (PX4Controller_B.Selector5[13] - frac_3[1]) / rtb_ske_0;

  /* MinMax: '<S14>/Min' incorporates:
   *  Product: '<S16>/Divide2'
   */
  rtb_ske_0 = rtb_Divide2[0];
  for (int32_T i{0}; i < 13; i++) {
    rtb_ske_0 = std::fmin(rtb_ske_0, rtb_Divide2[i + 1]);
  }

  /* End of MinMax: '<S14>/Min' */

  /* Saturate: '<S14>/Saturation4' */
  if (rtb_ske_0 > PX4Controller_P.Saturation4_UpperSat) {
    rtb_Saturation7_g = PX4Controller_P.Saturation4_UpperSat;
  } else if (rtb_ske_0 < PX4Controller_P.Saturation4_LowerSat) {
    rtb_Saturation7_g = PX4Controller_P.Saturation4_LowerSat;
  } else {
    rtb_Saturation7_g = rtb_ske_0;
  }

  /* End of Saturate: '<S14>/Saturation4' */
  if (rtmIsMajorTimeStep(PX4Controller_M)) {
    for (int32_T i{0}; i < 14; i++) {
      /* Selector: '<S14>/Selector4' incorporates:
       *  Constant: '<S14>/allLimits'
       */
      PX4Controller_B.Selector4[i] = PX4Controller_P.allLimits_Value[i << 1];
    }
  }

  /* Saturate: '<S14>/Saturation6' */
  if (rtb_Sum2_d[0] > PX4Controller_P.Saturation6_UpperSat) {
    rtb_ske_0 = PX4Controller_P.Saturation6_UpperSat;
  } else if (rtb_Sum2_d[0] < PX4Controller_P.Saturation6_LowerSat) {
    rtb_ske_0 = PX4Controller_P.Saturation6_LowerSat;
  } else {
    rtb_ske_0 = rtb_Sum2_d[0];
  }

  /* Product: '<S14>/Divide2' incorporates:
   *  Selector: '<S14>/Selector4'
   *  Sum: '<S14>/Sum1'
   *  Sum: '<S15>/Sum2'
   */
  rtb_Sum2_f[0] = 1.0 / rtb_ske_0 * (PX4Controller_B.Selector4[0] - rtb_spe);

  /* Saturate: '<S14>/Saturation6' */
  if (rtb_Sum2_d[1] > PX4Controller_P.Saturation6_UpperSat) {
    rtb_ske_0 = PX4Controller_P.Saturation6_UpperSat;
  } else if (rtb_Sum2_d[1] < PX4Controller_P.Saturation6_LowerSat) {
    rtb_ske_0 = PX4Controller_P.Saturation6_LowerSat;
  } else {
    rtb_ske_0 = rtb_Sum2_d[1];
  }

  /* Product: '<S14>/Divide2' incorporates:
   *  Selector: '<S14>/Selector4'
   *  Sum: '<S14>/Sum1'
   *  Sum: '<S15>/Sum2'
   */
  rtb_Sum2_f[1] = 1.0 / rtb_ske_0 * (PX4Controller_B.Selector4[1] - rtb_skeCmd);

  /* Saturate: '<S14>/Saturation6' */
  if (rtb_Sum2_d[2] > PX4Controller_P.Saturation6_UpperSat) {
    rtb_ske_0 = PX4Controller_P.Saturation6_UpperSat;
  } else if (rtb_Sum2_d[2] < PX4Controller_P.Saturation6_LowerSat) {
    rtb_ske_0 = PX4Controller_P.Saturation6_LowerSat;
  } else {
    rtb_ske_0 = rtb_Sum2_d[2];
  }

  /* Product: '<S14>/Divide2' incorporates:
   *  Selector: '<S14>/Selector4'
   *  Sum: '<S14>/Sum1'
   *  Sum: '<S15>/Sum2'
   */
  rtb_Sum2_f[2] = 1.0 / rtb_ske_0 * (PX4Controller_B.Selector4[2] - rtb_Gain2);

  /* Saturate: '<S14>/Saturation6' */
  if (rtb_Sum2_d[3] > PX4Controller_P.Saturation6_UpperSat) {
    rtb_ske_0 = PX4Controller_P.Saturation6_UpperSat;
  } else if (rtb_Sum2_d[3] < PX4Controller_P.Saturation6_LowerSat) {
    rtb_ske_0 = PX4Controller_P.Saturation6_LowerSat;
  } else {
    rtb_ske_0 = rtb_Sum2_d[3];
  }

  /* Product: '<S14>/Divide2' incorporates:
   *  Selector: '<S14>/Selector4'
   *  Sum: '<S14>/Sum1'
   *  Sum: '<S15>/Sum2'
   */
  rtb_Sum2_f[3] = 1.0 / rtb_ske_0 * (PX4Controller_B.Selector4[3] -
    rtb_Saturation1);
  for (int32_T i{0}; i < 6; i++) {
    /* Saturate: '<S14>/Saturation6' */
    rtb_f3 = rtb_Sum2_d[i + 4];
    if (rtb_f3 > PX4Controller_P.Saturation6_UpperSat) {
      rtb_f3 = PX4Controller_P.Saturation6_UpperSat;
    } else if (rtb_f3 < PX4Controller_P.Saturation6_LowerSat) {
      rtb_f3 = PX4Controller_P.Saturation6_LowerSat;
    }

    rtb_Sum2_f[i + 4] = (PX4Controller_B.Selector4[i + 4] -
                         rtb_MatrixMultiply1_k[i]) * (1.0 / rtb_f3);
  }

  /* Saturate: '<S14>/Saturation6' */
  if (rtb_Sum2_d[10] > PX4Controller_P.Saturation6_UpperSat) {
    rtb_ske_0 = PX4Controller_P.Saturation6_UpperSat;
  } else if (rtb_Sum2_d[10] < PX4Controller_P.Saturation6_LowerSat) {
    rtb_ske_0 = PX4Controller_P.Saturation6_LowerSat;
  } else {
    rtb_ske_0 = rtb_Sum2_d[10];
  }

  /* Product: '<S14>/Divide2' incorporates:
   *  Selector: '<S14>/Selector4'
   *  Sum: '<S14>/Sum1'
   *  Sum: '<S15>/Sum2'
   */
  rtb_Sum2_f[10] = 1.0 / rtb_ske_0 * (PX4Controller_B.Selector4[10] - frac_0[0]);

  /* Saturate: '<S14>/Saturation6' */
  if (rtb_Sum2_d[12] > PX4Controller_P.Saturation6_UpperSat) {
    rtb_ske_0 = PX4Controller_P.Saturation6_UpperSat;
  } else if (rtb_Sum2_d[12] < PX4Controller_P.Saturation6_LowerSat) {
    rtb_ske_0 = PX4Controller_P.Saturation6_LowerSat;
  } else {
    rtb_ske_0 = rtb_Sum2_d[12];
  }

  /* Product: '<S14>/Divide2' incorporates:
   *  Selector: '<S14>/Selector4'
   *  Sum: '<S14>/Sum1'
   *  Sum: '<S15>/Sum2'
   */
  rtb_Sum2_f[12] = 1.0 / rtb_ske_0 * (PX4Controller_B.Selector4[12] - frac_3[0]);

  /* Saturate: '<S14>/Saturation6' */
  if (rtb_Sum2_d[11] > PX4Controller_P.Saturation6_UpperSat) {
    rtb_ske_0 = PX4Controller_P.Saturation6_UpperSat;
  } else if (rtb_Sum2_d[11] < PX4Controller_P.Saturation6_LowerSat) {
    rtb_ske_0 = PX4Controller_P.Saturation6_LowerSat;
  } else {
    rtb_ske_0 = rtb_Sum2_d[11];
  }

  /* Product: '<S14>/Divide2' incorporates:
   *  Selector: '<S14>/Selector4'
   *  Sum: '<S14>/Sum1'
   *  Sum: '<S15>/Sum2'
   */
  rtb_Sum2_f[11] = 1.0 / rtb_ske_0 * (PX4Controller_B.Selector4[11] - frac_0[1]);

  /* Saturate: '<S14>/Saturation6' */
  if (rtb_Sum2_d[13] > PX4Controller_P.Saturation6_UpperSat) {
    rtb_ske_0 = PX4Controller_P.Saturation6_UpperSat;
  } else if (rtb_Sum2_d[13] < PX4Controller_P.Saturation6_LowerSat) {
    rtb_ske_0 = PX4Controller_P.Saturation6_LowerSat;
  } else {
    rtb_ske_0 = rtb_Sum2_d[13];
  }

  /* Product: '<S14>/Divide2' incorporates:
   *  Selector: '<S14>/Selector4'
   *  Sum: '<S14>/Sum1'
   *  Sum: '<S15>/Sum2'
   */
  rtb_Sum2_f[13] = 1.0 / rtb_ske_0 * (PX4Controller_B.Selector4[13] - frac_3[1]);

  /* MinMax: '<S14>/Min1' incorporates:
   *  Sum: '<S15>/Sum2'
   */
  rtb_ske_0 = rtb_Sum2_f[0];
  for (int32_T i{0}; i < 13; i++) {
    rtb_ske_0 = std::fmin(rtb_ske_0, rtb_Sum2_f[i + 1]);
  }

  /* End of MinMax: '<S14>/Min1' */

  /* Saturate: '<S14>/Saturation7' */
  if (rtb_ske_0 > PX4Controller_P.Saturation7_UpperSat) {
    rtb_Saturation7 = PX4Controller_P.Saturation7_UpperSat;
  } else if (rtb_ske_0 < PX4Controller_P.Saturation7_LowerSat) {
    rtb_Saturation7 = PX4Controller_P.Saturation7_LowerSat;
  } else {
    rtb_Saturation7 = rtb_ske_0;
  }

  /* End of Saturate: '<S14>/Saturation7' */

  /* Product: '<S14>/Divide1' */
  for (int32_T i{0}; i < 14; i++) {
    rtb_Sum2_d[i] = rtb_Saturation7_g * rtb_Sum2_d[i] * rtb_Saturation7;
  }

  /* End of Product: '<S14>/Divide1' */

  /* Product: '<S14>/Matrix Multiply1' incorporates:
   *  Interpolation_n-D: '<S25>/Interpolation Using Prelookup'
   */
  for (int32_T i{0}; i < 6; i++) {
    rtb_Sum2_a[i] = 0.0;
  }

  for (int32_T i_0{0}; i_0 < 14; i_0++) {
    for (int32_T i{0}; i < 6; i++) {
      rtb_Sum2_a[i] += rtb_InterpolationUsingPrelook_g[6 * i_0 + i] *
        rtb_Sum2_d[i_0];
    }
  }

  /* End of Product: '<S14>/Matrix Multiply1' */

  /* Interpolation_n-D: '<S24>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S24>/Constant'
   */
  frac_g[0] = rtb_f1;
  frac_g[1] = rtb_f2;
  bpIndex_g[0] = rtb_k1;
  bpIndex_g[1] = rtb_k2;
  for (int32_T i{0}; i < 84; i++) {
    if (PX4Controller_P.Constant_Value_nh[i] > 83.0) {
      bpIndex_g[2] = 83U;
    } else if (PX4Controller_P.Constant_Value_nh[i] >= 0.0) {
      bpIndex_g[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_nh[i]);
    } else {
      bpIndex_g[2] = 0U;
    }

    rtb_InterpolationUsingPrelook_l[i] = intrp2d_l_pw(bpIndex_g, frac_g,
      &PX4Controller_P.InterpolationUsingPrelookup__bc[bpIndex_g[2] << 2], 2U);
  }

  /* End of Interpolation_n-D: '<S24>/Interpolation Using Prelookup' */
  if (rtmIsMajorTimeStep(PX4Controller_M)) {
    for (int32_T i{0}; i < 14; i++) {
      /* Selector: '<S15>/Selector5' incorporates:
       *  Constant: '<S15>/allLimits'
       */
      PX4Controller_B.Selector5_l[i] = PX4Controller_P.allLimits_Value_e[(i << 1)
        + 1];
    }
  }

  /* Sum: '<S14>/Sum2' */
  rtb_Divide2[0] = rtb_Sum2_d[0] + rtb_spe;
  rtb_Divide2[1] = rtb_Sum2_d[1] + rtb_skeCmd;
  rtb_Divide2[2] = rtb_Sum2_d[2] + rtb_Gain2;
  rtb_Divide2[3] = rtb_Sum2_d[3] + rtb_Saturation1;
  for (int32_T i{0}; i < 6; i++) {
    rtb_Divide2[i + 4] = rtb_Sum2_d[i + 4] + rtb_MatrixMultiply1_k[i];
  }

  rtb_Divide2[10] = frac_0[0] + rtb_Sum2_d[10];
  rtb_Divide2[12] = frac_3[0] + rtb_Sum2_d[12];
  rtb_Divide2[11] = frac_0[1] + rtb_Sum2_d[11];
  rtb_Divide2[13] = frac_3[1] + rtb_Sum2_d[13];
  std::memcpy(&rtb_Sum2_d[0], &rtb_Divide2[0], 14U * sizeof(real_T));

  /* End of Sum: '<S14>/Sum2' */

  /* Interpolation_n-D: '<S23>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S23>/Constant'
   */
  frac_h[0] = rtb_f1;
  frac_h[1] = rtb_f2;
  bpIndex_h[0] = rtb_k1;
  bpIndex_h[1] = rtb_k2;
  for (int32_T i{0}; i < 84; i++) {
    if (PX4Controller_P.Constant_Value_c[i] > 83.0) {
      bpIndex_h[2] = 83U;
    } else if (PX4Controller_P.Constant_Value_c[i] >= 0.0) {
      bpIndex_h[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_c[i]);
    } else {
      bpIndex_h[2] = 0U;
    }

    rtb_InterpolationUsingPrelook_g[i] = intrp2d_l_pw(bpIndex_h, frac_h,
      &PX4Controller_P.InterpolationUsingPrelookup_T_a[bpIndex_h[2] << 2], 2U);
  }

  /* End of Interpolation_n-D: '<S23>/Interpolation Using Prelookup' */

  /* Product: '<S13>/Product1' incorporates:
   *  Constant: '<S13>/Constant1'
   *  Constant: '<S13>/Constant5'
   */
  rtb_InterpolationUsingPrelook_d[0] = rtb_InterpolationUsingPrelook_k *
    PX4Controller_P.Constant1_Value_k[0];
  rtb_InterpolationUsingPrelook_d[1] = PX4Controller_P.Constant5_Value[0] *
    PX4Controller_P.Constant1_Value_k[1];
  rtb_InterpolationUsingPrelook_d[2] = PX4Controller_P.Constant5_Value[1] *
    PX4Controller_P.Constant1_Value_k[2];
  rtb_InterpolationUsingPrelook_d[3] = rtb_Assignment1_idx_0 *
    PX4Controller_P.Constant1_Value_k[3];
  rtb_InterpolationUsingPrelook_d[4] = rtb_Assignment1_idx_1 *
    PX4Controller_P.Constant1_Value_k[4];
  rtb_InterpolationUsingPrelook_d[5] = rtb_Assignment1_idx_2 *
    PX4Controller_P.Constant1_Value_k[5];

  /* Product: '<S15>/Matrix Multiply' incorporates:
   *  Interpolation_n-D: '<S25>/Interpolation Using Prelookup'
   */
  std::memset(&rtb_Sum2_f[0], 0, 14U * sizeof(real_T));
  for (int32_T i_0{0}; i_0 < 6; i_0++) {
    for (int32_T i{0}; i < 14; i++) {
      rtb_Sum2_f[i] += rtb_InterpolationUsingPrelook_g[14 * i_0 + i] *
        rtb_InterpolationUsingPrelook_d[i_0];
    }
  }

  /* End of Product: '<S15>/Matrix Multiply' */

  /* Product: '<S15>/Divide' incorporates:
   *  Product: '<S16>/Divide2'
   *  Selector: '<S15>/Selector5'
   *  Sum: '<S15>/Sum'
   */
  for (int32_T i{0}; i < 14; i++) {
    /* Saturate: '<S15>/Saturation5' */
    rtb_f3 = rtb_Sum2_f[i];
    if (rtb_f3 > PX4Controller_P.Saturation5_UpperSat_a) {
      rtb_f3 = PX4Controller_P.Saturation5_UpperSat_a;
    } else if (rtb_f3 < PX4Controller_P.Saturation5_LowerSat_g) {
      rtb_f3 = PX4Controller_P.Saturation5_LowerSat_g;
    }

    /* End of Saturate: '<S15>/Saturation5' */
    rtb_Divide2[i] = (PX4Controller_B.Selector5_l[i] - rtb_Sum2_d[i]) / rtb_f3;
  }

  /* End of Product: '<S15>/Divide' */

  /* MinMax: '<S15>/Min' incorporates:
   *  Product: '<S16>/Divide2'
   */
  rtb_ske_0 = rtb_Divide2[0];
  for (int32_T i{0}; i < 13; i++) {
    rtb_ske_0 = std::fmin(rtb_ske_0, rtb_Divide2[i + 1]);
  }

  /* End of MinMax: '<S15>/Min' */

  /* Saturate: '<S15>/Saturation4' */
  if (rtb_ske_0 > PX4Controller_P.Saturation4_UpperSat_k) {
    rtb_Saturation7_g = PX4Controller_P.Saturation4_UpperSat_k;
  } else if (rtb_ske_0 < PX4Controller_P.Saturation4_LowerSat_k) {
    rtb_Saturation7_g = PX4Controller_P.Saturation4_LowerSat_k;
  } else {
    rtb_Saturation7_g = rtb_ske_0;
  }

  /* End of Saturate: '<S15>/Saturation4' */
  if (rtmIsMajorTimeStep(PX4Controller_M)) {
    for (int32_T i{0}; i < 14; i++) {
      /* Selector: '<S15>/Selector4' incorporates:
       *  Constant: '<S15>/allLimits'
       */
      PX4Controller_B.Selector4_b[i] = PX4Controller_P.allLimits_Value_e[i << 1];
    }
  }

  /* Product: '<S15>/Divide2' incorporates:
   *  Product: '<S16>/Divide2'
   *  Selector: '<S15>/Selector4'
   *  Sum: '<S15>/Sum1'
   */
  for (int32_T i{0}; i < 14; i++) {
    /* Saturate: '<S15>/Saturation6' */
    rtb_f3 = rtb_Sum2_f[i];
    if (rtb_f3 > PX4Controller_P.Saturation6_UpperSat_h) {
      rtb_f3 = PX4Controller_P.Saturation6_UpperSat_h;
    } else if (rtb_f3 < PX4Controller_P.Saturation6_LowerSat_h) {
      rtb_f3 = PX4Controller_P.Saturation6_LowerSat_h;
    }

    /* End of Saturate: '<S15>/Saturation6' */
    rtb_Divide2[i] = 1.0 / rtb_f3 * (PX4Controller_B.Selector4_b[i] -
      rtb_Sum2_d[i]);
  }

  /* End of Product: '<S15>/Divide2' */

  /* MinMax: '<S15>/Min1' incorporates:
   *  Product: '<S16>/Divide2'
   */
  rtb_ske_0 = rtb_Divide2[0];
  for (int32_T i{0}; i < 13; i++) {
    rtb_ske_0 = std::fmin(rtb_ske_0, rtb_Divide2[i + 1]);
  }

  /* End of MinMax: '<S15>/Min1' */

  /* Saturate: '<S15>/Saturation7' */
  if (rtb_ske_0 > PX4Controller_P.Saturation7_UpperSat_p) {
    rtb_Saturation7 = PX4Controller_P.Saturation7_UpperSat_p;
  } else if (rtb_ske_0 < PX4Controller_P.Saturation7_LowerSat_d) {
    rtb_Saturation7 = PX4Controller_P.Saturation7_LowerSat_d;
  } else {
    rtb_Saturation7 = rtb_ske_0;
  }

  /* End of Saturate: '<S15>/Saturation7' */

  /* Product: '<S15>/Divide1' */
  for (int32_T i{0}; i < 14; i++) {
    rtb_Sum2_f[i] = rtb_Saturation7_g * rtb_Sum2_f[i] * rtb_Saturation7;
  }

  /* End of Product: '<S15>/Divide1' */

  /* Product: '<S15>/Matrix Multiply1' incorporates:
   *  Interpolation_n-D: '<S26>/Interpolation Using Prelookup'
   */
  for (int32_T i{0}; i < 6; i++) {
    rtb_MatrixMultiply1_k[i] = 0.0;
  }

  for (int32_T i_0{0}; i_0 < 14; i_0++) {
    for (int32_T i{0}; i < 6; i++) {
      rtb_MatrixMultiply1_k[i] += rtb_InterpolationUsingPrelook_l[6 * i_0 + i] *
        rtb_Sum2_f[i_0];
    }
  }

  /* End of Product: '<S15>/Matrix Multiply1' */

  /* Interpolation_n-D: '<S26>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S26>/Constant'
   */
  frac_i[0] = rtb_f1;
  frac_i[1] = rtb_f2;
  bpIndex_i[0] = rtb_k1;
  bpIndex_i[1] = rtb_k2;
  for (int32_T i{0}; i < 84; i++) {
    if (PX4Controller_P.Constant_Value_p[i] > 83.0) {
      bpIndex_i[2] = 83U;
    } else if (PX4Controller_P.Constant_Value_p[i] >= 0.0) {
      bpIndex_i[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_p[i]);
    } else {
      bpIndex_i[2] = 0U;
    }

    rtb_InterpolationUsingPrelook_l[i] = intrp2d_l_pw(bpIndex_i, frac_i,
      &PX4Controller_P.InterpolationUsingPrelookup_T_g[bpIndex_i[2] << 2], 2U);
  }

  /* End of Interpolation_n-D: '<S26>/Interpolation Using Prelookup' */
  if (rtmIsMajorTimeStep(PX4Controller_M)) {
    for (int32_T i{0}; i < 14; i++) {
      /* Selector: '<S16>/Selector5' incorporates:
       *  Constant: '<S16>/allLimits'
       */
      PX4Controller_B.Selector5_k[i] = PX4Controller_P.allLimits_Value_i[(i << 1)
        + 1];
    }
  }

  /* Sum: '<S15>/Sum2' */
  for (int32_T i{0}; i < 14; i++) {
    rtb_Sum2_f[i] += rtb_Sum2_d[i];
  }

  /* End of Sum: '<S15>/Sum2' */

  /* Interpolation_n-D: '<S25>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S25>/Constant'
   */
  frac_j[0] = rtb_f1;
  frac_j[1] = rtb_f2;
  bpIndex_j[0] = rtb_k1;
  bpIndex_j[1] = rtb_k2;
  for (int32_T i{0}; i < 84; i++) {
    if (PX4Controller_P.Constant_Value_pv[i] > 83.0) {
      bpIndex_j[2] = 83U;
    } else if (PX4Controller_P.Constant_Value_pv[i] >= 0.0) {
      bpIndex_j[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_pv[i]);
    } else {
      bpIndex_j[2] = 0U;
    }

    rtb_InterpolationUsingPrelook_g[i] = intrp2d_l_pw(bpIndex_j, frac_j,
      &PX4Controller_P.InterpolationUsingPrelookup__is[bpIndex_j[2] << 2], 2U);
  }

  /* End of Interpolation_n-D: '<S25>/Interpolation Using Prelookup' */

  /* Product: '<S13>/Product2' incorporates:
   *  Constant: '<S13>/Constant2'
   *  Constant: '<S13>/Constant5'
   */
  rtb_InterpolationUsingPrelook_d[0] = rtb_InterpolationUsingPrelook_k *
    PX4Controller_P.Constant2_Value_e[0];
  rtb_InterpolationUsingPrelook_d[1] = PX4Controller_P.Constant5_Value[0] *
    PX4Controller_P.Constant2_Value_e[1];
  rtb_InterpolationUsingPrelook_d[2] = PX4Controller_P.Constant5_Value[1] *
    PX4Controller_P.Constant2_Value_e[2];
  rtb_InterpolationUsingPrelook_d[3] = rtb_Assignment1_idx_0 *
    PX4Controller_P.Constant2_Value_e[3];
  rtb_InterpolationUsingPrelook_d[4] = rtb_Assignment1_idx_1 *
    PX4Controller_P.Constant2_Value_e[4];
  rtb_InterpolationUsingPrelook_d[5] = rtb_Assignment1_idx_2 *
    PX4Controller_P.Constant2_Value_e[5];

  /* Product: '<S16>/Matrix Multiply' incorporates:
   *  Interpolation_n-D: '<S25>/Interpolation Using Prelookup'
   */
  std::memset(&rtb_Sum2_d[0], 0, 14U * sizeof(real_T));
  for (int32_T i_0{0}; i_0 < 6; i_0++) {
    for (int32_T i{0}; i < 14; i++) {
      rtb_Sum2_d[i] += rtb_InterpolationUsingPrelook_g[14 * i_0 + i] *
        rtb_InterpolationUsingPrelook_d[i_0];
    }
  }

  /* End of Product: '<S16>/Matrix Multiply' */

  /* Product: '<S16>/Divide' incorporates:
   *  Product: '<S16>/Divide2'
   *  Selector: '<S16>/Selector5'
   *  Sum: '<S16>/Sum'
   */
  for (int32_T i{0}; i < 14; i++) {
    /* Saturate: '<S16>/Saturation5' */
    rtb_f3 = rtb_Sum2_d[i];
    if (rtb_f3 > PX4Controller_P.Saturation5_UpperSat_f) {
      rtb_f3 = PX4Controller_P.Saturation5_UpperSat_f;
    } else if (rtb_f3 < PX4Controller_P.Saturation5_LowerSat_f) {
      rtb_f3 = PX4Controller_P.Saturation5_LowerSat_f;
    }

    /* End of Saturate: '<S16>/Saturation5' */
    rtb_Divide2[i] = (PX4Controller_B.Selector5_k[i] - rtb_Sum2_f[i]) / rtb_f3;
  }

  /* End of Product: '<S16>/Divide' */

  /* MinMax: '<S16>/Min' incorporates:
   *  Product: '<S16>/Divide2'
   */
  rtb_ske_0 = rtb_Divide2[0];
  for (int32_T i{0}; i < 13; i++) {
    rtb_ske_0 = std::fmin(rtb_ske_0, rtb_Divide2[i + 1]);
  }

  /* End of MinMax: '<S16>/Min' */

  /* Saturate: '<S16>/Saturation4' */
  if (rtb_ske_0 > PX4Controller_P.Saturation4_UpperSat_m) {
    rtb_InterpolationUsingPrelook_k = PX4Controller_P.Saturation4_UpperSat_m;
  } else if (rtb_ske_0 < PX4Controller_P.Saturation4_LowerSat_o) {
    rtb_InterpolationUsingPrelook_k = PX4Controller_P.Saturation4_LowerSat_o;
  } else {
    rtb_InterpolationUsingPrelook_k = rtb_ske_0;
  }

  /* End of Saturate: '<S16>/Saturation4' */
  if (rtmIsMajorTimeStep(PX4Controller_M)) {
    for (int32_T i{0}; i < 14; i++) {
      /* Selector: '<S16>/Selector4' incorporates:
       *  Constant: '<S16>/allLimits'
       */
      PX4Controller_B.Selector4_g[i] = PX4Controller_P.allLimits_Value_i[i << 1];
    }
  }

  /* Product: '<S16>/Divide2' incorporates:
   *  Selector: '<S16>/Selector4'
   *  Sum: '<S16>/Sum1'
   */
  for (int32_T i{0}; i < 14; i++) {
    /* Saturate: '<S16>/Saturation6' */
    rtb_f3 = rtb_Sum2_d[i];
    if (rtb_f3 > PX4Controller_P.Saturation6_UpperSat_p) {
      rtb_f3 = PX4Controller_P.Saturation6_UpperSat_p;
    } else if (rtb_f3 < PX4Controller_P.Saturation6_LowerSat_g) {
      rtb_f3 = PX4Controller_P.Saturation6_LowerSat_g;
    }

    /* End of Saturate: '<S16>/Saturation6' */
    rtb_Divide2[i] = 1.0 / rtb_f3 * (PX4Controller_B.Selector4_g[i] -
      rtb_Sum2_f[i]);
  }

  /* End of Product: '<S16>/Divide2' */

  /* MinMax: '<S16>/Min1' incorporates:
   *  Product: '<S16>/Divide2'
   */
  rtb_ske_0 = rtb_Divide2[0];
  for (int32_T i{0}; i < 13; i++) {
    rtb_ske_0 = std::fmin(rtb_ske_0, rtb_Divide2[i + 1]);
  }

  /* End of MinMax: '<S16>/Min1' */

  /* Saturate: '<S16>/Saturation7' */
  if (rtb_ske_0 > PX4Controller_P.Saturation7_UpperSat_i) {
    rtb_Saturation7_g = PX4Controller_P.Saturation7_UpperSat_i;
  } else if (rtb_ske_0 < PX4Controller_P.Saturation7_LowerSat_g) {
    rtb_Saturation7_g = PX4Controller_P.Saturation7_LowerSat_g;
  } else {
    rtb_Saturation7_g = rtb_ske_0;
  }

  /* End of Saturate: '<S16>/Saturation7' */

  /* Product: '<S16>/Divide1' */
  for (int32_T i{0}; i < 14; i++) {
    rtb_Sum2_d[i] = rtb_InterpolationUsingPrelook_k * rtb_Sum2_d[i] *
      rtb_Saturation7_g;
  }

  /* End of Product: '<S16>/Divide1' */

  /* Sum: '<S2>/Sum2' incorporates:
   *  Interpolation_n-D: '<S26>/Interpolation Using Prelookup'
   *  Product: '<S16>/Matrix Multiply1'
   */
  for (int32_T i{0}; i < 6; i++) {
    rtb_Assignment1_idx_2 = 0.0;
    for (int32_T i_0{0}; i_0 < 14; i_0++) {
      rtb_Assignment1_idx_2 += rtb_InterpolationUsingPrelook_l[6 * i_0 + i] *
        rtb_Sum2_d[i_0];
    }

    rtb_Sum2_a[i] = (rtb_Sum2_a[i] + rtb_MatrixMultiply1_k[i]) +
      rtb_Assignment1_idx_2;
  }

  /* End of Sum: '<S2>/Sum2' */

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
  PX4Controller_Y.CmdBusOut.FCmd = rtb_Switch1_b;
  PX4Controller_Y.CmdBusOut.manualActuation =
    PX4Controller_U.CmdBusIn.manualActuation;
  PX4Controller_Y.CmdBusOut.trackLine[0] = PX4Controller_U.CmdBusIn.trackLine[0];
  PX4Controller_Y.CmdBusOut.eulerCmd[0] = rtb_Switch7_idx_0;
  PX4Controller_Y.CmdBusOut.wCmd[0] = rtb_Switch7_d_idx_0;
  PX4Controller_Y.CmdBusOut.MCmd[0] = frac_2[0];
  PX4Controller_Y.CmdBusOut.MSat[0] = rtb_Sum2_a[3];
  PX4Controller_Y.CmdBusOut.trackLine[1] = PX4Controller_U.CmdBusIn.trackLine[1];
  PX4Controller_Y.CmdBusOut.eulerCmd[1] = rtb_Switch7_idx_1;
  PX4Controller_Y.CmdBusOut.wCmd[1] = rtb_Switch7_d_idx_1;
  PX4Controller_Y.CmdBusOut.MCmd[1] = frac_2[1];
  PX4Controller_Y.CmdBusOut.MSat[1] = rtb_Sum2_a[4];
  PX4Controller_Y.CmdBusOut.trackLine[2] = PX4Controller_U.CmdBusIn.trackLine[2];
  PX4Controller_Y.CmdBusOut.wCmd[2] = rtb_Switch7_d_idx_2;
  PX4Controller_Y.CmdBusOut.MCmd[2] = frac_2[2];
  PX4Controller_Y.CmdBusOut.MSat[2] = rtb_Sum2_a[5];
  PX4Controller_Y.CmdBusOut.FSat = rtb_Sum2_a[0];
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

  /* Sum: '<S16>/Sum2' */
  for (int32_T i{0}; i < 14; i++) {
    rtb_Sum2_d[i] += rtb_Sum2_f[i];
  }

  /* End of Sum: '<S16>/Sum2' */

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
    PX4Controller_B.Switch1[0] = rtb_Sum2_d[0];
    PX4Controller_B.Switch1[1] = rtb_Sum2_d[1];
    PX4Controller_B.Switch1[2] = rtb_Sum2_d[2];
    PX4Controller_B.Switch1[3] = rtb_Sum2_d[3];
  }

  /* End of Switch: '<S11>/Switch1' */

  /* Switch: '<S11>/Switch2' incorporates:
   *  Inport: '<Root>/CmdBusIn'
   */
  if (PX4Controller_U.CmdBusIn.manualActuation >
      PX4Controller_P.Switch2_Threshold) {
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
      PX4Controller_B.Switch2[i] = rtb_Sum2_d[i + 4];
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
    PX4Controller_B.Switch3[0] = rtb_Sum2_d[10];
    PX4Controller_B.Switch3[1] = rtb_Sum2_d[11];
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
    PX4Controller_B.Switch4[0] = rtb_Sum2_d[12];
    PX4Controller_B.Switch4[1] = rtb_Sum2_d[13];
  }

  /* End of Switch: '<S11>/Switch4' */

  /* Interpolation_n-D: '<S35>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S35>/Constant'
   */
  frac_k[0] = rtb_f1;
  frac_k[1] = rtb_f2;
  bpIndex_k[0] = rtb_k1;
  bpIndex_k[1] = rtb_k2;
  if (PX4Controller_P.Constant_Value_cn[0] > 1.0) {
    bpIndex_k[2] = 1U;
  } else if (PX4Controller_P.Constant_Value_cn[0] >= 0.0) {
    bpIndex_k[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_cn[0]);
  } else {
    bpIndex_k[2] = 0U;
  }

  rtb_Switch7_d_idx_0 = intrp2d_l_pw(bpIndex_k, frac_k,
    &PX4Controller_P.InterpolationUsingPrelookup__ho[bpIndex_k[2] << 2], 2U);

  /* Product: '<S35>/Product' */
  PX4Controller_B.Product[0] = frac[0] * rtb_Switch7_d_idx_0;

  /* Interpolation_n-D: '<S35>/Interpolation Using Prelookup' incorporates:
   *  Constant: '<S35>/Constant'
   */
  if (PX4Controller_P.Constant_Value_cn[1] > 1.0) {
    bpIndex_k[2] = 1U;
  } else if (PX4Controller_P.Constant_Value_cn[1] >= 0.0) {
    bpIndex_k[2] = static_cast<uint32_T>(PX4Controller_P.Constant_Value_cn[1]);
  } else {
    bpIndex_k[2] = 0U;
  }

  rtb_Switch7_d_idx_0 = intrp2d_l_pw(bpIndex_k, frac_k,
    &PX4Controller_P.InterpolationUsingPrelookup__ho[bpIndex_k[2] << 2], 2U);

  /* Product: '<S35>/Product' */
  PX4Controller_B.Product[1] = frac[1] * rtb_Switch7_d_idx_0;

  /* Interpolation_n-D: '<S80>/Interpolation Using Prelookup' */
  frac_l[0] = rtb_f1;
  frac_l[1] = rtb_f2;
  bpIndex_l[0] = rtb_k1;
  bpIndex_l[1] = rtb_k2;

  /* Product: '<S80>/Product' incorporates:
   *  Interpolation_n-D: '<S80>/Interpolation Using Prelookup'
   */
  PX4Controller_B.Product_g = rtb_sebErr * intrp2d_l_pw(bpIndex_l, frac_l,
    PX4Controller_P.InterpolationUsingPrelookup_T_m, 2U);

  /* Sum: '<S78>/Sum1' incorporates:
   *  Sum: '<S78>/steDot'
   *  Sum: '<S78>/steDotCmd '
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

  /* Derivatives for Integrator: '<S77>/Integrator1' */
  if (!PX4Controller_B.OR) {
    _rtXdot->Integrator1_CSTATE = PX4Controller_B.Product_g;
  } else {
    /* level reset is active */
    _rtXdot->Integrator1_CSTATE = 0.0;
  }

  /* End of Derivatives for Integrator: '<S77>/Integrator1' */

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

  /* Derivatives for TransferFcn: '<S78>/Low pass' */
  _rtXdot->Lowpass_CSTATE = PX4Controller_P.Lowpass_A *
    PX4Controller_X.Lowpass_CSTATE;
  _rtXdot->Lowpass_CSTATE += PX4Controller_B.steDotErr;

  /* Derivatives for Integrator: '<S78>/Integrator' */
  if (!PX4Controller_B.OR) {
    _rtXdot->Integrator_CSTATE_f = PX4Controller_B.steErr;
  } else {
    /* level reset is active */
    _rtXdot->Integrator_CSTATE_f = 0.0;
  }

  /* End of Derivatives for Integrator: '<S78>/Integrator' */
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
  PX4Controller_P.allLimits_Value[0] = rtMinusInf;
  PX4Controller_P.allLimits_Value[2] = rtMinusInf;
  PX4Controller_P.allLimits_Value[3] = rtInf;
  PX4Controller_P.allLimits_Value[4] = rtMinusInf;
  PX4Controller_P.allLimits_Value[5] = rtInf;
  PX4Controller_P.allLimits_Value[6] = rtMinusInf;
  PX4Controller_P.allLimits_Value[7] = rtInf;
  PX4Controller_P.allLimits_Value[12] = rtMinusInf;
  PX4Controller_P.allLimits_Value[13] = rtInf;
  PX4Controller_P.allLimits_Value[14] = rtMinusInf;
  PX4Controller_P.allLimits_Value[15] = rtInf;
  PX4Controller_P.allLimits_Value[22] = rtMinusInf;
  PX4Controller_P.allLimits_Value[23] = rtInf;
  PX4Controller_P.allLimits_Value[26] = rtMinusInf;
  PX4Controller_P.allLimits_Value[27] = rtInf;
  PX4Controller_P.Saturation5_UpperSat = rtInf;
  PX4Controller_P.Saturation6_LowerSat = rtMinusInf;
  PX4Controller_P.allLimits_Value_e[0] = rtMinusInf;
  PX4Controller_P.allLimits_Value_e[2] = rtMinusInf;
  PX4Controller_P.allLimits_Value_e[3] = rtInf;
  PX4Controller_P.allLimits_Value_e[4] = rtMinusInf;
  PX4Controller_P.allLimits_Value_e[5] = rtInf;
  PX4Controller_P.allLimits_Value_e[6] = rtMinusInf;
  PX4Controller_P.allLimits_Value_e[7] = rtInf;
  PX4Controller_P.allLimits_Value_e[12] = rtMinusInf;
  PX4Controller_P.allLimits_Value_e[13] = rtInf;
  PX4Controller_P.allLimits_Value_e[14] = rtMinusInf;
  PX4Controller_P.allLimits_Value_e[15] = rtInf;
  PX4Controller_P.allLimits_Value_e[22] = rtMinusInf;
  PX4Controller_P.allLimits_Value_e[23] = rtInf;
  PX4Controller_P.allLimits_Value_e[26] = rtMinusInf;
  PX4Controller_P.allLimits_Value_e[27] = rtInf;
  PX4Controller_P.Saturation5_UpperSat_a = rtInf;
  PX4Controller_P.Saturation6_LowerSat_h = rtMinusInf;
  PX4Controller_P.allLimits_Value_i[0] = rtMinusInf;
  PX4Controller_P.allLimits_Value_i[2] = rtMinusInf;
  PX4Controller_P.allLimits_Value_i[3] = rtInf;
  PX4Controller_P.allLimits_Value_i[4] = rtMinusInf;
  PX4Controller_P.allLimits_Value_i[5] = rtInf;
  PX4Controller_P.allLimits_Value_i[6] = rtMinusInf;
  PX4Controller_P.allLimits_Value_i[7] = rtInf;
  PX4Controller_P.allLimits_Value_i[12] = rtMinusInf;
  PX4Controller_P.allLimits_Value_i[13] = rtInf;
  PX4Controller_P.allLimits_Value_i[14] = rtMinusInf;
  PX4Controller_P.allLimits_Value_i[15] = rtInf;
  PX4Controller_P.allLimits_Value_i[22] = rtMinusInf;
  PX4Controller_P.allLimits_Value_i[23] = rtInf;
  PX4Controller_P.allLimits_Value_i[26] = rtMinusInf;
  PX4Controller_P.allLimits_Value_i[27] = rtInf;
  PX4Controller_P.Saturation5_UpperSat_f = rtInf;
  PX4Controller_P.Saturation6_LowerSat_g = rtMinusInf;

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

  /* InitializeConditions for Integrator: '<S77>/Integrator1' */
  PX4Controller_X.Integrator1_CSTATE = PX4Controller_P.Integrator1_IC;

  /* InitializeConditions for Integrator: '<S3>/Integrator' */
  PX4Controller_X.Integrator_CSTATE[0] = PX4Controller_P.Integrator_IC;
  PX4Controller_X.Integrator_CSTATE[1] = PX4Controller_P.Integrator_IC;

  /* InitializeConditions for TransferFcn: '<S78>/Low pass' */
  PX4Controller_X.Lowpass_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S78>/Integrator' */
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
