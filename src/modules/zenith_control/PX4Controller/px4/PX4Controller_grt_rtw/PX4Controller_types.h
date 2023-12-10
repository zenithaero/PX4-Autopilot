/*
 * PX4Controller_types.h
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

#ifndef RTW_HEADER_PX4Controller_types_h_
#define RTW_HEADER_PX4Controller_types_h_
#include "rtwtypes.h"

/* Model Code Variants */
#ifndef DEFINED_TYPEDEF_FOR_CmdBus_
#define DEFINED_TYPEDEF_FOR_CmdBus_

struct CmdBus
{
  real_T rc[4];
  real_T trackLine[3];
  real_T hCmd;
  real_T hDotCmd;
  real_T tasCmd;
  real_T tasDotCmd;
  real_T manualAttitude;
  real_T eulerCmd[3];
  real_T manualRate;
  real_T wCmd[3];
  real_T manualFM;
  real_T MCmd[3];
  real_T FCmd;
  real_T manualActuation;
  real_T eulerSat[3];
  real_T wSat[3];
  real_T MSat[3];
  real_T FSat;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_StateBus_
#define DEFINED_TYPEDEF_FOR_StateBus_

struct StateBus
{
  real_T vDotMeas[3];
  real_T vDotBody[3];
  real_T qBody[4];
  real_T wDotBody[3];
  real_T vNed[3];
  real_T eulerBody[3];
  real_T vBody[3];
  real_T wBody[3];
  real_T xNed[3];
  real_T lla[3];
  real_T alpha;
  real_T beta;
  real_T tas;
  real_T vAirBody[3];
  real_T wAirBody[3];
  real_T q;
  real_T rho;
  real_T gamma;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_ActBus_
#define DEFINED_TYPEDEF_FOR_ActBus_

struct ActBus
{
  real_T motors[4];
  real_T ailerons[6];
  real_T elevators[2];
  real_T rudders[2];
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_EstBus_
#define DEFINED_TYPEDEF_FOR_EstBus_

struct EstBus
{
  real_T vNedEst[3];
  real_T eulerBodyEst[3];
  real_T wBodyEst[3];
  real_T wDotBodyEst[3];
  real_T xNedEst[3];
  real_T tasDotEst;
  real_T hDotEst;
  real_T alphaEst;
  real_T betaEst;
  real_T tasEst;
  real_T hEst;
  real_T qEst;
};

#endif

/* Parameters (default storage) */
typedef struct P_PX4Controller_T_ P_PX4Controller_T;

/* Forward declaration for rtModel */
typedef struct tag_RTM_PX4Controller_T RT_MODEL_PX4Controller_T;

#endif                                 /* RTW_HEADER_PX4Controller_types_h_ */
