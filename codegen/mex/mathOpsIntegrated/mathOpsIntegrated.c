/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mathOpsIntegrated.c
 *
 * Code generation for function 'mathOpsIntegrated'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mathOpsIntegrated.h"
#include "adder.h"

/* Function Definitions */
void mathOpsIntegrated(const emlrtStack *sp, real_T in1, real_T in2, real_T
  *added, real_T *multed)
{
  (void)sp;
  *added = adder(in1, in2);
  *multed = in1 * in2;
}

/* End of code generation (mathOpsIntegrated.c) */
