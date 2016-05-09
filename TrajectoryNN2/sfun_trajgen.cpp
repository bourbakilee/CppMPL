/* 2015.12.20, LI Yunsheng */
/*
±‡“Î÷∏¡Ó£∫
mex -I"D:\Dev\boost_1_60_0" -I"D:\ProDocs\CppMPL\CppMPL\libsqlite3" -I"D:\Dev\eigen3" -I"D:\ProDocs\CppMPL\CppMPL\TrajectoryGeneration" -L"D:\ProDocs\CppMPL\CppMPL\x64\Release" -lsqlite3 -lTrajectoryGeneration s_fun_spiral3_calc.cpp

*/
// *******************************************************************
// **** To build this mex function use: mex s_fun_spiral3_calc.cpp ****
// *******************************************************************



#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME  sfun_trajgen

// Need to include simstruc.h for the definition of the SimStruct and
// its associated macro definitions.
#include <simstruc.h>

#include "TrajectoryNN2.h"
using namespace TrajectoryNN2;

#define IS_PARAM_DOUBLE(pVal) (mxIsNumeric(pVal) && !mxIsLogical(pVal) &&\
!mxIsEmpty(pVal) && !mxIsSparse(pVal) && !mxIsComplex(pVal) && mxIsDouble(pVal))

// Function: mdlInitializeSizes ===============================================
// Abstract:
//    The sizes information is used by Simulink to determine the S-function
//    block's characteristics (number of inputs, outputs, states, etc.).
static void mdlInitializeSizes(SimStruct *S)
{
	// No expected parameters
	ssSetNumSFcnParams(S, 0);

	// Parameter mismatch will be reported by Simulink
	if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
		return;
	}

	// Specify I/O
	if (!ssSetNumInputPorts(S, 4)) return;
	//In0
	ssSetInputPortWidth(S, 0, 5);
	ssSetInputPortDirectFeedThrough(S, 0, 1);
	ssSetInputPortDataType(S, 0, SS_DOUBLE);
	ssSetInputPortComplexSignal(S, 0, COMPLEX_NO);
	ssSetInputPortRequiredContiguous(S, 0, 1);
	//In1
	ssSetInputPortWidth(S, 1, 5);
	ssSetInputPortDirectFeedThrough(S, 1, 1);
	ssSetInputPortDataType(S, 1, SS_DOUBLE);
	ssSetInputPortComplexSignal(S, 1, COMPLEX_NO);
	ssSetInputPortRequiredContiguous(S, 1, 1);
	//In2: flag
	ssSetInputPortWidth(S, 2, 1);
	ssSetInputPortDirectFeedThrough(S, 2, 1);
	ssSetInputPortDataType(S, 2, SS_DOUBLE);
	ssSetInputPortComplexSignal(S, 2, COMPLEX_NO);
	ssSetInputPortRequiredContiguous(S, 2, 1);
	//In3: time
	ssSetInputPortWidth(S, 3, 1);
	ssSetInputPortDirectFeedThrough(S, 3, 1);
	ssSetInputPortDataType(S, 3, SS_DOUBLE);
	ssSetInputPortComplexSignal(S, 3, COMPLEX_NO);
	ssSetInputPortRequiredContiguous(S, 3, 1);


	if (!ssSetNumOutputPorts(S, 1)) return;

	//Out0
	ssSetOutputPortWidth(S, 0, 5);
	ssSetOutputPortDataType(S, 0, SS_DOUBLE);
	ssSetOutputPortComplexSignal(S, 0, COMPLEX_NO);

	ssSetNumSampleTimes(S, 1);
	ssSetNumRWork(S, 0);
	ssSetNumIWork(S, 0);
	ssSetNumPWork(S, 1); ////////////////////
	ssSetNumModes(S, 0);
	ssSetNumNonsampledZCs(S, 0);

	ssSetOptions(S,
		SS_OPTION_WORKS_WITH_CODE_REUSE |
		SS_OPTION_EXCEPTION_FREE_CODE);
}


// Function: mdlInitializeSampleTimes =========================================
// Abstract:
//   This function is used to specify the sample time(s) for your
//   S-function. You must register the same number of sample times as
//   specified in ssSetNumSampleTimes.
static void mdlInitializeSampleTimes(SimStruct *S)
{
	ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
	ssSetOffsetTime(S, 0, 0.0);
	ssSetModelReferenceSampleTimeDefaultInheritance(S);
}

// Function: mdlStart =======================================================
// Abstract:
//   This function is called once at start of model execution. If you
//   have states that should be initialized once, this is the place
//   to do it.
#define MDL_START
static void mdlStart(SimStruct *S)
{
	ssGetPWork(S)[0] = (void *) new Planner;
}

// Function: mdlOutputs =======================================================
// Abstract:
//   In this function, you compute the outputs of your S-function
//   block.
static void mdlOutputs(SimStruct *S, int_T tid)
{
	Planner *planner = (Planner *)ssGetPWork(S)[0];   // retrieve C++ object from
	const real_T   *u0 = (const real_T*)ssGetInputPortSignal(S, 0);
	const real_T   *u1 = (const real_T*)ssGetInputPortSignal(S, 1);
	const real_T   *u2 = (const real_T*)ssGetInputPortSignal(S, 2); // flag
	const real_T   *u3 = (const real_T*)ssGetInputPortSignal(S, 3); // time
	real_T        *y0 = (real_T *)ssGetOutputPortRealSignal(S, 0);
	//
	double state_i[5] = { u0[0],u0[1],u0[2],u0[3], u0[4] };
	double state_g[5] = { u1[0],u1[1],u1[2],u1[3], u1[4] };
	double flag = u2[0];
	double time = u3[0];
	double state[5] = { -1.,-1.,-1.,-1.,-1. };
	if (flag && (!planner->start_planning) && planner->end_planning)
	{
		planner->update(state_i, state_g, time);
	}
	if (time > planner->end_time)
	{
		planner->reset();
	}
	if (planner->r[4] > 0. && planner->u[3] > 0.)
	{
		double delta_t = time - planner->start_time;
		double v = planner->u[0] + planner->u[1] * delta_t + planner->u[2] * delta_t*delta_t;
		planner->interp(delta_t, state);
	}
	y0[0] = state[0];
	y0[1] = state[1];
	y0[2] = state[2];
	y0[3] = state[3];
	y0[4] = state[4];
	// UNUSED_ARG(tid);                             // object
}

// Function: mdlTerminate =====================================================
// Abstract:
//   In this function, you should perform any actions that are necessary
//   at the termination of a simulation.  For example, if memory was
//   allocated in mdlStart, this is the place to free it.
static void mdlTerminate(SimStruct *S)
{
	Planner *planner = (Planner *)ssGetPWork(S)[0]; // retrieve and destroy C++
	delete planner;                                  // object in the termination
}


// Required S-function trailer
#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
