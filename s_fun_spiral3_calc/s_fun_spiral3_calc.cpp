/* 2015.12.20, LI Yunsheng */
/*
±‡“Î÷∏¡Ó£∫
mex -I"D:\Dev\boost_1_60_0" -I"D:\ProDocs\CppMPL\CppMPL\libsqlite3" -I"D:\Dev\eigen3" -I"D:\ProDocs\CppMPL\CppMPL\TrajectoryGeneration" -L"D:\ProDocs\CppMPL\CppMPL\x64\Release" -lsqlite3 -lTrajectoryGeneration s_fun_spiral3_calc.cpp

*/
// *******************************************************************
// **** To build this mex function use: mex s_fun_spiral3_calc.cpp ****
// *******************************************************************

#include <spiral3.h>

#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME  s_fun_spiral3_calc

// Need to include simstruc.h for the definition of the SimStruct and
// its associated macro definitions.
#include <simstruc.h>

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
	if (!ssSetNumInputPorts(S, 2)) return;
	//In0
	ssSetInputPortWidth(S, 0, 4);
	ssSetInputPortDirectFeedThrough(S, 0, 1);
	ssSetInputPortDataType(S, 0, SS_DOUBLE);
	ssSetInputPortComplexSignal(S, 0, COMPLEX_NO);
	ssSetInputPortRequiredContiguous(S, 0, 1);
	//In1
	ssSetInputPortWidth(S, 1, 4);
	ssSetInputPortDirectFeedThrough(S, 1, 1);
	ssSetInputPortDataType(S, 1, SS_DOUBLE);
	ssSetInputPortComplexSignal(S, 1, COMPLEX_NO);
	ssSetInputPortRequiredContiguous(S, 1, 1);

	//Out0
	if (!ssSetNumOutputPorts(S, 1)) return;
	ssSetOutputPortWidth(S, 0, 5);
	ssSetOutputPortDataType(S, 0, SS_DOUBLE);
	ssSetOutputPortComplexSignal(S, 0, COMPLEX_NO);
	ssSetNumSampleTimes(S, 1);
	ssSetNumRWork(S, 0);
	ssSetNumIWork(S, 0);
	ssSetNumPWork(S, 0);
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

}

// Function: mdlOutputs =======================================================
// Abstract:
//   In this function, you compute the outputs of your S-function
//   block.
static void mdlOutputs(SimStruct *S, int_T tid)
{
	sqlite3* db = nullptr;
	int rc = sqlite3_open("InitialGuessTable.db", &db);

	// Get data addresses of I/O
	const real_T   *u0 = (const real_T*)ssGetInputPortSignal(S, 0);
	const real_T   *u1 = (const real_T*)ssGetInputPortSignal(S, 1);
	real_T        *y0 = (real_T *)ssGetOutputPortRealSignal(S, 0);
	//
	double q0[4] = { u0[0],u0[1],u0[2],u0[3] };
	double q1[4] = { u1[0],u1[1],u1[2],u1[3] };
	double r[5];
	spiral3::spiral3(r, q0, q1, db);
	y0[0] = r[0];
	y0[1] = r[1];
	y0[2] = r[2];
	y0[3] = r[3];
	y0[4] = r[4];
	sqlite3_close(db);
}

// Function: mdlTerminate =====================================================
// Abstract:
//   In this function, you should perform any actions that are necessary
//   at the termination of a simulation.  For example, if memory was
//   allocated in mdlStart, this is the place to free it.
static void mdlTerminate(SimStruct *S)
{

}


// Required S-function trailer
#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
