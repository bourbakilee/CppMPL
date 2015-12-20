#include <spiral3.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 1)//
	{
		mexErrMsgTxt("the number of inputs is must 1");
	}
	if (nlhs != 1)
	{
		mexErrMsgTxt("the number of outputs is must 1");
	}

	//----------------Inputs-------------------
	//double *qs, *qg;
	//qs = mxGetPr(prhs[0]); //q_s (x,y,theta,k)
	//qg = mxGetPr(prhs[1]); //q_g
	/*
	double *x0, *y0, *theta0, *k0, *x1, *y1, *theta1, *k1;
	x0 = mxGetPr(prhs[0]);
	y0 = mxGetPr(prhs[1]);
	theta0 = mxGetPr(prhs[2]);
	k0 = mxGetPr(prhs[3]);
	x1 = mxGetPr(prhs[0]);
	y1 = mxGetPr(prhs[1]);
	theta1 = mxGetPr(prhs[2]);
	k1 = mxGetPr(prhs[3]);
	*/
	double *q = mxGetPr(prhs[0]);
	double qs[] = { q[0],q[1],q[2],q[3] };
	double qg[] = { q[4],q[5],q[6],q[7] };
	//---------------End of Inputs--------------

	//---------------Open DataBase--------------
	sqlite3* db = nullptr;
	int rc = sqlite3_open("InitialGuessTable.db", &db);
	if (rc) {
		mexErrMsgTxt("Can't open database");
	}
	else {
		//-------------calculate path parameters-----
		plhs[0] = mxCreateDoubleMatrix(5, 1, mxREAL); //r=(a,b,c,d,sg)
		double *r = mxGetPr(plhs[0]);
		spiral3::spiral3(r, qs, qg, db);
		//------------Outputs-------------------------
	}
	sqlite3_close(db);
}