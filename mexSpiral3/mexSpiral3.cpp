#include <spiral3.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 2)//
	{
		mexErrMsgTxt("the number of inputs is must 2");
	}
	if (nlhs != 1)
	{
		mexErrMsgTxt("the number of outputs is must 1");
	}

	//----------------Inputs-------------------
	double *qs, *qg;
	qs = mxGetPr(prhs[0]); //q_s (x,y,theta,k)
	qg = mxGetPr(prhs[1]); //q_g
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