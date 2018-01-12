#include <math.h>
#include "initial_data_and_functions.h"

FunctionParams g_p_init(double x, double y)
/* initialise the gaussian parameters 
*/
{
	FunctionParams z;
	z.mu = x;
	z.sigma = y;
	return z;
}

InitialData initialdata_init(double a, double b, double c, double d)
{
	InitialData x;
	x.stepsize = a;
	x.initial = b;
	x.final_val = c;
	x.error = d;
	return x;
}

TestCall test_call_init(char* function_name, 
InitialData init_data, FunctionParams func_params,
double (*function)(double, FunctionParams))
{

	TestCall x;
	x.function_name = function_name;
	x.init_data = init_data;
	x.func_params = func_params;
	x.function = function;
	return x;
}


double gaussian(double x, FunctionParams params)
{
	return (1/(sqrt(2*M_PI)*params.sigma))*
	exp(-pow(x-params.mu, 2)/(2*pow(params.sigma, 2)));
}

double strange_cos(double x, FunctionParams params)
{
	return(x*pow(cos(2*M_PI*pow(x,2)),2));
}

double exp_minus_x_sq(double x, FunctionParams params)
{
	return (exp(-pow(x,2)));
}

double rev_sqrt(double x, FunctionParams params)
{
	return (1/sqrt(x));
}
