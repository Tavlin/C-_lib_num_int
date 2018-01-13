#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "integrals.h"

// COLOURS
#define RED     "\x1b[31m"
#define GREEN   "\x1b[32m"
#define YELLOW  "\x1b[33m"
#define BLUE    "\x1b[34m"
#define MAGENTA "\x1b[35m"
#define CYAN    "\x1b[36m"
#define RESET   "\x1b[0m"



//planed:
//changing all given parameters to pointers. Gonna happen... someday... I guess
// printf function to give it a better look and fix some double posting

// if lower boundary > upper boundary
void swap (double x, double y)
{
	if (x != y)
	{
		double b;
		x = b;
		x = y;
		y = b;
	}
}

// function to print solutions of Integrals
void printerplus(TestCall test_call,const char * integral_name,
double (*integral)(InitialData A, FunctionParams params,
double (*func)(double, FunctionParams), double eps))
{
	double integral_val = (*integral)(test_call.init_data, test_call.func_params,
	test_call.function, test_call.init_data.error);
	
	if(integral_val >= 1000)
	{
		integral_val = INFINITY;
	}
	
	printf("\n*******************************************************************\
***********\n");
	printf("Integraltype:\t%s\nFunction:\t%s\nIntegralvalue = %+16.10lf\n",
	integral_name, test_call.function_name, integral_val);
	printf("lower bound \t= %+16.10lf\tupper bound  = %+16.10lf or infinite\n\
#steps \t\t\t\t= %+16lf\terror margin = %+16.10lf\n",
	test_call.init_data.initial, test_call.init_data.final_val,
	test_call.init_data.stepsize, test_call.init_data.error);
	printf("mu \t\t\t\t\t\t= %+16.10lf\tsigma\t\t\t\t = %+16.10lf\n",
	test_call.func_params.mu, test_call.func_params.sigma);
	printf("*********************************************************************\
*********\n");
}


//left Riemann sum
double left_riemann_sum(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
//initial instead of init, cuz __init__ is SPECIAL, at least somewhere else :)

{
	int swaped;
	// checks that initial value is allways smaller then final value 
	if(A.final_val < A.initial)
	{	
		swaped = 1;
		swap(A.final_val, A.initial);
	}
	
	// start value for the sum
	double left_sum = 0;
	
	double stepsize = (A.final_val-A.initial)/A.stepsize; 
	
	for(double i = 1; i <= A.stepsize; i++)
	{
		left_sum += stepsize *
		(*func)((A.initial + (((i-1)*(A.final_val-A.initial))/A.stepsize)), params);
	}
	
	if(swaped ==1)
	{
		return (-left_sum);
	}
	
	else
	{
		return left_sum;
	}
}


//right Rieman sum
double right_riemann_sum(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
{

	int swaped;
	// checks that initial value is allways smaller then final value 
	if(A.final_val < A.initial)
	{	
		swaped = 1;
		swap(A.final_val, A.initial);
	}
	
	double right_sum = 0;
	
	double stepsize = (A.final_val-A.initial)/A.stepsize;
	
	for(double i = 1; i <= A.stepsize; i++)
	{
		right_sum += stepsize *
		(*func)(((A.initial+(i*(A.final_val-A.initial))/A.stepsize)), params);
	}
	
	if(swaped ==1)
	{
		return (-right_sum);
	}
	
	else
	{
		return right_sum;
	}
}

//Trapezodial Rule
double trapezodial_integral(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
{

	int swaped;
	// checks that initial value is allways smaller then final value 
	if(A.final_val < A.initial)
	{
		swaped = 1;
		swap(A.final_val, A.initial);
	}
	
	// upper function value with pointy pointer
	double right_step = 0;
	double * pright_step;
	pright_step = &right_step;
	
	// lower function value with pointy pointer
	double left_step = 0;
	double * pleft_step;
	pleft_step = &left_step;
	
	double integral_val = 0;
	
	double stepsize = (A.final_val-A.initial)/A.stepsize;
	
	*pleft_step = (*func)((A.initial), params);
	
	for(double i = 1; i <= A.stepsize; i++)
	{
		*pright_step = (*func)((A.initial+(i*(A.final_val-A.initial))/A.stepsize), params);
		
		integral_val += stepsize/2.0 * 
		((*pleft_step) + (*pright_step));
		
		(*pleft_step) = (*pright_step);
	}
	
	if(swaped ==1)
	{
		return (-integral_val);
	}
	
	else
	{
		return integral_val;
	}
	
}

//Simpson's Rule
double simpson_integral(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
{

	int swaped;
	// checks that initial value is allways smaller then final value 
	if(A.final_val < A.initial)
	{
		swaped = 1;
		swap(A.final_val, A.initial);
	}
	
	// upper function value with pointy pointer
	double right_step = 0;
	double * pright_step;
	pright_step = &right_step;
	
	// lower function value with pointy pointer
	double left_step = 0;
	double * pleft_step;
	pleft_step = &left_step;
	
	double integral_val = 0;
	
	double stepsize = (A.final_val-A.initial)/A.stepsize;
	
	// initial lower function value
	*pleft_step = (*func)((A.initial), params);
	
	for(double i = 1; i <= A.stepsize; i++)
	{
		*pright_step = (*func)((A.initial+(i*(A.final_val-A.initial))/A.stepsize), params);
		
		integral_val += (stepsize/6.0) *
		((*pleft_step) + (*pright_step) + 
		(4.0*(*func)(((A.initial+(i*(A.final_val-A.initial))/A.stepsize) +
		(A.initial + (((i-1)*(A.final_val-A.initial))/A.stepsize)))/2.0,params)));
		
		(*pleft_step) = (*pright_step);
		
	}
	
	if(swaped ==1)
	{
		return (-integral_val);
	}
	
	else
	{
		return integral_val;
	}

}


//Trapezodial Rule semi adaptive stepsize
double trapezodial_integral_sas(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
{
	int swaped;
	// checks that initial value is allways smaller then final value 
	if(A.final_val < A.initial)
	{
		swaped = 1;
		swap(A.final_val, A.initial);
	}
	
	// reset # steps to 1
	A.stepsize = 1.0;
	
	double integral_val = 0;
	double* p_integral_val = &integral_val;
	
	double previous_int_val = 0;
	double* p_previous_int_val = &previous_int_val;
	
	double stepsize = 1.0/A.stepsize;
	double h = (A.final_val - A.initial);
	
	double new_int_val = 0;
	double* p_new_int_val = &new_int_val;
	
	// initial int value N = 1
	(*p_previous_int_val) = (h/2.0) *
	((*func)(A.initial,params) + (*func)(A.final_val,params));

	// so that N is up to date with N = 2 BEFORE entering the while loop!
	A.stepsize *= 2.0;
	h /= 2.0;

	// initial int value N = 2
	(*p_integral_val) = (*p_previous_int_val)/2.0 + 
	h * ((*func)((A.final_val + A.initial)/2.0, params));
	

	while(fabs((*p_previous_int_val)-(*p_integral_val)) > eps || A.stepsize < 4.0)
	{
		
		A.stepsize *= 2.0;
		h /= 2.0;
		stepsize = 1.0/A.stepsize;
	
		// 2*stepsize so the calculation is only done once per stepsize
		// otherwise it would need to calc for every for iteration
		double tt_stepsize = 2.0*stepsize;
	
		// new maximum value for the new midpoints:
		double max_val_midpoints = (A.stepsize-1.0)/A.stepsize;
	
		(*p_previous_int_val) = (*p_integral_val);
	
		*p_new_int_val = 0;
	
		for(double j = 1; j <= A.stepsize; j += 2)
		{
			(*p_new_int_val) += h * 
			(*func)((A.initial + j*h), params);
			
		}
	
		(*p_integral_val) = (*p_new_int_val) + (*p_integral_val)/2.0;

	}
	
	printf("Number of steps used in last calculation = %.0lf\n\n", A.stepsize);
	
	if(swaped ==1)
	{
		return (-integral_val);
	}
	
	else
	{
		return integral_val;
	}
	
	
}


// semi adaptive midpoint rule
double midpoint_int(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
{
	int swaped;
	// checks that initial value is allways smaller then final value 
	if(A.final_val < A.initial)
	{
		swaped = 1;
		swap(A.final_val, A.initial);
	}
	
	// reset # steps to 1
	A.stepsize = 1.0;
	
	double integral_val = 0;
	double* p_integral_val = &integral_val;
	
	double previous_int_val = 0;
	double* p_previous_int_val = &previous_int_val;
	
	double stepsize = 1.0/A.stepsize;
	double h = (A.final_val - A.initial);
	
	double midpoint = (A.final_val + A.initial)/2.0;
	
	double new_int_val = 0;
	double* p_new_int_val = &new_int_val;
	
	// initial int value N = 1
	(*p_previous_int_val) = h * (*func)(midpoint, params);
	
	// so that N is up to date with N = 3 BEFORE entering the while loop!
	A.stepsize *= 3.0;
	h /= 3.0;
	
	// initial int value N = 3
	(*p_integral_val) = (*p_previous_int_val)/3.0 + h * 
	(((*func)(midpoint - h, params)) + ((*func)(midpoint + h, params)));
	
	while(fabs((*p_previous_int_val)-(*p_integral_val)) > eps || A.stepsize < 9.0)
	{
	
		if(A.stepsize >= 19683)
		{
			if(*p_integral_val >= 1000)
			{
				return *p_integral_val;
				break;
			}
		}
		
		
		A.stepsize *= 3.0;
		h /= 3.0;
	
		(*p_previous_int_val) = (*p_integral_val);
	
		*p_new_int_val = 0;
	
		for(double j = -((A.stepsize/2.0) - 0.5); j <= ((A.stepsize/2.0) - 0.5); j += 3)
		{
			(*p_new_int_val) += h * (*func)(midpoint + (j*h), params);
			
		}
	
		for(double i = (-((A.stepsize/2.0) - 0.5)) + 2; i <= ((A.stepsize/2.0) - 0.5); i += 3)
		{
			(*p_new_int_val) += h * (*func)(midpoint + (i*h), params);
			
		}
	
		(*p_integral_val) = (*p_new_int_val) + (*p_integral_val)/3.0;

	}
	
	//printf("Number of steps used in last calculation = %.0lf\n\n", A.stepsize);
	
	if(swaped ==1)
	{
		return (-integral_val);
	}
	
	else
	{
		return integral_val;
	}


	
	
}

// semi adaptive midpoint rule with [a,inf)
double midpoint_int_to_inf(InitialData  A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
{
	
	// check if lower bound is < 0, so split is needed in [a,0.25] nad [0.25,inf)
	// 0.25 so the calculation for both integrals does not take too long/too many
	// steps
	
	double f_int;
	
	if(A.initial < 0)
	{
		A.final_val = 0.25;
		f_int = midpoint_int(A, params,(*func), eps);
		A.initial = 0;
		A.final_val = 4;
	}
	
	if(A.initial == 0)
	{
		A.final_val = 1;
		f_int = midpoint_int(A, params,(*func), eps);
		A.initial = 0;
		A.final_val = 1;
	
	}
	
	if(A.initial > 0)
	{
		A.final_val = 1/A.initial;
		A.initial = 0;
	}
	
	// reset # steps to 1
	A.stepsize = 1.0;
	
	double integral_val = 0;
	double* p_integral_val = &integral_val;
	
	double previous_int_val = 0;
	double* p_previous_int_val = &previous_int_val;
	
	double stepsize = 1.0/A.stepsize;
	double h = (A.final_val - A.initial);
	
	double midpoint = (A.final_val + A.initial)/2.0;
	
	double new_int_val = 0;
	double* p_new_int_val = &new_int_val;
	
	// initial int value N = 1
	(*p_previous_int_val) = h/pow(midpoint,2) * (*func)(1.0/midpoint, params);
	
	// so that N is up to date with N = 3 BEFORE entering the while loop!
	A.stepsize *= 3.0;
	h /= 3.0;
	
	// initial int value N = 3
	(*p_integral_val) = (*p_previous_int_val)/3.0 +  h/pow(midpoint-h,2)*
	((*func)(midpoint - h, params)) + h/pow(midpoint+h,2) *
	((*func)(midpoint + h, params));
	
	while(fabs((*p_previous_int_val)-(*p_integral_val)) > eps || A.stepsize < 9.0)
	{
		if(A.stepsize >= 19683)
		{
			if(*p_integral_val >= 1000)
			{
				return *p_integral_val;
				break;
			}
		}
		
		A.stepsize *= 3.0;
		h /= 3.0;
	
		(*p_previous_int_val) = (*p_integral_val);
	
		*p_new_int_val = 0;
	
		for(double j = -((A.stepsize/2.0) - 0.5); j <= ((A.stepsize/2.0) - 0.5); j += 3)
		{
			(*p_new_int_val) += h/pow((midpoint + (j*h)),2) * 
			(*func)(1/(midpoint + (j*h)), params);
			
		}
	
		for(double i = (-((A.stepsize/2.0) - 0.5)) + 2; i <= ((A.stepsize/2.0) - 0.5); i += 3)
		{
			(*p_new_int_val) += h/pow((midpoint + (i*h)),2) * 
			(*func)(1/(midpoint + (i*h)), params);
			
		}
	
		(*p_integral_val) = (*p_new_int_val) + (*p_integral_val)/3.0;

	}
	integral_val += f_int;
  //printf("Number of steps used in last calculation = %.0lf\n\n", A.stepsize);
	
	return integral_val;
	
}

double monte_carlo_integral(InitialData A, FunctionParams params,
double(*func)(double, FunctionParams), double eps)
{
	InitialData* pA = &A;
	int swaped;
	// checks that initial value is allways smaller then final value 
	if(A.final_val < A.initial)
	{		
	
		swaped = 1;
		swap(A.final_val, A.initial);
	}
	
	double integral_value = 0;
	double h = (A.final_val-A.initial);
	
	integral_value = h * (mean_value)((*func), A, params);
	
	double error = (statistical_error)((*func), A, params);

	while(error > eps)
	{
		(*pA).stepsize *= 10;
		integral_value = h * (mean_value)((*func), (*pA), params);
		error = (statistical_error)((*func), (*pA), params);
	}
	
	if(swaped ==1)
	{
		return (-integral_value);
	}
	
	else
	{
		return integral_value;
	}
}


// <Y>
double mean_value(double(*func)(double, FunctionParams), InitialData A,
FunctionParams params)
{
	//seed rand()
	srand((unsigned)time(NULL)); 
	
	double sum = 0;
	
	for(int i = 0; i < A.stepsize; i++)
	{
		int r = (rand()%((int)A.final_val-(int)A.initial)) + (int)A.initial;
		sum += (*func)(r, params);
		
	}
	return (sum/A.stepsize);
}

// <Y²>
double standard_deviation(double(*func)(double, FunctionParams), InitialData A, 
FunctionParams params)
{
	//seed rand()
	srand((unsigned)time(NULL)); 
	
	double sum = 0;
	
	for(int i = 0; i < A.stepsize; i++)
	{
		int r = (rand()%((int)A.final_val-(int)A.initial)) + (int)A.initial;
		sum += pow((*func)(r, params),2);
		
	}
	return (sum/A.stepsize);

}

// statistical error function
double statistical_error(double(*func)(double, FunctionParams), InitialData A, 
FunctionParams params)
{
	double stat_error = 0;
	double h = (A.final_val-A.initial);
	
	return (h * 
	sqrt(((standard_deviation)((*func), A, params) - 
	pow((mean_value)((*func), A, params),2))/A.stepsize));
	
}




