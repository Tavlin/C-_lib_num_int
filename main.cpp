#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "integrals.h"




//definition of size_of_list
template <class T, const unsigned int size>
const unsigned int size_of_array(T(&)[size])
{
	return size;
}

int main (int argc, char *argv[])
{

	char* function_name_list[4] = 
	{"Gaussian", "x*cos²(2*pi*x²)", "exp^(-x²)", "1/sqrt(x)"};
		
	char* integral_name_list[7] = 
	{"Left Riemann Sum", "Right Riemann Sum", "Trapezodial Rule",
	"Simpson's Rule", "Trapezodial Rule with semi adaptive stepsize",
	"Midpoint Rule with semiadatipve stepsize",
	"Midpoint Rule with semiadatipve stepsize and open bondary"};
	
	double (*function_list[4])(double x, FunctionParams params) =
	{gaussian, strange_cos, exp_minus_x_sq, rev_sqrt};
	
	double (* integral_list[7])(InitialData A, FunctionParams params,
	double(*func)(double, FunctionParams), double eps) =
	{left_riemann_sum, right_riemann_sum, trapezodial_integral, simpson_integral,
	trapezodial_integral_sas, midpoint_int, midpoint_int_to_inf};
	
	
	double mus[1] = {0};
	double sigmas[1] = {1};
	double Ns[2] = {1000,100000};
	double lower_bounds[4] = {-100,-10,-2,10};
	double upper_bounds[4] = {0,2,50,1000};
	double error = 0.0001;
	
	//initialise a list of initial data
	unsigned int init_data_list_lenght = (size_of_array(Ns) * 
	size_of_array(lower_bounds));
	InitialData init_data_list[init_data_list_lenght];
	
	int num_1 = 0;
	for(int i = 0; i < size_of_array(Ns); i++)
	{
		for(int j = 0; j < size_of_array(lower_bounds); j++)
		{
			
			InitialData A;
			InitialData* pA;
			pA = &A;
			*pA = initialdata_init(Ns[i], lower_bounds[j], upper_bounds[j], error);
			init_data_list[num_1] = *pA;
			num_1++;			
		}
		
	}
	
	//initialise a list of function parameters
	//expecting that leght of different parameterlists is equal for all parameters
	FunctionParams function_params_list[size_of_array(mus)];
	
	for(int k = 0; k < size_of_array(mus); k++)
	{
		FunctionParams A;
		FunctionParams* pA;
		pA = &A;
		*pA = g_p_init(mus[k], sigmas[k]);
		function_params_list[k] = *pA;
	}
	
	unsigned int test_call_list_lenght =
	size_of_array(function_name_list) * init_data_list_lenght * size_of_array(mus);
	TestCall test_call_list[test_call_list_lenght];
	
	// big loop to generate a big variety of possible test calls
	TestCall TC;
	TestCall* pTC;
	pTC = &TC;
	int num_2 = 0;
	for(int l = 0; l < size_of_array(function_name_list); l++)
	{
		for(int m = 0; m < init_data_list_lenght; m++)
		{
			for(int n = 0; n < size_of_array(mus); n++)
			{
				

				*pTC = test_call_init(function_name_list[l],
				init_data_list[m], function_params_list[n], function_list[l]);
				test_call_list[num_2] = *pTC;
				num_2++;
				
			}
		}
	}
	
	
	int check = atoi(argv[1]);
	//printf("Enter mode (1 for [a,b] or 2 for [a,inf)): ");
	//scanf("%i", &check);
		
	if(check == 1)
	{
		for(int o = 0; o < test_call_list_lenght; o++)
		{
			for(int p = 0; p < (size_of_array(integral_list) - 1); p++)
			{
				printerplus(test_call_list[o], integral_name_list[p], integral_list[p]);
			}
		}
	}
	
	else if(check == 2)
	{		
		for(int o = 0; o < test_call_list_lenght; o++)
		{
			for(int p = 6; p < (size_of_array(integral_list)); p++)
			{
				printerplus(test_call_list[o], integral_name_list[p], integral_list[p]);
			}
		}
	}
	return 0;
}
