

// Definition of function parameters
struct function_parameters
{
	double mu;
	double sigma;
};

// Definition for the integral-bounds and the stepsize
struct initial_data
{
	double stepsize;
	double initial;
	double final_val;
	double error;
};

// structure for function parameters
typedef struct function_parameters FunctionParams;

// structure for the integral-bounds and the stepsize
typedef struct initial_data InitialData;


// Definition for test calling
struct test_call
{
	char* function_name;
	InitialData init_data;
	FunctionParams func_params;
	double (*function)(double, FunctionParams);
};

// struct for test calls
typedef struct test_call TestCall;

// declaration
FunctionParams g_p_init(double x, double y);

// declaration
InitialData initialdata_init(double a, double b, double c, double d);

// decleration
TestCall test_call_init(char* function_name, InitialData init_data,
FunctionParams func_params, double (*function)(double, FunctionParams));


// Function declerations:

double gaussian(double x, FunctionParams params);

double strange_cos(double x, FunctionParams params);

double exp_minus_x_sq(double x, FunctionParams params);

double rev_sqrt(double x, FunctionParams params);
