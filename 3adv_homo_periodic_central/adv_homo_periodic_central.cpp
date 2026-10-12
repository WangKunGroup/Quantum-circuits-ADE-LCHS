#include "QPanda.h"
#include <math.h>
#include <complex>
#include <string>
#include <chrono>

using namespace QPanda;
using namespace std;
using namespace chrono;

qcomplex_t get_dagger(qcomplex_t x)
{
    return qcomplex_t(x.real(), -x.imag());
}

int getPowerOf2(int n)
{
    return 1 << n;
}

string decimalToBinaryWithPadding(int length, int num)
{
    string binary = "";
    for(int i = length - 1; i >= 0; --i)
    {
        binary += ((num >> i) & 1) ? '1' : '0';
    }
    return binary;
}

double f(double x)
{
    return 1 / (PI * (1 + x * x));
}

qcomplex_t kernal(double x, double beta)
{
    qcomplex_t temp;
    temp = qcomplex_t(2.0 * PI, 0.0) * exp(qcomplex_t(-1.0 * pow(2, beta), 0)) * exp(pow(qcomplex_t(1, x), beta));
    return qcomplex_t(1.0, 0) / temp;
}

double get_norm1(qcomplex_t x)
{
    return sqrt(pow(x.real(), 2) + pow(x.imag(), 2));
}

double u_initial(double x)
{
    double temp = sin(2.0 * PI * x);
    return temp;
}

QCircuit U_gate(QVec qvec, int j, double lambda)
{
    int m;
    QCircuit U;
    U << H(qvec[j]);
    U << U1(qvec[j], lambda);
    for(m = 0; m < j; m++)
    {
        U << CNOT(qvec[j], qvec[m]);
    }
    return U;
}

QCircuit W_gate(QVec qvec, int j, double gamma, double tau, double lambda)
{
    int m;
    QCircuit U = U_gate(qvec, j, -1.0 * lambda);
    QVec control_qubit;
    for(m = 0; m < j; m++)
    {
        control_qubit.push_back(qvec[m]);
    }
    
    QCircuit W_j;
    W_j << U.dagger();
    W_j << RZ(qvec[j], 2.0 * gamma * tau).control(control_qubit);
    W_j << U;

    return W_j;
}

QCircuit V_gate(QVec qvec, int n_minus_1, double gamma, double tau, double lambda)
{
    int m;
    QCircuit U = U_gate(qvec, n_minus_1, lambda);
    QVec control_qubit;
    for(m = 0; m < n_minus_1; m++)
    {
        control_qubit.push_back(qvec[m]);
    }
    
    QCircuit V;
    V << U.dagger();
    for(m = 0; m < n_minus_1; m++)
    {
        V << X(qvec[m]);
    }
    V << RZ(qvec[n_minus_1], 2.0 * gamma * tau).control(control_qubit);
    for(m = 0; m < n_minus_1; m++)
    {
        V << X(qvec[m]);
    }
    V << U;

    return V;
}

int main()
{
    auto start = system_clock::now();
    //Periodic boundary condition

    int i, j;

    int variable_qubit_num = 9; // number of qubits for spatial discretization points, n
    int r = 8; // Trotter steps

    int variable_point_num = getPowerOf2(variable_qubit_num); // number of spatial discretization points, N

    double a = 1.0; // advection coefficient a

    double T = 0.01; // evolution time T
    double domain_left = 0.0; // computational domain
    double domain_right = 1.0;
    double l = domain_right - domain_left;
    double delta_x = l / variable_point_num;

    vector<double> u0(variable_point_num);
    double u0_norm2 = 0.0;
    for(i = 0; i < variable_point_num; i++)
    {
        u0[i] = u_initial(domain_left + i * delta_x + 0.5 * delta_x);
        u0_norm2 += u0[i] * u0[i];
    }
    u0_norm2 = sqrt(u0_norm2);
    for(i = 0; i < variable_point_num; i++)
    {
        u0[i] = u0[i] / u0_norm2;
    }

    printf("u0_norm2 = %.6f\n", u0_norm2);

    vector<double> u(variable_point_num);

    //0.Apply for the quantum machine and qubits
    QuantumMachine *machine = initQuantumMachine(CPU);
    auto qvec_u = machine->qAllocMany(variable_qubit_num);

    QProg prog;
    //1.Prepare the initial state

    prog << amplitude_encode(qvec_u, u0);

    //2.Select oracle
    double tau = T / r;
    printf("tau = %.6f\n", tau);
    for(int index = 0; index < r; index++)
    {
        prog << W_gate(qvec_u, 0, -0.5 * a / delta_x, 0.5 * tau, 0.5 * PI);

        for(i = 1; i < variable_qubit_num; i++)
        {
            prog << W_gate(qvec_u, i, -0.5 * a / delta_x, tau, 0.5 * PI);
        }

        prog << V_gate(qvec_u, variable_qubit_num - 1, -0.5 * a / delta_x, tau, 0.5 * PI);

        prog << W_gate(qvec_u, 0, -0.5 * a / delta_x, 0.5 * tau, 0.5 * PI);
    }

    //3.Run quantum program and get the result
    QVec qvec_all;
    for(i = 0; i < variable_qubit_num; i++)
    {
        qvec_all.push_back(qvec_u[i]);
    }

    const size_t all_qubit_num = variable_qubit_num;

    auto result = probRunDict(prog, qvec_all);

    double rate = 0.0;

    double xi;
    double ui;
    double lchs_ui;
    double sign;

    FILE* fp;
    string folder_path = "./3adv_homo_periodic_central/";

    string filename = "n=";
    filename += to_string(variable_qubit_num);
    filename += "r=";
    filename += to_string(r);
    filename +=".plt";

    string full_path = folder_path + filename;

    fp = fopen(full_path.c_str(), "w");
    fprintf(fp, "%s", filename.c_str());
    fprintf(fp, "\nVariables=x,      u_exact,      u_LCHS,      ,u_error\n");

    double error_infty = 0.0;
    double error_1 = 0.0;
    double error_2 = 0.0;
    double error_temp;

    double u_infty = 0.0;
    double u_1 = 0.0;
    double u_2 = 0.0;
    double u_temp;

    string str;
    for(i = 0; i < variable_point_num; i++)
    {

        xi = domain_left + i * delta_x + 0.5 * delta_x;
        ui = u_initial(xi - a * T);
        if(ui < 0.0) sign = -1.0;
        else sign = 1.0;

        str = decimalToBinaryWithPadding(all_qubit_num, i);
        rate += result[str];

        lchs_ui = sqrt(result[str]) * u0_norm2 * sign;

        u_temp = fabs(ui);
        u_infty = (u_infty > u_temp) ? u_infty : u_temp;
        u_1 += u_temp;
        u_2 += pow(u_temp, 2);

        error_temp = fabs(ui - lchs_ui);
        error_infty = (error_infty > error_temp) ? error_infty : error_temp;
        error_1 += error_temp;
        error_2 += pow(error_temp, 2);

        printf("i = %d, x = %.6f, u_0 = %.6f, u_exact = %.6f, u_LCHS = %.6f, str = %s\n", i, xi, u0[i] * u0_norm2, ui, lchs_ui, str.c_str());

        fprintf(fp, "%.10f  %.10f  %.10f  %.10f\n", xi, ui, lchs_ui, error_temp);
    }

    fclose(fp);

    string filename1 = "n=";
    filename1 += to_string(variable_qubit_num);
    filename1 += "r=";
    filename1 += to_string(r);
    filename1 +="error.dat";

    string full_path1 = folder_path + filename1;
    fp = fopen(full_path1.c_str(), "w");

    double relative_error_infty = error_infty / u_infty;
    double relative_error_1 = error_1 / u_1;
    double relative_error_2 = sqrt(error_2 / u_2);

    error_1 *= delta_x;
    error_2 = sqrt(error_2 * delta_x);

    printf("%s\n", full_path1.c_str());
    fprintf(fp, "%s\n", full_path1.c_str());

    printf("relative_error_infty = %.12f\n", relative_error_infty);
    fprintf(fp, "relative_error_infty = %.12f\n", relative_error_infty);

    printf("relative_error_1 = %.12f\n", relative_error_1);
    fprintf(fp, "relative_error_1 = %.12f\n", relative_error_1);

    printf("relative_error_2 = %.12f\n", relative_error_2);
    fprintf(fp, "relative_error_2 = %.12f\n", relative_error_2);

    printf("error_infty = %.12f\n", error_infty);
    fprintf(fp, "error_infty = %.12f\n", error_infty);

    printf("error_1 = %.12f\n", error_1);
    fprintf(fp, "error_1 = %.12f\n", error_1);

    printf("error_2 = %.12f\n", error_2);
    fprintf(fp, "error_2 = %.12f\n", error_2);

    printf("rate = %.6f\n", rate);
    fprintf(fp, "rate = %.6f\n", rate);

    printf("sqrt(rate) = %.6f\n", sqrt(rate));
    fprintf(fp, "sqrt(rate) = %.6f\n", sqrt(rate));

    auto end = system_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "Total time: " << double(duration.count()) * microseconds::period::num / microseconds::period::den << "s" << endl;

    fprintf(fp, "Total time: %.6fs\n", double(duration.count()) * microseconds::period::num / microseconds::period::den);

    fclose(fp);

    return 0;
}