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
    for (int i = length - 1; i >= 0; --i) {
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

qcomplex_t kernal_lgh(double k, double gamma, double c)
{
    qcomplex_t temp;
    temp = qcomplex_t(exp(c - (k * k + 1) / (4 * gamma * gamma)) / (1 + k * k), 0.0) * exp(qcomplex_t(0.0, -1.0 * k * c)) / PI;
    return temp;
}

double get_norm1(qcomplex_t x)
{
    return sqrt(pow(x.real(), 2) + pow(x.imag(), 2));
}

double u_initial(double x, double y)
{
    double l1 = 1.0;
    double l2 = 1.0;
    double temp = cos(PI * x / (2.0 * l1)) * cos(PI * y / (2.0 * l2));
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

QCircuit exp_siama11_gate(QVec qvec, int n, double gamma, double tau)
{
    int m;
    QVec control_qubit;
    for(m = 0; m < n - 1; m++)
    {
        control_qubit.push_back(qvec[m]);
    }

    QCircuit G;
    G << U1(qvec[n - 1], gamma * tau).control(control_qubit);

    return G;
}

QCircuit exp_siama00_gate(QVec qvec, int n, double gamma, double tau)
{
    int m;
    QVec control_qubit;
    for(m = 0; m < n - 1; m++)
    {
        control_qubit.push_back(qvec[m]);
    }

    QCircuit G;
    for(m = 0; m < n; m++)
    {
        G << X(qvec[m]);
    }
    G << U1(qvec[n - 1], gamma * tau).control(control_qubit);
    for(m = 0; m < n; m++)
    {
        G << X(qvec[m]);
    }

    return G;
}

int main()
{
    auto start = system_clock::now();

    //Robin boundary condition 2d

    int d = 2;

    int i, j, k;
    double temp;

    int m = 4; // number of qubits for LCHS, m
    
    vector<int> n(d);
    n[0] = 7;
    n[1] = 7;

    int r = 1; // Trotter steps
    int temp_K = 1774;
    double K = 1.0 * temp_K; // integration truncation interval is [-K, K]
    K = 17.74;

    int n_total = 0;
    vector<int> N(d);
    int N_total = 1;
    for(i = 0; i < d; i++)
    {
        n_total += n[i];

        N[i] = getPowerOf2(n[i]);
        N_total *= N[i];
    }

    int M = getPowerOf2(m); // number of LCHS coefficients, M

    vector<double> a(d);
    a[0] = 1.0;
    a[1] = 1.0;

    double T = 0.04; // evolution time T

    vector<double> left(d); // left side of the computational domain
    vector<double> right(d); // right side of the computational domain

    left[0] = 0.0;
    right[0] = 1.0;
    left[1] = 0.0;
    right[1] = 1.0;

    vector<double> length(d);
    for(i = 0; i < d; i++)
    {
        length[i] = right[i] - left[i];
    }

    vector<double> dx(d);
    for(i = 0; i < d; i++)
    {
        dx[i] = length[i] / N[i];
    }

    vector<double> u0(N_total);
    double u0_norm2 = 0.0;

    // assign 2D data: the x direction is the lowest layer, then the y direction, and so on
    for(j = 0; j < N[1]; j++)
    {
        for(i = 0; i < N[0]; i++)
        {
            temp = u_initial(left[0] + i*dx[0] + 0.5*dx[0], left[1] + j*dx[1] + 0.5*dx[1]);
            u0[j*N[0] + i] = temp;
            u0_norm2 += temp * temp;
        }
    }
    
    u0_norm2 = sqrt(u0_norm2);

    for(i = 0; i < N_total; i++)
    {
        u0[i] = u0[i] / u0_norm2;
    }

    double dk = 2 * K / M;
    
    printf("dk = %.6f\n", dk);

    vector<double> temp_k(M);
    vector<qcomplex_t> integrate_coeff(M);
    vector<qcomplex_t> sqrt_integrate_coeff(M);
    vector<qcomplex_t> sqrt_integrate_coeff_dagger(M);
    double integrate_coeff_norm1 = 0.0;

    double K_minus_dk = K - 0.5 * dk;

    double eps = 0.001;
    double c = 0.4;
    double gamma = sqrt(c + log((1 + 1/(2*PI)) / eps)) / c;

    for(i = 0; i < M; i++)
    {
        integrate_coeff[i] = qcomplex_t(dk, 0.0) * kernal_lgh(i * dk - K_minus_dk, gamma, c);

        integrate_coeff_norm1 += get_norm1(integrate_coeff[i]);
    }

    printf("integrate_coeff_norm1 = %.6f\n", integrate_coeff_norm1);
    printf("u0_norm2 = %.6f\n", u0_norm2);
    printf("integrate_coeff_norm1 * u0_norm2 = %.6f\n", integrate_coeff_norm1 * u0_norm2);

    for(i = 0; i < M; i++)
    {
        sqrt_integrate_coeff[i] = pow(integrate_coeff[i] / qcomplex_t(integrate_coeff_norm1, 0.0), 0.5);
        sqrt_integrate_coeff_dagger[i] = get_dagger(sqrt_integrate_coeff[i]);

    }

    //0.Prepare the quantum machine and allocate qubits
    auto qvm = CPUQVM();
    qvm.init();
    auto qvec_ancilla = qvm.qAllocMany(m);

    // d dimensions, n[i] qubits per dimension
    vector<QVec> qvec1(d);
    for(i = 0; i < d; i++)
    {
        qvec1[i] = qvm.qAllocMany(n[i]);
    }

    // system qubits
    QVec qvec_system;
    for(k = 0; k < d; k++)
    {
        for(i = 0; i < n[k]; i++)
        {
            qvec_system.push_back(qvec1[k][i]);
        }
    }

    QProg prog;
    //1.Prepare the initial state
    prog << amplitude_encode(qvec_ancilla, sqrt_integrate_coeff);
    prog << amplitude_encode(qvec_system, u0);

    //2.Select oracle
    double tau = T / r;
    printf("tau = %.6f\n", tau);
    for(int index = 0; index < r; index++)
    {
        for(k = 0; k < d; k++) // 2nd-order Trotter
        {
            prog << W_gate(qvec1[k], 0, K_minus_dk * a[k] * a[k] / (dx[k] * dx[k]), 0.5 * tau, 0.0);
            for(j = 0; j < m; j++)
            {
                prog << W_gate(qvec1[k], 0, -1.0 * dk * a[k] * a[k] * getPowerOf2(j) / (dx[k] * dx[k]), 0.5 * tau, 0.0).control(qvec_ancilla[j]);
            }

            for(i = 1; i < n[k]; i++)
            {
                prog << W_gate(qvec1[k], i, K_minus_dk * a[k] * a[k] / (dx[k] * dx[k]), tau, 0.0);
            }
            prog << exp_siama11_gate(qvec1[k], n[k], 1.0 * K_minus_dk * a[k] * a[k] / (dx[k] * dx[k]), tau);
            prog << exp_siama00_gate(qvec1[k], n[k], -1.0 * K_minus_dk * a[k] * a[k] / (dx[k] * dx[k]), tau);

            for(j = 0; j < m; j++)
            {
                for(i = 1; i < n[k]; i++)
                {
                    prog << W_gate(qvec1[k], i, -1.0 * dk * a[k] * a[k] * getPowerOf2(j) / (dx[k] * dx[k]), tau, 0.0).control(qvec_ancilla[j]);
                }
                prog << exp_siama11_gate(qvec1[k], n[k], -1.0 * dk * a[k] * a[k] * getPowerOf2(j) / (dx[k] * dx[k]), tau).control(qvec_ancilla[j]);
                prog << exp_siama00_gate(qvec1[k], n[k], dk * a[k] * a[k] * getPowerOf2(j) / (dx[k] * dx[k]), tau).control(qvec_ancilla[j]);
            }

            prog << W_gate(qvec1[k], 0, K_minus_dk * a[k] * a[k] / (dx[k] * dx[k]), 0.5 * tau, 0.0);
            for(j = 0; j < m; j++)
            {
                prog << W_gate(qvec1[k], 0, -1.0 * dk * a[k] * a[k] * getPowerOf2(j) / (dx[k] * dx[k]), 0.5 * tau, 0.0).control(qvec_ancilla[j]);
            }
        }
    }
    
    for(k = 0; k < d; k++)
    {
        for(j = 0; j < m; j++)
        {
            prog << U1(qvec_ancilla[j], -2.0 * dk * T * a[k] * a[k] * getPowerOf2(j) / (dx[k] * dx[k]));
        }
        prog << U4(qvec_ancilla[0], 2.0 * K_minus_dk * T * a[k] * a[k] / (dx[k] * dx[k]), 0.0, 0.0, 0.0);
    }

    //3.Dagger of preparation oracle
    prog << amplitude_encode(qvec_ancilla, sqrt_integrate_coeff_dagger).dagger();

    //Run quantum program
    QVec qvec_all;
    for(k = 0; k < d; k++)
    {
        for(i = 0; i < n[k]; i++)
        {
            qvec_all.push_back(qvec1[k][i]);
        }
    }
    for(i = 0; i < m; i++)
    {
        qvec_all.push_back(qvec_ancilla[i]);
    }

    const size_t all_qubit_num = m + n_total;

    auto result = qvm.probRunDict(prog, qvec_all);

    temp = exp(-(a[0] / pow(length[0], 2) + a[1] / pow(length[1], 2)) * pow(PI, 2) * T / 4);

    double rate = 0.0;

    double xi;
    double yi;
    double ui;
    double lchs_ui;
    double sign;

    FILE* fp;
    string folder_path = "./8diff_homo_robin_central_2d/";

    string filename = "n0=";
    filename += to_string(n[0]);
    filename += "n1=";
    filename += to_string(n[1]);
    filename += "m=";
    filename += to_string(m);
    filename += "r=";
    filename += to_string(r);
    filename += "K=";
    filename += to_string(temp_K);
    filename +=".plt";

    string full_path = folder_path + filename;

    fp = fopen(full_path.c_str(), "w");
    fprintf(fp, "%s", filename.c_str());
    fprintf(fp, "\nVariables=x,      y,      u_exact,      u_LCHS,      u_error\n");
    fprintf(fp, "zone t=\"Frame 0\"i=%d,j=%d,f=point\n", N[0], N[1]);

    double error_infty = 0.0;
    double error_1 = 0.0;
    double error_2 = 0.0;
    double error_temp;

    double u_infty = 0.0;
    double u_1 = 0.0;
    double u_2 = 0.0;
    double u_temp;

    string str;

    for(j = 0; j < N[1]; j++)
    {
        for(i = 0; i < N[0]; i++)
        {
            xi = left[0] + i * dx[0] + 0.5 * dx[0];
            yi = left[1] + j * dx[1] + 0.5 * dx[1];
            ui = temp * u_initial(xi, yi);
            if (ui < 0) sign = -1.0;
            else sign = 1.0;

            str = decimalToBinaryWithPadding(all_qubit_num, j * N[0] + i);
            rate += result[str];

            lchs_ui = sqrt(result[str]) * integrate_coeff_norm1 * u0_norm2 * sign;

            u_temp = fabs(ui);
            u_infty = (u_infty > u_temp) ? u_infty : u_temp;
            u_1 += u_temp;
            u_2 += pow(u_temp, 2);

            error_temp = fabs(ui - lchs_ui);
            error_infty = (error_infty > error_temp) ? error_infty : error_temp;
            error_1 += error_temp;
            error_2 += pow(error_temp, 2);

            fprintf(fp, "%.10f  %.10f  %.10f  %.10f  %.10f\n", xi, yi, ui, lchs_ui, error_temp);
        }
    }

    fclose(fp);

    string filename1 = "n0=";
    filename1 += to_string(n[0]);
    filename1 += "n1=";
    filename1 += to_string(n[1]);
    filename1 += "m=";
    filename1 += to_string(m);
    filename1 += "r=";
    filename1 += to_string(r);
    filename1 += "K=";
    filename1 += to_string(temp_K);
    filename1 +="error.dat";

    string full_path1 = folder_path + filename1;
    fp = fopen(full_path1.c_str(), "w");

    double relative_error_infty = error_infty / u_infty;
    double relative_error_1 = error_1 / u_1;
    double relative_error_2 = sqrt(error_2 / u_2);

    error_1 *= dx[0] * dx[1];
    error_2 = sqrt(error_2 * dx[0] * dx[1]);

    printf("%s\n", filename1.c_str());
    fprintf(fp, "%s\n", filename1.c_str());

    fprintf(fp, "integrate_coeff_norm1 = %.6f\n", integrate_coeff_norm1);
    fprintf(fp, "u0_norm2 = %.6f\n", u0_norm2);
    fprintf(fp, "integrate_coeff_norm1 * u0_norm2 = %.6f\n", integrate_coeff_norm1 * u0_norm2);

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

    printf("energy = exp(-(a1^2/l1^2+a2^2/l2^2)*pi^2*t) = %.6f\n", temp);
    fprintf(fp, "energy = exp(-(a1^2/l1^2+a2^2/l2^2)*pi^2*t) = %.6f\n", temp);

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