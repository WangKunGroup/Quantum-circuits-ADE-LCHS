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

qcomplex_t kernal_improve(double x, double beta)
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

double u_initial(double x)
{
    double temp = 0.0;
    return temp;
}

double u_exact(double x, double t)
{
    double a = 1.0;
    double l = 1.0;
    double temp = (sin(2 * PI * x) - sin(2 * PI * (x - a * t))) / (a * 2 * PI);
    return temp;
}

double right_b_f(double x)
{
    double l = 1.0;
    double temp = cos(2 * PI * x / l);
    return temp;
}

double right_b_g(double t)
{
    double temp = 1.0;
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
    for (m = 0; m < n_minus_1; m++)
    {
        control_qubit.push_back(qvec[m]);
    }

    QCircuit V;
    V << U.dagger();
    for (m = 0; m < n_minus_1; m++)
    {
        V << X(qvec[m]);
    }
    V << RZ(qvec[n_minus_1], 2.0 * gamma * tau).control(control_qubit);
    for (m = 0; m < n_minus_1; m++)
    {
        V << X(qvec[m]);
    }
    V << U;

    return V;
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

QCircuit Select_U_gate(QVec qvec_u, QVec qvec_ancilla2, int n, int m, int r, double delta_x, double delta_k, double K_minus_delta_k, double T, double a) // T = delta_T * 2^{k - 1}
{
    // 2-order TS decomposition
    int i, j;
    int index;
    QCircuit Sel;

    double tau = T / r;

    for(index = 0; index < r; index++)
    {
        // Circuit of L
        Sel << W_gate(qvec_u, 0, 0.5 * a * K_minus_delta_k / delta_x, 0.5 * tau, 0.0);
        for (j = 0; j < m; j++)
        {
            Sel << W_gate(qvec_u, 0, -0.5 * a * delta_k * getPowerOf2(j) / delta_x, 0.5 * tau, 0.0).control(qvec_ancilla2[j]);
        }

        for (i = 1; i < n; i++)
        {
            Sel << W_gate(qvec_u, i, 0.5 * a * K_minus_delta_k / delta_x, tau, 0.0);
        }
        for (j = 0; j < m; j++)
        {
            for (i = 1; i < n; i++)
            {
                Sel << W_gate(qvec_u, i, -0.5 * a * delta_k * getPowerOf2(j) / delta_x, tau, 0.0).control(qvec_ancilla2[j]);
            }
        }

        Sel << V_gate(qvec_u, n - 1, 0.5 * a * K_minus_delta_k / delta_x, tau, 0.0);
        for (j = 0; j < m; j++)
        {
            Sel << V_gate(qvec_u, n - 1, -0.5 * a * delta_k * getPowerOf2(j) / delta_x, tau, 0.0).control(qvec_ancilla2[j]);
        }

        Sel << W_gate(qvec_u, 0, 0.5 * a * K_minus_delta_k / delta_x, 0.5 * tau, 0.0);
        for (j = 0; j < m; j++)
        {
            Sel << W_gate(qvec_u, 0, -0.5 * a * delta_k * getPowerOf2(j) / delta_x, 0.5 * tau, 0.0).control(qvec_ancilla2[j]);
        }

        // Circuit of H
        Sel << W_gate(qvec_u, 0, -0.5 * a / delta_x, 0.5 * tau, 0.5 * PI);

        for (i = 1; i < n; i++)
        {
            Sel << W_gate(qvec_u, i, -0.5 * a / delta_x, tau, 0.5 * PI);
        }
        Sel << V_gate(qvec_u, n - 1, -0.5 * a / delta_x, tau, 0.5 * PI);

        Sel << W_gate(qvec_u, 0, -0.5 * a / delta_x, 0.5 * tau, 0.5 * PI);
    }

    for (j = 0; j < m; j++)
    {
        Sel << U1(qvec_ancilla2[j], -1.0 * a * T * delta_k * getPowerOf2(j) / delta_x);
    }
    Sel << U4(qvec_ancilla2[0], a * T * K_minus_delta_k / delta_x, 0.0, 0.0, 0.0);

    return Sel;
}

int main()
{
    auto start = system_clock::now();

    //Periodic boundary condition

    int i, j;

    int l = 5; // number of qubits for outer LCU, l
    int m = 4; // number of qubits for LCHS, m
    int n = 8; // number of qubits for spatial discretization points, n

    int L = getPowerOf2(l); // number of outer LCU coefficients, L
    int M = getPowerOf2(m); // number of LCHS coefficients, M
    int N = getPowerOf2(n); // number of spatial discretization points, N

    double T = 0.05; // evolution time T
    double delta_T = T / L; // time interval for the integration of the inhomogeneous term, delta_T

    int r = 8; // Trotter steps

    int temp_k = 8;
    double K = 1.0 * temp_k; // integration truncation interval is [-K, K]
    K = 8;

    double a = 1.0;
    double domain_left = 0.0; // computational domain
    double domain_right = 1.0;
    double domain_length = domain_right - domain_left;
    double delta_x = domain_length / N;

    // initialize the right-hand side b
    vector<double> b(N);
    double b_norm2 = 0.0;

    for(i = 0; i < N; i++)
    {
        b[i] = right_b_f(domain_left + i * delta_x + 0.5 * delta_x);
        b_norm2 += b[i] * b[i];
    }
    b_norm2 = sqrt(b_norm2);

    for(i = 0; i < N; i++)
    {
        b[i] = b[i] / b_norm2;
    }

    double delta_k = 2 * K / M;
    
    printf("delta_k = %.6f\n", delta_k);

    // LCHS coefficients
    vector<double> k(M);
    vector<qcomplex_t> integrate_coeff(M);
    vector<qcomplex_t> sqrt_integrate_coeff(M);
    vector<qcomplex_t> sqrt_integrate_coeff_dagger(M);
    double integrate_coeff_norm1 = 0.0;

    double K_minus_delta_k = K - 0.5 * delta_k;

    double eps = 0.001;
    double c = 0.4;
    double gamma = sqrt(c + log((1 + 1/(2*PI)) / eps)) / c;
    
    for(i = 0; i < M; i++)
    {
        k[i] = i * delta_k - K_minus_delta_k;
        integrate_coeff[i] = qcomplex_t(delta_k, 0.0) * kernal_lgh(k[i], gamma, c);

        integrate_coeff_norm1 += get_norm1(integrate_coeff[i]);
    }

    printf("integrate_coeff_norm1 = %.6f\n", integrate_coeff_norm1);
    printf("b_norm2 = %.6f\n", b_norm2);
    printf("integrate_coeff_norm1 * b_norm2 = %.6f\n", integrate_coeff_norm1 * b_norm2);

    for(i = 0; i < M; i++)
    {
        sqrt_integrate_coeff[i] = pow(integrate_coeff[i] / qcomplex_t(integrate_coeff_norm1, 0.0), 0.5);
        sqrt_integrate_coeff_dagger[i] = get_dagger(sqrt_integrate_coeff[i]);

        cout << "sqrt_integrate_coeff["<< i << "] = " << sqrt_integrate_coeff[i] << "  sqrt_integrate_coeff_dagger[" << i << "] = " << sqrt_integrate_coeff_dagger[i] << endl;
    }

    // outer LCU coefficients
    vector<double> combination_coeff(L);
    vector<double> sqrt_combination_coeff(L);
    double combination_coeff_norm1 = 0.0;

    for(i = 0; i < L; i++)
    {
        combination_coeff[i] = delta_T * right_b_g((L - i) * delta_T - 0.5 * delta_T);
        combination_coeff_norm1 += combination_coeff[i];

        printf("combination_coeff[%d] = %.6f\n", i, combination_coeff[i]);
    }

    printf("combination_coeff_norm1 = %.6f\n", combination_coeff_norm1);

    for(i = 0; i < L; i++)
    {
        sqrt_combination_coeff[i] = sqrt(combination_coeff[i] / combination_coeff_norm1);
        printf("sqrt_combination_coeff[%d] = %.6f\n", i, sqrt_combination_coeff[i]);
    }

    vector<double> u(N);

    QuantumMachine *machine = initQuantumMachine(CPU);
    auto qvec_ancilla1 = machine->qAllocMany(l);
    auto qvec_ancilla2 = machine->qAllocMany(m);
    auto qvec_u = machine->qAllocMany(n);

    QProg prog;
    // 1.Prepare the LCHS and outer LCU coefficients
    prog << amplitude_encode(qvec_ancilla1, sqrt_combination_coeff);
    prog << amplitude_encode(qvec_ancilla2, sqrt_integrate_coeff);

    // 2.Select-Preparation Oracle
    prog << amplitude_encode(qvec_u, b);

    // 3.Construct Select-U(l,j) Oracle

    // first layer: 0.5 * delta_T
    prog << Select_U_gate(qvec_u, qvec_ancilla2, n, m, r, delta_x, delta_k, K_minus_delta_k, 0.5 * delta_T, a);

    // select operator on the outer LCU coefficients
    for(i = 0; i < l; i++)
    {
        prog << Select_U_gate(qvec_u, qvec_ancilla2, n, m, r, delta_x, delta_k, K_minus_delta_k, delta_T * getPowerOf2(i), a).control(qvec_ancilla1[i]);
    }

    // 4.Dagger of preparation oracle
    prog << amplitude_encode(qvec_ancilla2, sqrt_integrate_coeff_dagger).dagger();
    prog << amplitude_encode(qvec_ancilla1, sqrt_combination_coeff).dagger();

    // 5.Run quantum program
    QVec qvec_all;
    for(i = 0; i < n; i++)
    {
        qvec_all.push_back(qvec_u[i]);
    }
    for(i = 0; i < m; i++)
    {
        qvec_all.push_back(qvec_ancilla2[i]);
    }
    for(i = 0; i < l; i++)
    {
        qvec_all.push_back(qvec_ancilla1[i]);
    }

    const size_t all_qubit_num = l + m + n;

    auto result = probRunDict(prog, qvec_all);

    double rate = 0.0;

    double xi;
    double ui;
    double lchs_ui;
    double sign;

    FILE* fp;
    string folder_path = "./6adv_inhomo_periodic_upwind/";

    string filename = "n=";
    filename += to_string(n);
    filename += "m=";
    filename += to_string(m);
    filename += "l=";
    filename += to_string(l);
    filename += "r=";
    filename += to_string(r);
    filename += "K=";
    filename += to_string(temp_k);
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
    for(i = 0; i < N; i++)
    {
        xi = domain_left + i * delta_x + 0.5 * delta_x;
        ui = u_exact(xi, T);
        if(ui < 0) sign = -1.0;
        else sign = 1.0;

        str = decimalToBinaryWithPadding(all_qubit_num, i);
        rate += result[str];

        lchs_ui = sqrt(result[str]) * combination_coeff_norm1 * integrate_coeff_norm1 * b_norm2 * sign;

        u_temp = fabs(ui);
        u_infty = (u_infty > u_temp) ? u_infty : u_temp;
        u_1 += u_temp;
        u_2 += pow(u_temp, 2);

        error_temp = fabs(ui - lchs_ui);
        error_infty = (error_infty > error_temp) ? error_infty : error_temp;
        error_1 += error_temp;
        error_2 += pow(error_temp, 2);

        printf("i = %d, x = %.6f, u_exact = %.6f, u_LCHS = %.6f, str = %s\n", i, xi, ui, lchs_ui, str.c_str());

        fprintf(fp, "%.10f  %.10f  %.10f  %.10f\n", xi, ui, lchs_ui, error_temp);
    }

    fclose(fp);

    string filename1 = "n=";
    filename1 += to_string(n);
    filename1 += "m=";
    filename1 += to_string(m);
    filename1 += "l=";
    filename1 += to_string(l);
    filename1 += "r=";
    filename1 += to_string(r);
    filename1 += "K=";
    filename1 += to_string(temp_k);
    filename1 +="error.dat";

    string full_path1 = folder_path + filename1;
    fp = fopen(full_path1.c_str(), "w");

    double relative_error_infty = error_infty / u_infty;
    double relative_error_1 = error_1 / u_1;
    double relative_error_2 = sqrt(error_2 / u_2);

    error_1 *= delta_x;
    error_2 = sqrt(error_2 * delta_x);

    printf("%s\n", filename.c_str());
    fprintf(fp, "%s\n", filename.c_str());

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