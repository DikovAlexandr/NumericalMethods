#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>

///////////////
// Constants //
///////////////

//  Condition constants
const double pho_l = 1.0; // Left boundary density
const double u_l = 0; // Left boundary gas velocity (along x)
const double p_l = 1.0; // Left boundary pressure

const double pho_r = 0.125; // Right boundary density
const double u_r = 0; // Right boundary gas velocity (along x)
const double p_r = 0.1; // Right boundary pressure

// Main part of programm
const double X = 1.0; // Coordinate boundary [meters]
const double T = 0.2; // Time boundary [seconds]

const double nx = 2000; // Number of steps on dimension coordinate
const double gamma = 1.4; // Polytropic index
const double cu = 0.99; // Courant number

const double dx = X / nx; // Coordinate step

/////////////////////////
// Condition functions //
/////////////////////////

// Function to set left boundary condition
std::vector<double> LeftBoundaryCondition(double t) {
    return {pho_l, pho_l * u_l, p_l / (gamma - 1) + 0.5 * pho_l * u_l * u_l};
}

// Function to set right boundary condition
std::vector<double> RightBoundaryCondition(double t) {
    return {pho_r, pho_r * u_r, p_r / (gamma - 1) + 0.5 * pho_r * u_r * u_r};
}

// Function to set initial condition
std::vector<double> InitialCondition(int i) {
    if ((i - 1) * dx < 0.5) {
        return LeftBoundaryCondition(0);
    } else {
        return RightBoundaryCondition(0);
    }
}

///////////////////////////
// Calculating functions //
///////////////////////////

// Function to calculate eigenvalues
std::vector<double> EigenValue(const std::vector<double>& Q) {
    double rho = Q[0];
    double u = Q[1] / rho;
    double etot = Q[2];
    double eint = etot - 0.5 * rho * u * u;
    double p = eint * (gamma - 1);
    double c = sqrt(gamma * p / rho);
    return {u, u + c, u - c};
}

// Function to calculate the flow values in Q
std::vector<double> f(const std::vector<double>& Q) {
    double rho = Q[0];
    double u = Q[1] / rho;
    double etot = Q[2];
    double eint = etot - 0.5 * rho * u * u;
    double p = eint * (gamma - 1);
    return {Q[1], rho * u * u + p, (etot + p) * u};
}

// Function to perform HLL method
std::vector<double> HLL(const std::vector<std::vector<double>>& Q, double dt) {
    std::vector<double> leigenvalues = EigenValue(Q[0]);
    std::vector<double> reigenvalues = EigenValue(Q[1]);
    double S_l = *std::min_element(leigenvalues.begin(), leigenvalues.end());
    double S_r = *std::max_element(reigenvalues.begin(), reigenvalues.end());
    if (S_r < 0) {
        return f(Q[1]);
    } else if (S_l > 0) {
        return f(Q[0]);
    } else {
        std::vector<double> Qhll(3, 0.0);
        std::vector<double> Fhll(3, 0.0);
        for (int i = 0; i < 3; ++i) {
            Qhll[i] = (f(Q[0])[i] - f(Q[1])[i] + Q[1][i] * S_r - Q[0][i] * S_l) / (S_r - S_l);
            Fhll[i] = f(Q[1])[i] + S_r * (Qhll[i] - Q[1][i]);
        }
        return Fhll;
    }
}

// Function to calculate the value of quantity in the cell to the next time layer
std::vector<std::vector<double>> CalculateCell(const std::vector<std::vector<double>>& Q, const std::function<std::vector<double>(const std::vector<std::vector<double>>&, double)>& F, double dt) {
    std::vector<std::vector<double>> result(Q.size(), std::vector<double>(Q[0].size(), 0.0));
    for (size_t i = 1; i < Q[0].size() - 1; ++i) {
        for (size_t j = 0; j < Q.size(); ++j) {
            result[j][i] = Q[j][i + 1] - dt / dx * (F(Q, dt)[j] - F(Q, dt)[j]);
        }
    }
    return result;
}

////////////
// Output //
////////////

// Function to calculate variables rho, u, total energy and p
std::vector<std::vector<double>> CalculateVariables(const std::vector<std::vector<double>>& Q) {
    std::vector<double> x(nx, 0.0);
    for (int i = 0; i < nx; ++i) {
        x[i] = dx * i;
    }

    std::vector<double> rho = Q[0];
    std::vector<double> u = Q[1] / rho;
    std::vector<double> etot = Q[2];
    std::vector<double> eint(rho.size());
    std::vector<double> p(rho.size());

    for (size_t i = 0; i < rho.size(); ++i) {
        eint[i] = etot[i] - 0.5 * rho[i] * u[i] * u[i];
        p[i] = eint[i] * (gamma - 1);
    }

    return {x, rho, u, etot, eint, p};
}

// Function to save variables to CSV
void SaveVariablesToCSV(const std::vector<std::vector<double>>& variables, const std::string& filename) {
    std::ofstream file(filename);

    if (file.is_open()) {
        for (size_t i = 0; i < variables[0].size(); ++i) {
            for (size_t j = 0; j < variables.size(); ++j) {
                file << variables[j][i];
                if (j < variables.size() - 1) {
                    file << ",";
                }
            }
            file << "\n";
        }

        file.close();
    } else {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
    }
}

int main() {
    // Creating a 2-dimensional calculation field
    std::vector<std::vector<double>> Q0 = InitialCondition(nx);

    double t = 0;
    while (t < T) {
        // Finding the time step
        double lambda = std::abs(Q0[1][0] / Q0[0][0]);
        for (int i = 0; i < nx; ++i) {
            std::vector<double> eigenvalues = EigenValue(Q0[i]);
            lambda = std::max({std::abs(eigenvalues[0]), std::abs(eigenvalues[1]), std::abs(eigenvalues[2]), lambda});
        }
        double dt = cu * dx / lambda;

        t += dt;

        // Creating the next time layer
        std::vector<std::vector<double>> Q1 = Q0;
        // Apply boundary conditions
        Q1[0].front() = LeftBoundaryCondition(t)[0];
        Q1[1].front() = LeftBoundaryCondition(t)[1];
        Q1[2].front() = LeftBoundaryCondition(t)[2];

        Q1[0].back() = RightBoundaryCondition(t)[0];
        Q1[1].back() = RightBoundaryCondition(t)[1];
        Q1[2].back() = RightBoundaryCondition(t)[2];

        for (int i = 2; i < nx - 1; ++i) {
            Q1[0][i] = CalculateCell(Q0[0][i - 1:i + 1], HLL, dt);
            Q1[1][i] = CalculateCell(Q0[1][i - 1:i + 1], HLL, dt);
            Q1[2][i] = CalculateCell(Q0[2][i - 1:i + 1], HLL, dt);
        }

        std::vector<std::vector<double>> variables = CalculateVariables(Q1);

        // Save variables to CSV
        SaveVariablesToCSV(variables, "rho.csv"); // Save rho values to CSV
        SaveVariablesToCSV(variables, "u.csv");   // Save u values to CSV
        SaveVariablesToCSV(variables, "etot.csv"); // Save etot values to CSV
        SaveVariablesToCSV(variables, "p.csv");    // Save p values to CSV

        Q0 = Q1;
    }

    // Print or use the final results as needed
    return 0;
}