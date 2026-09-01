#include <stdio.h>
#include <math.h>

#define N 31 

int main() {
    double x[N], A[N];
    double dx = 3.0 / (N - 1); 
    
    double rho[N], T[N], V[N];
    double U1[N], U2[N], U3[N];
double F1[N], F2[N], F3[N], J2[N];
double F1_bar[N], F2_bar[N], F3_bar[N], J2_bar[N];
    double U1_bar[N], U2_bar[N], U3_bar[N];
    double dt = 0.0267; // Stable time step based on the Courant number 
    double gamma = 1.4;

    for (int i = 0; i < N; i++) {
        // Spatial and Area setup
        x[i] = i * dx;
        A[i] = 1.0 + 2.2 * pow((x[i] - 1.5), 2);
        
        // Piecewise initial conditions
        if (x[i] >= 0.0 && x[i] <= 0.5) {
            rho[i] = 1.0;
            T[i] = 1.0;
        } else if (x[i] > 0.5 && x[i] <= 1.5) {
            rho[i] = 1.0 - 0.366 * (x[i] - 0.5);
            T[i] = 1.0 - 0.167 * (x[i] - 0.5);
        } else {
            rho[i] = 0.634 - 0.3879 * (x[i] - 1.5);
            T[i] = 0.833 - 0.3507 * (x[i] - 1.5);
        }
        
        V[i] = 0.59 / (rho[i] * A[i]);
        
        // Conservation variables
        U1[i] = rho[i] * A[i];
        U2[i] = rho[i] * A[i] * V[i];
        U3[i] = U1[i] * (T[i] / (gamma - 1.0) + (gamma / 2.0) * pow(V[i], 2));
        
        
    }
// Time-marching loop (1400 steps allows the transient flow to reach a steady state)
    for (int k = 0; k < 1400; k++) {
        
        // 1. Calculate Fluxes (F) across the entire nozzle
        for (int i = 0; i < N; i++) {
            F1[i] = U2[i];
            F2[i] = (pow(U2[i], 2) / U1[i]) + ((gamma - 1.0) / gamma) * (U3[i] - (gamma / 2.0) * pow(U2[i], 2) / U1[i]);
            F3[i] = gamma * (U2[i] * U3[i] / U1[i]) - (gamma * (gamma - 1.0) / 2.0) * (pow(U2[i], 3) / pow(U1[i], 2));
        }

        // 2. PREDICTOR STEP (Forward Spatial Difference)
        for (int i = 0; i < N - 1; i++) {
            // Calculate the source term J2 (pressure gradient through the changing area)
            J2[i] = ((gamma - 1.0) / gamma) * (U3[i] - (gamma / 2.0) * pow(U2[i], 2) / U1[i]) * ((A[i+1] - A[i]) / dx) / A[i];
            
            // Predict the intermediate solution vectors (U_bar)
            U1_bar[i] = U1[i] - (dt / dx) * (F1[i+1] - F1[i]);
            U2_bar[i] = U2[i] - (dt / dx) * (F2[i+1] - F2[i]) + dt * J2[i];
            U3_bar[i] = U3[i] - (dt / dx) * (F3[i+1] - F3[i]);
        }
        
        // 3. Calculate Predicted Fluxes (F_bar) from U_bar
        for (int i = 0; i < N - 1; i++) {
            F1_bar[i] = U2_bar[i];
            F2_bar[i] = (pow(U2_bar[i], 2) / U1_bar[i]) + ((gamma - 1.0) / gamma) * (U3_bar[i] - (gamma / 2.0) * pow(U2_bar[i], 2) / U1_bar[i]);
            F3_bar[i] = gamma * (U2_bar[i] * U3_bar[i] / U1_bar[i]) - (gamma * (gamma - 1.0) / 2.0) * (pow(U2_bar[i], 3) / pow(U1_bar[i], 2));
        }

        // 4. CORRECTOR STEP (Backward Spatial Difference)
        for (int i = 1; i < N - 1; i++) {
            J2_bar[i] = ((gamma - 1.0) / gamma) * (U3_bar[i] - (gamma / 2.0) * pow(U2_bar[i], 2) / U1_bar[i]) * ((A[i] - A[i-1]) / dx) / A[i];
            
            U1[i] = 0.5 * (U1[i] + U1_bar[i] - (dt / dx) * (F1_bar[i] - F1_bar[i-1]));
            U2[i] = 0.5 * (U2[i] + U2_bar[i] - (dt / dx) * (F2_bar[i] - F2_bar[i-1]) + dt * J2_bar[i]);
            U3[i] = 0.5 * (U3[i] + U3_bar[i] - (dt / dx) * (F3_bar[i] - F3_bar[i-1]));
        }

        // 5. BOUNDARY CONDITIONS
        // Inlet: Density and Temperature fixed at stagnation, Momentum extrapolated
        U1[0] = A[0]; 
        U2[0] = 2.0 * U2[1] - U2[2];
        double V0 = U2[0] / U1[0];
        U3[0] = U1[0] * (1.0 / (gamma - 1.0) + (gamma / 2.0) * pow(V0, 2));

        // Outlet: All variables extrapolated for supersonic flow
        U1[N-1] = 2.0 * U1[N-2] - U1[N-3];
        U2[N-1] = 2.0 * U2[N-2] - U2[N-3];
        U3[N-1] = 2.0 * U3[N-2] - U3[N-3];
    }
// Decode final solution vectors back to primitive variables
    // Create and open a CSV file in write mode
    FILE *fp = fopen("nozzle_data.csv", "w");
    
    // Write the column headers
    fprintf(fp, "Node,x,Mach,Density\n");

    // Decode final solution vectors back to primitive variables
    for (int i = 0; i < N; i++) {
        rho[i] = U1[i] / A[i];
        V[i] = U2[i] / U1[i];
        T[i] = (gamma - 1.0) * ((U3[i] / U1[i]) - (gamma / 2.0) * pow(V[i], 2));
        double Mach = V[i] / sqrt(T[i]); 
        
        // Print to the terminal (so you can still see it)
        printf("Node %02d | x: %.2f | Mach: %.4f | Density: %.4f\n", i, x[i], Mach, rho[i]);
        
        // Save the exact same data to the CSV file, separated by commas
        fprintf(fp, "%d,%.2f,%.4f,%.4f\n", i, x[i], Mach, rho[i]);
    }
    
    // Close the file to save it safely
    fclose(fp);
    printf("\nData successfully exported to nozzle_data.csv\n");
    return 0;
}
