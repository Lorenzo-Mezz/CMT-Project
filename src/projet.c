#include <stdio.h>
#include <math.h> 
#include <stdlib.h>
#include <string.h>

// Definitions and constants
#define MAX_LINE_LENGTH 90
#define MAX_ROWS 301

// Structure to save the data of each line of the Csv file
typedef struct {
    char time[11];
    float surface_temp; // surface temperature (C)
    float atm_co2;      // atmospheric CO2 (ppm)
    float bottom_temp;    // bottom temperature(C)
    float surface_biomass;  // biomass concentration at depth 0 (mg/m3)
} DataRow;

// --- Function for CSV lecture ---

/**
 * Read a csv file that contain temporal data 
 * @param co2_temp_for_c Name of csv file
 * @param data Table of structures DataRow.
 * @param max_rows Maximal numer of lines to read.
 * @return Number of read lines or -1 if there is an error.
 */
int read_csv(const char *co2_temp_for_c, DataRow data[], int max_rows) {
    FILE *file = fopen(co2_temp_for_c, "r");
    if (file == NULL) {
        printf("Error when opening the file '%s'.\n", co2_temp_for_c );
        return -1;
    }

    char line[MAX_LINE_LENGTH];
    int row_count = 0;

    while (fgets(line, sizeof(line), file) != NULL) {
        if (row_count >= max_rows) {
            break; 
        }

        // Ignore the first line
        if (row_count == 0) {
            row_count++;
            continue;
        }

        if (sscanf(line, "%11s,%f,%f,%f,%f", 
                   data[row_count].time, 
                   &data[row_count].surface_temp, 
                   &data[row_count].bottom_temp,
                   &data[row_count].atm_co2, 
                   &data[row_count].surface_biomass) == 5) {
             row_count++;
        }
    }

    fclose(file);
    return row_count;
}

// Constants for the biomass and reaction model (global variables)
double R = 8.314;           // Constant for ideal gases (J/(mol*K))
double A = 0.0077;          // Constant for biomass growth (day^-1)
double E = 2537.0;          // Activation energy (cal/mol)
double I0 = 1084.93;        // Surface Irradiance (cal/(cm2*day))
double k = 0.547;           // Coefficient of light attenuation (m-1)
const int n_depth = 6;                //number of depth we're interested in
const int n_csv = 301;                //number of lines in our csv
// Constant for thermal diffusion
double alpha = 6.6 * pow(10, -5); // Thermal diffusivity (m2/s)

// --- Modelisation Functions ---

//Solves the 1D heat diffusion equation (explicit Euler scheme). Updates the T table with the final temperature profile after n_steps.

void diffusion_1d(double* T, double alpha, double dt, double dz, int n_points, int n_steps) {
    double* T_new = (double*)malloc(n_points * sizeof(double));
    if (T_new == NULL) {
        printf("Memory allocation error for T_new.\n");
        return; 
    }

    // Calculation of the stability factor (Must be r <= 0.5 for stability)
    double r = alpha * dt / (dz * dz); 

    printf("--- Thermal Diffusion Simulation (r=%.6f) ---\n", r);
    printf("Step 0 (Initial): ");
    for (int i = 0; i < n_points; i++) {
        printf("%.2f ", T[i]);
    }
    printf("\n");
    if (r > 0.5) {
        printf("WARNING: The stability factor r (%.6f) is > 0.5. UNSTABLE results.\n", r);
    }
    
    // Simulation loop
    for (int t = 0; t < n_steps; t++) {
        
        // Calculation of the new profile for interior points
        for (int i = 1; i < n_points - 1; i++) {
            T_new[i] = T[i] + r * (T[i + 1] - 2 * T[i] + T[i - 1]);
        }
        
        // Dirichlet boundary conditions (fixed temperatures T[0] and T[n_points-1])
        T_new[0] = T[0]; 
        T_new[n_points - 1] = T[n_points - 1]; 

        // Update T for the next iteration
        for (int i = 0; i < n_points; i++) {
            T[i] = T_new[i];
        }

        // Regular display
        if ((t + 1) % 1000 == 0) { 
            printf("Step %d: ", t + 1);
            for (int i = 0; i < n_points; i++) {
                printf("%.2f ", T[i]);
            }
            printf("\n");
        }
    }
    
    // Displaying the final status
    printf("Step %d (Final): ", n_steps);
    for (int i = 0; i < n_points; i++) {
        printf("%.2f ", T[i]);
    }
    printf("\n");
    
    printf("--- Diffusion simulation completed ---\n");
    free(T_new); 
}
int extract_columns(DataRow data[], 
                    int count, 
                    double** surface_temp_array, 
                    double** bottom_temp_array,
                    double** atm_co2_array,  
                    double** surface_biomass_array) {

    *surface_temp_array = (double*)malloc(count * sizeof(double));
    *bottom_temp_array = (double*)malloc(count * sizeof(double));
    *atm_co2_array = (double*)malloc(count * sizeof(double));
    *surface_biomass_array = (double*)malloc(count * sizeof(double));
    for (int i = 0; i < count; i++) {
        int data_index = i + 1; 
        
        (*surface_temp_array)[i] = data[data_index].surface_temp;
        (*bottom_temp_array)[i] = data[data_index].bottom_temp;
        (*atm_co2_array)[i] = data[data_index].atm_co2;
        (*surface_biomass_array)[i] = data[data_index].surface_biomass;
    }

    return 0;
}
double irradiance(int z_index, double dz) {
    double z_meters = z_index * dz;
    double Iz = I0 * exp(-k * z_meters);
    return Iz;
}

//Convert Celsius into Kelvin.
double temp(double t_celsius) {
    double T_kelvin = 273.15 + t_celsius;
    return T_kelvin;
}

//Calculate the Henry constant Kh (mg/L) as function of T.
double Kh(double temperature_z) {
    double kh = 1450 * exp(-2400 * (1/temp(temperature_z) - 1/298)); 
    return kh; // mg/L 
}

//Converts atmospheric CO2 (ppm) to dissolved CO2 (mg/L).
double convertCO2(double atm_co2, double temperature_z) { 
    double co2_atm = atm_co2 * pow(10, -6); // conversion from ppm to atm
    double diss_CO2 = Kh(temperature_z) * co2_atm;
    return diss_CO2; // mg/L
}

//Calculates the half-saturation constant Ks (mg/L) as a function of T (Celsius).

double Ks(double temperature_z) {
    double T_kelvin = temp(temperature_z);
    double ks = -0.2067 * pow(T_kelvin, 2) + 125.6 * T_kelvin -19029; 
    return ks; // mg/L 
}

void biomass_growth(double* T_array, double ppm_CO2, int *z_index, double dz, double* growth_rate_output){
    double* temperature_z = (double*)malloc(n_depth * sizeof(double));
    double *iz_val = (double*)malloc(n_depth * sizeof(double));
    double *ks_val = (double*)malloc(n_depth * sizeof(double));
    double* diss_co2 = (double*)malloc(n_depth* sizeof(double));
    double* growth_rate = (double*)malloc(n_depth * sizeof(double));

    for (int i = 0; i < n_depth; i++) {
        ks_val[i] = Ks(temperature_z[i]);
        iz_val[i] = irradiance(z_index[i], dz); //Term dependant of Irradiance
        temperature_z[i] = T_array[z_index[i]]; //Temperature at depth z (Celsius)
        diss_co2[i] = convertCO2(ppm_CO2, temperature_z[i]); //Terms dependant of T and CO2

        // Formule complète :
        growth_rate[i] = A * exp(-E/(R * temp(temperature_z[i]))) * (diss_co2[i] / (diss_co2[i] + ks_val[i])) * (iz_val[i] / I0); 
        growth_rate_output[i] = growth_rate[i]; 
        
    }              

    free(diss_co2);
    free(ks_val);
    free(iz_val);
    free(growth_rate);
    free(temperature_z);

}
double biomass_depth(double B, int z, double dz){
    double biomass_z = B * irradiance(z, dz) / I0;
    return biomass_z;
}
double biomass(double mu, double dz, int z, double B0){
    double b = biomass_depth(B0, z, dz) + 30*mu*biomass_depth(B0, z, dz);
    return b;
}


int main(int argc, char * argv[]) {
    DataRow data[MAX_ROWS]; 
    int row_count = read_csv("./data/co2_temp_for_c.csv", data, MAX_ROWS);
    
    if (row_count <= 1) { 
        if (row_count == -1) {
            // Message already displayed by read_csv
        } else {
            printf("Error: The CSV file is empty or contains only the header.\n");
        }
        return 1;
    }
    // values for thermal diffusion
    double dz = 1.0;  // Space step (m)
    double dt = 500; // Time step (seconds) - Choose dt for r <= 0.5 
    int n_points = 195; // Number of points in space, because bottom temperature was taken at 195m
    int n_steps = 10000;  // Number of time steps (simulates approximately 58 days)

    double r_value = alpha * dt / (dz * dz); //calculation of the stability factor for diffusion

    // Extracting data from the initial CSV file
    double* surfacetemp = (double*)malloc(n_csv * sizeof(double));
    double* bottomtemp = (double*)malloc(n_csv * sizeof(double));
    double* atmco2 = (double*)malloc(n_csv * sizeof(double));
    double* surfacebiomass = (double*)malloc(n_csv * sizeof(double));
    double* T = (double*)malloc(n_points * sizeof(double));

    //creation of z_index, depth tables for the biomass_growth function
    int *z_index = (int*)malloc(n_depth * sizeof(int));
    int z_values[] = {0, 1, 2, 5, 10, 30};  // Values to assign at z_index
    for (int i = 0; i < n_depth; i++) {
        z_index[i] = z_values[i]; // Initialisation with values
    }

    // creation of the CSV file that stores the biomass_growth values
    FILE* output_csv = fopen("./data/biomass_growth_results.csv", "w");
    if (output_csv == NULL) {
        printf("Error opening output file.\n");
        return 1;
    }
    fprintf(output_csv, "Index, Time, Depth_0,Depth_1,Depth_2,Depth_5,Depth_10,Depth_30\n");

    // creation of the CSV file that stores the biomass values
    FILE* outputcsv = fopen("./data/biomass_results.csv", "w");
    if (outputcsv == NULL) {
        printf("Error opening output file.\n");
        return 1;
    }
    fprintf(outputcsv, "Index, Time, Depth_0,Depth_1,Depth_2,Depth_5,Depth_10,Depth_30\n");

    for (int j = 1; j < n_csv; j++) {
        surfacetemp[j] = data[j].surface_temp;
        bottomtemp[j] = data[j].bottom_temp;
        atmco2[j] = data[j].atm_co2;
        surfacebiomass[j] = data[j].surface_biomass;
        
        // Initialisation of table of temperatures
        if (T == NULL) {
            printf("Memory allocation error for T.\n");
        return 1;
        }

        double surface_temp_j = surfacetemp[j];
        double bottom_temp_j = bottomtemp[j];

        // Initialise the temperature array (T) for this iteration
        for (int i = 0; i < n_points; i++) {
            // Dirichlet boundary conditions
            if (i == 0) {
                T[i] = surface_temp_j;  // Surface temperature
            } else if (i == n_points - 1) {
                T[i] = bottom_temp_j;   // Bottom temperature
            } else {
                T[i] = (surface_temp_j + bottom_temp_j) / 2;  // Initial temperature for other points
            }
        }

        // Calculate thermal diffusion for this moment
        diffusion_1d(T, alpha, dt, dz, n_points, n_steps);
        printf("\n--- Thermal diffusion simulation for time %d ---\n", j);
        printf("Surface temperature (T[0]): %.2f C\n", T[0]);
        printf("Bottom temperature (T[n_points-1]): %.2f C\n", T[n_points-1]);

        double growth_rate_output[n_depth];
        biomass_growth(T, atmco2[j], z_index, dz, growth_rate_output); // Calling the biomass function and retrieving results

        // Writing results to the first CSV file
        fprintf(output_csv, "%d,%s", j + 1, data[j].time);  // Time
        for (int i = 0; i < n_depth; i++) {
            fprintf(output_csv, ",%.30f", growth_rate_output[i]);  // Growth value at each depth
        }
        fprintf(output_csv, "\n");

        // Writing results to the second CSV file
        fprintf(outputcsv, "%d,%s", j + 1, data[j].time);  // Time

        // Loop for each column 
        for (int i = 0; i < n_depth; i++) {
            double biom = biomass(growth_rate_output[i], dz, z_index[i], surfacebiomass[j]);
            fprintf(outputcsv, ",%.30f", biom);  
        }

        fprintf(outputcsv, "\n");  
    }
    
    fclose(output_csv);
    fclose(outputcsv);

    //Displaying the parameters of the diffusion simulation
    printf("Simulation parameters:\n");
    printf("  dt (Time step): %.1f s\n", dt);
    printf("  n_steps (Total): %d\n", n_steps);
    printf("  r (Stability factor): %.6f\n", r_value);
    printf("  Total simulated time: %.2f hours (%.0f s)\n", (double)n_steps * dt / 3600.0, (double)n_steps * dt);
    


    printf("\n--- Biomass Results (After Diffusion) ---\n");

    double* depth_m = (double*)malloc(n_depth * sizeof(double));
    double* T_final = (double*)malloc(n_depth * sizeof(double));
    double* Iz_final = (double*)malloc(n_depth * sizeof(double));
    double* Iz_ratio = (double*)malloc(n_depth * sizeof(double));


    // Calculation and assignment of values in preallocated tables
    for (int i = 0; i < n_depth; i++) {
        depth_m[i] = dz * z_index[i];            // Depth
        T_final[i] = T[z_index[i]];              // Temperature at corresponding depth
        Iz_final[i] = irradiance(z_index[i], dz); // Irradiance at corresponding depth
        Iz_ratio[i] = Iz_final[i] / I0;           // Irradiance fraction
    }

    //Display of irradiance fraction and depth results
    for (int i = 0; i < n_depth; i++) {
        printf("Calculation depth  (Index T[%d]): %.1f m\n", i, depth_m[i]);

        printf("Irradiance (Iz/I0): %f\n", Iz_ratio[i]);
    }

    // Freeing allocated memory
    free(depth_m);
    free(T_final);
    free(Iz_final);
    free(Iz_ratio);
    free(z_index);
    free(atmco2);
    free(surfacetemp);
    free(bottomtemp);
    free(surfacebiomass);
    free(T);

    return 0;
} 
