#include <stdio.h>
#include <math.h> 
#include <stdlib.h>
#include <string.h>

// Définitions de constantes
#define MAX_LINE_LENGTH 90
#define MAX_ROWS 301

// Structure pour stocker les données de chaque ligne du CSV
typedef struct {
    char time[11];
    float surface_temp; // Température de surface (C)
    float atm_co2;      // CO2 atmosphérique (ppm)
    float bottom_temp;
    float surface_biomass;  // Température du fond (C)
} DataRow;

// --- Fonctions de Lecture CSV ---

/**
 * Lit un fichier CSV contenant des données temporelles.
 * @param co2_temp_for_c Nom du fichier CSV.
 * @param data Tableau de structures DataRow.
 * @param max_rows Nombre maximal de lignes à lire.
 * @return Le nombre de lignes lues (sans compter l'en-tête) ou -1 en cas d'erreur.
 */
int read_csv(const char *co2_temp_for_c, DataRow data[], int max_rows) {
    FILE *file = fopen(co2_temp_for_c, "r");
    if (file == NULL) {
        printf("Erreur lors de l'ouverture du fichier '%s'.\n", co2_temp_for_c );
        return -1;
    }

    char line[MAX_LINE_LENGTH];
    int row_count = 0;

    while (fgets(line, sizeof(line), file) != NULL) {
        if (row_count >= max_rows) {
            break; 
        }

        // Ignorer l'en-tête (première ligne)
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

// Constantes pour le modèle de biomasse et de réaction (variables globales)
double R = 8.314;           // Constante des gaz parfaits (J/(mol*K))
double A = 0.0077;          // Constante du taux de croissance max (jour^-1)
double E = 2537.0;          // Énergie d'activation (cal/mol)
double I0 = 1084.93;        // Irradiance de surface (cal/(cm2*jour))
double k = 0.547;           // Coefficient d'atténuation de la lumière (m-1)
const int n_depth = 6;                //number of depth we're interested in
const int n_csv = 301;                //number of lines in our csv
// Constante pour la diffusion thermique
double alpha = 6.6 * pow(10, -5); // Diffusivité thermique (m2/s)

// --- Fonctions de Modélisation ---

//Résout l'équation de diffusion thermique 1D (schéma explicite de Euler). Met à jour le tableau T avec le profil de température final après n_steps.

void diffusion_1d(double* T, double alpha, double dt, double dz, int n_points, int n_steps) {
    double* T_new = (double*)malloc(n_points * sizeof(double));
    if (T_new == NULL) {
        printf("Erreur d'allocation mémoire pour T_new.\n");
        return; 
    }

    // Calcul du facteur de stabilité (Doit être r <= 0.5 pour la stabilité)
    double r = alpha * dt / (dz * dz); 

    printf("--- Simulation de Diffusion Thermique (r=%.6f) ---\n", r);
    printf("Step 0 (Initial): ");
    for (int i = 0; i < n_points; i++) {
        printf("%.2f ", T[i]);
    }
    printf("\n");
    if (r > 0.5) {
        printf("AVERTISSEMENT: Le facteur de stabilite r (%.6f) est > 0.5. Resultats INSTABLES.\n", r);
    }
    
    // Boucle de simulation
    for (int t = 0; t < n_steps; t++) {
        
        // Calcul du nouveau profil pour les points intérieurs
        for (int i = 1; i < n_points - 1; i++) {
            T_new[i] = T[i] + r * (T[i + 1] - 2 * T[i] + T[i - 1]);
        }
        
        // Conditions aux limites de Dirichlet (Températures T[0] et T[n_points-1] fixes)
        T_new[0] = T[0]; 
        T_new[n_points - 1] = T[n_points - 1]; 

        // Mise à jour de T pour la prochaine itération
        for (int i = 0; i < n_points; i++) {
            T[i] = T_new[i];
        }

        // Affichage régulier
        if ((t + 1) % 1000 == 0) { 
            printf("Step %d: ", t + 1);
            for (int i = 0; i < n_points; i++) {
                printf("%.2f ", T[i]);
            }
            printf("\n");
        }
    }
    
    // Affichage de l'état final
    printf("Step %d (Final): ", n_steps);
    for (int i = 0; i < n_points; i++) {
        printf("%.2f ", T[i]);
    }
    printf("\n");
    
    printf("--- Simulation de diffusion terminee ---\n");
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

//Convertit Celsius en Kelvin.
double temp(double t_celsius) {
    double T_kelvin = 273.15 + t_celsius;
    return T_kelvin;
}

//Calcule la constante de Henry Kh (mg/L) en fonction de T (Celsius).
double Kh(double temperature_z) {
    double kh = 1450 * exp(-2400 * (1/temp(temperature_z) - 1/298)); 
    return kh; // mg/L 
}

//Convertit le CO2 atmosphérique (ppm) en CO2 dissous (mg/L).
double convertCO2(double atm_co2, double temperature_z) { 
    double co2_atm = atm_co2 * pow(10, -6); // conversion de ppm à atm
    double diss_CO2 = Kh(temperature_z) * co2_atm;
    return diss_CO2; // mg/L
}

//Calcule la constante de semi-saturation Ks (mg/L) en fonction de T (Celsius).

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
        iz_val[i] = irradiance(z_index[i], dz); //Terme dépendant de l'Irradiance
        temperature_z[i] = T_array[z_index[i]]; //Température à la profondeur z_index (en Celsius)
        diss_co2[i] = convertCO2(ppm_CO2, temperature_z[i]); //Termes dépendant de T et CO2

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
            // Message déjà affiché par read_csv
        } else {
            printf("Erreur: Le fichier CSV est vide ou ne contient que l'en-tête.\n");
        }
        return 1;
    }

    double dz = 1.0;  // Pas d'espace (m)
    double dt = 500; // Pas de temps (secondes) - Choisir dt pour r <= 0.5 (ici r ~ 0.4)
    int n_points = 195; // Nombre de points dans l'espace, beacause bottom temperature was taken at 195m
    int n_steps = 10000;  // Nombre de pas de temps (simule environ 11 jours)

    double r_value = alpha * dt / (dz * dz); //calcul du facteur de stabilité pour la diffusion

    // Extraction des données du CSV initial
    double* surfacetemp = (double*)malloc(n_csv * sizeof(double));
    double* bottomtemp = (double*)malloc(n_csv * sizeof(double));
    double* atmco2 = (double*)malloc(n_csv * sizeof(double));
    double* surfacebiomass = (double*)malloc(n_csv * sizeof(double));
    double* T = (double*)malloc(n_points * sizeof(double));

    //création de z_index, tableaux des profondeurs pour la fonction biomass_growth
    int *z_index = (int*)malloc(n_depth * sizeof(int));
    int z_values[] = {0, 1, 2, 5, 10, 30};  // Valeurs à assigner à z_index
    for (int i = 0; i < n_depth; i++) {
        z_index[i] = z_values[i]; // Initialisation avec les valeurs
    }

    // création du CSV qui stocke les valeurs de biomass_growth
    FILE* output_csv = fopen("./data/biomass_growth_results.csv", "w");
    if (output_csv == NULL) {
        printf("Erreur lors de l'ouverture du fichier de sortie.\n");
        return 1;
    }
    fprintf(output_csv, "Index, Time, Depth_0,Depth_1,Depth_2,Depth_5,Depth_10,Depth_30\n");

    // création du CSV qui stocke les valeurs de biomasse
    FILE* outputcsv = fopen("./data/biomass_results.csv", "w");
    if (outputcsv == NULL) {
        printf("Erreur lors de l'ouverture du fichier de sortie.\n");
        return 1;
    }
    fprintf(outputcsv, "Index, Time, Depth_0,Depth_1,Depth_2,Depth_5,Depth_10,Depth_30\n");

    for (int j = 1; j < n_csv; j++) {
        surfacetemp[j] = data[j].surface_temp;
        bottomtemp[j] = data[j].bottom_temp;
        atmco2[j] = data[j].atm_co2;
        surfacebiomass[j] = data[j].surface_biomass;
        
        // Initialisation du tableau de températures
        if (T == NULL) {
            printf("Erreur d'allocation mémoire pour T.\n");
        return 1;
        }

        double surface_temp_j = surfacetemp[j];
        double bottom_temp_j = bottomtemp[j];

        // Initialiser le tableau des températures (T) pour cette itération
        for (int i = 0; i < n_points; i++) {
            // Conditions aux limites de Dirichlet
            if (i == 0) {
                T[i] = surface_temp_j;  // Température de surface
            } else if (i == n_points - 1) {
                T[i] = bottom_temp_j;   // Température de fond
            } else {
                T[i] = (surface_temp_j + bottom_temp_j) / 2;  // Température initiale pour les autres points
            }
        }

        // Calcul de la diffusion thermique pour cet instant
        diffusion_1d(T, alpha, dt, dz, n_points, n_steps);
        printf("\n--- Simulation de diffusion thermique pour le temps %d ---\n", j);
        printf("Température à la surface (T[0]): %.2f C\n", T[0]);
        printf("Température au fond (T[n_points-1]): %.2f C\n", T[n_points-1]);

        double growth_rate_output[n_depth];
        biomass_growth(T, atmco2[j], z_index, dz, growth_rate_output); // Appel de la fonction biomasse et récupération des résultats

        // Écriture des résultats dans le 1 er fichier CSV
        fprintf(output_csv, "%d,%s", j + 1, data[j].time);  // Temps
        for (int i = 0; i < n_depth; i++) {
            fprintf(output_csv, ",%.30f", growth_rate_output[i]);  // Valeur de croissance pour chaque profondeur
        }
        fprintf(output_csv, "\n");

        // Écriture des résultats dans le 2 eme fichier CSV
        fprintf(outputcsv, "%d,%s", j + 1, data[j].time);  // Temps

        // Boucle pour chaque colonne (n_csv)
        for (int i = 0; i < n_csv; i++) {
            double biom = biomass(growth_rate_output[i], dz, z_index[i], surfacebiomass[j]);
            fprintf(outputcsv, ",%.30f", biom);  // Remplir avec la biomasse calculée
        }

        fprintf(outputcsv, "\n");  // Nouvelle ligne après chaque ligne de profondeur
    }
    
    fclose(output_csv);
    fclose(outputcsv);

    //Affichage des paramètres de la simulation de la diffusion
    printf("Parametres de simulation:\n");
    printf("  dt (Pas de temps): %.1f s\n", dt);
    printf("  n_steps (Total): %d\n", n_steps);
    printf("  r (Facteur de Stabilite): %.6f\n", r_value);
    printf("  Temps simule total: %.2f heures (%.0f s)\n", (double)n_steps * dt / 3600.0, (double)n_steps * dt);
    


    printf("\n--- Resultats de Biomasse (Apres Diffusion) ---\n");

    double* depth_m = (double*)malloc(n_depth * sizeof(double));
    double* T_final = (double*)malloc(n_depth * sizeof(double));
    double* Iz_final = (double*)malloc(n_depth * sizeof(double));
    double* Iz_ratio = (double*)malloc(n_depth * sizeof(double));


    // Calcul et assignation des valeurs dans les tableaux préalloués
    for (int i = 0; i < n_depth; i++) {
        depth_m[i] = dz * z_index[i];            // Profondeur
        T_final[i] = T[z_index[i]];              // Température à la profondeur correspondante
        Iz_final[i] = irradiance(z_index[i], dz); // Irradiance à la profondeur correspondante
        Iz_ratio[i] = Iz_final[i] / I0;           // Fraction d'irradiance
    }

    // Affichage des résultats
    for (int i = 0; i < n_depth; i++) {
        printf("Profondeur de calcul (Index T[%d]): %.1f m\n", i, depth_m[i]);

        printf("Irradiance (Iz/I0): %f\n", Iz_ratio[i]);
    }

    // Libération de la mémoire allouée
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
