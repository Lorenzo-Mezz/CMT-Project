#include <stdio.h>
#include <math.h> 
void diffusion_1d(double* T, double alpha, double dt, double dz, int n_points, int n_steps) {
    double* T_new = (double*)malloc(n_points * sizeof(double));  // Tableau pour les nouvelles valeurs de T
    for (int t = 0; t < n_steps; t++) {
        // Appliquer la méthode de différences finies pour la dérivée seconde
        for (int i = 1; i < n_points - 1; i++) {
            T_new[i] = T[i] + (alpha * dt / (dz * dz)) * (T[i + 1] - 2 * T[i] + T[i - 1]);
        }

        // Mettre à jour les températures pour le prochain pas de temps
        for (int i = 1; i < n_points - 1; i++) {
            T[i] = T_new[i];
        }

        // Affichage des résultats à chaque étape
        printf("Step %d: ", t);
        for (int i = 0; i < n_points; i++) {
            printf("%.2f ", T[i]);
        }
        printf("\n");
    }
    free(T_new);  // Libérer la mémoire allouée pour T_new
}

int main() {

    return 0;
}
double irradiance(z) {
    double Iz = I0 * exp(-k * z)
    return Iz
}
double temp(t) {
    double T = 273.15 + t;
    return T;
}
double Kh(double T) {
    double kh = 1450*exp(-2400*(1/temp(T) - 1/298)); // mg/L équation de van't hoff pour le CO2
    return kh; //mg/L 
}
double convertCO2(double ppm_CO2, double T) { //convert from atmospheric CO2 to dissolved CO2
    double atm_CO2 = ppm_CO2 * pow(10, -6); //conversion de ppm à atm pour pouvoir utiliser la loi de Henry
    double diss_CO2 = Kh(temp(T)) * atm_CO2;
    return diss_CO2; //mg/L
}
double Ks(double T) {
    double ks = -0.2067 * pow(temp(T), 2) + 125.6 * temp(T) -19029; // mg/L
    return ks; //mg/L 
}
double biomass_growth(double T, double ppm_CO2, double z){
    double E_J = E * 4.184; // conersions en Joule/mol
    growth_rate = A * exp(-E_J/(R * temp(T))) * convertCO2(ppm_CO2, temp(T)) / (convertCO2(ppm_CO2, temp(T)) + Ks(temp(T))) * Irradiance(z) / I0 * ; // co2 en mg/L et T en Kelvin
    return growth_rate;
}
int main(int argc, char * argv[]) {
    double dz = 0.1;  // Pas de position (distance entre les points)
    double dt = 0.01; // Pas de temps
    int n_points = 20; // Nombre de points dans l'espace
    int n_steps = 50;  // Nombre de pas de temps

    // Initialisation du tableau de températures (T)
    double* T = (double*)malloc(n_points * sizeof(double));
    for (int i = 0; i < n_points; i++) {
        if (i == n_points / 2) {
            T[i] = 100.0;  // Température initiale élevée au centre
        } else {
            T[i] = 20.0;  // Température de fond plus faible
        }
    }

    // Affichage des températures initiales
    printf("Initial temperatures:\n");
    for (int i = 0; i < n_points; i++) {
        printf("%.2f ", T[i]);
    }
    printf("\n");

    // Résolution de l'équation de diffusion thermique
    diffusion_1d(T, D, dt, dz, n_points, n_steps);

    free(T);  // Libérer la mémoire allouée pour Tdouble deltaH = - 93.6; // kj/mol enthalpie standard
    double ks = 0;
    double kh = 0;
    double Ks_298K = 3.3 * pow(10, -2); // mol/L
    double growth_rate = 0;
    double T = 0;
    double Iz = 0;
    double atm_CO2 = 0;
    double diss_CO2 = 0;
    double R = 8.314; // constante des gaz parfait
    double A = 0.0077; //day^-1 constante
    double E = 2537.0; //cal/mol énergie d'activation
    double v = 24.45; //volume occupé par une mole de gaz parfait (L/mol)
    double M = 44.01; //g/mol masse molaire du CO2
    double alpha = pow(10, -7); //thermal molecular diffusivity (m2/s)
    double I0 = 1084.93; // Surface irradiance (cal/cm2day)
    double k = 1.97; // coefficient of light attenuation due to water itself (m-1)
    double gr = biomass_growth(309.15, 315, 15);
    printf("le taux d'agrandissement est %f", gr);
    return 0;
}
