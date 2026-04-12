#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <iomanip>
/*
7-01-2025
working 
need to check for higher particles
*/
namespace fs = std::filesystem;
using namespace std;

int main() {
    // Folder and parameter file path
    string folder_path = fs::current_path().string();
    ifstream param_file(folder_path + "/ccf-vcfparameters.txt");
    if (!param_file) {
        cerr << "Error opening parameter file" << endl;
        return 1;
    }

    // Reading parameters
    string line;
    vector<string> data;
    while (getline(param_file, line)) {
        data.push_back(line);
    }
    param_file.close();

    // Validate parameter count
    if (data.size() < 13) {
        cerr << "Insufficient parameters in the input file." << endl;
        return 1;
    }

    // Extract parameters
    int number_of_agents = stoi(data[0]);
    double Lx = stod(data[1]);
    double Ly = stod(data[2]);
    int number_of_timesteps = stoi(data[3]);
    double velocity = stod(data[4]);
    double dt = stod(data[5]);
    double half_angle = stod(data[6]);
    double noise = stod(data[7]);
    double density = stod(data[8]);
    double radius = stod(data[9]);
    int number_of_trials = stoi(data[10]);
    double dr = stod(data[11]);
    int time_interval = stoi(data[12]);

    cout << "Length of system = " << Lx << " * " << Ly << endl;
    cout << "Number of time steps = " << number_of_timesteps << endl;
    cout << "Number of agents = " << number_of_agents << endl;
    cout << "Density = " << density << endl;
    cout << "Noise = " << noise << endl;
    cout << "Total number of trials = " << number_of_trials << endl;

    // Derived variables
    double rmax = Lx / 2;
    int bins = static_cast<int>(rmax / dr);
    vector<double> correlation_length_ccf_zero(number_of_timesteps, 0);
    vector<double> correlation_length_vcf_zero(number_of_timesteps, 0);    
    vector<double> correlation_time(number_of_timesteps, 0);
    vector<vector<double>> connected(bins, vector<double>(number_of_timesteps, 0));
    vector<vector<double>> velocitycorr(bins, vector<double>(number_of_timesteps, 0));
    int time_con_correlation;

    // Start time tracking
    time_t start_time = time(NULL);

    // Main trial loop
    for (int trial = 0; trial < number_of_trials; trial++) {
        cout << "Trial number " << trial << endl;
        int time_step_counter = 0;
        time_t start_time_trial = time(NULL);
        for (int t = 0; t < number_of_timesteps; t++) {
            if(t<10) time_interval=1;
            if(t>10) time_interval=10;
            if(t>100) time_interval=50;
            if(t>1000)time_interval=100;
                        
            if (t % time_interval == 0) {
                cout << t << " >> ";

                // Initialize position and velocity vectors
                vector<double> position_x(number_of_agents, 0.0);
                vector<double> position_y(number_of_agents, 0.0);
                vector<double> theta(number_of_agents, 0.0);
                vector<vector<double>> deviation(2, vector<double>(number_of_agents, 0.0));

                // File paths
                string theta_filename = folder_path + "/flockingdata/theta_" + to_string(trial) + "_" + to_string(static_cast<int>(t * dt)) + "_.dat";
                string posix_filename = folder_path + "/flockingdata/positionx_" + to_string(trial) + "_" + to_string(static_cast<int>(t * dt)) + "_.dat";
                string posiy_filename = folder_path + "/flockingdata/positiony_" + to_string(trial) + "_" + to_string(static_cast<int>(t * dt)) + "_.dat";

                // Open files
                ifstream theta_file(theta_filename), posix_file(posix_filename), posiy_file(posiy_filename);
                if (!theta_file.is_open() || !posix_file.is_open() || !posiy_file.is_open()) {
                    cerr << "Failed to open input files at timestep " << t << endl;
                    return 1;
                }

                // Read data
                for (int i = 0; i < number_of_agents; i++) {
                    if (!(theta_file >> theta[i]) || !(posix_file >> position_x[i]) || !(posiy_file >> position_y[i])) {
                        cerr << "Error reading data from files at timestep " << t << endl;
                        return 1;
                    }
                }

                // Calculate average velocity
                double vy_avg = 0.0, vx_avg = 0.0;
                for (int i = 0; i < number_of_agents; i++) {
                    vy_avg += velocity * sin(theta[i]);
                    vx_avg += velocity * cos(theta[i]);
                }
                vx_avg /= number_of_agents;
                vy_avg /= number_of_agents;

                // Calculate velocity deviations
                for (int i = 0; i < number_of_agents; i++) {
                    deviation[0][i] = velocity * cos(theta[i]) - vx_avg;
                    deviation[1][i] = velocity * sin(theta[i]) - vy_avg;
                }

                // Average deviation
                double dvx_avg = 0.0, dvy_avg = 0.0;
                for (int i = 0; i < number_of_agents; i++) {
                    dvx_avg += deviation[0][i];
                    dvy_avg += deviation[1][i];
                }
                dvx_avg /= number_of_agents;
                dvy_avg /= number_of_agents;
                double dv_avg_sq = (dvx_avg*dvx_avg) + (dvy_avg*dvy_avg);
                double v_avg_sq=(vx_avg*vy_avg)+(vy_avg*vy_avg);  
                // Correlation calculation
                double r = dr;
                double vcf0=0.0,ccf0 = 0.0;
                int zero_crossing_counter_ccf = 0;
                int zero_crossing_counter_vcf = 0;
                //import rij data 
                ifstream rijdata(folder_path+"/rij_data/rij_"+to_string(trial)+"_"+to_string(int(dt*t))+"_.txt");
                if(!rijdata.is_open()) {cout<<folder_path+"/rij_data/rij_"+to_string(trial)+"_"+to_string(int(dt*t))+"_.txt not open";}
                vector<vector <double>> rij;
                while(getline(rijdata,line)){
                    istringstream r(line);
                    double value;
                    vector <double>row;
                    while(r>>value){
                        row.push_back(value);
                    }
                    rij.push_back(row);
                }


                for (int r_bin = 0; r_bin < bins; r_bin++) {
                    double ccf = 0.0;
                    double vcf = 0.0;
                    int particles_in_bin = 0;

                    for (int i = 0; i < number_of_agents; i++) {
                        for (int j = i+1; j <= i; j++) {
                            
                            if (rij[i][j] >= r && rij[i][j] < r + dr) {
                                ccf += (deviation[0][i] * deviation[0][j] + deviation[1][i] * deviation[1][j]);
                                vcf+= velocity*velocity*( (cos(theta[i]) *cos(theta[j])) +  (sin(theta[i]) *sin(theta[j])) );
                                particles_in_bin++;
             time_t finish_time_trial = time(NULL);               }
                        }
                    }

                    if (particles_in_bin != 0) {
                        ccf /= particles_in_bin;
                        vcf/= particles_in_bin;
                    }
                    ccf -= dv_avg_sq;
                    vcf -= v_avg_sq;

                    // Update connected and velocitycorr matrix
                    connected[r_bin][time_step_counter] += ccf;
                    velocitycorr[r_bin][time_step_counter]+=vcf;

                    // Find correlation length
                    if (ccf * ccf0 < 0 && zero_crossing_counter_ccf == 0) {
                        correlation_length_ccf_zero[time_step_counter] += r - (ccf * dr / (ccf - ccf0));
                        zero_crossing_counter_ccf++;
                    }
                    ccf0 = ccf;
                    if (vcf * vcf0 < 0 && zero_crossing_counter_vcf == 0) {
                        correlation_length_vcf_zero[time_step_counter] += r - (vcf * dr / (vcf - vcf0));
                        correlation_time[time_step_counter] = t * dt;
                        zero_crossing_counter_vcf++;
                    
                    }
                    vcf0 = vcf;                    
                    r += dr;
                }
                time_step_counter++;
            }
        }
        cout<<number_of_timesteps<<endl;
        time_con_correlation = time_step_counter;

        time_t finish_time_trial = time(NULL);
        cout << "\nTime taken to calculate trial : " << finish_time_trial - start_time_trial << " seconds" << endl;
    }

    // Average velocity correlations
    for (int t = 0; t < time_con_correlation; t++) {
        correlation_length_vcf_zero[t] /= number_of_trials;
        correlation_length_ccf_zero[t] /= number_of_trials;
        for (int r_bin = 0; r_bin < bins; r_bin++) {
            velocitycorr[r_bin][t] /= number_of_trials;
            connected[r_bin][t] /= number_of_trials;
        }
    }
    vector <double> correlation_length_vcf_e,correlation_length_ccf_e;
    for (int t = 0; t < time_con_correlation; t++) {
        double v0=0.0,c0=0.0;
        int one_over_e_counter_ccf=0,one_over_e_counter_vcf=0;
        for (int r_bin = 0; r_bin < bins; r_bin++) {
            double vcf=velocitycorr[r_bin][t] - (exp(-1)*velocitycorr[0][t]) ;
            double ccf=connected[r_bin][t] - (exp(-1)*connected[0][t]) ;
            if(vcf*v0<0 && one_over_e_counter_vcf==0){correlation_length_vcf_e.push_back((r_bin*dr) - (vcf*dr)/(vcf-v0));  one_over_e_counter_vcf++;}
            if(c0*ccf<0 && one_over_e_counter_ccf==0){correlation_length_ccf_e.push_back((r_bin*dr) - (ccf*dr)/(ccf-c0));  one_over_e_counter_ccf++;}
            v0=vcf;
            c0=ccf;
        }
    }


    // Write outputs
    ofstream velocorrzero_file(folder_path + "/correlation_data/velocorrlength_vs_time_zerocrossing.dat");
    if (!velocorrzero_file) {
        cerr << "Error writing velocorrlength_vs_time_zerocrossing.dat" << endl;
        return 1;
    }
    for (int t = 0; t < time_con_correlation; t++) {
        velocorrzero_file << correlation_time[t] << " " << correlation_length_vcf_zero[t] << "\n";
    }
    velocorrzero_file.close();

    ofstream velocorr_file(folder_path + "/correlation_data/velocorrlength_vs_time_one_over_e.dat");
    if (!velocorr_file) {
        cerr << "Error writing velocorrlength_vs_time_oneovere.dat" << endl;
        return 1;
    }
    for (int t = 0; t < time_con_correlation; t++) {
        velocorr_file << correlation_time[t] << " " << correlation_length_vcf_e[t] << "\n";
    }
    velocorr_file.close();    

    for (int t = 0; t < time_con_correlation; t++) {
        ofstream velocorr_vs_r_file(folder_path + "/correlation_data/velocitycorrelation_vs_r_" + to_string(static_cast<int>(dt * correlation_time[t])) + "_.dat");
        if (!velocorr_vs_r_file) {
            cerr << "Error writing velocorrelation_vs_r file" << endl;
            return 1;
        }
        for (int r_bin = 0; r_bin < bins; r_bin++) {
            velocorr_vs_r_file << dr * (r_bin) << " " << velocitycorr[r_bin][t] / velocitycorr[0][t] << "\n";
        }
        velocorr_vs_r_file.close();
    }


    // Write outputs
    ofstream corrz_file(folder_path + "/correlation_data/connectedcorrlength_vs_time_zerocrossing.dat");
    if (!corrz_file) {
        cerr << "Error writing corrlength_vs_time_zero.dat" << endl;
        return 1;
    }
    for (int t = 0; t < time_con_correlation; t++) {
        corrz_file << correlation_time[t] << " " << correlation_length_ccf_zero[t] << "\n";
    }
    corrz_file.close();

    ofstream corr_file(folder_path + "/correlation_data/connectedcorrlength_vs_time_one_over_e.dat");
    if (!corr_file) {
        cerr << "Error writing corrlength_vs_time_oneovere.dat" << endl;
        return 1;
    }
    for (int t = 0; t < time_con_correlation; t++) {
        corr_file << correlation_time[t] << " " << correlation_length_ccf_e[t] << "\n";
    }
    corr_file.close();

    for (int t = 0; t < time_con_correlation; t++) {
        ofstream corr_vs_r_file(folder_path + "/correlation_data/connectedcorrelation_vs_r_" + to_string(static_cast<int>(dt * correlation_time[t])) + "_.dat");
        if (!corr_vs_r_file) {
            cerr << "Error writing correlation_vs_r file" << endl;
            return 1;
        }
        for (int r_bin = 0; r_bin < bins; r_bin++) {
            corr_vs_r_file << dr * (r_bin) << " " << connected[r_bin][t] / connected[0][t] << "\n";
        }
        corr_vs_r_file.close();
    }


    // Finish
    time_t finish_time = time(NULL);
    cout << "\nTotal time elapsed: " << finish_time - start_time << " seconds" << endl;

    return 0;
}
