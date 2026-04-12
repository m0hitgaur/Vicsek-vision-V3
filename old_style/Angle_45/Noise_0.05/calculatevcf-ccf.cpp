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

    

    // Start time tracking
    time_t start_time = time(NULL);

    // Main trial loop
    for (int trial = 0; trial < number_of_trials; trial++) {
        cout << "Trial number " << trial << endl;
        int time_step_counter = 0;
        time_t start_time_trial = time(NULL);
        vector<double> correlation_length_ccf_zero(number_of_timesteps, 0),correlation_length_ccf_one_over_e(number_of_timesteps, 0);
        vector<double> correlation_length_vcf_zero(number_of_timesteps, 0),correlation_length_vcf_one_over_e(number_of_timesteps, 0);    
        vector<double> correlation_time(number_of_timesteps, 0);
        vector<vector<double>> connected(bins, vector<double>(number_of_timesteps, 0));
        vector<vector<double>> velocitycorr(bins, vector<double>(number_of_timesteps, 0));    
        
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
                double v_avg_sq=(vx_avg*vx_avg)+(vy_avg*vy_avg);  
                // Correlation calculation
                double r = dr;
                double vcf0=0.0,ccf0=0.0,vcf0_e=0.0,ccf0_e=0.0;
                int zero_crossing_counter_ccf = 0;
                int zero_crossing_counter_vcf = 0;
                int one_over_e_counter_ccf = 0;
                int one_over_e_counter_vcf = 0;
                double vcf1,ccf1;
                for (int r_bin = 0; r_bin < bins; r_bin++) {
                    double ccf = 0.0;
                    double vcf = 0.0;
                    double ccf_e = 0.0;
                    double vcf_e = 0.0;                    
                    int particles_in_bin = 0;

                    for (int i = 0; i < number_of_agents; i++) {
                        for (int j = 0; j <= i; j++) {
                            double dx = position_x[i] - position_x[j];
                            double dy = position_y[i] - position_y[j];

                            dx -= Lx * nearbyint(dx / Lx); // Periodic boundary
                            dy -= Ly * nearbyint(dy / Ly);

                            double dist = sqrt(dx * dx + dy * dy);                            
                            if (dist >= r && dist < r + dr) {
                                ccf += (deviation[0][i] * deviation[0][j] + deviation[1][i] * deviation[1][j]);
                                vcf+= velocity*velocity*( (cos(theta[i]) *cos(theta[j])) +  (sin(theta[i]) *sin(theta[j])) );
                                particles_in_bin++;
                            }
                        }
                    }

                    if (particles_in_bin != 0) {
                        ccf /= particles_in_bin;
                        vcf/= particles_in_bin;
                    }
                    ccf -= dv_avg_sq;
                    vcf -= v_avg_sq;

                    if(r_bin==0){vcf1=vcf;ccf1=ccf;vcf/=vcf1;ccf/=ccf1;}
                    
                    if(r_bin!=0){
                        vcf/=vcf1;
                        ccf/=ccf1;
                        vcf_e=vcf-exp(-1);
                        ccf_e=ccf-exp(-1);
                    }

                    // Update connected and velocitycorr matrix
                    connected[r_bin][time_step_counter] = ccf;
                    velocitycorr[r_bin][time_step_counter]=vcf;
                    //if(trial==0&& time_step_counter==20)cout<<t<<" "<<r_bin<<" "<<time_step_counter<<" "<<connected[r_bin][time_step_counter]<<" "<<velocitycorr[r_bin][time_step_counter]<<"\n";

                    // Find correlation length
                    if (ccf * ccf0 < 0 && zero_crossing_counter_ccf == 0) {
                        correlation_length_ccf_zero[time_step_counter] = r - (ccf * dr / (ccf - ccf0));
                        zero_crossing_counter_ccf++;
                    }
                    ccf0 = ccf;

                    if (vcf * vcf0 < 0 && zero_crossing_counter_vcf == 0) {
                        correlation_length_vcf_zero[time_step_counter] = r - (vcf * dr / (vcf - vcf0));
                        correlation_time[time_step_counter] = t * dt;
                        zero_crossing_counter_vcf++;
                    
                    }
                    vcf0 = vcf;
                    
                    
                    if (ccf_e * ccf0_e < 0 && one_over_e_counter_ccf == 0) {
                        correlation_length_ccf_one_over_e[time_step_counter] = r - (ccf_e * dr / (ccf_e - ccf0_e));
                        one_over_e_counter_ccf++;
                    }
                    ccf0_e = ccf_e;
                    


                    if (vcf_e * vcf0_e < 0 && one_over_e_counter_vcf == 0) {
                        correlation_length_vcf_one_over_e[time_step_counter] = r - (vcf_e * dr / (vcf_e - vcf0_e));
                        one_over_e_counter_vcf++;
                    }
                    vcf0_e = vcf_e;
        
                    r += dr;
                }
                time_step_counter++;
                
                ofstream file1(folder_path + "/correlation_data/connectedcorrelation_vs_r_"+to_string(trial)+"_"+to_string(t)+"_.dat");
                if (!file1) {cerr << "Error writing connectedcorrelation_vs_r_"+to_string(trial)+"_"+to_string(t)+"_.dat" << endl;return 1;}
                for (int r_bin = 0; r_bin < bins; r_bin++) {file1 << dr * (r_bin) << " " << connected[r_bin][time_step_counter-1]  << "\n";}
                file1.close();
                ofstream file2(folder_path + "/correlation_data/velocitycorrelation_vs_r_"+to_string(trial)+"_"+to_string(t)+"_.dat");
                if (!file2) {cerr << "Error writing velocitycorrela100tion_vs_r_"+to_string(trial)+"_"+to_string(t)+"_.dat" << endl;return 1;}
                for (int r_bin = 0; r_bin < bins; r_bin++) {file2 << dr * (r_bin) << " " << velocitycorr[r_bin][time_step_counter-1]  << "\n";}
                file2.close();

            }

        }
        
        cout<<number_of_timesteps<<endl;
        

        ofstream file3(folder_path + "/correlation_data/velocorrlength_vs_time_zerocrossing_"+to_string(trial)+"_.dat");
        if (!file3) {cerr << "Error writing velocorrlength_vs_time_zerocrossing_"+to_string(trial)+"_.dat" << endl;return 1;}
        for (int t = 0; t < time_step_counter; t++) {file3 << correlation_time[t] << " " << correlation_length_vcf_zero[t] << "\n";}
        file3.close();

        ofstream file4(folder_path + "/correlation_data/velocorrlength_vs_time_one_over_e_"+to_string(trial)+"_.dat");
        if (!file4) {cerr << "Error writing velocorrlength_vs_time_zerocrossing_"+to_string(trial)+"_.dat" << endl;return 1;}
        for (int t = 0; t < time_step_counter; t++) {file4 << correlation_time[t] << " " << correlation_length_vcf_one_over_e[t] << "\n";}
        file4.close();


        ofstream file5(folder_path + "/correlation_data/connectedcorrlength_vs_time_zerocrossing_"+to_string(trial)+"_.dat");
        if (!file5) {cerr << "Error writing connectedcorrlength_vs_time_zerocrossing_"+to_string(trial)+"_.dat" << endl;return 1;}
        for (int t = 0; t < time_step_counter; t++) {file5 << correlation_time[t] << " " << correlation_length_ccf_zero[t] << "\n";}
        file5.close();

        ofstream file6(folder_path + "/correlation_data/connectedcorrlength_vs_time_one_over_e_"+to_string(trial)+"_.dat");
        if (!file6) {cerr << "Error writing connectedcorrlength_vs_time_zerocrossing_"+to_string(trial)+"_.dat" << endl;return 1;}
        for (int t = 0; t < time_step_counter; t++) {file6 << correlation_time[t] << " " << correlation_length_ccf_one_over_e[t] << "\n";}
        file6.close();

        



        time_t finish_time_trial = time(NULL);
        cout << "\nTime taken to calculate trial : " << finish_time_trial - start_time_trial << " seconds" << endl;
    }




    // Finish
    time_t finish_time = time(NULL);
    cout << "\nTotal time elapsed: " << finish_time - start_time << " seconds" << endl;

    return 0;
}
