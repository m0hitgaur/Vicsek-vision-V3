#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <sstream>
#include <thread>
#include <filesystem>
#include <iomanip>
#include "locker.h"
#include <algorithm> 
using namespace std;
namespace fs = filesystem;


struct Particle {
    double x, y;
    double vx, vy;
    double dvx=0.0, dvy=0.0;    
    vector <int> neighbours;
    double sigma=0.5;
};

class calculation {
private:
    const int N=640;                     // Number of particles
    const double Lx=32;                 // Box size
    const double Ly=32.0;                 //  Box size
    const double dt=1.0;                 // Timestep
    const double v0=0.05;                 // Magnitude of v0
    const int long_neighbours=4;     // Number of long range neighbours
    const double noise;              // Noise strength
    const double half_angle;         // Half of vision angle
    const double rc=1;                 // Cutoff radius
    const double beta=1;
    const int numberoftrials=2;              
    const int trial_start=0;                 // Trial number    
    const int tmax=2.0e3;                  // Total nujmber of time steps  
    string folder_path;              // Directory for saving data 
    vector<bool> time_record;        // Time points to record
    vector<int> time_array;               // Time points recorded 
    vector<Particle> particles;
    const int cluster_cutoff=10;
    const double r_c=0.5;
public:

    calculation(double noise_input,double angle):noise(noise_input),half_angle(angle*(M_PI/180.0)){
        particles.resize(N);
        ostringstream angle_str,Noise_str,L_str;
        angle_str<< fixed << setprecision(0) <<angle;
        Noise_str<< fixed << setprecision(2) <<noise;
        folder_path="data_vicsek_scalar/Angle_" + angle_str.str() +"/Noise_" + Noise_str.str()+"/"+"Avg_nearest_neighbour_data/";
        create_directory(folder_path);  
    }

    void load_snapshot(int step, int trial) {
        ifstream file(folder_path + "config_data/trial_" + to_string(trial) + "/config_" +to_string(step) + ".csv");
        if (!file.is_open()) {throw runtime_error("Could not open snapshot file config_data/trial_" + to_string(trial) + "/config_" +to_string(step) + ".csv");}
        string line;  
        // Skip header: "x,y,vx,vy"
        if (!getline(file, line)) {throw runtime_error("Snapshot file is empty");}
        vector<Particle> loaded_particles;
        while (getline(file, line)) {
            if (line.empty()) continue;
            stringstream ss(line);
            string item;
            Particle p;
            if (!getline(ss, item, ',')) break;
            p.x = stod(item);
            if (!getline(ss, item, ',')) break;
            p.y = stod(item);
            if (!getline(ss, item, ',')) break;
            p.vx = stod(item);
            if (!getline(ss, item, ',')) break;
            p.vy = stod(item);
            loaded_particles.push_back(p);
        }
        
        particles = move(loaded_particles);
    }
    void initialize_time_to_record(){
        for (int t = 0; t < tmax; t++) {
            bool should_record = false;
                if (t <= 10) should_record = true;
                else if(t>=10 && t<100&& t%10==0) should_record = true;                   
                else if (t<1000 && t >= 100 && t % 50 == 0) should_record = true;   
                else if (t>=1000 && t % 100 == 0) should_record=true;
                if (should_record)time_array.push_back(t);                        
                }
    }
    double minimum_image(double dx,double L){
        if(dx>L/2) dx=dx-L;
        if(dx<-L/2) dx=dx+L ;  
        return dx;
    }
    void start_calculation() {
        initialize_time_to_record();
        for (int trial = trial_start; trial < numberoftrials; trial++){
            time_t start_time_trial = time(NULL);
            cout<<"\n"<<"Trial number "<<trial<<fixed<<setprecision(2)<<" || Packing Fraction : "<<(N/(Lx*Ly))<<" | Noise = "<<noise<<" | Angle = "<<180*(half_angle/M_PI)<<" | N = "<<N<<" || "<<endl ;
            vector<double> correlation_length_ccf_zero(time_array.size(), 0),correlation_length_ccf_one_over_e(time_array.size(), 0);
            vector<double> correlation_length_vcf_zero(time_array.size(), 0),correlation_length_vcf_one_over_e(time_array.size(), 0);    
            vector<double> correlation_time(time_array.size(), 0),connected(bins, 0),velocitycorr(bins, 0);
            double timestarted=static_cast <double>(time(NULL));    
            for (int t=0;t<time_array.size();t++){
                vector<vector<double>> deviation(2, vector<double>(N, 0.0));
                load_snapshot(time_array[t],trial);
                 // Calculate average v0
                double vy_avg = 0.0, vx_avg = 0.0;
                for (Particle p:particles) {
                    vy_avg += p.vx;
                    vx_avg += p.vy;
                }
                vx_avg /= N;
                vy_avg /= N;

                // Calculate v0 deviations
                for (Particle p:particles) {
                    p.dvx = p.vx - vx_avg;
                    p.dvy = p.vy - vy_avg;
                }

                // Average deviation
                double dvx_avg = 0.0, dvy_avg = 0.0;
                for (Particle p:particles) {
                    dvx_avg += p.dvx;
                    dvy_avg += p.dvy;
                }
                dvx_avg /= N;
                dvy_avg /= N;

                double dv_avg_sq = hypot(dvx_avg,dvy_avg);
                double v_avg_sq= hypot(vx_avg,vy_avg);  
                // Correlation calculation
                double r = -dr;
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

                    for (int i = 0; i < particles.size(); i++) {
                        for (int j = 0; j <= i; j++) {
                            double dx = particles[i].x - particles[j].x;
                            double dy = particles[i].y - particles[j].y;

                            dx -= Lx * nearbyint(dx / Lx); // Periodic boundary
                            dy -= Ly * nearbyint(dy / Ly);

                            double dist = sqrt(dx * dx + dy * dy);                            
                            if (dist > r && dist <= r + dr) {
                                ccf += (particles[i].dvx*particles[j].dvx) + (particles[i].dvy*particles[j].dvy);
                                vcf+= (particles[i].vx*particles[j].vx) + (particles[i].vy*particles[j].vy);
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
                    connected[r_bin]= ccf;
                    velocitycorr[r_bin]=vcf;
                    // Find correlation length
                    if (ccf * ccf0 < 0 && zero_crossing_counter_ccf == 0) {
                        correlation_length_ccf_zero[t] = r - (ccf * dr / (ccf - ccf0));
                        zero_crossing_counter_ccf++;
                    }
                    ccf0 = ccf;

                    if (vcf * vcf0 < 0 && zero_crossing_counter_vcf == 0) {
                        correlation_length_vcf_zero[t] = r - (vcf * dr / (vcf - vcf0));
                        correlation_time[t] = t * dt;
                        zero_crossing_counter_vcf++;
                    
                    }
                    vcf0 = vcf;
                      
                    if (ccf_e * ccf0_e < 0 && one_over_e_counter_ccf == 0) {
                        correlation_length_ccf_one_over_e[t] = r - (ccf_e * dr / (ccf_e - ccf0_e));
                        one_over_e_counter_ccf++;
                    }
                    ccf0_e = ccf_e;
                    


                    if (vcf_e * vcf0_e < 0 && one_over_e_counter_vcf == 0) {
                        correlation_length_vcf_one_over_e[t] = r - (vcf_e * dr / (vcf_e - vcf0_e));
                        one_over_e_counter_vcf++;
                    }
                    vcf0_e = vcf_e;
        
                    r += dr;
                }
                
                ofstream file1(folder_path + "/correlation_data/connectedcorrelation_vs_r_"+to_string(trial)+"_"+to_string(t)+"_.dat");
                if (!file1) {cerr << "Error writing connectedcorrelation_vs_r_"+to_string(trial)+"_"+to_string(t)+"_.dat" << endl;}
                for (int r_bin = 0; r_bin < bins; r_bin++) {file1 << dr * (r_bin) << " " << connected[r_bin]  << "\n";}
                file1.close();
                ofstream file2(folder_path + "/correlation_data/v0correlation_vs_r_"+to_string(trial)+"_"+to_string(t)+"_.dat");
                if (!file2) {cerr << "Error writing v0correla100tion_vs_r_"+to_string(trial)+"_"+to_string(t)+"_.dat" << endl;}
                for (int r_bin = 0; r_bin < bins; r_bin++) {file2 << dr * (r_bin) << " " << velocitycorr[r_bin]  << "\n";}
                file2.close();




                
                if (time_array[t] % 500 == 0){
                    double frac = static_cast<double>(time_array[t]) / tmax;
                    double now  = static_cast<double>(time(NULL));
                    double elapsed = now - timestarted;

                    std::ostringstream label;
                    label << "[trial " << trial<< ", angle=" << 180.0 * (half_angle / M_PI)<< ", noise=" << noise << "]";
                    print_progress_threadsafe(label.str(), frac, elapsed);
                }
            }    
                
            ofstream file3(folder_path + "/correlation_data/velocorrlength_vs_time_zerocrossing_"+to_string(trial)+"_.dat");
            if (!file3) {cerr << "Error writing velocorrlength_vs_time_zerocrossing_"+to_string(trial)+"_.dat" << endl;}
            for (int t = 0; t < time_array.size(); t++) {file3 << correlation_time[t] << " " << correlation_length_vcf_zero[t] << "\n";}
            file3.close();

            ofstream file4(folder_path + "/correlation_data/velocorrlength_vs_time_one_over_e_"+to_string(trial)+"_.dat");
            if (!file4) {cerr << "Error writing velocorrlength_vs_time_zerocrossing_"+to_string(trial)+"_.dat" << endl;}
            for (int t = 0; t < time_array.size(); t++) {file4 << correlation_time[t] << " " << correlation_length_vcf_one_over_e[t] << "\n";}
            file4.close();


            ofstream file5(folder_path + "/correlation_data/connectedcorrlength_vs_time_zerocrossing_"+to_string(trial)+"_.dat");
            if (!file5) {cerr << "Error writing connectedcorrlength_vs_time_zerocrossing_"+to_string(trial)+"_.dat" << endl;}
            for (int t = 0; t < time_array.size(); t++) {file5 << correlation_time[t] << " " << correlation_length_ccf_zero[t] << "\n";}
            file5.close();

            ofstream file6(folder_path + "/correlation_data/connectedcorrlength_vs_time_one_over_e_"+to_string(trial)+"_.dat");
            if (!file6) {cerr << "Error writing connectedcorrlength_vs_time_zerocrossing_"+to_string(trial)+"_.dat" << endl;}
            for (int t = 0; t < time_array.size(); t++) {file6 << correlation_time[t] << " " << correlation_length_ccf_one_over_e[t] << "\n";}
            file6.close();

            time_t finish_time_trial = time(NULL);
            cout << "\nTime taken to calculate trial : " << finish_time_trial - start_time_trial << " seconds" << endl;
                
        } 

    }

};




int main() { 
    
    
    //double noise=1.0e-1 ;     // strength of noise
    //double half_angle=M_PI;   // Half of the vision angle  
    
    vector<double> angles = {180};
    vector<double> noises = { 0.05 }; 
    int seed=12345;           // random seed
    time_t trial_time,start_time=time(NULL) , finish_time;
    vector<thread> threads;
    unsigned int max_threads = thread::hardware_concurrency();
    if (max_threads == 0) max_threads = 4; 
 
    for (double angle : angles) {
        for (double noise : noises) {

            // If too many threads are running, wait for some to finish
            while (threads.size() >= max_threads) {
                threads.back().join();
                threads.pop_back();
            }

            // Launch a thread for this (angle, noise)
            threads.emplace_back([=]() {
                calculation calc(noise, angle);
                calc.start_calculation();
                
            });
        }
    }

    // Join remaining threads
    for (auto &th : threads) {
        th.join();
    }

    cout << "\nTotal time elapsed : " << time(NULL) - start_time << " seconds ";
    
    return 0;
}
