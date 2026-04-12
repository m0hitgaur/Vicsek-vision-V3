#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <numeric>
/*
3-02-2025

*/
namespace fs = std::filesystem;
using namespace std;

int main() {
    // Folder and parameter file path
    string folder_path = fs::current_path().string();
    ifstream param_file(folder_path + "/fluctuationparameters.txt");
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
    int numberofgrids = stoi(data[11]);
    
    cout << "Length of system = " << Lx << " * " << Ly << endl;
    cout << "Number of time steps = " << number_of_timesteps << endl;
    cout << "Number of agents = " << number_of_agents << endl;
    cout << "Density = " << density << endl;
    cout << "Noise = " << noise << endl;
    cout << "Total number of trials = " << number_of_trials << endl;

    int time_interval;
    // Start time tracking
    time_t start_time = time(NULL);

    // Main trial loop
    
    vector <int> timearray;
    for (int t = 0; t < number_of_timesteps; t++) {
        if(t<10) time_interval=1;
        if(t>10) time_interval=10;
        if(t>100) time_interval=50;
        if(t>1000)time_interval=100;
        if(t%time_interval==0){timearray.push_back(t);
        }
    }

    vector <vector<double>> avgnneighbour(number_of_trials,vector <double> (timearray.size(),0));
 
    for (int trial = 0; trial < number_of_trials; trial++) {
        
        cout << "Trial number " << trial << endl;
        time_t start_time_trial = time(NULL);
        
        ofstream file(folder_path+"/Avg_nearest_neighbour_data/avg_near_neigh_"+to_string(trial)+".dat");
        for (int t = 0; t < timearray.size(); t++) {
            // Initialize position and velocity vectors
            vector<double> position_x(number_of_agents, 0.0);
            vector<double> position_y(number_of_agents, 0.0);
            vector<double> theta(number_of_agents, 0.0);


            // File paths
            string theta_filename = folder_path + "/flockingdata/theta_" + to_string(trial) + "_" + to_string(static_cast<int>(timearray[t])) + "_.dat";
            string posix_filename = folder_path + "/flockingdata/positionx_" + to_string(trial) + "_" + to_string(static_cast<int>(timearray[t])) + "_.dat";
            string posiy_filename = folder_path + "/flockingdata/positiony_" + to_string(trial) + "_" + to_string(static_cast<int>(timearray[t])) + "_.dat";

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

            vector<int> count(number_of_agents,1); // starting from 1 as the particle itself is always counted            
            for (int i = 0; i < number_of_agents; i++) {
                
                for (int j = i+1; j < number_of_agents; j++) 
                {   double dx=position_x[j] - position_x[i];
                    double dy=position_y[j] - position_y[i];
                    if(dx>Lx/2) dx=dx-Lx;
                    else if(dx<-Lx/2) dx=dx+Lx ;  
                    if(dy>Ly/2) dy=dy-Ly;
                    else if(dy<-Ly/2) dy=dy+Ly ;                 

                    double rij = sqrt(pow(dx, 2) + pow(dy, 2));
                    
                    if(rij <= radius ){count[j]++;count[i]++;} 
                                            
                }    
            }
            for(int i=0;i<number_of_agents;i++)avgnneighbour[trial][t]+=count[i];
            avgnneighbour[trial][t]/=number_of_agents;
            file<<timearray[t]<<" "<<avgnneighbour[trial][t]<<"\n";            
        }
        file.close();


        time_t finish_time_trial = time(NULL);
        cout << "\nTime taken to calculate trial : " << finish_time_trial - start_time_trial << " seconds" << endl;
    }
    vector <double> tr_avg(timearray.size(),0);
    for(int t=0;t<timearray.size();t++){
        for (int tr=0;tr<number_of_trials;tr++)tr_avg[t]+=avgnneighbour[tr][t];
    tr_avg[t]/=number_of_trials;
    }
    

    ofstream file2(folder_path+"/Avg_nearest_neighbour_data/avg_near_neigh.dat");
    for(int t=0;t<timearray.size();t++)file2<<timearray[t]<<" "<< tr_avg[t]<<"\n";
    file2.close();

    // Finish
    time_t finish_time = time(NULL);
    cout << "\nTotal time elapsed: " << finish_time - start_time << " seconds" << endl;

    return 0;
}
