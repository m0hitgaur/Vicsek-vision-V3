#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <numeric>
#include <algorithm> 

/*
26-06-2025

*/
namespace fs = std::filesystem;
using namespace std;






int main() {
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
    int cluster_cutoff=25;
    cout << "Length of system = " << Lx << " * " << Ly << endl;
    cout << "Number of time steps = " << number_of_timesteps << endl;
    cout << "Number of agents = " << number_of_agents << endl;
    cout << "Density = " << density << endl;
    cout<<"Angle = "<<(half_angle*180)/3.14<<endl;
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
        if(t%time_interval==0){timearray.push_back(t);}
    }


    vector <double>r_t(timearray.size(),0);
    for (int trial = 0; trial < number_of_trials; trial++) {
        
        cout << "Trial number " << trial << endl;
        time_t start_time_trial = time(NULL);
        
        
        for (int t = 0; t < timearray.size() ; t++) {
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

        
            vector <vector <int>> adjacent(number_of_agents);
            for (int i = 0; i < number_of_agents; i++) {
                for (int j = i+1; j <number_of_agents ; j++) 
                {   double dx=position_x[j] - position_x[i];
                    double dy=position_y[j] - position_y[i];
                    if(dx>Lx/2) dx=dx-Lx;
                    else if(dx<-Lx/2) dx=dx+Lx ;  
                    if(dy>Ly/2) dy=dy-Ly;
                    else if(dy<-Ly/2) dy=dy+Ly ;     
                    double r_ = sqrt(pow(dx, 2) + pow(dy, 2));
                    if(r_ <= 0.5){adjacent[i].push_back(j);adjacent[j].push_back(i);}             
        
                }
            }
            vector<vector<int>> clusters;   
            vector<bool>visited(number_of_agents,false);
            for(int i=0;i<number_of_agents;i++){
                if(!visited[i]){
                    vector<int> q;
                    vector<int> current_cluster;
                    q.push_back(i);
                    visited[i]=true;
                    current_cluster.push_back(i);
                    while(!q.empty()){
                        int u=q[0];
                        q.erase(q.begin());  
                        for(int element:adjacent[u]){
                            if(!visited[element]){
                                visited[element]=true;
                                q.push_back(element);
                                current_cluster.push_back(element);
                            }
                        }
                    }
                    sort(current_cluster.begin(),current_cluster.end());
                    clusters.push_back(current_cluster);
                }
            }

            vector<int> cluster_index(number_of_agents,0);
            int totalclusters=1;
            for(auto row:clusters){
                for(int element:row){
                    cluster_index[element]=totalclusters;
                }
                totalclusters++;
            }

            ofstream f(folder_path+"/cluster_index/clusterindex_"+to_string(trial)+"_"+to_string(timearray[t])+"_.txt");
            for(int a:cluster_index){f<<a<<"\n";}
            f.close();
           
            
            vector<double> radius;
            for(int j=0;j<clusters.size();j++){
                if(clusters[j].size()>=cluster_cutoff){

                    double alpha_x=0.0,alpha_y=0.0,beta_x=0.0,beta_y=0.0;
                    for(int i:clusters[j]){
                        
                        double theta_x=(2*M_PI)*position_x[i]/Lx;
                        double theta_y=(2*M_PI)*position_y[i]/Ly;
                        
                        alpha_x+=cos(theta_x);beta_x+=sin(theta_x);
                        alpha_y+=cos(theta_y);beta_y+=sin(theta_y);
                    }
                    alpha_x/=clusters[j].size();
                    alpha_y/=clusters[j].size();
                    beta_x/=clusters[j].size();
                    beta_y/=clusters[j].size();    

                    double x_cm=(Lx/(2*M_PI))*(atan2(-beta_x,-alpha_x)+M_PI);
                    double y_cm=(Ly/(2*M_PI))*(atan2(-beta_y,-alpha_y)+M_PI);
                    
                    double rg_sq=0.0;
                    for(int i:clusters[j]){

                        double dx=position_x[i]-x_cm;
                        double dy=position_y[i]-y_cm;

                        if(dx>Lx/2) dx=dx-Lx;
                        else if(dx<-Lx/2) dx=dx+Lx ;  
                        if(dy>Ly/2) dy=dy-Ly;
                        else if(dy<-Ly/2) dy=dy+Ly ;                        
                        rg_sq+=(pow(dx,2)+pow(dy,2));
                    }
                
                    rg_sq/=clusters[j].size();
                    
                    radius.push_back(sqrt(rg_sq));
                }
            }

 
            double r__=0;
            for(double i:radius)r__+=i;
            r_t[t]+=r__/radius.size();

            ofstream fi(folder_path+"/cluster_index/radius_"+to_string(trial)+"_"+to_string(timearray[t])+"_.txt");
            for(auto a:radius){fi<<a<<"\n";}
            fi.close();

            ofstream fil(folder_path+"/cluster_index/mass_"+to_string(trial)+"_"+to_string(timearray[t])+"_.txt");
            
            vector <double> mass;
            for(auto row:clusters){if(row.size()>=cluster_cutoff){mass.push_back(row.size());}}
            for(auto element:mass){fil<<element<<"\n";}
            fil.close();            


    
    }
    time_t finish_time_trial = time(NULL);
    cout << "\nTime taken to calculate trial : " << finish_time_trial - start_time_trial << " seconds" << endl;

    }
    
    // Finish
    time_t finish_time = time(NULL);
    cout << "\nTotal time elapsed: " << finish_time - start_time << " seconds" << endl;

    ofstream fill(folder_path+"/cluster_index/radius_vs_t_.txt");
    for(int i=0;i<timearray.size();i++){fill<<timearray[i]<<" "<<r_t[i]/number_of_trials<<"\n";}
    fill.close();



    return 0;
}
