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

// standard Vicsek model
struct Particle {
    double x, y;
    double vx, vy;
    vector <int> neighbours;
    double sigma=0.5;
};

class calculation {
private:
    const int N=640;                     // Number of particles
    const double Lx=32;                 // Box size
    const double Ly=32.0;                 //  Box size
    const double dt=1.0;                 // Timestep
    const double v0=0.05;                 // Magnitude of velocity
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
public:
    calculation(double noise_input,double angle)
               :noise(noise_input),half_angle(angle*(M_PI/180.0)){
        particles.resize(N);
        ostringstream angle_str,Noise_str,L_str;
        angle_str<< fixed << setprecision(0) <<angle;
        Noise_str<< fixed << setprecision(2) <<noise;
        folder_path="data_vicsek_scalar/Angle_" + angle_str.str() +"/Noise_" + Noise_str.str()+"/";
        create_directory(folder_path+"cluster_data");
        
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
    void calculate_neighbours(){
        for(Particle p:particles)p.neighbours.clear();
            for (int i = 0; i < N; i++) {
                for (int j = i+1; j <N ; j++) 
                {   double dx=minimum_image(particles[j].x - particles[i].x,Lx);
                    double dy=minimum_image(particles[j].y - particles[i].y,Ly);     
                    double r = hypot(dx,dy);
                    if(r <= 0.5){particles[i].neighbours.push_back(j);particles[j].neighbours.push_back(i);}            
                }
        }
    }
    double minimum_image(double dx,double L){
        if(dx>L/2) dx=dx-L;
        if(dx<-L/2) dx=dx+L ;  
        return dx;
    }
    vector<vector<int>>  BFS(){
        vector<vector<int>> cluster;   
        vector<bool>visited(N,false);
            for(int i=0;i<N;i++){
                if(!visited[i]){
                    vector<int> q;
                    vector<int> current_cluster;
                    q.push_back(i);
                    visited[i]=true;
                    current_cluster.push_back(i);
                    while(!q.empty()){
                        int u=q[0];
                        q.erase(q.begin());  
                        for(int element:particles[u].neighbours){
                            if(!visited[element]){
                                visited[element]=true;
                                q.push_back(element);
                                current_cluster.push_back(element);
                            }
                        }
                    }
                    sort(current_cluster.begin(),current_cluster.end());
                    cluster.push_back(current_cluster);
                }
            }
        return cluster;
    }
    vector<double> calculate_rg(vector<vector<int>> clusters){
        vector<double> radius;
        for(int j=0;j<clusters.size();j++){
            if(clusters[j].size()>=cluster_cutoff){

                double alpha_x=0.0,alpha_y=0.0,beta_x=0.0,beta_y=0.0;
                for(int i:clusters[j]){
                    
                    double theta_x=(2*M_PI)*particles[i].x/Lx;
                    double theta_y=(2*M_PI)*particles[i].y/Ly;
                    
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
                    double dx=minimum_image(particles[i].x - x_cm,Lx);
                    double dy=minimum_image(particles[i].y - y_cm,Ly); 
                    rg_sq+=hypot(dx,dy);
                }
            
                rg_sq/=clusters[j].size();
                
                radius.push_back(sqrt(rg_sq));
            }
        }
        return radius;
    }
    
    void save_cluster_index(vector<vector<int>> clusters,int trial,int step){
        vector<int> cluster_index(N,0);
        int totalclusters=1;
        for(auto row:clusters){
            for(int element:row){
                cluster_index[element]=totalclusters;
            }
            totalclusters++;
        }
        ofstream f(folder_path+"cluster_data/clusterindex_"+to_string(trial)+"_"+to_string(time_array[step])+"_.txt");
        for(int a:cluster_index){f<<a<<"\n";}
        f.close();
    }


    void save_rg(vector<double>radius,int trial,int t){
        ofstream fi(folder_path+"cluster_data/radius_"+to_string(trial)+"_"+to_string(time_array[t])+"_.txt");
        for(auto a:radius){fi<<a<<"\n";}
        fi.close();
    }
    void save_mass(vector<vector<int>> clusters,int trial,int t){
        ofstream fil(folder_path+"cluster_data/mass_"+to_string(trial)+"_"+to_string(time_array[t])+"_.txt");
        vector <double> mass;
        for(auto row:clusters){if(row.size()>=cluster_cutoff){mass.push_back(row.size());}}
        for(auto element:mass){fil<<element<<"\n";}
        fil.close();
    }
    void start_calculation() {
        initialize_time_to_record();
        vector <double>r_t(time_array.size(),0);
        for (int trial = trial_start; trial < numberoftrials; trial++){
            cout<<"\n"<<"Trial number "<<trial<<fixed<<setprecision(2)<<" || Packing Fraction : "<<(N/(Lx*Ly))<<" | Noise = "<<noise<<" | Angle = "<<180*(half_angle/M_PI)<<" | N = "<<N<<" || "<<endl ;
            double timestarted=static_cast <double>(time(NULL));    
            for (int t=0;t<time_array.size();t++){
                load_snapshot(time_array[t],trial);
                calculate_neighbours();
                vector<vector<int>> clusters=BFS();   
                save_cluster_index(clusters,trial,t);
                save_mass(clusters,trial,t);
                
                vector<double> radius =calculate_rg(clusters);
                save_rg(radius,trial,t);
                double r__=0;
                for(double i:radius)r__+=i;
                r_t[t]+=r__/radius.size();

                if (time_array[t] % 500 == 0){print_progress(static_cast<double>(time_array[t])/static_cast<double>(tmax),timestarted);}                
                } 
                
            
            print_progress(1.0,static_cast<double>(timestarted));
        } 

        ofstream fill(folder_path+"/cluster_index/radius_vs_t_.txt");
        for(int i=0;i<time_array.size();i++){fill<<time_array[i]<<" "<<r_t[i]/numberoftrials<<"\n";}
        fill.close();
        cout << "\nSimulation complete. Recorded " << time_array.size() << " snapshots." << endl;
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
