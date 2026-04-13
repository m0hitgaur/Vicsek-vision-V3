#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <format>
#include <thread>
#include <filesystem>
#include <iomanip>
#include "locker.h"
using namespace std;
namespace fs = filesystem;

// standard Vicsek model
bool create_directory(const string& path) {
    try {
        fs::create_directories(path);
        return fs::is_directory(path);
    } 
    catch (const fs::filesystem_error& e) {
        cerr << "Error creating directory '" << path << "': " << e.what() << '\n';
        return false;
    }
}
struct Particle {
    double x, y;
    double vx, vy;
    double ax, ay;
    double theta_avg;
    double x_new,y_new;
    double vx_new,vy_new;
    vector <int> neighbours;
    double sigma=0.5;
};

class Simulation {
private:
    const int N=640;                     // Number of particles
    const double Lx=16.0;                 // Box size
    const double Ly=16.0;                 //  Box size
    const double dt=1.0;                 // Timestep
    const double v0=1.0e-2;                 // Magnitude of velocity
    const int long_neighbours=4;     // Number of long range neighbours
    const double noise;              // Noise strength
    const double half_angle;         // Half of vision angle
    const double rc=1;                 // Cutoff radius
    const double beta=1.0;              // Aligning strength   
    const int trial;                 // Trial number    
    const int tmax;                  // Total nujmber of time steps  
    string folder_path;              // Directory for saving data 
    vector<bool> time_record;        // Time points to record
    vector<int> times;               // Time points recorded 
    vector<double> order;            // Order parameter data
    vector<Particle> particles;
    mt19937 gen; 
    uniform_real_distribution<double> rng_uniform_symm{-1.0, 1.0};
    uniform_real_distribution<double> rng_uniform_one{0, 1};
    normal_distribution<double> rng_normal_unity{0.0, 1.0};

public:
    Simulation(double noise_input,double angle,int trial_input,int tmax_input,int seed)
               :noise(noise_input),half_angle(angle),trial(trial_input),tmax(tmax_input){
        particles.resize(N);
        initialize_particles();
        initialize_time_to_record();
        folder_path="data/";
        create_directory(folder_path + "order_data");
        create_directory(folder_path+"config_data/trial_"+ to_string(trial)+"/");
        
    }

    vector<double> get_order_data(){
        return order;}
    vector<int> get_time_data(){
        return times;}

    void update_neigbours(){      
        for(Particle &p:particles){
            p.neighbours.clear();            
            for(int i=0;i<particles.size();i++){
                double r=rij(p,particles[i]);
                if(r<=2*rc ) {p.neighbours.push_back(i);}
            }
        }
    }
    void initialize_particles() {
        gen.seed(12345 + 10 * trial);
        int grid_size = static_cast<int>(ceil(sqrt(N)));
        double spacing_x = Lx / grid_size;
        double spacing_y = Ly / grid_size;
        int i=0;
        for(int j=0;j<static_cast<int>(grid_size) &&i<N;j++){
            for(int k=0;k<static_cast<int>(grid_size)&&i<N;k++){

                particles[i].x = (j ) * spacing_x;
                particles[i].y = (k ) * spacing_y;
                particles[i].x_new = particles[i].x;
                particles[i].y_new = particles[i].y;
                
                double theta=M_PI*rng_uniform_symm(gen);
                particles[i].vx = v0*cos(theta);
                particles[i].vy = v0*sin(theta);
                particles[i].vx_new = 0;
                particles[i].vy_new = 0;
                
                particles[i].ax = 0;
                particles[i].ay = 0;
                
                particles[i].theta_avg = 0.0;
                i++;                
            } 
        }
        
        for(int t=0;t<1500;t++){
            if(t%50==0)update_neigbours();
            velocity_alignment(M_PI,5);
            velocity_update();
            position_update();
            EndTimeStep();
            }
        
    }          
    double dot_product(double theta,double dx,double dy,double rij){
        return ( (cos(theta) * (dx)) + (sin(theta) * (dy)) )/(rij);
    }
    void pbc_position(Particle & p){
        // Periodic boundary conditions
        if (p.x_new < 0) p.x_new = fmod(p.x_new,Lx) +Lx;
        if (p.x_new > Lx) p.x_new = fmod(p.x_new,Lx);
        if (p.y_new < 0) p.y_new = fmod(p.y_new,Ly) +Ly;
        if (p.y_new > Ly) p.y_new = fmod(p.y_new,Ly);
    }
    double minimum_image(double dx,double L){
        if(dx>L/2) dx=dx-L;
        if(dx<-L/2) dx=dx+L ;  
        return dx;
    }

     double rij( Particle& p_i,  Particle& p_j) {
        double dx = p_i.x - p_j.x;
        double dy = p_i.y - p_j.y;
        dx=minimum_image(dx,Lx);
        dy=minimum_image(dy,Ly);
        double r=hypot(dx,dy);
        if (r < 1e-5) r = 1e-5;
        return r;
    }   

    void velocity_alignment(double alpha,double eta) {    
        vector<double> avg_vx(N,0);
        vector<double> avg_vy(N,0);
        double newtheta;
        vector<int> count_v(N,1); // starting from 1 as the particle itself is always counted
        
        for (int i=0;i<particles.size();i++) {
            double theta_i=atan2(particles[i].vy,particles[i].vx);
            avg_vx[i]+=cos(theta_i);   // just so that we can include the particle itself in average
            avg_vy[i]+=sin(theta_i);   // just so that we can include the particle itself in average
            
            for (int j:particles[i].neighbours) 
            {   
                double theta_j=atan2(particles[j].vy,particles[j].vx);
                double dx=minimum_image(particles[j].x - particles[i].x,Lx);
                double dy=minimum_image(particles[j].y - particles[i].y,Ly);                 
                double rij = hypot(dx,dy);
                double innerproduct_i=  dot_product(theta_i,dx,dy,rij); 
                double innerproduct_j=  dot_product(theta_j,-dx,-dy,rij); 
                        
                if( rij<=rc && innerproduct_i >= cos(alpha) )      
                    {avg_vy[i]+= sin(theta_j);
                    avg_vx[i]+=cos(theta_j);
                    count_v[i]++;}
                
                if( rij<=rc && innerproduct_j >= cos(alpha) ) 
                    {avg_vy[j]+= sin(theta_i);
                    avg_vx[j]+=cos(theta_i);
                    count_v[j]++;}
            }
            
        if(count_v[i]!=0)
        {avg_vx[i] /= static_cast<double>(count_v[i]);
        avg_vy[i]/=static_cast<double>(count_v[i]);}
        
        double totalforce_x=beta*avg_vx[i] ;
        double totalforce_y=beta*avg_vy[i] ;
        
        newtheta=atan2(totalforce_y,totalforce_x) + (rng_uniform_symm(gen))*(eta/2);
        
        if(newtheta<-M_PI)newtheta= fmod(newtheta , M_PI)+M_PI;
        else if(newtheta>M_PI)newtheta= fmod(newtheta , M_PI)-M_PI;
        
        particles[i].theta_avg = newtheta;
        }
            
        }
    void velocity_update(){
        for (Particle & p : particles) {    
            p.vx_new = v0* cos(p.theta_avg) ;
            p.vy_new = v0* sin(p.theta_avg) ; 
        }
    }
    void position_update(){ 
        for (Particle & p : particles) {    
            p.x_new += p.vx * dt ;
            p.y_new += p.vy * dt ;      
            pbc_position(p);
        }
    }
    void EndTimeStep(){
       for (Particle &p : particles) {
            p.vx = p.vx_new;
            p.vy = p.vy_new;
            p.x = p.x_new;
            p.y = p.y_new;      
            }         
    } 

    double velocity_order_parameter() {
        double sum_vx = 0, sum_vy = 0;
        
        for ( Particle & p : particles) {
            sum_vx += p.vx;         
            sum_vy += p.vy;    
        }
        
        double va = hypot(sum_vx,sum_vy) /(v0* N);
        return va;
    }   
    void save_snapshot(int step,int trial) {
        ofstream file(folder_path+"config_data/trial_"+to_string(trial)+"/config_" +to_string(step) + ".csv");
        
        file<<"x,y,vx,vy\n";
        for ( Particle & p : particles) {
            file << p.x << "," << p.y << ","
                 << p.vx << "," << p.vy << "\n";
               
        }
        
        file.close();
    }  
    
    void save_order_data(){
        ofstream order_file(folder_path+"order_data/order_parameter_"+to_string(trial)+"_.csv");
        order_file<<"va,t\n";
        for (int i=0;i<order.size();i++)order_file<<order[i]<<","<<times[i]<<"\n";      
        order_file.close();

    }

    void initialize_time_to_record(){
        for (int t = 0; t < tmax; t++) {
            bool should_record = false;
                if (t <= 10) should_record = true;
                if(t>=10 && t<100&& t%10==0) should_record = true;                   
                if (t<1000 && t >= 100 && t % 50 == 0) should_record = true;   
                if (t>=1000 && t % 100 == 0) should_record = true;                        
                time_record.push_back(should_record);
                }
    }

    
    void start_run() {
        ofstream f(folder_path+"parameters.csv");
        string head="N,Lx,Ly,half_angle,noise,v0,dt,maxiter,trial\n";
        f<< head;
        f<<N<<","<<Lx<<","<<Ly<<","<<half_angle<<","<<noise<<","<<v0
        <<","<<dt<<","<<tmax<<","<<trial; 
        f.close();
        cout<<"\n"<<"Trial number "<<trial<<fixed<<setprecision(2)<<" || Packing Fraction : "<<(N/(Lx*Ly))<<" | Noise = "<<noise<<" | Angle = "<<180*(half_angle/M_PI)<<" | N = "<<N<<" || "<<endl ;
        
        double timestarted=static_cast <double>(time(NULL));
        int time_counter=0;
        int time_update_neighbours=50;//static_cast<int>(1/(2*v0));

        for (int t = 0; t < tmax; t++) {             
            if(t%time_update_neighbours==0) update_neigbours();
            if (time_record[t]){
                order.push_back(velocity_order_parameter());    
                save_snapshot(t,trial);
                times.push_back(t);
                } 
            if (t % 500 == 0){print_progress(static_cast<double>(t)/static_cast<double>(tmax),timestarted);}
            
            velocity_alignment(half_angle,noise);
            velocity_update();
            position_update(); 
            EndTimeStep();

        } 
        
        print_progress(1.0,static_cast<double>(timestarted));
        save_order_data();
        cout << "\nSimulation complete. Recorded " << times.size() << " snapshots." << endl;
    }

};




int main() { 
    double noise=1.0e-1 ;     // strength of noise
    double half_angle=M_PI;   // Half of the vision angle  
    int tmax = 1.0e3;         // Maximum time
    int numberoftrials=1;     // Number of trials
    int trialstart=0;         // Starting trial number 
    int seed=12345;           // random seed
    time_t trial_time,start_time=time(NULL) , finish_time;
    for(int trial=trialstart;trial<numberoftrials;trial++){ 
        trial_time=time(NULL);
        Simulation sim(noise,half_angle,trial,tmax,seed);
        sim.start_run();
        cout<<"Time to calculate trial = "  <<time(NULL)-trial_time<<" seconds ";  
        
    }   
    
    cout<<"\n"<<"Total time elapsed : "<< time(NULL) - start_time <<" seconds ";
    
    return 0;
}