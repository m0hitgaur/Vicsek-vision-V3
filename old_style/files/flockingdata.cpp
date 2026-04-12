#define _USE_MATH_DEFINES
//11-12-2024
#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <string>      
#include <filesystem>
#include <vector>
using namespace std;

namespace fs = std::filesystem;
time_t trial_time,start_time=time(NULL) , finish_time;


int seedvalue;
int    Number_of_active_agents ; // number of particles
double Lx ;                        // length of the system
double Ly ;                        // length of the system
int    maxiter ;                  // total number of time steps 
double velocity ;                // magnitude of velocity 
double dt ;                       // time step 
double half_angle;
double noise;
double density ;  // density
double radius ;  // radius of interaction 
int numberoftrials; 
int trialstart;

random_device rd; 
mt19937 gen(rd());
uniform_real_distribution<double> angle_dis(-noise/2,noise/2);


void visualize( vector <vector<double>>& position,vector<double>& theta,int Number_of_active_agents,double half_angle,double radius,double noise,double velocity,double dt,double Lx,double Ly,int time ,int trial)
{
    string folder_path = fs::current_path().string();
       
    string theta_filename = folder_path+"/flockingdata/theta_"+ to_string(trial) +"_"+to_string(int(time*dt))+ "_.dat";
    string posix_filename = folder_path+"/flockingdata/positionx_" + to_string(trial) +"_"+to_string(int(time*dt))+"_.dat";
    string posiy_filename = folder_path+"/flockingdata/positiony_" + to_string(trial) +"_"+to_string(int(time*dt))+ "_.dat";

    ofstream theta_file(theta_filename);
    ofstream posix_file(posix_filename);
    ofstream posiy_file(posiy_filename);

    if (!theta_file.is_open() || !posix_file.is_open() || !posiy_file.is_open()) 
        cerr << "Failed to open files for writing" << endl;



    for (int i = 0; i < Number_of_active_agents; i++) 
    {   theta_file << theta[i] << "\n";
        posix_file << position[0][i] << "\n";
        posiy_file << position[1][i] << "\n";
    }

    theta_file.close();
    posix_file.close();
    posiy_file.close();



}


double orderparameter(vector <vector<double>>& position,vector<double>& theta,int Number_of_active_agents,double half_angle,double radius,double noise,double velocity,double dt,double Lx,double Ly,int time ,int trial)
  {

    double avg_vx=0.0,avg_vy=0.0,v=0.0 ;
    for(int j=0 ; j<Number_of_active_agents ; j++)
    {
    avg_vx+=cos(theta[j]);
    avg_vy+=sin(theta[j]);

    }
    avg_vx/=Number_of_active_agents ;
    avg_vy/=Number_of_active_agents ;
    v=sqrt((avg_vx*avg_vx) + (avg_vy*avg_vy));
    return v;
    
  }


int main() {
    string folder_path = fs::current_path().string();
    ifstream f(folder_path+"/parameters.txt");
    if (!f) {
        cerr << "Error opening file" << endl;
        return 1;}
    string line;
    int i=0;
    vector<string> data;
    while (getline(f, line)) {
        data.push_back(line) ;
    i++;
    }


    int    Number_of_active_agents = stoi(data[0]) ; // number of particles
    double Lx = stod(data[1]);                        // length of the system
    double Ly = stod(data[2]);                        // length of the system
    int    maxiter = stoi(data[3]);                  // total number of time steps 
    double velocity = stod(data[4]);                // magnitude of velocity 
    double dt = stod(data[5]);                       // time step 
    double half_angle=stod(data[6]);
    double noise=stod(data[7]);
    double density = stod(data[8]);  // density
    double radius = stod(data[9]);  // radius of interaction 
    int numberoftrials=stoi(data[10]); 
    int trialstart=stoi(data[11]);
    int seedvalue=stoi(data[12]);
    int timeinterval_flocking=stoi(data[13]);
    double w=0.5;
    cout<<"Initializing variables.........."<<"\n";
    cout<<"Density"<<"\n";
    cout<<"Noise= "<<noise<<"\n";
    cout<<"Velocity= "<<velocity<<"\n";
    cout<<"Time step= "<<dt<<"\n";
    cout<<"Length of system = "<<Lx<<"*"<<Ly<<"\n";
    cout<<"Number of trials = "<<numberoftrials<<"\n";
    cout<<"Total number of timesteps = "<<maxiter<<"\n";


    
  
    for(int trial=trialstart;trial<numberoftrials;trial++)
    {   trial_time=time(NULL);
        cout<<"\n"<<"Trial number : "<<trial<< " Out of "<<numberoftrials<<"    ";
        
        vector<vector<double>> position(2, vector<double>(Number_of_active_agents,0));vector<vector<double>> newposition(2, vector<double>(Number_of_active_agents,0));
        vector<double> theta(Number_of_active_agents,0); vector<double> newtheta(Number_of_active_agents,0); 
        vector<vector<double>> rij(Number_of_active_agents, vector<double>(Number_of_active_agents,0));
        vector<double> order,ordertime;
        //initialize
        srand(seedvalue);   
        int max=10000;
        vector <int> checkx(max,0),checky(max,0);
        if(max<Number_of_active_agents)cout<<"Number of active agents more than initializable";
        for (int j = 0; j < Number_of_active_agents; j++)
        {   
            int flagx=0,flagy=0;      
            while(flagx==0)
            {
                int randx=rand()%max;
                if(checkx[randx]==0) 
                {   position[0][j] = (double(randx)/max)*Lx;
                    checkx[randx]++;
                    flagx++;
                }
            } 
                
            while(flagy==0)
            {
                int randy=rand()%max;
                if(checky[randy]==0)
                { 
                    position[1][j] = (double(randy)/max)*Ly;
                    checky[randy]++;
                    flagy++;
                }
            }

            theta[j] = ((double(rand()%100)/100)*2*M_PI)-M_PI;
        }
            seedvalue+=10;
        //end initialize

        
        for (int time = 0; time < maxiter; time++)
        {   
            if(time<10) timeinterval_flocking=1;
            if(time>10) timeinterval_flocking=10;
            if(time>100) timeinterval_flocking=50;
            if(time>1000) timeinterval_flocking=100;

            if(time%100==0)cout<<time<<">>";
            if(time%timeinterval_flocking==0)
            {   
                double avg_vx=0.0,avg_vy=0.0,v=0.0 ;
                for(int j=0 ; j<Number_of_active_agents ; j++)
                {
                avg_vx+=cos(theta[j]);
                avg_vy+=sin(theta[j]);

                }
                avg_vx/=Number_of_active_agents ;
                avg_vy/=Number_of_active_agents ;
                v=sqrt((avg_vx*avg_vx) + (avg_vy*avg_vy));
                
                order.push_back(v);
                visualize(position,theta,Number_of_active_agents,half_angle,radius,noise,velocity,dt,Lx,Ly,time, trial);
                ordertime.push_back(time);

            }
               
        //update theta implemented from paper 
            vector<double> avg_vx(Number_of_active_agents,0);
            vector<double> avg_vy(Number_of_active_agents,0);
            vector<int> count_v(Number_of_active_agents,0); // starting from 1 as the particle itself is always counted
            vector<double> avg_sx(Number_of_active_agents,0);
            vector<double> avg_sy(Number_of_active_agents,0);            
            vector<int> count_s(Number_of_active_agents,0);
            for (int i = 0; i < Number_of_active_agents; i++) {
                //avg_vx[i]+=cos(theta[i]);   // just so that we can include the particle itself in average
                //avg_vy[i]+=sin(theta[i]);   // just so that we can include the particle itself in average
        
                for (int j = i+1; j < Number_of_active_agents; j++) 
                {   double dx=position[0][j] - position[0][i];
                    double dy=position[1][j] - position[1][i];
                    if(dx>Lx/2) dx=dx-Lx;
                    if(dx<-Lx/2) dx=dx+Lx ;  
                    if(dy>Ly/2) dy=dy-Ly;
                    if(dy<-Ly/2) dy=dy+Ly ;                 

                    rij[i][j] = sqrt(pow(dx, 2) + pow(dy, 2));
                    rij[j][i] = rij[i][j];
                    double innerproduct_i=((cos(theta[i]) * (dx)) + (sin(theta[i]) * (dy)) )/(rij[i][j]); 
                    double innerproduct_j= -1 * ( (cos(theta[j])*(dx) )+(sin(theta[j]) * (dy)))/(rij[j][i]); 
                            
                    if(rij[i][j] <= radius && innerproduct_i >= cos(half_angle) ) 
                        {avg_vy[i]+= sin(theta[j]);
                        avg_vx[i]+=cos(theta[j]);
                        count_v[i]++;}
                    
                    if(rij[j][i] <= radius && innerproduct_j >= cos(half_angle) ) 
                        {avg_vy[j]+= sin(theta[i]);
                        avg_vx[j]+=cos(theta[i]);
                        count_v[j]++;}

                                        
                    if(count_s[i]<7 && rij[i][j] > radius && rij[i][j] <= 2*radius && innerproduct_i >= cos(half_angle) ) 
                        {avg_sy[i]+= sin(theta[j]);
                        avg_sx[i]+=cos(theta[j]);
                        count_s[i]++;}
                    
                    if(count_s[j]<7 && rij[i][j] > radius &&rij[j][i] <= 2*radius && innerproduct_j >= cos(half_angle) ) 
                        {avg_sy[j]+= sin(theta[i]);
                        avg_sx[j]+=cos(theta[i]);
                        count_s[j]++;}                    

                }
                
            if(count_v[i]!=0)
            {avg_vx[i] /= static_cast<double>(count_v[i]);
            avg_vy[i]/=static_cast<double>(count_v[i]);}

            if(count_s[i]!=0)
            {avg_sx[i] /= static_cast<double>(count_s[i]);
            avg_sy[i]/=static_cast<double>(count_s[i]);}
            
            double Sx=w*avg_vx + (1-w)avg_sx;
            double Sy=w*avg_vy + (1-w)avg_sy;

            newtheta[i] = acos(((Sx*cos(theta[i]))+(Sy*sin(theta[i])))/(velocity)) + ((double(rand()%1000)/1000)*(noise))-(noise/2);
            if(newtheta[i]<-M_PI)newtheta[i]= fmod(newtheta[i] , M_PI)+M_PI;
            else if(newtheta[i]>M_PI)newtheta[i]= fmod(newtheta[i] , M_PI)-M_PI;

            }
            //end update theta
            
            //update position
            for (int i = 0; i < Number_of_active_agents; i++) {
                newposition[0][i] = position[0][i] + velocity * dt * cos(theta[i]);
                newposition[1][i] = position[1][i] + velocity * dt * sin(theta[i]);
                if(newposition[0][i]<0)newposition[0][i]= Lx+fmod(newposition[0][i] , Lx);
                else if(newposition[0][i]>Lx)newposition[0][i]= fmod(newposition[0][i] , Lx);
                if(newposition[1][i]<0) newposition[1][i]= Ly+fmod(newposition[1][i] , Ly);
                else if(newposition[1][i]>Ly) newposition[1][i]= fmod(newposition[1][i] , Ly);

            }
            //end update position
        //reinitializing 
            
            theta = newtheta;
            position = newposition;


        
        }
        cout<<maxiter<<"\n";    
        cout<<"Time to calculate trial = "  <<time(NULL)-trial_time<<" seconds  for   Angle : "+to_string(half_angle*180/M_PI)+" | Noise : "+to_string(noise)+" | Density : "+to_string(density)+" | N = "+to_string(Number_of_active_agents)<<endl;  
         //pushing out order data

        ofstream orderpara(folder_path+"/orderdata/order_vs_time_"+to_string(trial)+"_.txt");
        for(int j=0;j<int(order.size());j++) orderpara<<dt*ordertime[j]<<" "<<order[j]<<"\n"; 
            
        orderpara.close();
        order.clear();
        ordertime.clear();
    }
    ofstream detail(folder_path+"/flockingdata/details.txt");
    detail<<Number_of_active_agents<<"\n"<<Lx<<"\n"<<Ly<<"\n"<<half_angle<<"\n"<<noise<<"\n"<<density<<"\n"<<dt<<"\n"<<maxiter<<"\n"<<numberoftrials; 
    detail.close();
    

        
    
    finish_time=time(NULL);
    cout<<"\n"<<"Total time elapsed : "<< finish_time - start_time <<" seconds ";
    
}


