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



    

    // Start time tracking
    time_t start_time = time(NULL);
    int numberoftimes=0;
    for (int t = 0; t < number_of_timesteps; t++) {
        if(t<10) time_interval=1;
        if(t>10) time_interval=10;
        if(t>100) time_interval=50;
        if(t>1000)time_interval=100;
        if(t%time_interval==0)numberoftimes++;  
    }
    
    vector<vector<double>> bindercum_2(number_of_trials, vector<double>(numberoftimes, 0.0)),bindercum_4(number_of_trials, vector<double>(numberoftimes, 0.0)),firstcum(number_of_trials, vector<double>(numberoftimes, 0.0)),secondcum(number_of_trials, vector<double>(numberoftimes, 0.0)),thirdcum(number_of_trials, vector<double>(numberoftimes, 0.0)),fourthcum(number_of_trials, vector<double>(numberoftimes, 0.0));
    vector<double> time_U2;
    // Main trial loop
    for (int trial = 0; trial < number_of_trials; trial++) {
        cout << "Trial number " << trial << endl;
        int time_step_counter = 0;
        time_t start_time_trial = time(NULL);
        ofstream file1(folder_path+"/Binder_Cumulant/U2_vs_t_"+to_string(trial)+".dat"),file2(folder_path+"/Binder_Cumulant/U4_vs_t_"+to_string(trial)+".dat");
        int time_array=0;
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
                double vy_firstcum = 0.0,vy_secondmom=0.0,vy_thirdmom=0.0,vy_fourthmom=0.0;
                for (int i = 0; i < number_of_agents; i++) {
                    double vy_temp= pow(sin(theta[i]),1)  ;
                    vy_firstcum+=vy_temp;
                    vy_secondmom+=vy_temp*vy_temp;
                    vy_thirdmom+=vy_temp*vy_temp*vy_temp;
                    vy_fourthmom+=vy_temp*vy_temp*vy_temp*vy_temp;
    
                }
                vy_firstcum/=number_of_agents;
                vy_secondmom/=number_of_agents;
                vy_thirdmom/=number_of_agents;
                vy_fourthmom/=number_of_agents;
                
                double vy_secondcum= (vy_secondmom) - pow(vy_firstcum,2);
                double vy_thirdcum= (vy_thirdmom) - (3*(vy_secondmom*vy_firstcum)) + (2*pow(vy_firstcum,3));
                double vy_fourthcum= (vy_fourthmom) - (4*(vy_thirdmom*vy_firstcum)) - (3*pow(vy_secondmom,2)) + (12*vy_secondmom*pow(vy_firstcum,2)) - (6*pow(vy_firstcum,4));                
                
                double vx_firstcum = 0.0,vx_secondmom=0.0,vx_thirdmom=0.0,vx_fourthmom=0.0;
                for (int i = 0; i < number_of_agents; i++) {
                    double vx_temp= pow(cos(theta[i]),1)  ;
                    vx_firstcum+=vx_temp;
                    vx_secondmom+=vx_temp*vx_temp;
                    vx_thirdmom+=vx_temp*vx_temp*vx_temp;
                    vx_fourthmom+=vx_temp*vx_temp*vx_temp*vx_temp;
    
                }
                vx_firstcum/=number_of_agents;
                vx_secondmom/=number_of_agents;
                vx_thirdmom/=number_of_agents;
                vx_fourthmom/=number_of_agents;
                if(trial==99)cout<<"\n"<<"firstmom= "<<vx_firstcum<<" "<<"secmom= "<<vx_secondmom<<"3rdmom= "<<vx_thirdmom<<" "<<"4thmom= "<<vx_fourthmom<<"\n";
                double vx_secondcum= (vx_secondmom) - pow(vx_firstcum,2);
                double vx_thirdcum= (vx_thirdmom) - (3*(vx_secondmom*vx_firstcum)) + (2*pow(vx_firstcum,3));
                double vx_fourthcum= (vx_fourthmom) - (4*(vx_thirdmom*vx_firstcum)) - (3*pow(vx_secondmom,2)) + (12*vx_secondmom*pow(vx_firstcum,2)) - (6*pow(vx_firstcum,4));                
                if(trial==99)cout<<"firstcumm= "<<vx_firstcum<<" "<<"seccum= "<<vx_secondcum<<"3rdcum= "<<vx_thirdcum<<" "<<"4thcum= "<<vx_fourthcum<<"\n";
 
                firstcum[trial][time_array]=vx_firstcum;
                secondcum[trial][time_array]=vx_secondcum;
                thirdcum[trial][time_array]=vx_thirdcum;
                fourthcum[trial][time_array]=vx_fourthcum;
                //calculating fourth order binder cumulant U4
                double bindercum_4_temp= 1- ((vx_fourthmom)/(3*vx_secondmom*vx_secondmom));
                bindercum_4[trial][time_array]=bindercum_4_temp;                

                //calculating second order binder cumulant U2
                
                double bindercum_2_temp= vx_secondcum;
                if(trial==0)time_U2.push_back(dt*t);
                bindercum_2[trial][time_array]=bindercum_2_temp;
                
               
                file1<<t*dt<<" "<<bindercum_2_temp<<"\n";
                file2<<t*dt<<" "<<bindercum_4_temp<<"\n";
                

            time_array++;    
            }

            
           
        }
        file1.close();
        file2.close();
        time_t finish_time_trial = time(NULL);
        cout << "\nTime taken to calculate trial : " << finish_time_trial - start_time_trial << " seconds" << endl;
    }



    vector <double> U2(time_U2.size(),0);
    ofstream file_o(folder_path+"/Binder_Cumulant/U2_vs_t.dat");
    for(int t=0;t<time_U2.size();t++){
        for(int trial=0;trial<number_of_trials;trial++){
            U2[t]+=bindercum_2[trial][t];
            
        }
        U2[t]/=number_of_trials;
    }

    for(int t=0;t<time_U2.size();t++){
        file_o<<time_U2[t]<<" "<<U2[t]<<"\n";
    }
    file_o.close();
 
    vector <double> U1(time_U2.size(),0);
    ofstream file_o1(folder_path+"/Binder_Cumulant/cum1_vs_t.dat");
    for(int t=0;t<time_U2.size();t++){
        for(int trial=0;trial<number_of_trials;trial++){
            U1[t]+=firstcum[trial][t];
            
        }
        U1[t]/=number_of_trials;
    }

    for(int t=0;t<time_U2.size();t++){
        file_o1<<time_U2[t]<<" "<<U1[t]<<"\n";
    }
    file_o1.close();    


    vector <double> cum2(time_U2.size(),0);
    ofstream file_o2(folder_path+"/Binder_Cumulant/cum2_vs_t.dat");
    for(int t=0;t<time_U2.size();t++){
        for(int trial=0;trial<number_of_trials;trial++){
            cum2[t]+=secondcum[trial][t];
            
        }
        U1[t]/=number_of_trials;
    }

    for(int t=0;t<time_U2.size();t++){
        file_o2<<time_U2[t]<<" "<<cum2[t]<<"\n";
    }
    file_o2.close();   


    vector <double> cum3(time_U2.size(),0);
    ofstream file_o3(folder_path+"/Binder_Cumulant/cum3_vs_t.dat");
    for(int t=0;t<time_U2.size();t++){
        for(int trial=0;trial<number_of_trials;trial++){
            cum3[t]+=thirdcum[trial][t];
            
        }
        cum3[t]/=number_of_trials;
    }

    for(int t=0;t<time_U2.size();t++){
        file_o3<<time_U2[t]<<" "<<cum3[t]<<"\n";
    }
    file_o3.close();       
    
    vector <double> cum4(time_U2.size(),0);
    ofstream file_o4(folder_path+"/Binder_Cumulant/cum4_vs_t.dat");
    for(int t=0;t<time_U2.size();t++){
        for(int trial=0;trial<number_of_trials;trial++){
            cum4[t]+=fourthcum[trial][t];
            
        }
        cum4[t]/=number_of_trials;
    }

    for(int t=0;t<time_U2.size();t++){
        file_o4<<time_U2[t]<<" "<<cum4[t]<<"\n";
    }
    file_o4.close();       


    
    
    vector <double> U4(time_U2.size(),0);
    ofstream file_oo(folder_path+"/Binder_Cumulant/U4_vs_t.dat");
    for(int t=0;t<time_U2.size();t++){
        for(int trial=0;trial<number_of_trials;trial++){
            U4[t]+=bindercum_4[trial][t];
            
        }
        U4[t]/=number_of_trials;
    }

    for(int t=0;t<time_U2.size();t++){
        file_oo<<time_U2[t]<<" "<<U4[t]<<"\n";
    }
    
    
    // Finish
    time_t finish_time = time(NULL);
    cout << "\nTotal time elapsed: " << finish_time - start_time << " seconds" << endl;

    return 0;
}


