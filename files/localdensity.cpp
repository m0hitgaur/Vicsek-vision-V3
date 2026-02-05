#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <filesystem>
#include <iomanip>

namespace fs = std::filesystem;
using namespace std;

int main() {
    string folder_path = fs::current_path().string();

    vector<string> data;
    ifstream details(folder_path + "/densityparameters.txt");
    string line;
    while (getline(details, line)) {
        data.push_back(line);
    }
    time_t start_time=time(NULL) , finish_time;

    int    numberofagents = stoi(data[0]) ;          // number of particles
    double Lx = stod(data[1]);                        // length of the system
    double Ly = stod(data[2]);                        // length of the system
    int    Numberoftimesteps = stoi(data[3]);                  // total number of time steps 
    double velocity = stod(data[4]);                // magnitude of velocity 
    double dt = stod(data[5]);                       // time step 
    double Halfangle=stod(data[6]);
    double noise=stod(data[7]);
    double density = stod(data[8]);  // density
    double radius = stod(data[9]);  // radius of interaction 
    double numberoftrials=stoi(data[10]); 
    int    gridsize=stoi(data[11]); // bin size for r
    int    timeinterval=stoi(data[12]);



    
    float binsize = static_cast<float>(Lx) / gridsize;
   

    cout << "Length of system = " << Lx << " * " << Ly << endl;
    cout << "Number of time steps = " << Numberoftimesteps << endl;
    cout << "Number of agents = " << numberofagents << endl;
    cout << "Density = " << density << endl;
    cout << "Noise = " << noise << endl;
    cout << "Total number of trials = " << numberoftrials << endl;

    for (int trial = 0; trial < numberoftrials; trial++) {
        cout <<"Trial number "<<trial<<endl<<"Creating empty lists to store file names.........." << endl;

        cout << "Reading files.........." << endl;
        for (int time = 0; time < Numberoftimesteps; time++) {
            if(time<10)timeinterval=1;
            if(time>10)timeinterval=10;
            if(time>100)timeinterval=100;
            if(time>1000)timeinterval=200;
            if(time>1500)timeinterval=500;            
            if (time % timeinterval == 0) {    
                vector<double> positionx;
                vector<double> positiony;
                vector<double> theta;
                string theta_filename = folder_path+"/flockingdata/theta_" + to_string(trial) + "_" + to_string(static_cast<int>(time * dt)) + "_.dat";
                string posix_filename = folder_path+"/flockingdata/positionx_" + to_string(trial) + "_" + to_string(static_cast<int>(time * dt)) + "_.dat";
                string posiy_filename = folder_path+"/flockingdata/positiony_" + to_string(trial) + "_" + to_string(static_cast<int>(time * dt)) + "_.dat";

                ifstream theta_file(theta_filename);
                ifstream posix_file(posix_filename);
                ifstream posiy_file(posiy_filename);

                if (!theta_file.is_open() || !posix_file.is_open() || !posiy_file.is_open()) {
                    cerr << "Failed to open files for reading" << endl;
                    return 1;
                }
                // Reinitialize variables
                vector<vector<double>> localorderx(gridsize, vector<double>(gridsize, 0));
                vector<vector<double>> localordery(gridsize, vector<double>(gridsize, 0));
                vector<vector<double>> localorder(gridsize, vector<double>(gridsize, 0));
                vector<vector<int>> localordercount(gridsize, vector<int>(gridsize, 0));
                vector<vector<double>> localdensity(gridsize, vector<double>(gridsize, 0));
                
                for (int i = 0; i < numberofagents; i++) 
                {
                    string t, x, y;
                    theta_file >> t;
                    posix_file >> x;
                    posiy_file >> y;
                    positionx.push_back(stod(x));
                    positiony.push_back(stod(y));
                    theta.push_back(stod(t));
                }


                for (int i = 0; i < numberofagents; i++) 
                {   int xi = static_cast<int>(positionx[i] / binsize)%gridsize;
                    int yi = static_cast<int>(positiony[i] / binsize)%gridsize;
                    localorderx[xi][yi] += cos(theta[i]);
                    localordery[xi][yi] += sin(theta[i]);
                    localordercount[xi][yi]++;
                }
                for (int i = 0; i < gridsize; i++) 
                {
                    for (int j = 0; j < gridsize; j++) 
                    {
                        if (localordercount[i][j] != 0) 
                        {   localorderx[i][j] /= localordercount[i][j];
                            localordery[i][j] /= localordercount[i][j];
                            localorder[i][j] = sqrt(pow(localorderx[i][j], 2) + pow(localordery[i][j], 2));
                            localdensity[i][j] = localordercount[i][j] / (binsize * binsize);
                        }
                    }
                }

                string filename1 = folder_path + "/density_heatmaps/localorderimage_" +to_string(trial)+"_"+ to_string(static_cast<int>(dt * time)) + "_.dat";
                ofstream file1(filename1);
                if (!file1.is_open()) {
                    cerr << "Failed to open file for writing: " << filename1 << endl;
                    return 1;
                }
                for (int i = 0; i < gridsize; i++) {
                    for (int j = 0; j < gridsize; j++) {
                        file1 << i + 1 << " " << j + 1 << " " << localorder[i][j] << "\n";
                    }
                    file1 << "\n";
                }

                string filename2 = folder_path + "/density_heatmaps/localorder_" + to_string(trial) + "_" + to_string(static_cast<int>(dt * time)) + "_.dat";
                ofstream file2(filename2);
                if (!file2.is_open()) {
                    cerr << "Failed to open file for writing: " << filename2 << endl;
                    return 1;
                }

                for (int i = 0; i < gridsize; i++) {
                    for (int j = 0; j < gridsize; j++) {
                        file2 << localorder[i][j] << "\n";
                    }
                }
                string filename3 = folder_path + "/density_heatmaps/localdensityimage_" +to_string(trial)+"_"+ to_string(static_cast<int>(dt * time)) + "_.dat";
                ofstream file3(filename3);
                if (!file3.is_open()) {
                    cerr << "Failed to open file for writing: " << filename3 << endl;
                    return 1;
                }

                for (int i = 0; i < gridsize; i++) {
                    for (int j = 0; j < gridsize; j++) {
                        file3 << i + 1 << " " << j + 1 << " " << localdensity[i][j] << "\n";
                    }
                    file3 << "\n";
                }

                string filename4= folder_path + "/density_heatmaps/localdensity_" + to_string(trial) + "_" + to_string(static_cast<int>(dt * time)) + "_.dat";
                ofstream file4(filename4);
                if (!file4.is_open()) {
                    cerr << "Failed to open file for writing: " << filename4 << endl;
                    return 1;
                }

                for (int i = 0; i < gridsize; i++) {
                    for (int j = 0; j < gridsize; j++) {
                        file4 << localdensity[i][j] << "\n";
                    }
                }           
        
            }
        }
    }



    finish_time=time(NULL);
    cout<<"\n"<<"Total time elapsed : "<< finish_time - start_time <<" seconds ";
}
