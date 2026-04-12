#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <filesystem>
#include <vector>
using namespace std;
namespace fs = std::filesystem;
string folder_path = fs::current_path().string();







int main()
{
int    Number_of_active_agents =2560  ; // number of particles
double Lx = 32;                        // length of the system
double Ly = 32 ;                        // length of the system
int    maxiter = 10000;                  // total number of time steps 
double velocity = 0.01;                // magnitude of velocity 
double dt = 1;         
int numberoftrials=100;              // time step
int trialstart=0  ;          //starting value of trial
int timeinterval_flocking=10;

// for connected correlation function /velocity correlation function
double binsize_ccf=0.499;
int timeinterval_ccf=100; // interval after which ccf/vcf is calculated

int gridsize=10;
int timeinterval_density=100; // interval after which density and order is calculated

// for fluctuation
int numberofgrids=10; //number of grid points


vector <double> half_angle={M_PI_4,M_PI};
vector <string> noise={"0.05"};
double density = Number_of_active_agents/(Lx*Ly);  // density
double radius = 1;  // radius of interaction 
int seedvalue=19809876;
for(int i=0;i<trialstart;i++)
    seedvalue+=10;

int numberofangles=half_angle.size();
int numberofnoises=noise.size();
for(int i=0;i<numberofangles;i++)
{
    for (int j=0;j<numberofnoises;j++)
    {
        stringstream s1;
        s1<< fixed <<setprecision(0)<< half_angle[i]*(180/3.14159);
        string anglefilename=s1.str();

        ofstream parameter(folder_path+"/Angle_"+anglefilename+"/Noise_"+ noise[j]+"/parameters.txt");
        if(!parameter.is_open())
            cout<<"/Angle_"+anglefilename+"/Noise_"+ noise[j]+"parameters.txt not open"<<"\n";

        parameter<<Number_of_active_agents<<"\n";
        parameter<<Lx<<"\n";
        parameter<<Ly<<"\n";
        parameter<<maxiter<<"\n";
        parameter<<velocity<<"\n";
        parameter<<dt<<"\n";
        parameter<<half_angle[i]<<"\n";
        parameter<<noise[j]<<"\n";
        parameter<<density<<"\n";
        parameter<<radius<<"\n";
        parameter<<numberoftrials<<"\n";
        parameter<<trialstart<<"\n";
        parameter<<seedvalue<<"\n";
        parameter<<timeinterval_flocking<<"\n";
    }
}

for(int i=0;i<numberofangles;i++)
{
    for (int j=0;j<numberofnoises;j++)
    {
        stringstream s1;
        s1<< fixed <<setprecision(0)<< half_angle[i]*(180/3.14159);
        string anglefilename=s1.str();

        ofstream parameter(folder_path+"/Angle_"+anglefilename+"/Noise_"+ noise[j]+"/ccf-vcfparameters.txt");
        if(!parameter.is_open())
            cout<<"/Angle_"+anglefilename+"/Noise_"+ noise[j]+"ccf-vcfparameters.txt not open"<<"\n";

        parameter<<Number_of_active_agents<<"\n";
        parameter<<Lx<<"\n";
        parameter<<Ly<<"\n";
        parameter<<maxiter<<"\n";
        parameter<<velocity<<"\n";
        parameter<<dt<<"\n";
        parameter<<half_angle[i]<<"\n";
        parameter<<noise[j]<<"\n";
        parameter<<density<<"\n";
        parameter<<radius<<"\n";
        parameter<<numberoftrials<<"\n";
        parameter<<binsize_ccf<<"\n";
        parameter<<timeinterval_ccf<<"\n";
    }
}


for(int i=0;i<numberofangles;i++)
{
    for (int j=0;j<numberofnoises;j++)
    {
        stringstream s1;
        s1<< fixed <<setprecision(0)<< half_angle[i]*(180/3.14159);
        string anglefilename=s1.str();

        ofstream parameter(folder_path+"/Angle_"+anglefilename+"/Noise_"+ noise[j]+"/densityparameters.txt");
        if(!parameter.is_open())
            cout<<"/Angle_"+anglefilename+"/Noise_"+ noise[j]+"densityparameters.txt not open"<<"\n";

        parameter<<Number_of_active_agents<<"\n";
        parameter<<Lx<<"\n";
        parameter<<Ly<<"\n";
        parameter<<maxiter<<"\n";
        parameter<<velocity<<"\n";
        parameter<<dt<<"\n";
        parameter<<half_angle[i]<<"\n";
        parameter<<noise[j]<<"\n";
        parameter<<density<<"\n";
        parameter<<radius<<"\n";
        parameter<<numberoftrials<<"\n";
        parameter<<gridsize<<"\n";
        parameter<<timeinterval_density<<"\n";
    }
}


for(int i=0;i<numberofangles;i++)
{
    for (int j=0;j<numberofnoises;j++)
    {
        stringstream s1;
        s1<< fixed <<setprecision(0)<< half_angle[i]*(180/3.14159);
        string anglefilename=s1.str();

        ofstream parameter(folder_path+"/Angle_"+anglefilename+"/Noise_"+ noise[j]+"/fluctuationparameters.txt");
        if(!parameter.is_open())
            cout<<"/Angle_"+anglefilename+"/Noise_"+ noise[j]+"fluctuationparameters.txt not open"<<"\n";

        parameter<<Number_of_active_agents<<"\n";
        parameter<<Lx<<"\n";
        parameter<<Ly<<"\n";
        parameter<<maxiter<<"\n";
        parameter<<velocity<<"\n";
        parameter<<dt<<"\n";
        parameter<<half_angle[i]<<"\n";
        parameter<<noise[j]<<"\n";
        parameter<<density<<"\n";
        parameter<<radius<<"\n";
        parameter<<numberoftrials<<"\n";
        parameter<<numberofgrids<<"\n";

    }
}



}
