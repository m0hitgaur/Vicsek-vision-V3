#pragma once
#include <iostream>
#include <time.h>
#include <mutex>
#include <filesystem>
using namespace std;
std::mutex cout_mutex;
namespace fs = filesystem;
void print_progress(double progress,double timestarted) 

{   double current_time=static_cast <double>(time(NULL));
    double elapsed_time=current_time - timestarted;
    if(progress>0.0)
        {double avg_time_remaining= elapsed_time*(1.0/progress -1.0);
        const int bar_width = 50;
        cout << "\rProgress: [";
        int pos = bar_width * progress;
        for (int i = 0; i < bar_width; i++) {
            if (i < pos) cout << "=";
            else if (i == pos) cout << ">";
            else cout << " ";
        }
        cout << "] " << int(progress * 100.0) << " %"<<"  Estimated time remaining: "<< avg_time_remaining<<" secs ("<<avg_time_remaining/60<<" mins) ("<<avg_time_remaining/3600<<" hrs)"<<flush;
         }

    else  cout << "] " << int(progress * 100.0) << " %"<<"  Estimated time remaining: N/A secs"<<flush;
     
}


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

void print_progress_threadsafe(const std::string& label,
                               double fraction,
                               double elapsed_seconds)
{
    std::lock_guard<std::mutex> lock(cout_mutex);

    std::cout << label << " " << std::fixed << std::setprecision(1)<< (fraction * 100.0) << "% "<< "(elapsed: " << elapsed_seconds << " s)"<< std::endl;
}