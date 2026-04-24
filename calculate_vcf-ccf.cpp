#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <map>
#include <algorithm>

using namespace std;
namespace fs = std::filesystem;

struct ParticleSnapshot {
    double x, y;
    double vx, vy;
};

struct Parameters {
    int N;
    double Lx, Ly;
    double half_angle;
    double noise;
    double v0;
    double dt;
    int tmax;
    int trial;
    int seed;
};

bool read_parameters(const fs::path &param_file, Parameters &p) {
    ifstream f(param_file);
    if (!f.is_open()) {
        cerr << "Could not open parameters file: " << param_file << endl;
        return false;
    }

    string line;
    // read header
    if (!getline(f, line)) {
        cerr << "Empty parameters file: " << param_file << endl;
        return false;
    }

    // read data line (only first line is used)
    if (!getline(f, line)) {
        cerr << "No data line in parameters file: " << param_file << endl;
        return false;
    }

    stringstream ss(line);
    char comma;
    // N,Lx,Ly,half_angle,noise,v0,dt,maxiter,trial,seed
    if (!(ss >> p.N >> comma >> p.Lx >> comma >> p.Ly >> comma
          >> p.half_angle >> comma >> p.noise >> comma >> p.v0
          >> comma >> p.dt >> comma >> p.tmax >> comma >> p.trial
          >> comma >> p.seed)) {
        cerr << "Error parsing parameters in file: " << param_file << endl;
        return false;
    }

    return true;
}

bool read_snapshot(const fs::path &file_path, vector<ParticleSnapshot> &snap) {
    ifstream f(file_path);
    if (!f.is_open()) {
        cerr << "Could not open snapshot file: " << file_path << endl;
        return false;
    }

    string line;
    // Skip header line: x,y,vx,vy
    if (!getline(f, line)) {
        cerr << "Empty snapshot file: " << file_path << endl;
        return false;
    }

    snap.clear();
    while (getline(f, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        ParticleSnapshot p;
        char comma;
        if (!(ss >> p.x >> comma >> p.y >> comma >> p.vx >> comma >> p.vy)) {
            cerr << "Error parsing snapshot line in file: " << file_path << endl;
            return false;
        }
        snap.push_back(p);
    }

    return true;
}

int extract_step_from_filename(const string &filename) {
    // expects: config_<step>.csv
    // find between "config_" and ".csv"
    const string prefix = "config_";
    const string suffix = ".csv";

    size_t pos1 = filename.find(prefix);
    size_t pos2 = filename.rfind(suffix);
    if (pos1 == string::npos || pos2 == string::npos || pos2 <= pos1 + prefix.size()) {
        return -1;
    }
    string num_str = filename.substr(pos1 + prefix.size(), pos2 - (pos1 + prefix.size()));
    try {
        return stoi(num_str);
    } catch (...) {
        return -1;
    }
}

int main() {
    fs::path base_dir = "data_";

    if (!fs::exists(base_dir) || !fs::is_directory(base_dir)) {
        cerr << "Base directory 'data_' does not exist." << endl;
        return 1;
    }

    // Iterate over Angle_* directories
    for (const auto &angle_entry : fs::directory_iterator(base_dir)) {
        if (!angle_entry.is_directory()) continue;
        string angle_dir_name = angle_entry.path().filename().string();
        if (angle_dir_name.rfind("Angle_", 0) != 0) continue;

        // Iterate over Noise_* directories inside each Angle_*
        for (const auto &noise_entry : fs::directory_iterator(angle_entry.path())) {
            if (!noise_entry.is_directory()) continue;
            string noise_dir_name = noise_entry.path().filename().string();
            if (noise_dir_name.rfind("Noise_", 0) != 0) continue;

            fs::path param_file = noise_entry.path() / "parameters.csv";
            Parameters params;
            if (!read_parameters(param_file, params)) {
                cerr << "Skipping directory (no valid parameters): " << noise_entry.path() << endl;
                continue;
            }

            // Prepare correlation_data directory
            fs::path corr_dir = noise_entry.path() / "correlation_data";
            if (!fs::exists(corr_dir)) {
                fs::create_directories(corr_dir);
            }

            cout << "Processing: " << angle_dir_name << " / " << noise_dir_name
                 << " (trial " << params.trial << ")" << endl;

            // config_data/trial_<trial>/
            fs::path config_root = noise_entry.path() / "config_data";
            if (!fs::exists(config_root) || !fs::is_directory(config_root)) {
                cerr << "config_data directory missing at: " << config_root << endl;
                continue;
            }

            // There may be multiple trial_* directories if you later run more trials.
            for (const auto &trial_entry : fs::directory_iterator(config_root)) {
                if (!trial_entry.is_directory()) continue;
                string trial_dir_name = trial_entry.path().filename().string();
                if (trial_dir_name.rfind("trial_", 0) != 0) continue;

                // Collect all snapshot files config_*.csv, and their time index
                vector<pair<int, fs::path>> snapshots;  // (step, path)

                for (const auto &cfg_file : fs::directory_iterator(trial_entry.path())) {
                    if (!cfg_file.is_regular_file()) continue;
                    string fname = cfg_file.path().filename().string();
                    if (fname.rfind("config_", 0) != 0) continue;
                    if (fname.size() < 5 || fname.substr(fname.size() - 4) != ".csv") continue;

                    int step = extract_step_from_filename(fname);
                    if (step >= 0) {
                        snapshots.emplace_back(step, cfg_file.path());
                    }
                }

                if (snapshots.empty()) {
                    cerr << "No snapshot files in " << trial_entry.path() << endl;
                    continue;
                }

                // Sort by time step
                sort(snapshots.begin(), snapshots.end(),
                     [](const auto &a, const auto &b){ return a.first < b.first; });

                // Take the earliest snapshot as reference t0
                int t0 = snapshots.front().first;
                fs::path t0_file = snapshots.front().second;

                vector<ParticleSnapshot> ref_snap;
                if (!read_snapshot(t0_file, ref_snap)) {
                    cerr << "Failed to read reference snapshot: " << t0_file << endl;
                    continue;
                }

                if ((int)ref_snap.size() != params.N) {
                    cerr << "Warning: N from parameters (" << params.N
                         << ") != number of particles in snapshot (" << ref_snap.size()
                         << ") for " << t0_file << endl;
                }

                double v0 = params.v0;
                double v0_sq = v0 * v0;
                if (v0_sq == 0.0) {
                    cerr << "v0 = 0 in parameters, cannot normalise VCF by v0^2." << endl;
                    continue;
                }

                // Compute VCF for each snapshot relative to t0
                vector<int> times;
                vector<double> vcf_values;

                for (const auto &entry : snapshots) {
                    int t = entry.first;
                    fs::path file_path = entry.second;

                    vector<ParticleSnapshot> snap;
                    if (!read_snapshot(file_path, snap)) {
                        cerr << "Skipping snapshot due to read error: " << file_path << endl;
                        continue;
                    }
                    if (snap.size() != ref_snap.size()) {
                        cerr << "Skipping snapshot " << file_path
                             << " because size mismatch with t0." << endl;
                        continue;
                    }

                    double dot_sum = 0.0;
                    int N_here = static_cast<int>(snap.size());
                    for (int i = 0; i < N_here; ++i) {
                        double dot = ref_snap[i].vx * snap[i].vx + ref_snap[i].vy * snap[i].vy;
                        dot_sum += dot;
                    }

                    double C_t = (dot_sum / (static_cast<double>(N_here))) / v0_sq;
                    times.push_back(t - t0);  // time difference
                    vcf_values.push_back(C_t);
                }

                // Save VCF to file
                // vcf file name: vcf_<trial>.csv (extract trial number from directory name)
                string trial_id = trial_dir_name.substr(string("trial_").size());
                fs::path out_file = corr_dir / ("vcf_trial_" + trial_id + ".csv");

                ofstream out(out_file);
                if (!out.is_open()) {
                    cerr << "Could not open output file for writing VCF: " << out_file << endl;
                    continue;
                }

                out << "t,vcf\n";
                for (size_t i = 0; i < times.size(); ++i) {
                    out << times[i] << "," << vcf_values[i] << "\n";
                }
                out.close();

                cout << "   Saved VCF for " << trial_dir_name
                     << " to " << out_file << endl;
            } // end loop over trial_* directories
        }     // end loop over Noise_*
    }         // end loop over Angle_*

    cout << "Finished computing VCF for all available parameter sets." << endl;
    return 0;
}