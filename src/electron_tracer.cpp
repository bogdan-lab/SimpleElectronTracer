#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
#include "electron_tracer.h"
#include <sstream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <random>
#include <cmath>
#include <future>
#include <omp.h>

//#include <thread>


using namespace std;



PtStatistics::PtStatistics(const Point& pp, const vector<double>& pV,
                           const int p_vol_count, const int p_surf_count){
    p = pp;
    V = pV;
    vol_count = p_vol_count;
    surf_count = p_surf_count;
}

void Surface::SaveInfo(const PtStatistics& pt_stat, const vector<Point>& pt_tragectory){
     #pragma omp critical (saving)
    {
        stat.push_back(pt_stat);
        //tragectories.push_back(pt_tragectory);
    }
}


void Surface::WriteInfo(const size_t number_of_simulated_particles, const double pressure){
    string name = label + "_stat.txt";
    const char* st_file_name = name.c_str();
    FILE* st_f;
    st_f = fopen(st_file_name, "a");
    fprintf(st_f, "#Number of particles simulated in this run %lu\n", number_of_simulated_particles);
    fprintf(st_f, "#Surface reflection coefficient %lf\n", this->R);
    fprintf(st_f, "#Pressure in simulateion %lf Pa\n", pressure);
    fprintf(st_f, "#x\ty\tz\tVx\tVy\tVz\tvolume count\tsurface_count\n");
    for(size_t i=0; i<stat.size(); i++){
        fprintf(st_f, "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%i\t%i\n",
                stat[i].p.x, stat[i].p.y, stat[i].p.z, stat[i].V[0], stat[i].V[1], stat[i].V[2],
                stat[i].vol_count, stat[i].surf_count);
    }
    fclose(st_f);
    name = label + "_tragectories.txt";
    const char* tg_file_name = name.c_str();
    FILE* tg_f;
    tg_f = fopen(tg_file_name, "a");
    fprintf(tg_f, "#x\ty\tz\n");
    for(size_t i=0; i<tragectories.size(); i++){
        for(size_t j=0; j<tragectories[i].size(); j++){
              fprintf(tg_f, "%.8e\t%.8e\t%.8e\n", tragectories[i][j].x, tragectories[i][j].y, tragectories[i][j].z);
        }
        fprintf(tg_f, "NEXT\n");
    }
    fclose(tg_f);
}


Point Surface::GetRandomPoint(default_random_engine& rnd_gen) const {
    Point p;
    uniform_real_distribution<double> rnd(0.0,1.0);
    if (coor_flag==0){
        p.x = coor_val;
        p.y = rnd(rnd_gen)*(y_bnd[1]-y_bnd[0]) + y_bnd[0];
        p.z = rnd(rnd_gen)*(z_bnd[1]-z_bnd[0]) + z_bnd[0];
    } else if (coor_flag==1){
        p.y = coor_val;
        p.x = rnd(rnd_gen)*(x_bnd[1]-x_bnd[0]) + x_bnd[0];
        p.z = rnd(rnd_gen)*(z_bnd[1]-z_bnd[0]) + z_bnd[0];
    } else{
        p.z = coor_val;
        p.x = rnd(rnd_gen)*(x_bnd[1]-x_bnd[0]) + x_bnd[0];
        p.y = rnd(rnd_gen)*(y_bnd[1]-y_bnd[0]) + y_bnd[0];
    }
    return p;
}

Surface::Surface(string &line){
    stringstream ss;
    ss << line;
    ss >> coor_flag;
    ss >> coor_val;
    ss >> R;
    vector<double> tmp_x(2, 0.0);
    vector<double> tmp_y(2, 0.0);
    vector<double> tmp_z(2, 0.0);
    ss >> tmp_x[0];
    ss >> tmp_x[1];
    x_bnd = tmp_x;
    ss >> tmp_y[0];
    ss >> tmp_y[1];
    y_bnd = tmp_y;
    ss >> tmp_z[0];
    ss >> tmp_z[1];
    z_bnd = tmp_z;
    ss >> label;
    int flag;
    ss >> flag;
    if (flag==1){
        save_stat = true;
    } else{
        save_stat = false;
    }
}

void Surface::PrintSurface(){
    cout << coor_flag << " " << coor_val << " " << R << " " << x_bnd[0] << " " << x_bnd[1] << " " <<
            y_bnd[0] << " " << y_bnd[1] << " " << z_bnd[0] << " " << z_bnd[1] << " " << label << " " <<
            save_stat << "\n";
}



void RunParticleGroup(vector<Surface>& walls, const double pressure, default_random_engine& rnd_gen, const size_t group_len){
    for(size_t i=0; i<group_len; i++){
        Particle pt(walls[0], rnd_gen);
        bool alive_flag = true;
        while (alive_flag){
            //get_collision distance
            double gas_dist = pt.GetDistanceInGas(pressure, rnd_gen);
            //compare with wall distance
            int ref_wall_idx = pt.GetReflectionSurfaceID(walls);
            double wall_dist = pt.GetDistanceToSurface(walls[ref_wall_idx]);
            if (wall_dist > gas_dist){
                pt.MakeGasCollision(gas_dist, rnd_gen);
            }
            else{
                alive_flag = pt.ReflectSurface(walls[ref_wall_idx], rnd_gen);
            }
        }
    }
}



int main(){
    //READ GEOMETRY
    string line;
    ifstream f("geometry.dat");
    vector<Surface> walls;
    while(getline(f, line)){
        if (line[0]=='#'){
            continue;
        }
        Surface s(line);
        walls.push_back(s);
    }
    f.close();
    const double pressure = 5.0;   //Pa
    const size_t group_len = 70000000;
    const size_t simulated_group_number = 7;
    printf("Pressure = %.3lf Pa\n" , pressure);

    #pragma omp parallel shared(walls) num_threads(7)
    {
        #pragma omp for
        for(size_t i=1; i<simulated_group_number+1; i++){
            unsigned seed = chrono::system_clock::now().time_since_epoch().count();
            seed*=i;
            default_random_engine rnd_gen(seed);
            RunParticleGroup(walls, pressure, rnd_gen, group_len);
        }
    }

    size_t simulated_particle_number = group_len*simulated_group_number;
    for(size_t i=0; i<walls.size(); i++){
        if (walls[i].GetSaveStatFlag()){
            walls[i].WriteInfo(simulated_particle_number, pressure);
        }
    }
    return 0;
}







