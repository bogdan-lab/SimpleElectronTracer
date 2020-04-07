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



PtStatistics::PtStatistics(const Point& pp, const vector<double>& pV, const int p_vol_count, const int p_surf_count){
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


void Surface::WriteInfo(){
    string name = label + "_stat.txt";
    const char* st_file_name = name.c_str();
    FILE* st_f;
    st_f = fopen(st_file_name, "a");
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


Point::Point(double xx, double yy, double zz){
    x = xx;
    y = yy;
    z = zz;
}

Point::Point(){
    x = 0.0;
    y = 0.0;
    z = 0.0;
}


int Surface::GetCoorFlag() const {return coor_flag;}
double Surface::GetCoorVal() const {return coor_val;}
vector<double> Surface::GetXbnd() const {return x_bnd;}
vector<double> Surface::GetYbnd() const {return y_bnd;}
vector<double> Surface::GetZbnd() const {return z_bnd;}
double Surface::GetRefl() const {return R;}
bool Surface::GetSaveStatFlag(){return save_stat;}


vector<double> Particle::GetRandVel(int direction, default_random_engine& rnd_gen) const {
    /*0 - for x, 1 - for y, 2 - for z, 3 - for 3D*/
    vector<double> v(3, 0.0);
    uniform_real_distribution<double> rnd(0.0, 1.0);
    double phi = rnd(rnd_gen)*2*M_PI;
    double costheta = rnd(rnd_gen);
    double theta = acos(costheta);
    if (direction==0){
        v[0] = costheta;
        v[1] = sin(theta)*cos(phi);
        v[2] = sin(theta)*sin(phi);
    }
    else if (direction==1){
        v[1] = costheta;
        v[0] = sin(theta)*cos(phi);
        v[2] = sin(theta)*sin(phi);
    }
    else if (direction==2){
        v[2] = costheta;
        v[1] = sin(theta)*cos(phi);
        v[0] = sin(theta)*sin(phi);
    }
    else {
        costheta = rnd(rnd_gen)*2 - 1;    //variate from -1 to 1
        theta = acos(costheta);
        v[2] = costheta;
        v[1] = sin(theta)*cos(phi);
        v[0] = sin(theta)*sin(phi);
    }
    return v;
}



double Particle::GetDistanceInGas(const double pressure, default_random_engine& rnd_gen) const{
    if (pressure==0.0){
        return 1e20;        //never collide with atom
    }
    uniform_real_distribution<double> rnd(0.0, 1.0);
    double sigma = 2e-16;       //cm^2
    double T = 300;         //K
    double k = 1.38e-23;     //J/K
    double N = pressure/(k*T)*1e-6;    //cm^-3
    double mfp = 1.0/(N*sigma);         //cm
    return mfp*log(1.0/(1.0 - rnd(rnd_gen)));
}



Particle::Particle(const Surface &s, default_random_engine& rnd_gen){
    p = s.GetRandomPoint(rnd_gen);
    //p = Point(50.0, 0.0, 50.0);
    vol_count = 0;
    surf_count = 0;
    vector<Point> tmp;
    tmp.push_back(p);
    tragectory = tmp;
    //rnd =  uniform_real_distribution<double>(0.0, 1.0);
    V = this->GetRandVel(1, rnd_gen);     //directed along positive y
    //V = {0.0, 1.0, 0.0};
}

Point Particle::GetCrossPoint(const Surface &s, bool& cross_flag){
    Point p_new(0.0, 0.0, 0.0);
    cross_flag = false;
    int coor_flag = s.GetCoorFlag();
    double coor_val = s.GetCoorVal();
    vector<double> x_bnd = s.GetXbnd();
    vector<double> y_bnd = s.GetYbnd();
    vector<double> z_bnd = s.GetZbnd();
    double t = 0.0;
    if (coor_flag==0){
        if((p.x<coor_val && V[0]>0.0) || (p.x>coor_val && V[0]<0.0)){
            //direction is correct
            p_new.x = coor_val;
            t = (coor_val - p.x)/V[0];
            p_new.y = p.y + V[1]*t;
            p_new.z = p.z + V[2]*t;
            if ((p_new.z>z_bnd[0] && p_new.z<z_bnd[1]) && (p_new.y>y_bnd[0] && p_new.y<y_bnd[1])){
                cross_flag = true;
            }
            else {
                cross_flag = false;
            }
        }
    }
    else if (coor_flag==1){
        if((p.y<coor_val && V[1]>0.0) || (p.y>coor_val && V[1]<0.0)){
            //direction is correct
            p_new.y = coor_val;
            t = (coor_val - p.y)/V[1];
            p_new.x = p.x + V[0]*t;
            p_new.z = p.z + V[2]*t;
            if ((p_new.z>z_bnd[0] && p_new.z<z_bnd[1]) && (p_new.x>x_bnd[0] && p_new.x<x_bnd[1])){
                cross_flag = true;
            }
            else {
                cross_flag = false;
            }
        }
    }
    else{
        if((p.z<coor_val && V[2]>0.0) || (p.z>coor_val && V[2]<0.0)){
            //direction is correct
            p_new.z = coor_val;
            t = (coor_val - p.z)/V[2];
            p_new.y = p.y + V[1]*t;
            p_new.x = p.x + V[0]*t;
            if ((p_new.x>x_bnd[0] && p_new.x<x_bnd[1]) && (p_new.y>y_bnd[0] && p_new.y<y_bnd[1])){
                cross_flag = true;
            }
            else {
                cross_flag = false;
            }
        }
    }
    return p_new;
}


double Particle::GetDistanceToSurface(const Surface &s){
    bool cross_flag;
    Point cross_point;
    cross_point = this->GetCrossPoint(s, cross_flag);
    double dist;
    if (cross_flag) {
        dist = sqrt((cross_point.x - this->p.x)*(cross_point.x - this->p.x) +
                    (cross_point.y - this->p.y)*(cross_point.y - this->p.y) +
                    (cross_point.z - this->p.z)*(cross_point.z - this->p.z));
    } else {
        dist = 1e20;
    }
    return dist;
}


int Particle::GetReflectionSurfaceID(const vector<Surface>& walls){
    double min_dist = 1e20;
    int idx = -1;
    for(int i=0; i<walls.size(); i++){
        double dist = this->GetDistanceToSurface(walls[i]);
        if (dist<min_dist){
            min_dist = dist;
            idx = i;
        }
    }
    return idx;
}


bool Particle::ReflectSurface(Surface &s){
    /*returns true/false for alive/dead particle*/
    bool cross_flag;
    Point p_new = this->GetCrossPoint(s, cross_flag);
    if (!cross_flag){
        printf("\nReflecting from surface I did not cross\n");
        printf("COOR_FLAG = %i\n", s.GetCoorFlag());
        printf("COOR_VAL = %lf\n", s.GetCoorVal());
        exit(1);
    }
    tragectory.push_back(p_new);
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    uniform_real_distribution<double> rnd(0.0,1.0);
    double roll = rnd(generator);
    if (roll>s.GetRefl()){
        //kill particle
        PtStatistics pt_stat(p_new, V, vol_count, surf_count);
        if (s.GetSaveStatFlag()){
            s.SaveInfo(pt_stat, tragectory);
        }
        return false;
    }
    else{
        //Do reflection
        int coor_flag = s.GetCoorFlag();
        V[coor_flag]*=-1;
        surf_count++;
        p = p_new;
        return true;
    }
}


void Particle::MakeGasCollision(double pt_dist, default_random_engine& rnd_gen){
    vol_count++;
    double t = pt_dist;     //since velocity is 1
    p.x+=V[0]*t;
    p.y+=V[1]*t;
    p.z+=V[2]*t;
    tragectory.push_back(p);
    V = this->GetRandVel(3, rnd_gen);
}


Point Particle::GetPosition() const{return p;}


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
                alive_flag = pt.ReflectSurface(walls[ref_wall_idx]);
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
    const double pressure = 0.0;   //Pa
    const size_t group_len = 70000000;
    printf("Pressure = %.3lf Pa\n" , pressure);

    #pragma omp parallel shared(walls) num_threads(7)
    {
        #pragma omp for
        for(size_t i=1; i<8; i++){
            unsigned seed = chrono::system_clock::now().time_since_epoch().count();
            seed*=i;
            default_random_engine rnd_gen(seed);
            RunParticleGroup(walls, pressure, rnd_gen, group_len);
        }
    }


    for(size_t i=0; i<walls.size(); i++){
        if (walls[i].GetSaveStatFlag()){
            walls[i].WriteInfo();
        }
    }
    return 0;
}







