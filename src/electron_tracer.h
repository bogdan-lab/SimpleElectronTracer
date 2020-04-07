#pragma once
#include <vector>
#include <string>
#include <random>
using namespace std;

struct Point{
Point(double xx, double yy, double zz);
Point();
double x;
double y;
double z;
};


struct PtStatistics
{
PtStatistics(const Point& pp, const vector<double>& pV, const int p_vol_count, const int p_surf_count);
Point p;
vector<double> V;  //velocity, v[0] - vx, v[1] = vy, v[2] = vz
int vol_count;
int surf_count;
};

class Surface{
public:
    Surface(string& line);
    void SaveInfo(const PtStatistics& pt_stat, const vector<Point>& pt_tragectory);
    void WriteInfo();
    Point GetRandomPoint(default_random_engine& rnd_gen) const;
    void PrintSurface();

    int GetCoorFlag() const;
    double GetCoorVal() const;
    vector<double> GetXbnd() const ;
    vector<double> GetYbnd() const ;
    vector<double> GetZbnd() const ;
    double GetRefl() const ;
    bool GetSaveStatFlag();
private:
    vector<PtStatistics> stat;      //vector of particle statistics
    bool save_stat;             //flag for saving statistics
    vector<vector<Point>> tragectories;     //saved particles trajectories
    int coor_flag;      //0 - for x, 1 - for y, 2 for z
    double coor_val;            //coordinate of surface position
    vector<double> x_bnd;       //surface boundaries not all required
    vector<double> y_bnd;
    vector<double> z_bnd;
    string label;           //surface name
    double R;           //reflection
};


class Particle{
public:
    Particle(const Surface& s, default_random_engine& rnd_gen);
    Point GetCrossPoint(const Surface& s, bool& cross_flag);
    double GetDistanceToSurface(const Surface& s);
    int GetReflectionSurfaceID(const vector<Surface>& walls);
    bool ReflectSurface(Surface& s);
    void MakeGasCollision(double pt_dist, default_random_engine& rnd_gen);
    Point GetPosition() const;
    double GetDistanceInGas(const double pressure, default_random_engine& rnd_gen) const;
    vector<double>GetRandVel(int direction, default_random_engine& rnd_gen) const;
private:
    Point p;
    vector<double> V;
    int vol_count;
    int surf_count;
    vector<Point> tragectory;
};





void RunParticleGroup(vector<Surface>& walls, const double mfp, default_random_engine& rnd_gen, const size_t group_len);
