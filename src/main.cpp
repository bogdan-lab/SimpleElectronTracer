#include <time.h>

#include "particle.hpp"
#include "surface.hpp"
#include "utils.hpp"
#include "reflector.hpp"



int main(){

    //plot geometry
    std::vector<Surface> walls;
    walls.reserve(6);
    //SURFACE 0 - YZ
    {
        std::string name = "surface_0_YZ";
        std::vector<Vector> contour {Vector(0.0, 0.0, 0.0),
                    Vector(0.0, 1.0, 0.0),
                    Vector(0.0, 1.0, 1.0),
                    Vector(0.0, 0.0, 1.0)};
        walls.push_back(Surface(name, contour,
                                std::make_unique<MirrorReflector>(0.9), true));
    }
    //SURFACE 1 - XZ
    {
        std::string name ="surface_1_XZ";
        std::vector<Vector> contour {Vector(0.0, 0.0, 0.0),
                                     Vector(0.0, 0.0, 1.0),
                                     Vector(1.0, 0.0, 1.0),
                                     Vector(1.0, 0.0, 0.0)};
        walls.push_back(Surface(name, contour,
                                std::make_unique<LambertianReflector>(0.9), true));

    }
    //SURFACE 2 - XY
    {
        std::string name = "surface_2_XY";
        std::vector<Vector> contour {Vector(0.0, 0.0, 0.0),
                                     Vector(1.0, 0.0, 0.0),
                                     Vector(1.0, 1.0, 0.0),
                                     Vector(0.0, 1.0, 0.0)};
        walls.push_back(Surface(name, contour,
                                std::make_unique<LambertianReflector>(0.9), true));

    }
    //SURFACE 3 - YZ
    {
        std::string name = "surface_3_YZ";
        std::vector<Vector> contour {Vector(1.0, 0.0, 0.0),
                                     Vector(1.0, 0.0, 1.0),
                                     Vector(1.0, 1.0, 1.0),
                                     Vector(1.0, 1.0, 0.0)};
        walls.push_back(Surface(name, contour,
                                std::make_unique<MirrorReflector>(0.9), true));

    }
    //SURFACE 4 - XZ
    {
        std::string name ="surface_4_XZ";
        std::vector<Vector> contour {Vector(0.0, 1.0, 0.0),
                                     Vector(1.0, 1.0, 0.0),
                                     Vector(1.0, 1.0, 1.0),
                                     Vector(0.0, 1.0, 1.0)};
        walls.push_back(Surface(name, contour,
                                std::make_unique<LambertianReflector>(0.9), true));

    }
    //SURFACE 5 - XY
    {
        std::string name = "surface_5_XY";
        std::vector<Vector> contour {Vector(0.0, 0.0, 1.0),
                                     Vector(0.0, 1.0, 1.0),
                                     Vector(1.0, 1.0, 1.0),
                                     Vector(1.0, 0.0, 1.0)};
        walls.push_back(Surface(name, contour,
                                std::make_unique<LambertianReflector>(0.9), true));

    }

    Background gas = {2e-16, 300.0, 5.0};

    //RUNNING PARTICLES
    std::mt19937 rnd_gen;
    rnd_gen.seed(static_cast<uint>(time(0)));
    for(size_t i=0; i<1000; i++){
        Particle pt(Vector(0.5, 0.5, 0.5), Vector(0.5, 0.5, 0.5), rnd_gen);
        while(pt.MakeStep(walls, gas, rnd_gen)){}
    }

    for(size_t i=0; i<walls.size(); i++){
        if(walls[i].IsSaveStat()){
            std::ofstream file;
            file.open(walls[i].GetName());
            if(!file){
                std::cerr << "Unable to open file\n";
                exit(1);
            }
            walls[i].SaveSurfaceParticles(file);
            file.close();
        }
    }
    return 0;
}
