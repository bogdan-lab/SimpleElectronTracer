#include <benchmark/benchmark.h>
#include <random>
#include <vector>
#include <string>

#include "particle.hpp"
#include "surface.hpp"
#include "utils.hpp"

std::vector<Surface> PrepareGeometry(){
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
}

static void MakeStepBenchmark(benchmark::State& state){
    //Preparing envirounment for running
    std::vector<Surface> walls = PrepareGeometry();
    Background gas = {2e-16, 300.0, 5.0};
    std::mt19937 rnd_gen;
    rnd_gen.seed(static_cast<uint>(time(0)));
    for(auto _ : state){
        Particle pt(Vector(0.5, 0.5, 0.0), Vector(0.0, 0.0, 1.0), rnd_gen);
        benchmark::DoNotOptimize(pt);
        pt.MakeStep(walls, gas, rnd_gen);
    }
}




BENCHMARK(MakeStepBenchmark);

BENCHMARK_MAIN();
