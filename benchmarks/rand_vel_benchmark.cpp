#include <benchmark/benchmark.h>
#include <random>
#include <cmath>

#include "utils.hpp"



Vector ACOSRandVel(std::mt19937& rnd_gen){
    std::uniform_real_distribution<double> rnd(0.0, 1.0);
    double cos_theta = rnd(rnd_gen)*2 - 1;
    double theta = std::acos(cos_theta);
    double phi = rnd(rnd_gen)*2*M_PI;
    return Vector(sin(theta)*sin(phi), sin(theta)*cos(phi), cos(theta));
}

Vector SQRTRandVel(std::mt19937& rnd_gen){
    std::uniform_real_distribution<double> rnd(0.0, 1.0);
    double cos_theta = rnd(rnd_gen)*2 - 1;
    double sin_theta = sqrt(1-cos_theta*cos_theta);
    double phi = rnd(rnd_gen)*2*M_PI;
    return Vector(sin_theta*sin(phi), sin_theta*cos(phi), cos_theta);
}



static void ACOSgenerator(benchmark::State& state){
    std::mt19937 rnd_gen;
    rnd_gen.seed(static_cast<uint>(time(0)));
    while(state.KeepRunning()){
        Vector x = ACOSRandVel(rnd_gen);
        benchmark::DoNotOptimize(x);
    }
}

static void SQRTgenerator(benchmark::State& state){
    std::mt19937 rnd_gen;
    rnd_gen.seed(static_cast<uint>(time(0)));
    while(state.KeepRunning()){
        Vector x = SQRTRandVel(rnd_gen);
        benchmark::DoNotOptimize(x);
    }
}


BENCHMARK(ACOSgenerator);
BENCHMARK(SQRTgenerator);

BENCHMARK_MAIN();
