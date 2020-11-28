#include <time.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <lyra/lyra.hpp>
#include<fmt/core.h>
#include <omp.h>

#include "particle.hpp"
#include "surface.hpp"
#include "math.hpp"
#include "reflector.hpp"
#include "loader.hpp"

int main(int argc, const char ** argv){
    std::string config_file = "NO_FILE_WAS_GIVEN";
    bool show_help = false;
    auto cli = lyra::cli()
            | lyra::opt(config_file, "config")["-c"]["--config"]
                ("Path to the json config file [no default value!]")
            | lyra::help(show_help);
    auto cmd_parse = cli.parse({argc, argv});
    if(show_help){
        std::cout << cli;
        return 0;
    }
    if (!cmd_parse){
        std::cerr << cmd_parse.errorMessage() << "\n";
    }

    json json_data = load_json_config(config_file);
    Background gas = load_background(json_data);
    std::vector<std::unique_ptr<Surface>> walls = load_geometry(json_data);
    size_t pt_num = json_data["particles"]["number"].get<size_t>();
    Vec3 source_point(json_data["particles"]["source_point"].get<std::vector<double>>());
    Vec3 direction(json_data["particles"]["direction"].get<std::vector<double>>());
    bool is_dir_random = json_data["particles"]["is_dir_random"].get<bool>();


    auto pt_generator = Particle::GetGenerator(is_dir_random);
    //************MAIN CYLE******************
    size_t thread_num = json_data["general"]["number_of_threads"].get<size_t>();
    if(thread_num<1) {std::cerr << "Wrong thread number\n"; exit(1);}
    omp_set_dynamic(0);
    omp_set_num_threads(static_cast<int>(thread_num));
    std::vector<size_t> thread_load(thread_num, pt_num/thread_num);
    thread_load.back() += pt_num % thread_num;
    #pragma omp parallel
    {
        size_t traced_pt_num = 0;
        size_t tid = static_cast<size_t>(omp_get_thread_num());

        std::ofstream output_file(fmt::format("thread_{:d}", tid),
                                  std::ios_base::app);
        if(!output_file.is_open()){
            std::cerr << fmt::format("Cannot open file for thread {:d}\n", tid);
            exit(1);
        }
        output_file << "#POS_X\tPOS_Y\tPOS_Z"
                     << "\tVX\tVY\tVZ\tVolumeCount\tSurfaceCount\n";
        std::mt19937 rnd_gen;
        rnd_gen.seed(static_cast<uint>(time(0))+tid);
        while(traced_pt_num<thread_load[tid]){
            traced_pt_num += pt_generator(source_point, direction, rnd_gen)
                    .Trace(walls, gas, rnd_gen, output_file);
            #pragma omp master
            {
                if((traced_pt_num+1)%(thread_load[tid]/10)==0){
                    std::cout << fmt::format("{:d} %\n",
                    (100*(traced_pt_num+1))/thread_load[tid]);
                }
            }
        }
    }
    //***********CYCLE END*******************
    return 0;
}

