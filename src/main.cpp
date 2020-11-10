#include <time.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <lyra/lyra.hpp>

#include "particle.hpp"
#include "surface.hpp"
#include "utils.hpp"
#include "reflector.hpp"

using json = nlohmann::json;

json load_json_config(const std::string& file_name){
    std::ifstream file(file_name);
    if(!file.is_open()){
        fprintf(stderr, "File %s cannot be open", file_name.c_str());
        exit(1);
    }
    json json_data;
    file >> json_data;
    return json_data;
}

Background load_background(const json& json_data){
    return {json_data["gas"]["sigma"].get<double>(),
            json_data["gas"]["temperature"].get<double>(),
            json_data["gas"]["pressure"].get<double>()};
}

Surface read_surface_parameters(const json& this_surf_data){
    std::string name = this_surf_data["name"].get<std::string>();
    std::vector<Vec3> contour;
    for(const auto& el : this_surf_data["contour"]){
        contour.push_back(Vec3(el.get<std::vector<double>>()));
    }
    std::string ref_type = this_surf_data["reflector_type"].get<std::string>();
    double R = this_surf_data["reflection_coefficient"].get<double>();
    bool stat_flag = this_surf_data["collect_statistics"].get<bool>();
    if(ref_type == "mirror"){
        return {name, contour, std::make_unique<MirrorReflector>(R), stat_flag};
    }
    else if (ref_type == "cosine"){
        return {name, contour, std::make_unique<LambertianReflector>(R), stat_flag};
    }
    else {
        fprintf(stderr, "unknown reflector type %s", ref_type.c_str());
        exit(1);
    }
}

std::vector<Surface> load_geometry(const json& json_data){
    std::vector<Surface> walls;
    for(const auto& el : json_data["geometry"]){
        walls.push_back(read_surface_parameters(el));
    }
    return walls;
}

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
    std::vector<Surface> walls = load_geometry(json_data);
    size_t pt_num = json_data["particles"]["number"].get<size_t>();
    Vec3 source_point(json_data["particles"]["source_point"].get<std::vector<double>>());
    Vec3 direction(json_data["particles"]["direction"].get<std::vector<double>>());
    bool is_dir_random = json_data["particles"]["is_dir_random"].get<bool>();

    std::mt19937 rnd_gen;
    rnd_gen.seed(static_cast<uint>(time(0)));

    auto pt_generator = [&](std::mt19937& rnd_gen){
        if (is_dir_random){
            return Particle(source_point, direction, rnd_gen);
        } else {
            return Particle(source_point, direction);
        }
    };
    //************MAIN CYLE******************
    for(size_t i=0; i<pt_num; i++){
        pt_generator(rnd_gen)
                .Trace(walls, gas, rnd_gen);
        if((i+1)%(pt_num/10)==0){
            printf("%.2lf %%\n" , static_cast<double>(100.0*i/pt_num));
        }
    }
    //***********CYCLE END*******************
    std::for_each(walls.cbegin(), walls.cend(), [](const Surface& s){
        s.SaveSurfaceParticles();
    });

    return 0;
}

