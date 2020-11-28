#include "loader.hpp"

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

std::unique_ptr<Surface> read_surface_parameters(const json& this_surf_data){
    std::string name = this_surf_data["name"].get<std::string>();
    std::vector<Vec3> contour;
    for(const auto& el : this_surf_data["contour"]){
        contour.push_back(Vec3(el.get<std::vector<double>>()));
    }
    std::string ref_type = this_surf_data["reflector_type"].get<std::string>();
    double R = this_surf_data["reflection_coefficient"].get<double>();
    bool stat_flag = this_surf_data["collect_statistics"].get<bool>();
    std::ofstream out_file;
    if(stat_flag){
        out_file.open(name, std::ios_base::app);
        if(!out_file.is_open()){fprintf(stderr, "could not open file\n"); exit(1);}
    }
    if(ref_type == "mirror"){
        return std::make_unique<Surface>(std::move(contour),
                   std::make_unique<MirrorReflector>(R),
                   std::move(out_file));
    }
    else if (ref_type == "cosine"){
        return std::make_unique<Surface>(std::move(contour),
              std::make_unique<LambertianReflector>(R),
                  std::move(out_file));
    }
    else {
        fprintf(stderr, "unknown reflector type %s", ref_type.c_str());
        exit(1);
    }
}

std::vector<std::unique_ptr<Surface>> load_geometry(const json& json_data){
    std::vector<std::unique_ptr<Surface>> walls;
    for(const auto& el : json_data["geometry"]){
        walls.push_back(read_surface_parameters(el));
    }
    if(!check_surface_orientations(walls)){
        fprintf(stderr, "Some surfaces has bad orientation. check contour numeration\n");
        exit(1);
    }
    return walls;
}

bool check_surface_orientations(const std::vector<std::unique_ptr<Surface>>& geo){
    for(const auto& s : geo){
        auto point = s->GetMassCenter();
        auto direction = s->GetNormal();
        auto found_crossetion = false;
        for(const auto& other_s : geo){
            if(other_s->GetCrossPoint(point, direction))
                found_crossetion=true;
        }
        if(!found_crossetion)
            return false;
    }
    return true;
}

