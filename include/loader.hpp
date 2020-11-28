#ifndef LOADER_HPP
#define LOADER_HPP

#include <nlohmann/json.hpp>
#include <memory>
#include <vector>

#include "particle.hpp"
#include "surface.hpp"

using json = nlohmann::json;

json load_json_config(const std::string& file_name);
Background load_background(const json& json_data);
std::unique_ptr<Surface> read_surface_parameters(const json& this_surf_data);
std::vector<std::unique_ptr<Surface>> load_geometry(const json& json_data);
bool check_surface_orientations(const std::vector<std::unique_ptr<Surface>>& geo);

#endif //LOADER_HPP
