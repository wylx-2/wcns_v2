#include "wcns_v2/io/config_reader.h"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

// ============================================================================
// Public interface
// ============================================================================

Config ConfigReader::read(const std::string& filename) {
    Config cfg;  // start with defaults

    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("ConfigReader: cannot open file \"" +
                                 filename + "\"");
    }

    std::string section;  // current section (lowercase)
    int line_no = 0;

    std::string line;
    while (std::getline(file, line)) {
        ++line_no;

        // Strip comment and trim
        line = strip_comment(line);
        line = trim(line);

        // Skip blank lines
        if (line.empty()) continue;

        // Check for section header: [section]
        if (line.front() == '[' && line.back() == ']') {
            section = trim(line.substr(1, line.size() - 2));
            // Lowercase for case-insensitive matching
            std::transform(section.begin(), section.end(),
                           section.begin(),
                           [](unsigned char c) { return std::tolower(c); });
            continue;
        }

        // Parse key = value
        std::string key, value;
        if (!parse_key_value(line, key, value)) {
            std::cerr << "ConfigReader: warning: skipping malformed line "
                      << line_no << ": \"" << line << "\"\n";
            continue;
        }

        // Lowercase key and section for matching
        std::string key_lower = key;
        std::transform(key_lower.begin(), key_lower.end(),
                       key_lower.begin(),
                       [](unsigned char c) { return std::tolower(c); });

        // Dispatch based on section
        if (section == "physical") {
            set_physical(cfg, key_lower, value);
        } else if (section == "reference") {
            set_reference(cfg, key_lower, value);
        } else if (section == "control") {
            set_control(cfg, key_lower, value);
        } else if (section == "scheme") {
            set_scheme(cfg, key_lower, value);
        } else if (section == "initialization") {
            set_initialization(cfg, key_lower, value);
        } else if (section == "output") {
            set_output(cfg, key_lower, value);
        } else {
            std::cerr << "ConfigReader: warning: unknown section \"[" << section
                      << "]\" at line " << line_no << ", ignoring\n";
        }
    }

    file.close();

    // Compute derived quantities
    cfg.finalize();

    std::cout << "Config file \"" << filename << "\" loaded successfully.\n";
    return cfg;
}

// ============================================================================
// Private helpers — string utilities
// ============================================================================

std::string ConfigReader::trim(const std::string& s) {
    auto start = s.begin();
    while (start != s.end() && std::isspace(static_cast<unsigned char>(*start)))
        ++start;

    auto end = s.end();
    do {
        --end;
    } while (std::distance(start, end) > 0 &&
             std::isspace(static_cast<unsigned char>(*end)));

    return std::string(start, end + 1);
}

std::string ConfigReader::strip_comment(const std::string& s) {
    // Find first '#' or ';' outside of quotes
    bool in_quote = false;
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '"') in_quote = !in_quote;
        if (!in_quote && (s[i] == '#' || s[i] == ';')) {
            return s.substr(0, i);
        }
    }
    return s;
}

bool ConfigReader::parse_key_value(const std::string& line,
                                     std::string& key, std::string& value) {
    auto eq_pos = line.find('=');
    if (eq_pos == std::string::npos) return false;

    key   = trim(line.substr(0, eq_pos));
    value = trim(line.substr(eq_pos + 1));

    return !key.empty();
}

// ============================================================================
// Section dispatchers
// ============================================================================

void ConfigReader::set_physical(Config& cfg, const std::string& key,
                                 const std::string& value) {
    if (key == "gamma") {
        cfg.gamma = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "prandtl") {
        cfg.Prandtl = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "reynolds" || key == "re") {
        cfg.Re = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "mach") {
        cfg.Mach = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "aoa" || key == "angle_of_attack") {
        cfg.AoA = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "beta" || key == "sideslip") {
        cfg.beta = static_cast<Real>(std::atof(value.c_str()));
    } else {
        std::cerr << "ConfigReader: warning: unknown key \""
                  << key << "\" in [physical] section\n";
    }
}

void ConfigReader::set_reference(Config& cfg, const std::string& key,
                                  const std::string& value) {
    if (key == "length" || key == "l_ref") {
        cfg.L_ref = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "rho" || key == "rho_ref") {
        cfg.rho_ref = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "temp" || key == "t_ref") {
        cfg.T_ref = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "gas_constant" || key == "r_gas") {
        cfg.R_gas = static_cast<Real>(std::atof(value.c_str()));
    } else {
        std::cerr << "ConfigReader: warning: unknown key \""
                  << key << "\" in [reference] section\n";
    }
}

void ConfigReader::set_control(Config& cfg, const std::string& key,
                                const std::string& value) {
    if (key == "cfl") {
        cfg.cfl = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "fixed_dt") {
        cfg.fixed_dt = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "max_time") {
        cfg.max_time = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "max_iter") {
        cfg.max_iter = static_cast<Int>(std::atoi(value.c_str()));
    } else if (key == "output_freq") {
        cfg.output_freq = static_cast<Int>(std::atoi(value.c_str()));
    } else if (key == "restart_freq") {
        cfg.restart_freq = static_cast<Int>(std::atoi(value.c_str()));
    } else if (key == "residual_freq") {
        cfg.residual_freq = static_cast<Int>(std::atoi(value.c_str()));
    } else if (key == "converge_tol") {
        cfg.converge_tol = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "time_scheme") {
        cfg.time_scheme = value;  // keep as-is
    } else if (key == "lu_sgs_kappa") {
        cfg.lu_sgs_kappa = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "ng" || key == "ghost_layers") {
        cfg.ng = static_cast<Int>(std::atoi(value.c_str()));
    } else if (key == "restart_file") {
        cfg.restart_file = value;  // basename for restart (e.g. "restart_000500")
    } else {
        std::cerr << "ConfigReader: warning: unknown key \""
                  << key << "\" in [control] section\n";
    }
}

void ConfigReader::set_scheme(Config& cfg, const std::string& key,
                               const std::string& value) {
    if (key == "grid_metrics") {
        cfg.grid_metrics = value;
    } else if (key == "interp_vars") {
        cfg.interp_vars = value;
    } else if (key == "interp_type") {
        cfg.interp_type = value;
    } else if (key == "mdcd_diss") {
        cfg.mdcd_diss = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "mdcd_disp") {
        cfg.mdcd_disp = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "mdcd_sai_ref") {
        cfg.mdcd_sai_ref = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "riemann_type") {
        cfg.riemann_type = value;
    } else if (key == "entropy_fix_eps") {
        cfg.entropy_fix_eps = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "viscous_type") {
        cfg.viscous_type = value;
    } else if (key == "mu_const") {
        cfg.mu_const = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "sutherland_s") {
        cfg.sutherland_S = static_cast<Real>(std::atof(value.c_str()));
    } else {
        std::cerr << "ConfigReader: warning: unknown key \""
                  << key << "\" in [scheme] section\n";
    }
}

void ConfigReader::set_initialization(Config& cfg, const std::string& key,
                                       const std::string& value) {
    if (key == "init_type") {
        cfg.init_type = value;
    } else if (key == "poiseuille_umax") {
        cfg.poiseuille_umax = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "poiseuille_y_min") {
        cfg.poiseuille_y_min = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "poiseuille_y_max") {
        cfg.poiseuille_y_max = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "body_force_type") {
        cfg.body_force_type = value;
    } else if (key == "body_force_x") {
        cfg.body_force_x = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "body_force_y") {
        cfg.body_force_y = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "body_force_z") {
        cfg.body_force_z = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "body_force_energy") {
        cfg.body_force_energy = (value == "true" || value == "1" || value == "yes");
    } else if (key == "wall_type") {
        cfg.wall_type = value;
    } else if (key == "wall_temperature") {
        cfg.wall_temperature = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "riemann_config") {
        cfg.riemann_config = static_cast<Int>(std::atoi(value.c_str()));
    } else if (key == "riemann_x_split") {
        cfg.riemann_x_split = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "riemann_y_split") {
        cfg.riemann_y_split = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "isentropic_vortex_strength") {
        cfg.isentropic_vortex_strength = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "isentropic_vortex_radius") {
        cfg.isentropic_vortex_radius = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "isentropic_vortex_xc") {
        cfg.isentropic_vortex_xc = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "isentropic_vortex_yc") {
        cfg.isentropic_vortex_yc = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "isentropic_vortex_u_inf") {
        cfg.isentropic_vortex_u_inf = static_cast<Real>(std::atof(value.c_str()));
    } else if (key == "isentropic_vortex_v_inf") {
        cfg.isentropic_vortex_v_inf = static_cast<Real>(std::atof(value.c_str()));
    } else {
        std::cerr << "ConfigReader: warning: unknown key \""
                  << key << "\" in [initialization] section\n";
    }
}

void ConfigReader::set_output(Config& cfg, const std::string& key,
                               const std::string& value) {
    if (key == "output_format") {
        cfg.output_format = value;
    } else if (key == "output_dir") {
        cfg.output_dir = value;
    } else if (key == "output_time_interval") {
        cfg.output_time_interval = static_cast<Real>(std::atof(value.c_str()));
    } else {
        std::cerr << "ConfigReader: warning: unknown key \""
                  << key << "\" in [output] section\n";
    }
}
