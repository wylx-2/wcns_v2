#pragma once

#include "wcns_v2/core/config.h"
#include <string>

/// @file config_reader.h
/// @brief Reader for INI-format configuration files.
///
/// Parses a simple INI-style file with [section] headers and key = value pairs.
/// Supports comments starting with '#' or ';', blank lines, and inline comments.
///
/// Expected sections:
///   [physical]       — gamma, prandtl, reynolds, mach, aoa, beta
///   [reference]      — length, rho, temp, gas_constant
///   [control]        — cfl, max_iter, output_freq, restart_freq, time_scheme,
///                       converge_tol, ng
///   [initialization] — init_type, poiseuille_umax, poiseuille_y_min,
///                       poiseuille_y_max, body_force_x/y/z, wall_type
///
/// All values have sensible defaults; only non-default values need to be
/// specified in the config file.

class ConfigReader {
public:
    /// Read a config file and return a populated Config struct.
    /// Calls config.finalize() before returning.
    /// Throws std::runtime_error on parse errors or missing file.
    static Config read(const std::string& filename);

private:
    /// Trim leading and trailing whitespace from a string.
    static std::string trim(const std::string& s);

    /// Remove inline comments (text after '#' or ';' that is not inside quotes).
    static std::string strip_comment(const std::string& s);

    /// Parse a "key = value" line. Returns true on success.
    static bool parse_key_value(const std::string& line,
                                 std::string& key, std::string& value);

    /// Populate the physical section from a key/value pair.
    static void set_physical(Config& cfg, const std::string& key,
                              const std::string& value);

    /// Populate the reference section from a key/value pair.
    static void set_reference(Config& cfg, const std::string& key,
                               const std::string& value);

    /// Populate the control section from a key/value pair.
    static void set_control(Config& cfg, const std::string& key,
                             const std::string& value);

    /// Populate the initialization section from a key/value pair.
    static void set_initialization(Config& cfg, const std::string& key,
                                    const std::string& value);
};
