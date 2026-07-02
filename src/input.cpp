#include "input.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace {

std::string trim(const std::string& s) {
    const auto a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    const auto b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

double toDouble(const std::string& s, const std::string& key) {
    try {
        size_t pos;
        double v = std::stod(s, &pos);
        if (trim(s.substr(pos)).empty()) return v;
    } catch (...) {}
    throw std::runtime_error("input: bad number for '" + key + "': " + s);
}

int toInt(const std::string& s, const std::string& key) {
    return static_cast<int>(toDouble(s, key));
}

bool toBool(const std::string& s, const std::string& key) {
    std::string v = s;
    std::transform(v.begin(), v.end(), v.begin(), ::tolower);
    if (v == "true" || v == "yes" || v == "on" || v == "1") return true;
    if (v == "false" || v == "no" || v == "off" || v == "0") return false;
    throw std::runtime_error("input: bad boolean for '" + key + "': " + s);
}

std::vector<double> toDoubleList(const std::string& s, const std::string& key) {
    std::istringstream ss(s);
    std::vector<double> out;
    double v;
    while (ss >> v) out.push_back(v);
    std::string rest;
    if (ss.clear(), ss >> rest)
        throw std::runtime_error("input: bad list for '" + key + "': " + s);
    return out;
}

}  // namespace

double DriveSpec::pressure(double t) const {
    if (table.empty()) return 0.0;
    if (t <= table.front().first) return table.front().second;
    if (t >= table.back().first) return table.back().second;
    for (size_t k = 1; k < table.size(); ++k) {
        if (t <= table[k].first) {
            const auto& [t0, p0] = table[k - 1];
            const auto& [t1, p1] = table[k];
            return p0 + (p1 - p0) * (t - t0) / (t1 - t0);
        }
    }
    return table.back().second;
}

InputDeck parseDeck(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("input: cannot open deck " + path);

    InputDeck deck;
    std::string section;      // current section name ("control", "material", ...)
    std::string sectionArg;   // e.g. material name
    std::string line;
    int lineno = 0;

    auto fail = [&](const std::string& msg) {
        throw std::runtime_error(path + ":" + std::to_string(lineno) + ": " + msg);
    };

    while (std::getline(in, line)) {
        ++lineno;
        auto hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        line = trim(line);
        if (line.empty()) continue;

        if (line.front() == '[') {
            if (line.back() != ']') fail("malformed section header");
            std::istringstream ss(line.substr(1, line.size() - 2));
            ss >> section;
            sectionArg.clear();
            std::getline(ss, sectionArg);
            sectionArg = trim(sectionArg);
            if (section == "material") {
                if (sectionArg.empty()) fail("[material] needs a name: [material NAME]");
                MaterialSpec m;
                m.name = sectionArg;
                deck.materials[sectionArg] = m;
            } else if (section == "layer") {
                deck.layers.emplace_back();
            } else if (section != "control" && section != "conduction" &&
                       section != "radiation" && section != "output" &&
                       section != "drive") {
                fail("unknown section [" + section + "]");
            }
            continue;
        }

        auto eq = line.find('=');
        if (eq == std::string::npos) fail("expected key = value");
        const std::string key = trim(line.substr(0, eq));
        const std::string val = trim(line.substr(eq + 1));
        if (key.empty() || val.empty()) fail("expected key = value");

        if (section == "control") {
            auto& c = deck.control;
            if      (key == "t_end")     c.t_end = toDouble(val, key);
            else if (key == "dt_init")   c.dt_init = toDouble(val, key);
            else if (key == "dt_max")    c.dt_max = toDouble(val, key);
            else if (key == "cfl")       c.cfl = toDouble(val, key);
            else if (key == "dt_growth") c.dt_growth = toDouble(val, key);
            else if (key == "max_steps") c.max_steps = static_cast<long>(toDouble(val, key));
            else if (key == "r_min")     c.r_min = toDouble(val, key);
            else if (key == "bc_outer")  c.bc_outer = val;
            else if (key == "c_quad")    c.c_quad = toDouble(val, key);
            else if (key == "c_lin")     c.c_lin = toDouble(val, key);
            else if (key == "T_floor")   c.T_floor = toDouble(val, key);
            else if (key == "temperatures") c.temperatures = toInt(val, key);
            else if (key == "geometry") {
                if      (val == "planar")      c.geometry = 1;
                else if (val == "cylindrical") c.geometry = 2;
                else if (val == "spherical")   c.geometry = 3;
                else fail("geometry must be planar|cylindrical|spherical");
            }
            else fail("unknown control key '" + key + "'");
        } else if (section == "conduction") {
            auto& c = deck.conduction;
            if      (key == "enabled")      c.enabled = toBool(val, key);
            else if (key == "flux_limiter") c.flux_limiter = toDouble(val, key);
            else if (key == "ln_lambda")    c.ln_lambda = (val == "auto") ? -1.0 : toDouble(val, key);
            else if (key == "ion_conduction")   c.ion_conduction = toBool(val, key);
            else if (key == "ion_flux_limiter") c.ion_flux_limiter = toDouble(val, key);
            else fail("unknown conduction key '" + key + "'");
        } else if (section == "radiation") {
            auto& rd = deck.radiation;
            if      (key == "enabled")  rd.enabled = toBool(val, key);
            else if (key == "bc_outer") rd.bc_outer = val;
            else fail("unknown radiation key '" + key + "'");
        } else if (section == "output") {
            auto& o = deck.output;
            if      (key == "directory")      o.directory = val;
            else if (key == "dt_dump")        o.dt_dump = toDouble(val, key);
            else if (key == "history_stride") o.history_stride = toInt(val, key);
            else fail("unknown output key '" + key + "'");
        } else if (section == "drive") {
            if (key == "table") {
                auto nums = toDoubleList(val, key);
                if (nums.empty() || nums.size() % 2 != 0)
                    fail("drive table must be pairs: t0 p0 t1 p1 ...");
                for (size_t k = 0; k + 1 < nums.size(); k += 2)
                    deck.drive.table.emplace_back(nums[k], nums[k + 1]);
                for (size_t k = 1; k < deck.drive.table.size(); ++k)
                    if (deck.drive.table[k].first <= deck.drive.table[k - 1].first)
                        fail("drive table times must be ascending");
            } else fail("unknown drive key '" + key + "'");
        } else if (section == "material") {
            auto& m = deck.materials[sectionArg];
            if      (key == "eos")   m.eos = val;
            else if (key == "gamma") m.gamma = toDouble(val, key);
            else if (key == "A")     m.A = toDouble(val, key);
            else if (key == "Z")     m.Z = toDouble(val, key);
            else if (key == "table") m.table_file = val;
            else if (key == "table_ion")      m.table_ion = val;
            else if (key == "table_electron") m.table_electron = val;
            else if (key == "ionization")     m.ionization = val;
            else if (key == "zbar_table")     m.zbar_table = val;
            else if (key == "opacity_table")  m.opacity_table = val;
            else if (key == "kappa_R")        m.kappa_R = toDouble(val, key);
            else if (key == "kappa_P")        m.kappa_P = toDouble(val, key);
            else fail("unknown material key '" + key + "'");
        } else if (section == "layer") {
            auto& l = deck.layers.back();
            if      (key == "material")  l.material = val;
            else if (key == "thickness") l.thickness = toDouble(val, key);
            else if (key == "zones")     l.zones = toInt(val, key);
            else if (key == "rho0")      l.rho0 = toDouble(val, key);
            else if (key == "T0")        l.T0 = toDouble(val, key);
            else if (key == "P0")        l.P0 = toDouble(val, key);
            else if (key == "Ti0")       l.Ti0 = toDouble(val, key);
            else if (key == "Te0")       l.Te0 = toDouble(val, key);
            else if (key == "Tr0")       l.Tr0 = toDouble(val, key);
            else if (key == "ratio")     l.ratio = toDouble(val, key);
            else fail("unknown layer key '" + key + "'");
        } else {
            fail("key '" + key + "' outside of any section");
        }
    }

    // --- validation -----------------------------------------------------
    if (deck.control.t_end <= 0.0) throw std::runtime_error("input: control.t_end must be > 0");
    if (deck.layers.empty()) throw std::runtime_error("input: at least one [layer] is required");
    for (size_t i = 0; i < deck.layers.size(); ++i) {
        const auto& l = deck.layers[i];
        const std::string tag = "layer " + std::to_string(i + 1);
        if (l.material.empty()) throw std::runtime_error("input: " + tag + " missing material");
        if (!deck.materials.count(l.material))
            throw std::runtime_error("input: " + tag + " references unknown material '" + l.material + "'");
        if (l.thickness <= 0.0) throw std::runtime_error("input: " + tag + " thickness must be > 0");
        if (l.zones <= 0) throw std::runtime_error("input: " + tag + " zones must be > 0");
        if (l.rho0 <= 0.0) throw std::runtime_error("input: " + tag + " rho0 must be > 0");
        if (l.T0 <= 0.0 && l.P0 <= 0.0 && !(l.Ti0 > 0.0 && l.Te0 > 0.0))
            throw std::runtime_error("input: " + tag +
                                     " needs T0 [eV], P0 [dyn/cm^2], or Ti0 and Te0");
        if (l.ratio <= 0.0) throw std::runtime_error("input: " + tag + " ratio must be > 0");
    }
    const bool twoT = (deck.control.temperatures == 2);
    if (deck.control.temperatures != 1 && deck.control.temperatures != 2)
        throw std::runtime_error("input: control.temperatures must be 1 or 2");
    for (const auto& [name, m] : deck.materials) {
        const std::string tag = "input: material " + name + ": ";
        if (m.eos != "ideal" && m.eos != "table")
            throw std::runtime_error(tag + "eos must be ideal|table");
        if (m.eos == "table" && !twoT && m.table_file.empty())
            throw std::runtime_error(tag + "table eos needs 'table = FILE'");
        if (m.eos == "table" && twoT && (m.table_ion.empty() || m.table_electron.empty()))
            throw std::runtime_error(tag + "table eos in 2T mode needs 'table_ion' "
                                     "and 'table_electron' files");
        if (m.A <= 0.0 || m.Z <= 0.0)
            throw std::runtime_error(tag + "A and Z must be > 0");
        if (m.ionization != "fixed" && m.ionization != "tf" && m.ionization != "table")
            throw std::runtime_error(tag + "ionization must be fixed|tf|table");
        if (m.ionization == "table" && m.zbar_table.empty())
            throw std::runtime_error(tag + "ionization = table needs 'zbar_table = FILE'");
        if (deck.radiation.enabled && m.opacity_table.empty() &&
            (m.kappa_R <= 0.0 || m.kappa_P <= 0.0))
            throw std::runtime_error(tag + "radiation is enabled; give 'opacity_table' "
                                     "or constant 'kappa_R' and 'kappa_P' [cm^2/g]");
    }
    if (deck.control.bc_outer != "wall" && deck.control.bc_outer != "pressure")
        throw std::runtime_error("input: control.bc_outer must be wall|pressure");
    if (deck.control.bc_outer == "pressure" && deck.drive.empty())
        throw std::runtime_error("input: bc_outer = pressure requires a [drive] table");
    if (deck.radiation.bc_outer != "insulated" && deck.radiation.bc_outer != "vacuum")
        throw std::runtime_error("input: radiation.bc_outer must be insulated|vacuum");

    return deck;
}
