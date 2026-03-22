#include "iqtree_interface.h"

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_map>

namespace clnj {

static const std::unordered_map<std::string, std::string> IQTREE_MODEL_MAP = {
    {"JC",  "JC69"}, {"JC+I",  "JC69"}, {"JC+G",  "JC69"}, {"JC+I+G",  "JC69"},
    {"K2P", "K2P"},  {"K2P+I", "K2P"},  {"K2P+G", "K2P"},  {"K2P+I+G", "K2P"},
    {"K80", "K2P"},  {"K80+I", "K2P"},  {"K80+G", "K2P"},  {"K80+I+G", "K2P"},
    {"HKY", "TN93"}, {"HKY+I", "TN93"}, {"HKY+G", "TN93"}, {"HKY+I+G", "TN93"},
    {"TN",  "TN93"}, {"TN+I",  "TN93"}, {"TN+G",  "TN93"}, {"TN+I+G",  "TN93"},
    {"TRN", "TN93"}, {"TRN+I", "TN93"}, {"TRN+G", "TN93"}, {"TRN+I+G", "TN93"},
    {"TNE", "TN93"}, {"TNE+I", "TN93"}, {"TNE+G", "TN93"}, {"TNE+I+G", "TN93"},
    {"SYM", "TN93"}, {"SYM+I", "TN93"}, {"SYM+G", "TN93"}, {"SYM+I+G", "TN93"},
    {"GTR", "TN93"}, {"GTR+I", "TN93"}, {"GTR+G", "TN93"}, {"GTR+I+G", "TN93"},
};

IqtreeResult run_iqtree_model_selection(const std::string& fasta_path, bool verbose) {
    // Create temporary directory
    char tmpdir_template[] = "/tmp/clnj_iqtree_XXXXXX";
    char* tmpdir = mkdtemp(tmpdir_template);
    if (!tmpdir) {
        std::cerr << "  ERROR: Failed to create temp directory\n";
        return {};
    }
    std::string prefix = std::string(tmpdir) + "/iqtree_mf";

    std::string cmd = "iqtree2 -s " + fasta_path
        + " -m MFP --mset JC,K2P,HKY,TN,GTR --mrate E,I,G,I+G"
        + " -n 0 --prefix " + prefix + " -T AUTO --safe"
        + " > /dev/null 2>&1";

    if (verbose)
        std::cout << "  IQ-TREE cmd: iqtree2 -s " << fasta_path
                  << " -m MFP --mset JC,K2P,HKY,TN,GTR --mrate E,I,G,I+G"
                  << " -n 0 --prefix " << prefix << " -T AUTO --safe\n";

    int ret = std::system(cmd.c_str());
    if (ret != 0) {
        std::cerr << "  WARNING: iqtree2 returned exit code " << ret << "\n";
    }

    // Parse the .iqtree output file
    IqtreeResult result;
    std::string iqtree_file = prefix + ".iqtree";
    std::ifstream in(iqtree_file);
    if (!in.is_open()) {
        std::cerr << "  ERROR: Cannot open " << iqtree_file << "\n";
        return result;
    }

    std::string line;
    std::regex re_model(R"(Best-fit model according to BIC:\s+(.+))");
    std::regex re_bic(R"(Bayesian information criterion \(BIC\) score:\s+([\d.]+))");
    std::regex re_gamma(R"(Gamma shape alpha:\s+([\d.]+))");
    std::smatch match;

    while (std::getline(in, line)) {
        if (std::regex_search(line, match, re_model))
            result.iqtree_model = match[1].str();
        if (std::regex_search(line, match, re_bic))
            result.bic = std::stod(match[1].str());
        if (std::regex_search(line, match, re_gamma))
            result.gamma_alpha = std::stod(match[1].str());
    }

    // Map IQ-TREE model to pipeline model
    std::string base_model = result.iqtree_model;
    // Remove {xxx} rate class notation
    auto brace = base_model.find('{');
    if (brace != std::string::npos) {
        auto brace_end = base_model.find('}', brace);
        if (brace_end != std::string::npos)
            base_model = base_model.substr(0, brace) + base_model.substr(brace_end + 1);
    }
    // Remove +F
    auto fpos = base_model.find("+F");
    if (fpos != std::string::npos)
        base_model = base_model.substr(0, fpos) + base_model.substr(fpos + 2);
    // Remove +G4, +I segments for mapping but check presence
    bool has_gamma = (result.iqtree_model.find("+G") != std::string::npos);
    if (!has_gamma) result.gamma_alpha = -1.0;

    // Uppercase for matching
    std::string upper_base = base_model;
    for (auto& c : upper_base) c = (char)std::toupper(c);

    result.pipeline_model = "TN93";  // default fallback
    // Match longest key first
    std::vector<std::pair<std::string, std::string>> sorted_map(IQTREE_MODEL_MAP.begin(), IQTREE_MODEL_MAP.end());
    std::sort(sorted_map.begin(), sorted_map.end(),
              [](auto& a, auto& b) { return a.first.size() > b.first.size(); });
    for (auto& [key, val] : sorted_map) {
        if (upper_base.find(key) == 0) {
            result.pipeline_model = val;
            break;
        }
    }

    if (verbose) {
        std::cout << "  IQ-TREE best model: " << result.iqtree_model
                  << "  (BIC=" << result.bic << ")\n";
        std::cout << "  Mapped to pipeline: " << result.pipeline_model;
        if (result.has_gamma())
            std::cout << "  +Gamma(alpha=" << result.gamma_alpha << ")";
        std::cout << "\n";
    }

    // Cleanup temp directory
    std::string cleanup = "rm -rf " + std::string(tmpdir);
    std::system(cleanup.c_str());

    return result;
}

}  // namespace clnj
