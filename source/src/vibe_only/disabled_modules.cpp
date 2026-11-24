// TITAN-VIBE (VIBE_ONLY=1): the members of the modules this build leaves out
// that the engine calls on every run. The engine holds one object of each
// analysis module and asks it whether it is enabled; here every module is
// built disabled and stays so, because option_allowlist.cpp refuses the options
// that would enable it. All other entry points of these modules are in
// absent_entry_points.cpp and stop the program if called (see absent.h).

#include "vibe_only/absent.h"

#include "cas_analysis.h"
#include "ml_extractor.h"
#include "tda.h"
#include "titan_energy_decomp_engine.h"
#include "titan_pbe.h"
#include "water_analysis.h"

// The CAS, feature-extraction, TDA and water analyses: objects with no work.
CASAnalyzer::CASAnalyzer() : global_total_frames(0), has_reference(false) {}
MLFeatureExtractor::MLFeatureExtractor() = default;
TDAAnalyzer::~TDAAnalyzer() = default;
WaterAnalyzer::WaterAnalyzer() = default;

// The Poisson-Boltzmann analyzer: TITAN-VIBE has no solver behind it.
namespace TitanPBE {

class AnalyzerImpl {};

Analyzer::Analyzer() = default;
Analyzer::~Analyzer() = default;
void Analyzer::configure(const arguments &) {}
bool Analyzer::enabled() const noexcept { return false; }

} // namespace TitanPBE

// The energy decompositions: the engine holds them through std::unique_ptr and
// never creates one here. Their constructors and destructors are defined here
// rather than in absent_entry_points.cpp because both need the accumulators as
// complete (here empty) types.
namespace titan_energy_decomp {

struct EnergyDecompositionEngine::Accumulators {};
struct EndpointEntityDecompositionEngine::Accumulators {};

EnergyDecompositionEngine::EnergyDecompositionEngine(EntityAxes, EntityAxes, EngineConfig) {
  titan_vibe_only::absent("titan_energy_decomp::EnergyDecompositionEngine::EnergyDecompositionEngine");
}
EndpointEntityDecompositionEngine::EndpointEntityDecompositionEngine(EntityAxes, EntityAxes, EngineConfig) {
  titan_vibe_only::absent(
      "titan_energy_decomp::EndpointEntityDecompositionEngine::EndpointEntityDecompositionEngine");
}
EnergyDecompositionEngine::~EnergyDecompositionEngine() = default;
EndpointEntityDecompositionEngine::~EndpointEntityDecompositionEngine() = default;

} // namespace titan_energy_decomp
