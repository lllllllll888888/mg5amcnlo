// Driver for Pythia 8. Reads an input file dynamically created on
// the basis of the inputs specified in MCatNLO_MadFKS_PY8.Script
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include "Pythia8Plugins/aMCatNLOHooks.h"
#include "Pythia8Plugins/CombineMatchingInput.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include <sstream>
#include <vector>

using namespace Pythia8;

int main() {
  Pythia pythia;

  // Register the q-cut variation list before reading the Pythia card so user
  // input "JetMatching:qCutList = 15 20 30 40" is accepted by the parser.
  // Empty default ⇒ legacy single-q-cut mode.
  pythia.settings.addPVec("JetMatching:qCutList", std::vector<double>(),
                          false, false, 0., 0.);

  string inputname="Pythia8.cmd",outputname="Pythia8.hep";

  pythia.readFile(inputname.c_str());

  //Create UserHooks pointer for the FxFX matching. Stop if it failed. Pass pointer to Pythia.
  CombineMatchingInput combined;
  //UserHooks* matching            = NULL;

  int nAbort=10;
  int nPrintLHA=1;
  int iAbort=0;
  int iPrintLHA=0;
  int iEventtot=pythia.mode("Main:numberOfEvents");
  int iEventshower=pythia.mode("Main:spareMode1");
  string evt_norm=pythia.word("Main:spareWord1");

  //FxFx merging
  bool isFxFx=pythia.flag("JetMatching:doFxFx");
  if (isFxFx) {
    combined.setHook(pythia);
    //matching = combined->getHook(pythia);
    //if (!matching) {
    //  std::cout << " Failed to initialise jet matching structures.\n"
    //            << " Program stopped.";
    //  return 1;
    //}
    //pythia.setUserHooksPtr(matching);
    int nJmax=pythia.mode("JetMatching:nJetMax");
    double Qcut=pythia.parm("JetMatching:qCut");
    double PTcut=pythia.parm("JetMatching:qCutME");
    if (Qcut <= PTcut || Qcut <= 0.) {
      std::cout << " \n";
      std::cout << "Merging scale (shower_card.dat) smaller than pTcut (run_card.dat)"
		<< Qcut << " " << PTcut << "\n";
      return 0;
    }
  }

  // q-cut variation: detect non-empty qCutList and grab a typed pointer to
  // the matching hook. Per-event accept flags (1 = kept, 0 = vetoed) for
  // each variation are written into the HepMC2 weight container under
  // names "FxFx_qCutAccept_<value>". These are NOT cross-section reweight
  // factors -- they are pure 0/1 flags. Cross-section at qcut X is
  //     sigma(X) = sum_evt  w_nominal[evt] * FxFx_qCutAccept_<X>[evt]
  // HepMC2 has no per-event attribute container, so the weight slot is the
  // only flexible per-event holder; the FxFx_qCutAccept_ prefix signals
  // that these are flags, not continuous reweights.
  std::vector<double> qCutList = pythia.settings.pvec("JetMatching:qCutList");
  bool doQCutVariation = isFxFx && !qCutList.empty();
  JetMatchingMadgraph* jmHook = nullptr;
  std::vector<std::string> qCutWeightNames;
  if (doQCutVariation) {
    jmHook = dynamic_cast<JetMatchingMadgraph*>(combined.hook.get());
    if (!jmHook) {
      std::cout << "Error: JetMatching:qCutList requires "
                << "JetMatchingMadgraph hook (FxFx scheme = 1)." << std::endl;
      return 1;
    }
    qCutWeightNames.reserve(qCutList.size());
    for (double q : qCutList) {
      std::ostringstream name;
      name << "FxFx_qCutAccept_" << q;
      qCutWeightNames.push_back(name.str());
    }
  }

  // Initialise Pythia.
  if (!pythia.init()) {
    cout << "Error: could not initialise Pythia" << endl;
    return 0;
  };

  HepMC::Pythia8ToHepMC ToHepMC;
  HepMC::IO_GenEvent ascii_io(outputname.c_str(), std::ios::out);
  // Do not store cross section information, as this will be done manually.
  ToHepMC.set_store_pdf(false);
  ToHepMC.set_store_proc(false);
  ToHepMC.set_store_xsec(false);

  // Cross section an error.
  double sigmaTotal  = 0.;
  double errorTotal  = 0.;

  for (int iEvent = 0; ; ++iEvent) {
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      break;
    }
    // the number of events read by Pythia so far
    int nSelected=pythia.info.nSelected();

    if (nSelected > iEventshower) break;
    if (pythia.info.isLHA() && iPrintLHA < nPrintLHA) {
      pythia.LHAeventList();
      pythia.info.list();
      pythia.process.list();
      pythia.event.list();
      ++iPrintLHA;
    }

    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    double evtweight = pythia.info.weight();
    double normhepmc;
    // ALWAYS NORMALISE HEPMC WEIGHTS TO SUM TO THE CROSS SECTION
    if (evt_norm != "sum") {
      normhepmc = 1. / double(iEventshower);
    } else {
      normhepmc = double(iEventtot) / double(iEventshower);
    }
    sigmaTotal += evtweight*normhepmc;
    hepmcevt->weights().push_back(evtweight*normhepmc);
    ToHepMC.fill_next_event( pythia, hepmcevt );
    // Add the weight of the current event to the cross section.
    // Report cross section to hepmc
    HepMC::GenCrossSection xsec;
    xsec.set_cross_section( sigmaTotal, pythia.info.sigmaErr() );
    hepmcevt->set_cross_section( xsec );

    // q-cut variation: push per-event accept flags as named HepMC weights
    // FxFx_qCutAccept_<q>. Pure 0/1 flags (1 = kept, 0 = vetoed at that qCut).
    // Cross-section at qcut X recovered as
    //     sigma(X) = sum_evt  w_nominal[evt] * FxFx_qCutAccept_<X>[evt]
    //
    // IMPORTANT: doVetoPartonLevelEarly is NOT called by Pythia for events
    // with no ME jets (e.g. Z+0j, where there's nothing to merge). For those
    // events, jmHook->getVetoVector() returns either empty (first event) or
    // STALE data from the previous hook call. To handle this safely, we
    // default every flag to 1 (accept) — Z+0j is unconditionally accepted at
    // any qCut — and only overwrite when the hook produced fresh decisions.
    //
    // Freshness signal: the JetMatching.h refactor pushes exactly
    // qCutListVec.size() entries when it runs. If getVetoVector().size()
    // doesn't match, the hook didn't fire for THIS event ⇒ keep defaults.
    if (doQCutVariation) {
      // Default: accept (Z+0j and any "hook didn't fire" case).
      for (size_t i = 0; i < qCutList.size(); ++i)
        hepmcevt->weights()[qCutWeightNames[i]] = 1.0;
      const std::vector<int>& flags = jmHook->getVetoVector();
      if (flags.size() == qCutList.size()) {
        for (size_t i = 0; i < flags.size(); ++i)
          hepmcevt->weights()[qCutWeightNames[i]] =
            static_cast<double>(1 - flags[i]);
      }
    }

    // Write the HepMC event to file. Done with it.
    ascii_io << hepmcevt;

    delete hepmcevt;
  }

  pythia.stat();
  if (isFxFx){
    std::cout << " \n";
    std::cout << "*********************************************************************** \n";
    std::cout << "*********************************************************************** \n";
    std::cout << "Cross section, including FxFx merging is: "
	      << sigmaTotal << "\n";
    std::cout << "*********************************************************************** \n";
    std::cout << "*********************************************************************** \n";
  }

  return 0;
}
