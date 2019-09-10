#include "SelectionToolBase.hh"

namespace selection {

  SelectionToolBase::SelectionToolBase() : core::PostProcessorBase() {}

  SelectionToolBase::~SelectionToolBase() {}

  bool SelectionToolBase::ProcessEvent(const gallery::Event &ev){
    selection::Event e = SelectionToolBase::GetSelectionToolEvent(ev);
    SelectionToolBase::ProcessEvent(e);
  } // ProcessEvent

  selection::Event SelectionToolBase::GetSelectionToolEvent(const gallery::Event &ev) {

    // Define the sbncode event details and event list to output
    std::vector<core::Event::Interaction> &truth    = ev.truth;
    std::vector<core::Event::RecoInteraction> &reco = ev.reco;

    // Call functions from LoadEvents
    // This is replacing 'EventSelectionHelper::LoadEventList
    // and will process a single event at a time
    //
    // Take an SBNCode processed file with (eventually) 
    // a reconstructed sample of events selected from 
    // external backgrounds
    //
    // Fill selection::Events with core::Event information

    // Fill event-based variables
    ParticleList mc_particles;       ///< vector of Monte Carlo particles
    ParticleList reco_particles;     ///< vector of reconstructed particles
    unsigned int interaction;        ///< interaction type of the event
    unsigned int scatter;            ///< scatter code for the event: physical process
    int          nu_pdg;             ///< Neutrino pdg code of the event
    int          init_pdg;           ///< Initial neutrino pdg code of the event
    bool         is_cc;              ///< whether the event contains and CC or NC interaction
    TVector3     reco_vertex;        ///< reconstructed neutrino vertex
    TVector3     mc_vertex;          ///< reconstructed neutrino vertex
    float        neutrino_energy;    ///< true neutrino energy
    float        neutrino_qsqr;      ///< true neutrino qsqr

    interaction     = truth.neutrino.gint;
    scatter         = truth.neutrino.gscatter;
    nu_pdg          = truth.neutrino.pdg;
    init_pdg        = truth.neutrino.initpdg;
    is_cc           = truth.neutrino.iscc;
    mc_vertex       = truth.neutrino.position;
    neutrino_energy = truth.neutrino.energy;
    neutrino_qsqr   = truth.neutrino.Q2;
    reco_vertex     = reco.reco_vertex;

    // Now deal with the individual particles
    
    return selection::Event(mc_particles,   
                            reco_particles, 
                            interaction,    
                            scatter,        
                            nu_pdg,         
                            init_pdg,       
                            is_cc,          
                            reco_vertex,    
                            mc_vertex,      
                            neutrino_energy,
                            neutrino_qsqr);
  
  } // GetSelectionToolEvent

}  // namespace core


