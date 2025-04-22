////////////////////////////////////////////////////////////////////////
//
// ChannelROIToWire class - Converts from ChannelROI to Wire data
//
// usher@slac.stanford.edu
//
////////////////////////////////////////////////////////////////////////

// C/C++ standard libraries
#include <string>
#include <vector>
#include <utility> // std::pair<>
#include <memory> // std::unique_ptr<>

// framework libraries
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/WireReadout.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/CoreUtils/zip.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/ArtDataHelper/WireCreator.h"

#include "sbnobj/ICARUS/TPC/ChannelROI.h"

///creation of calibrated signals on wires
namespace caldata {

class ChannelROIToWire : public art::EDProducer
{
public:
// create calibrated signals on wires. this class runs 
// an fft to remove the electronics shaping.     
    explicit ChannelROIToWire(fhicl::ParameterSet const& pset);
    void     produce(art::Event& evt); 
    void     beginJob(); 
    void     endJob();                 
    void     reconfigure(fhicl::ParameterSet const& p);
private:

    std::vector<art::InputTag>                                 fWireModuleLabelVec;         ///< vector of modules that made digits
    std::vector<std::string>                                   fOutInstanceLabelVec;        ///< The output instance labels to apply
    bool                                                       fDiagnosticOutput;           ///< secret diagnostics flag
    size_t                                                     fEventCount;                 ///< count of event processed

  const geo::WireReadoutGeom* fChannelMapAlg = &art::ServiceHandle<geo::WireReadout const>()->Get();
    
}; // class ChannelROIToWire

DEFINE_ART_MODULE(ChannelROIToWire)

//-------------------------------------------------
ChannelROIToWire::ChannelROIToWire(fhicl::ParameterSet const& pset) : EDProducer{pset}
{
    this->reconfigure(pset);

    for(const auto& wireLabel : fOutInstanceLabelVec)
    {
        produces<std::vector<recob::Wire>>(wireLabel);
    }
}

//////////////////////////////////////////////////////
void ChannelROIToWire::reconfigure(fhicl::ParameterSet const& pset)
{
    // Recover the parameters
    fWireModuleLabelVec    = pset.get<std::vector<art::InputTag>>("WireModuleLabelVec",   std::vector<art::InputTag>()={"decon1droi"});
    fOutInstanceLabelVec   = pset.get<std::vector<std::string>>  ("OutInstanceLabelVec",                            {"PHYSCRATEDATA"});
    fDiagnosticOutput      = pset.get< bool                     >("DiagnosticOutput",                                           false);

    if (fWireModuleLabelVec.size() != fOutInstanceLabelVec.size()) 
    {
        throw art::Exception(art::errors::Configuration) << " Configured " << fOutInstanceLabelVec.size()
          << " instance names (`OutInstanceLabelVec`) for " << fWireModuleLabelVec.size()
          << " input products (`WireModuleLabelVec`)\n";
    }

    return;
}

//-------------------------------------------------
void ChannelROIToWire::beginJob()
{
    fEventCount = 0;
} // beginJob

//////////////////////////////////////////////////////
void ChannelROIToWire::endJob()
{
}

//////////////////////////////////////////////////////
void ChannelROIToWire::produce(art::Event& evt)
{
    // We need to loop through the list of Wire data we have been given
    // This construct from Gianluca Petrillo who invented it and should be given all credit for it! 
    for(auto const& [channelLabel, instanceName] : util::zip(fWireModuleLabelVec, fOutInstanceLabelVec))
    {
        // make a collection of Wires
        std::unique_ptr<std::vector<recob::Wire>> wireCol = std::make_unique<std::vector<recob::Wire>>();

        mf::LogInfo("ChannelROIToWire") << "ChannelROIToWire, looking for ChannelROI data at " << channelLabel.encode();
    
        // Read in the collection of full length deconvolved waveforms
       const std::vector<recob::ChannelROI>& channelVec = evt.getProduct<std::vector<recob::ChannelROI>>(channelLabel);

        mf::LogInfo("ChannelROIToWire") << "--> Recovered ChannelROI data, size: " << channelVec.size();
    
        if (!channelVec.empty())
        {
            // Reserve the room for the output
            wireCol->reserve(channelVec.size());

            // Loop through the input ChannelROI collection
            for(const auto& channelROI : channelVec)
            {
                // Recover the channel and the view
                raw::ChannelID_t channel        = channelROI.Channel();
                geo::View_t      view           = fChannelMapAlg->View(channel);
                float            ADCScaleFactor = channelROI.ADCScaleFactor();

                // Create an ROI vector for output
                recob::Wire::RegionsOfInterest_t ROIVec;
    
                // Loop through the ROIs for this channel
                const recob::ChannelROI::RegionsOfInterest_t& channelROIs = channelROI.SignalROI();

                for(const auto& range : channelROIs.get_ranges())
                {
                    size_t startTick = range.begin_index();

                    std::vector<float> dataVec(range.data().size());

                    for(size_t binIdx = 0; binIdx < range.data().size(); binIdx++) dataVec[binIdx] = std::round(range.data()[binIdx] * ADCScaleFactor);

                    ROIVec.add_range(startTick, std::move(dataVec));
                }

                wireCol->push_back(recob::WireCreator(std::move(ROIVec),channel,view).move());
            }

            mf::LogInfo("ChannelROIToWire") << "--> Outputting Wire data, size: " << wireCol->size() << " with instance name: " << instanceName;

            // Time to stroe everything
            if(wireCol->empty()) mf::LogWarning("ChannelROIToWire") << "No wires made for this event.";
        }

        evt.put(std::move(wireCol), instanceName);
    }

    fEventCount++;

    return;
} // produce

} // end namespace caldata
