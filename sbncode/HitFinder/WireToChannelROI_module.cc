////////////////////////////////////////////////////////////////////////
//
// WireToChannelROI class - Convert recob::Wire to recob::ChannelROI
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

#include "sbnobj/ICARUS/TPC/ChannelROI.h"
#include "sbnobj/Common/Utilities/ChannelROICreator.h"

///creation of calibrated signals on wires
namespace caldata {

class WireToChannelROI : public art::EDProducer
{
public:
// create calibrated signals on wires. this class runs 
// an fft to remove the electronics shaping.     
    explicit WireToChannelROI(fhicl::ParameterSet const& pset);
    void     produce(art::Event& evt); 
    void     beginJob(); 
    void     endJob();                 
    void     reconfigure(fhicl::ParameterSet const& p);
private:

    std::vector<art::InputTag>  fWireModuleLabelVec;         ///< vector of modules that made digits
    std::vector<std::string>    fOutInstanceLabelVec;        ///< The output instance labels to apply
    bool                        fDiagnosticOutput;           ///< secret diagnostics flag
    size_t                      fEventCount;                 ///< count of event processed

  const geo::WireReadoutGeom* fChannelMapAlg = &art::ServiceHandle<geo::WireReadout const>()->Get();
    
}; // class WireToChannelROI

DEFINE_ART_MODULE(WireToChannelROI)

//-------------------------------------------------
WireToChannelROI::WireToChannelROI(fhicl::ParameterSet const& pset) : EDProducer{pset}
{
    this->reconfigure(pset);

    for(const auto& wireLabel : fOutInstanceLabelVec)
    {
        produces<std::vector<recob::ChannelROI>>(wireLabel);
    }
}

//////////////////////////////////////////////////////
void WireToChannelROI::reconfigure(fhicl::ParameterSet const& pset)
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
void WireToChannelROI::beginJob()
{
    fEventCount = 0;
} // beginJob

//////////////////////////////////////////////////////
void WireToChannelROI::endJob()
{
}

//////////////////////////////////////////////////////
void WireToChannelROI::produce(art::Event& evt)
{
    float ADCScaleFactor(recob::ChannelROI::defADCScaleFactor);

    // We need to loop through the list of Wire data we have been given
    // This construct from Gianluca Petrillo who invented it and should be given all credit for it! 
    for(auto const& [wireLabel, instanceName] : util::zip(fWireModuleLabelVec, fOutInstanceLabelVec))
    {
        // make a collection of Wires
        std::unique_ptr<std::vector<recob::ChannelROI>> channelROICol = std::make_unique<std::vector<recob::ChannelROI>>();

        mf::LogInfo("WireToChannelROI") << "WireToChannelROI, looking for Wire data at " << wireLabel.encode();
    
        // Read in the collection of full length deconvolved waveforms
       const std::vector<recob::Wire>& wireVec = evt.getProduct<std::vector<recob::Wire>>(wireLabel);

        mf::LogInfo("WireToChannelROI") << "--> Recovered Wire data, size: " << wireVec.size();
    
        if (!wireVec.empty())
        {
            // Reserve the room for the output
            channelROICol->reserve(wireVec.size());

            // Loop through the input ChannelROI collection
            for(const auto& wire : wireVec)
            {
                // Recover the channel and the view
                raw::ChannelID_t channel = wire.Channel();

                // Create an ROI vector for output
                recob::ChannelROI::RegionsOfInterest_t ROIVec;

                // Return the entire waveform (zero padded) so we can find maximum range
                std::vector<float> fullWaveform = wire.Signal();

                std::pair<std::vector<float>::const_iterator,std::vector<float>::const_iterator> minMaxPair = std::minmax_element(fullWaveform.begin(),fullWaveform.end());

                float maxValue = std::max(-*minMaxPair.first,*minMaxPair.second);

                ADCScaleFactor = std::max(std::min(float(512.),static_cast<float>(std::numeric_limits<short>::max())/maxValue),float(8.));
    
                // Loop through the ROIs for this channel
                const recob::Wire::RegionsOfInterest_t& wireROIs = wire.SignalROI();

                for(const auto& range : wireROIs.get_ranges())
                {
                    size_t startTick = range.begin_index();

                    std::vector<short int> dataVec(range.data().size());

                    for(size_t binIdx = 0; binIdx < range.data().size(); binIdx++) 
                    {
                        short int scaledIntADCVal = 
                            std::min(std::max(std::round(range.data()[binIdx] * ADCScaleFactor),static_cast<float>(std::numeric_limits<short>::min())),
                                                                                                static_cast<float>(std::numeric_limits<short>::max()));

                        dataVec[binIdx] = scaledIntADCVal;
                    }

                    ROIVec.add_range(startTick, std::move(dataVec));
                }

                channelROICol->push_back(recob::ChannelROICreator(std::move(ROIVec),channel,ADCScaleFactor).move());
            }

            mf::LogInfo("WireToChannelROI") << "--> Outputting ChannelROI data, size: " << channelROICol->size() << " with instance name: " << instanceName;

            // Time to stroe everything
            if(channelROICol->empty()) mf::LogWarning("WireToChannelROI") << "No wires made for this event.";
        }

        evt.put(std::move(channelROICol), instanceName);
    }

    fEventCount++;

    return;
} // produce

} // end namespace caldata
