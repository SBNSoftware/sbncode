#ifndef __sbnanalysis_selection_SelectionToolBase__
#define __sbnanalysis_selection_SelectionToolBase__

/**
 * Takes a gallery::Event and allows you to process the selection::Event format instead
 *
 * Author: R. Jones <rjones@hep.ph.liv.ac.uk>, 2019/05/30
 */

#include "GeneralAnalysisHelper.hh"
#include "EventSelectionHelper.hh"
#include "LoadEvents.hh"

#include "../core/ProcessorBase.hh"

namespace selection {

/**
 * \class core::SelectionToolBase
 * \brief Base class for translating the ProcessEvent function definition from 
 *        taking gallery events to taking selection (SelectionTool format) events in analysis modules
 *
 * See core::ProcessorBase for more details.
 */
class SelectionToolBase : public core::ProcessorBase {
public:
  /** Constructor. */
  SelectionToolBase();

  /** Destructor. */
  virtual ~SelectionToolBase();

  /**
   * \brief ProcessEvent using the SBNCode definition of the Event object
   *
   * \param  ev A single event, as a gallery::Event
   *
   * \return True to keep event
   */
  bool core::ProcessEvent(const gallery::Event& ev);
  
  /**
   * \brief  GetSelectionToolEvent to translate the definition of Event from SBNCode to the SelectionTool
   *
   * \param  gallery::Event the gallery Event as read by SBNCode
   *
   * \return selection::Event event format for analysis of the selection 
   */
  selection::Event GetSelectionToolEvent(const gallery::Event &ev);

  /**
   * \brief  ProcessEvent virtual function to process selection::Events
   *
   * \param  selection::Event selection tool event definition for analyses
   *
   * \return true if the event is valid for analysis
   *
   */
  virtual bool ProcessEvent(const selection::Event &ev) = 0; 

};

}  // namespace core

#endif  // __sbnanalysis_core_SelectionBase__

