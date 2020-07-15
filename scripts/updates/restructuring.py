#!/usr/bin/env python
#
# This script changes C++ code, CMake files and FHiCL configuration.
# 
# The fundamental patterns are set into a table.
# Additional substitutions are possible (it's a Python script after all...)
# 
# 
# Change log:
# 20200714 (petrillo@slac.stanford.edu)
#   original version
#

import sys, os, re

import SerialSubstitution # from 'larsoft'
from SerialSubstitution import AddProcessor, RunSubstitutor


# main substitution table
# (in Python 3.7 this is supposed to be the stored order as well):
Substitutions = {
  
  # ----------------------------------------------------------------------------
  'TriggerDataObjects': {
    'source': 'icaruscode/PMT/Trigger/Data',
    'dest':   'sbnobj/ICARUS/PMT/Trigger/Data',
    }, # 'TriggerDataObjects'
  # ----------------------------------------------------------------------------
  
  # ----------------------------------------------------------------------------
  'SBNDCRTObjects': {
    'source':  'sbndcode/CRT/CRTProducts',
    'dest':    'sbnobj/SBND/CRT',
    'headers': [ 'CRTData.hh', ],
    'namespaces': [ ( 'crt', 'sbnd::crt', ( 'CRTData', ) ), ],
    'plugins': [ ], # plugins (including the plain library) will require manual intervention
    }, # 'SBNDCRTObjects'
  # ----------------------------------------------------------------------------
  
  # ----------------------------------------------------------------------------
  'CommonCRTObjectsForSBND': {
    'source':  'sbndcode/CRT/CRTProducts',
    'dest':    'sbnobj/Common/CRT',
    'headers': [ 'CRTHit.hh', 'CRTTrack.hh', 'CRTTzero.hh', ],
    'namespaces': [
      ( 'crt', 'sbn::crt', ( 'CRTHit', 'CRTTrack', 'CRTTzero', )  ),
      ( 'sbnd::crt', 'sbn::crt', ( 'CRTHit', 'CRTTrack', 'CRTTzero', )  ),
      ],
    }, # 'CommonCRTObjectsForSBND'
  # ----------------------------------------------------------------------------
  
  # ----------------------------------------------------------------------------
  'ICARUSCRTObjects': {
    'source':  'icaruscode/CRT/CRTProducts',
    'dest':    'sbnobj/ICARUS/CRT',
    'headers': [ 'CRTData.hh', ],
    'plugins': [ ], # plugins (including the plain library) will require manual intervention
    }, # 'ICARUSCRTObjects'
  # ----------------------------------------------------------------------------
  
  # ----------------------------------------------------------------------------
  'CommonCRTObjectsForICARUS': {
    'source':  'icaruscode/CRT/CRTProducts',
    'dest':    'sbnobj/Common/CRT',
    'headers': [ 'CRTHit.hh', 'CRTTrack.hh', 'CRTTzero.hh', ],
    'namespaces': [
      ( 'icarus::crt', 'sbn::crt', ( 'CRTHit', 'CRTTrack', 'CRTTzero', )  ),
      ],
    }, # 'CommonCRTObjectsForICARUS'
  # ----------------------------------------------------------------------------
  
  # ----------------------------------------------------------------------------
  'ICARUSPurityObjects': {
    'source':  'icaruscode/IcarusObj',
    'dest':    'sbnobj/Common/Analysis',
    }, # 'ICARUSPurityObjects'
  # ----------------------------------------------------------------------------
  
} # Substitutions


################################################################################
def intoLibraryPrefix(path): return path.replace('/', '_')

def intoHeaderGuardPrefix(path): return path.replace('/', '_').replace('.', '_').upper()


################################################################################
if __name__ == "__main__":

  #############################################################################
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # FHiCL configuration
  #
  Subst = AddProcessor(SerialSubstitution.ProcessorClass("FHiCL"))

  Subst.AddFileType("fcl")
  
  for rule in Substitutions.values():
    
    # this is for the rare cases where a plugin is specified by full path, e.g.
    # module_type: "icaruscode/Analysis/AnalysisTree";

    # if we have a list of plugins, we translate only those
    plugins = rule.get('plugins', [ None, ])
    
    baseSrc = rule['source']
    baseDest = rule['dest']

    for plugin in plugins:
      Subst.AddSimplePattern(
        baseSrc + "_" + plugin if plugin else baseSrc,
        baseDest + "_" + plugin if plugin else baseDest,
        )
    # for
    
  # for
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # CMakeLists.txt
  #
  Subst = AddProcessor(SerialSubstitution.ProcessorClass("cmake"))

  Subst.AddFileNamePattern("CMakeLists.txt")
  
  # here we parse the substitution table
  for rule in Substitutions.values():
    
    # this is for the modules and libraries; for example, a replacement like
    # 'larsim/LArG4' into 'larsim/LegacyLArG4' may yield
    # changes like 'larsim_LArG4_LArG4_module' to 'larsim_LegacyLArG4_LArG4_module'
    # and 'larsim_LArG4' to 'larsim_LegacyLArG4':
    if 'plugins' in rule:
      # if we have a list of plugins, we translate only those
      baseSrc = intoLibraryPrefix(rule['source'])
      baseDest = intoLibraryPrefix(rule['dest'])

      for plugin in rule['plugins']:
        Subst.AddSimplePattern(
          baseSrc + "_" + plugin if plugin else baseSrc,
          baseDest + "_" + plugin if plugin else baseDest,
          )
      # for
    else:
      # no plugin list: we translate all the possible libraries
      Subst.AddSimplePattern(
        intoLibraryPrefix(rule['source']),
        intoLibraryPrefix(rule['dest'])
        )
    # if ... else
    
  # for substitution table rules
  
  # special rule: CRTProducts libraries had non-standard name sbndcode_CRTProducts
  Subst.AddWord("sbndcode_CRTProducts", "sbnobj_Common_CRT")
  Subst.AddWord("icaruscode_CRTProducts", "sbnobj_Common_CRT")


  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # C++ source code (including modules and services)
  #
  Subst = AddProcessor(SerialSubstitution.ProcessorClass("code"))

  Subst.AddFileType("h", "hh", "cc", "cpp", "cxx", "icc", "tcc" )
  
  for rule in Substitutions.values():
    
    headers = rule.get('headers', [ '' ])
      
    for header in headers:
      
      baseSrc = os.path.join(rule['source'], header)
      baseDest = os.path.join(rule['dest'], header)
      
      # this is for the paths, including the ones in the '#include' directives
      # and the '@file' Doxygen documentation
      (Subst.AddWord if header else Subst.AddSimplePattern)(
        baseSrc,
        baseDest
        )
      
      # this is for standard header guard macros; for example, a replacement like
      # 'larsim/LArG4' into 'larsim/LegacyLArG4' may yield changes like
      # '#define LARSIM_LARG4_PARTICLELIST_H' into '#define LARSIM_LEGACYLARG4_PARTICLELIST_H'
      Subst.AddSimplePattern(
        intoHeaderGuardPrefix(baseSrc),
        intoHeaderGuardPrefix(baseDest)
        )
    # for header
    
    
    for nsSpec in rule.get('namespaces', []):
      
      try:               (srcNS, destNS), classes = nsSpec, [ None, ]
      except ValueError: srcNS, destNS, classes = nsSpec
      
      for className in classes:
        srcClass = "::".join((srcNS, className)) if className else (srcNS + "::")
        destClass = "::".join((destNS, className)) if className else (destNS + "::")
        Subst.AddRegExPattern(r'([^:]|^)\b{}\b'.format(srcClass), r'\1{}'.format(destClass))
      # for class name
      
    # for namespaces
    
    
  # for substitution table rules

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #############################################################################

  sys.exit(RunSubstitutor())
# 
