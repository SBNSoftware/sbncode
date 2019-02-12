#!/usr/bin/env python
#----------------------------------------------------------------------
#
# Name: project_utilities.py
#
# Purpose: A python module containing various experiment-specific
#          python utility functions.
#
# Created: 28-Oct-2013  H. Greenlee
# Modified: 10-Jul-2017 T. Brooks (tbrooks@fnal.gov) - now works with 
#           LArSoft v06
# Modified: 28-Jan-2018 D.Barker (dominic.barker@sheffield.ac.uk) 
#           Altered sbnd version to work as an sbncode version
#----------------------------------------------------------------------

import os, pycurl
from StringIO import StringIO

# Don't fail (on import) if samweb is not available.

try:
    import samweb_cli
except ImportError:
    pass

def get_dropbox(filename):

    # Get metadata.

    md = {}
    exp = 'sbnd'
    if os.environ.has_key('SAM_EXPERIMENT'):
        exp = os.environ['SAM_EXPERIMENT']
    samweb = samweb_cli.SAMWebClient(experiment=exp)
    try:
        md = samweb.getMetadata(filenameorid=filename)
    except:
        pass

    # Extract the metadata fields that we need.
    
    file_type = ''
    group = ''
    data_tier = ''

    if md.has_key('file_type'):
        file_type = md['file_type']
    if md.has_key('group'):
        group = md['group']
    if md.has_key('data_tier'):
        data_tier = md['data_tier']

    if not file_type or not group or not data_tier:
        raise RuntimeError, 'Missing or invalid metadata for file %s.' % filename

    # Construct dropbox path.

    #path = '/sbnd/data/sbndsoft/dropbox/%s/%s/%s' % (file_type, group, data_tier)
    if os.environ.has_key('FTS_DROPBOX'):
        dropbox_root = os.environ['FTS_DROPBOX']
    else:
        dropbox_root = '/pnfs/sbnd/scratch/sbndpro/dropbox'
    path = '%s/%s/%s/%s' % (dropbox_root, file_type, group, data_tier)
    return path

# Return fcl configuration for experiment-specific sam metadata.

def get_sam_metadata(project, stage):
    result = 'services.FileCatalogMetadataSBND: {\n'
    if type(stage.fclname) == type('') or type(stage.fclname) == type(u''):
        result = result + '  FCLName: "%s"\n' % os.path.basename(stage.fclname)
    else:
        result = result + '  FCLName: "'
        for fcl in stage.fclname:
            result = result + '%s/' % os.path.basename(fcl)
        result = result[:-1]
        result = result + '"\n' 
    result = result + '  FCLVersion: "%s"\n' % project.release_tag
    result = result + '  ProjectName: "%s"\n' % project.name
    result = result + '  ProjectStage: "%s"\n' % stage.name
    result = result + '  ProjectVersion: "%s"\n' % project.release_tag
    result = result + '}\n'

    return result

# Function to return path to the setup_sbnd.sh script

def get_setup_script_path():

    #Got to cheat here and look for the setup script locally because this is not setup yet 
    SRCS = os.environ["MRB_SOURCE"]
    SBNUTIL  = "/sbncode/sbnutil/"
    SCRIPT_DIR = SRCS + SBNUTIL

    if os.path.isfile(SCRIPT_DIR+"setup_sbncode_grid-v08_03_00.sh"):
        setup_script = SCRIPT_DIR+"setup_sbncode_grid-v08_03_00.sh"
    else:
        raise RuntimeError, "Could not find setup script at "+SCRIPT_DIR
    
    return setup_script

# Construct dimension string for project, stage.

def dimensions(project, stage):

    data_tier = ''
    if ana:
        data_tier = stage.ana_data_tier
    else:
        data_tier = stage.data_tier
    dim = 'file_type %s' % project.file_type
    dim = dim + ' and data_tier %s' % data_tier
    dim = dim + ' and ub_project.name %s' % project.name
    dim = dim + ' and ub_project.stage %s' % stage.name
    dim = dim + ' and ub_project.version %s' % project.release_tag
    if stage.pubs_output:
        first_subrun = True
        for subrun in stage.output_subruns:
            if first_subrun:
                dim = dim + ' and run_number %d.%d' % (stage.output_run, subrun)
                first_subrun = False

                # Kluge to pick up mc files with wrong run number.

                if stage.output_run > 1 and stage.output_run < 100:
                    dim = dim + ',1.%d' % subrun
            else:
                dim = dim + ',%d.%d' % (stage.output_run, subrun)
    elif project.run_number != 0:
        dim = dim + ' and run_number %d' % project.run_number
    dim = dim + ' and availability: anylocation'
    return dim

def get_ups_products():
    return 'sbndcode'

