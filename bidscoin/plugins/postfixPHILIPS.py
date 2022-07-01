"""
This plugin runs after the dcm2niix2bids plugin when organizing data acquired through the Philips scanner. This involves the following fields being added to the json side cars:
 - EffectiveEchoSpacing
 - EstimatedTotalReadoutTime
 - PhaseEncodingDirection
 - SliceTiming
 - TaskName
The SeriesDescription field is also edited to include the acquisitions direction in AP/PA terms.

NOTE
 - SliceTiming is hard coded in this plugin. This is calculated using the sliceThickness and echoTiming header fields.
 - Intended for bidscoiner 3.7.0.

@author: Michael Mondaldi (mmonaldi@ufl.edu)
@created: 2021-11-30 by Michael Monaldi
@updated: 2022-03-13 by Tikahari Khanal (tikaharikhanal@ufl.edu)
"""

import logging
import dateutil.parser
import pandas as pd
import json
from typing import Union
from pathlib import Path
import shutil
import re
import nibabel as nib
try:
    from bidscoin import bidscoin, bids, physio
except ImportError:
    import bidscoin, bids, physio     # This should work if bidscoin was not pip-installed

LOGGER = logging.getLogger(__name__)



def test(options: dict) -> bool:
    """
    This plugin function tests the working of the plugin + its bidsmap options

    :param options: A dictionary with the plugin options, e.g. taken from the bidsmap['Options']['plugins']['README']
    :return:        True if the test was successful
    """

    LOGGER.warning(f'This is a demo-plugin test routine, validating its working with options: {options}')

    return True


def is_sourcefile(file: Path) -> str:
    """
    This plugin function assesses whether a sourcefile is of a supported dataformat

    :param file:    The sourcefile that is assessed
    :return:        The valid / supported dataformat of the sourcefile
    """

    if file.is_file():

        LOGGER.warning(f'This is a demo-plugin is_sourcefile routine, assessing whether "{file}" has a valid dataformat')
        return 'dataformat' if 'valid' == True else ''

    return ''



def get_attribute(dataformat: str, sourcefile: Path, attribute: str, options: dict) -> str:
    """
    This plugin function reads attributes from the supported sourcefile

    :param dataformat:  The dataformat of the sourcefile, e.g. DICOM of PAR
    :param sourcefile:  The sourcefile from which key-value data needs to be read
    :param attribute:   The attribute key for which the value needs to be retrieved
    :param options:     A dictionary with the plugin options, e.g. taken from the bidsmap['Options']
    :return:            The retrieved attribute value
    """
    if dataformat == 'DICOM':
        return bids.get_dicomfield(attribute, sourcefile)
    if dataformat == 'PAR':
        return bids.get_parfield(attribute, sourcefile)
    if dataformat in ('DICOM','PAR'):
        LOGGER.warning(f'This is a demo-plugin get_attribute routine, reading the {dataformat} "{attribute}" attribute value from "{sourcefile}"')

    return ''


def bidsmapper_plugin(session: Path, bidsmap_new: dict, bidsmap_old: dict, template: dict, store: dict) -> None:
    """
    All the logic to map the Philips PAR/XML fields onto bids labels go into this plugin function. The function is
    expecte to update / append new runs to the bidsmap_new data structure. The bidsmap options for this plugin can
    be found in:

    bidsmap_new/old['Options']['plugins']['README']

    See also the dcm2niix2bids plugin for reference implementation

    :param session:     The full-path name of the subject/session raw data source folder
    :param bidsmap_new: The study bidsmap that we are building
    :param bidsmap_old: Full BIDS heuristics data structure, with all options, BIDS labels and attributes, etc
    :param template:    The template bidsmap with the default heuristics
    :param store:       The paths of the source- and target-folder
    :return:
    """

    #LOGGER.info(f'This is a bidsmapper demo-plugin working on: {session}')


def bidscoiner_plugin(session: Path, bidsmap: dict, bidsfolder: Path) -> None:
    """
    The plugin to convert the runs in the source folder and save them in the bids folder. Each saved datafile should be
    accompanied with a json sidecar file. The bidsmap options for this plugin can be found in:

    bidsmap_new/old['Options']['plugins']['README']

    See also the dcm2niix2bids plugin for reference implementation

    :param session:     The full-path name of the subject/session raw data source folder
    :param bidsmap:     The full mapping heuristics from the bidsmap YAML-file
    :param bidsfolder:  The full-path name of the BIDS root-folder
    :return:            Nothing
    """
    plugin     = {'postfixPHILIPS': bidsmap['Options']['plugins']['postfixPHILIPS']}
    #LOGGER.warning(f'This is a bidscoiner plugin {plugin} working on: {session} -> {bidsfolder}')


# Get started and see what dataformat we have
    plugin2     = {'dcm2niix2bids': bidsmap['Options']['plugins']['dcm2niix2bids']}
    datasource = bids.get_datasource(session, plugin2)
    dataformat = datasource.dataformat
    if not dataformat:
        LOGGER.warning(f"No {__name__} sourcedata found in: {session}")
        return

    # Make a list of all the data sources / runs
    manufacturer = 'UNKNOWN'
    sources      = []
    if dataformat == 'DICOM':
        sources      = bidscoin.lsdirs(session)
        manufacturer = datasource.attributes('Manufacturer')
    elif dataformat == 'PAR':
        sources      = bids.get_parfiles(session)
        manufacturer = 'Philips Medical Systems'
    else:
        LOGGER.exception(f"Unsupported dataformat '{dataformat}'")

    # Get valid BIDS subject/session identifiers from the (first) DICOM- or PAR/XML source file
    subid, sesid = datasource.subid_sesid(bidsmap[dataformat]['subject'], bidsmap[dataformat]['session'])
    if not subid:
        LOGGER.warning(f'no SUBID found')
        return

    # Create the BIDS session-folder and a scans.tsv file
    bidsses = bidsfolder/subid/sesid


    if bidsses.is_dir():
        #LOGGER.warning(f'bids dir found --- bidsses is {bidsses}')
        bidsses_dwi = Path(bidsses / 'dwi')
        bidsses_fmap = Path(bidsses / 'fmap')
        bidsses_func = Path(bidsses / 'func')
        #LOGGER.warning(f'dwi dir found --- bidsses_dwi is {bidsses_dwi}')
        bidsFOLDERS = [bidsses_dwi, bidsses_fmap, bidsses_func]
        for mod_dir in bidsFOLDERS:
            for jsonfile in sorted((mod_dir).glob('sub-*.json')):
                #LOGGER.warning(f'JSONS found')
                # Load the existing meta-data
                with jsonfile.open('r') as json_fid:
                    jsondata = json.load(json_fid)
                # Modify Estimated* fields and set PhaseEncodingDirection based on SeriesDescription
                if 'EstimatedEffectiveEchoSpacing' in jsondata:
                    jsondata['EffectiveEchoSpacing']=jsondata['EstimatedEffectiveEchoSpacing']
                    jsondata['EstimatedEffectiveEchoSpacing']="GeneratedBy: dcm2niix + bidscoin"
                if 'EstimatedTotalReadoutTime' in jsondata:
                    jsondata['TotalReadoutTime']=jsondata['EstimatedTotalReadoutTime']
                    jsondata['EstimatedTotalReadoutTime']="GeneratedBy: dcm2niix + bidscoin"
                if 'PhaseEncodingDirection' not in jsondata:
                    if 'SeriesDescription' in jsondata:
                        phase_encoding_regex = re.search(r'.*(?i)(AP|PA).*', jsondata['SeriesDescription'])
                        try:
                            found_phase_encoding = phase_encoding_regex.groups()[0]
                            LOGGER.warning(f"found_phase_encoding: {found_phase_encoding} for:{jsonfile}")
                            #jsondata['PhaseEncodingDirection'] = found_phase_encoding
                        except Exception as e:
                            print(f"no regex match (e): {e}")
                            LOGGER.warning(f'Assumed PhaseEncodingDirection PA > j for {jsonfile}')
                            jsondata['PhaseEncodingDirection'] = 'j'
                
                if mod_dir == bidsses_func:
                    if 'SliceTiming' not in jsondata:
                        niifile = jsonfile.with_suffix('.nii')
                        niifile_gz = jsonfile.with_suffix('.nii.gz')
                        if niifile.is_file():
                            nii_img  = nib.load(niifile)
                        elif niifile_gz.is_file():
                            nii_img  = nib.load(niifile_gz)
                        else:
                            LOGGER.warning(f'NIFTI not found for SliceTiming: {niifile}')
                        nii_hdr = nii_img.header
                        try:
                            nrslices = nii_hdr.get_n_slices()
                        except Exception as niierror:
                            nrslices = nii_hdr.get_data_shape()[2]

                        LOGGER.warning(f'nrslices: {nrslices} , adding SliceTiming to {jsonfile}')
                        if nrslices == 69:
                            jsondata['SliceTiming'] = [0.000,0.065,0.130,0.196,0.261,0.326,0.391,0.457,0.522,0.587,0.652,0.717,0.783,0.848,0.913,0.978,1.043,1.109,1.174,1.239,1.304,1.370,1.435,0.000,0.065,0.130,0.196,0.261,0.326,0.391,0.457,0.522,0.587,0.652,0.717,0.783,0.848,0.913,0.978,1.043,1.109,1.174,1.239,1.304,1.370,1.435,0.000,0.065,0.130,0.196,0.261,0.326,0.391,0.457,0.522,0.587,0.652,0.717,0.783,0.848,0.913,0.978,1.043,1.109,1.174,1.239,1.304,1.370,1.435]
                        elif nrslices == 66:
                            jsondata['SliceTiming'] = [0.000,0.068,0.136,0.205,0.273,0.341,0.409,0.477,0.545,0.614,0.682,0.750,0.818,0.886,0.955,1.023,1.091,1.159,1.227,1.295,1.364,1.432,0.000,0.068,0.136,0.205,0.273,0.341,0.409,0.477,0.545,0.614,0.682,0.750,0.818,0.886,0.955,1.023,1.091,1.159,1.227,1.295,1.364,1.432,0.000,0.068,0.136,0.205,0.273,0.341,0.409,0.477,0.545,0.614,0.682,0.750,0.818,0.886,0.955,1.023,1.091,1.159,1.227,1.295,1.364,1.432]
                        else:
                            LOGGER.warning(f'{plugin} cannot determine SliceTiming nrslices:{nrslices} from NIFTI header for: {niifile}')
                            jsondata.pop('SliceTiming')
                    if 'TaskName' not in jsondata:
                        taskname = 'UNKNOWN'
                        LOGGER.warning(f'Empty TaskName for {jsonfile}, adding:{taskname}')
                        jsondata['TaskName'] = taskname


                # Save the collected meta-data to disk
                with jsonfile.open('w') as json_fid:
                    json.dump(jsondata, json_fid, indent=4)
    else:
        LOGGER.warning(f'NOTHING found -----145-----')
        LOGGER.warning(f'plugin is {plugin}')

# TODO: comment out
# Plot nifti by slice
#   import matplotlib.pyplot as plt
#
#   # Change the path to your path
#   path = 'path to img.nii.gz'
#   path)
#   nii_data = my_img.get_fdata()
#   nii_aff  = my_img.affine
#   nii_hdr  = my_img.header
#   print(nii_aff ,'\n',nii_hdr)
#   print(nii_data.shape)
#   if(len(nii_data.shape)==3):
#       for slice_Number in range(nii_data.shape[2]):
#           plt.imshow(nii_data[:,:,slice_Number ])
#           plt.show()
#   if(len(nii_data.shape)==4):
#       for frame in range(nii_data.shape[3]):
#           for slice_Number in range(nii_data.shape[2]):
#               plt.imshow(nii_data[:,:,slice_Number,frame])
#               plt.show()
