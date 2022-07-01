"""
This plugin replaces the dcm2nix2bids plugin in organizing data for the TrialNet Pancreas Volume Study (PIs Haller and Campbell-Thompson). Data was organized according to a custom standard which included the addition of a 'vp' key-value pair in the file name. 

NOTE
 - This plugin was used on data acquired at the following locations:
    - Australia
    - Barbara Davis
    - Colorado
    - Miami
    - Pitt
    - UF
 - A different run of this plugin was required over SKYRA and PRISMA data (separate bidsmaps were generated as well).
 - The data of interest is stored at USF. Contact Damon Lamb (dlamb@ufl.edu) for more information.
 - Intended for bidscoiner 3.6.3.

@author: Tikahari Khanal (tikaharikhanal@ufl.edu) Biomarker and Neuroimaging Core, Brain Rehabilitation Research Center (Director of Biomedical Informatics: Dr. Lamb) 
@created: 2022-11-30 by Tikahari Khanal
@updated: 2022-04-22 by Tikahari Khanal
"""

import logging
import dateutil.parser
import pandas as pd
import json
from pathlib import Path
from ruamel.yaml import YAML
from bidscoin.bids import *
from bidscoin.bids import entities as default_entities
try:
    from bidscoin import bidscoin, bids, physio
except ImportError:
    import bidscoin, bids, physio     # This should work if bidscoin was not pip-installed

LOGGER = logging.getLogger(__name__)
yaml = YAML()


# Define BIDScoin datatypes
bidscoindatatypes = ('fmap', 'anat', 'func', 'perf', 'dwi', 'pet', 'meg', 'eeg', 'ieeg', 'beh')           # NB: get_matching_run() uses this order to search for a match. TODO: sync with the modalities.yaml schema
ignoredatatype    = 'exclude'
unknowndatatype   = 'extra_data'

# Define the default paths
schema_folder     = Path(__file__).parent/'schema'
heuristics_folder = Path(__file__).parent/'heuristics'
bidsmap_template  = heuristics_folder/'bidsmap_template.yaml'

# Read the BIDS schema datatypes and entities
bidsdatatypes = {
    "dwi": [
        {
            "suffixes": [
                "dwi"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json",
                ".bvec",
                ".bval"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "acquisition": "optional",
                "direction": "optional",
                "run": "optional",
                "part": "optional"
            }
        },
        {
            "suffixes": [
                "sbref"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "acquisition": "optional",
                "direction": "optional",
                "run": "optional",
                "part": "optional"
            }
        }
    ],
    "meg": [
        {
            "suffixes": [
                "meg"
            ],
            "extensions": [
                "/",
                ".ds/",
                ".json",
                ".fif",
                ".sqd",
                ".con",
                ".raw",
                ".ave",
                ".mrk",
                ".kdf",
                ".mhd"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "task": "required",
                "acquisition": "optional",
                "run": "optional",
                "processing": "optional",
                "split": "optional"
            }
        },
        {
            "suffixes": [
                "headshape"
            ],
            "extensions": [
                "*"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "acquisition": "optional"
            }
        },
        {
            "suffixes": [
                "markers"
            ],
            "extensions": [
                ".sqd",
                ".mrk"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "task": "optional",
                "acquisition": "optional",
                "space": "optional"
            }
        },
        {
            "suffixes": [
                "coordsystem"
            ],
            "extensions": [
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "acquisition": "optional"
            }
        },
        {
            "suffixes": [
                "channels"
            ],
            "extensions": [
                ".json",
                ".tsv"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "task": "required",
                "acquisition": "optional",
                "run": "optional",
                "processing": "optional"
            }
        },
        {
            "suffixes": [
                "events"
            ],
            "extensions": [
                ".json",
                ".tsv"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "task": "required",
                "acquisition": "optional",
                "run": "optional"
            }
        },
        {
            "suffixes": [
                "photo"
            ],
            "extensions": [
                ".jpg"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "acquisition": "optional"
            }
        }
    ],
    "beh": [
        {
            "suffixes": [
                "stim",
                "physio"
            ],
            "extensions": [
                ".tsv.gz",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "task": "required",
                "acquisition": "optional",
                "run": "optional",
                "recording": "optional"
            }
        },
        {
            "suffixes": [
                "events",
                "beh"
            ],
            "extensions": [
                ".tsv",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "task": "required",
                "acquisition": "optional",
                "run": "optional"
            }
        }
    ],
    "fmap": [
        {
            "suffixes": [
                "phasediff",
                "phase1",
                "phase2",
                "magnitude1",
                "magnitude2",
                "magnitude",
                "fieldmap"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "acquisition": "optional",
                "run": "optional"
            }
        },
        {
            "suffixes": [
                "epi",
                "m0scan"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "acquisition": "optional",
                "ceagent": "optional",
                "direction": "required",
                "run": "optional"
            }
        },
        {
            "suffixes": [
                "TB1DAM"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "run": "optional",
                "acquisition": "optional",
                "ceagent": "optional",
                "reconstruction": "optional",
                "flip": "required",
                "inversion": "optional",
                "part": "optional"
            }
        },
        {
            "suffixes": [
                "TB1EPI"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "run": "optional",
                "acquisition": "optional",
                "ceagent": "optional",
                "reconstruction": "optional",
                "echo": "required",
                "flip": "required",
                "inversion": "optional",
                "part": "optional"
            }
        },
        {
            "suffixes": [
                "TB1AFI",
                "TB1TFL",
                "TB1RFM",
                "RB1COR"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "run": "optional",
                "acquisition": "optional",
                "ceagent": "optional",
                "reconstruction": "optional",
                "echo": "optional",
                "flip": "optional",
                "inversion": "optional",
                "part": "optional"
            }
        },
        {
            "suffixes": [
                "TB1SRGE"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "run": "optional",
                "acquisition": "optional",
                "ceagent": "optional",
                "reconstruction": "optional",
                "echo": "optional",
                "flip": "required",
                "inversion": "required",
                "part": "optional"
            }
        },
        {
            "suffixes": [
                "TB1map",
                "RB1map"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "run": "optional",
                "acquisition": "optional",
                "ceagent": "optional",
                "reconstruction": "optional"
            }
        }
    ],
    "pet": [
        {
            "suffixes": [
                "pet"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "task": "optional",
                "tracer": "optional",
                "reconstruction": "optional",
                "run": "optional"
            }
        },
        {
            "suffixes": [
                "blood"
            ],
            "extensions": [
                ".tsv",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "task": "optional",
                "tracer": "optional",
                "reconstruction": "optional",
                "run": "optional",
                "recording": "required"
            }
        },
        {
            "suffixes": [
                "events"
            ],
            "extensions": [
                ".tsv",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "task": "required",
                "tracer": "optional",
                "reconstruction": "optional",
                "run": "optional"
            }
        }
    ],
    "ieeg": [
        {
            "suffixes": [
                "ieeg"
            ],
            "extensions": [
                ".mefd/",
                ".json",
                ".edf",
                ".vhdr",
                ".eeg",
                ".vmrk",
                ".set",
                ".fdt",
                ".nwb"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "task": "required",
                "acquisition": "optional",
                "run": "optional"
            }
        },
        {
            "suffixes": [
                "channels"
            ],
            "extensions": [
                ".json",
                ".tsv"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "task": "required",
                "acquisition": "optional",
                "run": "optional"
            }
        },
        {
            "suffixes": [
                "coordsystem"
            ],
            "extensions": [
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "acquisition": "optional",
                "space": "optional"
            }
        },
        {
            "suffixes": [
                "electrodes"
            ],
            "extensions": [
                ".json",
                ".tsv"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "acquisition": "optional",
                "space": "optional"
            }
        },
        {
            "suffixes": [
                "events"
            ],
            "extensions": [
                ".json",
                ".tsv"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "task": "required",
                "acquisition": "optional",
                "run": "optional"
            }
        },
        {
            "suffixes": [
                "photo"
            ],
            "extensions": [
                ".jpg"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "acquisition": "optional"
            }
        }
    ],
    "perf": [
        {
            "suffixes": [
                "asl",
                "m0scan"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "acquisition": "optional",
                "reconstruction": "optional",
                "direction": "optional",
                "run": "optional"
            }
        },
        {
            "suffixes": [
                "aslcontext"
            ],
            "extensions": [
                ".tsv",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "acquisition": "optional",
                "reconstruction": "optional",
                "direction": "optional",
                "run": "optional"
            }
        },
        {
            "suffixes": [
                "asllabeling"
            ],
            "extensions": [
                ".jpg"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "acquisition": "optional",
                "reconstruction": "optional",
                "run": "optional"
            }
        }
    ],
    "anat": [
        {
            "suffixes": [
                "T1w",
                "T2w",
                "PDw",
                "T2starw",
                "FLAIR",
                "inplaneT1",
                "inplaneT2",
                "PDT2",
                "angio",
                "T2star",
                "FLASH",
                "PD"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "vp": "optional",
                "run": "optional",
                "acquisition": "optional",
                "ceagent": "optional",
                "reconstruction": "optional",
                "part": "optional"
            }
        },
        {
            "suffixes": [
                "T1map",
                "T2map",
                "T2starmap",
                "R1map",
                "R2map",
                "R2starmap",
                "PDmap",
                "MTRmap",
                "MTsat",
                "UNIT1",
                "T1rho",
                "MWFmap",
                "MTVmap",
                "PDT2map",
                "Chimap",
                "S0map",
                "M0map"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "run": "optional",
                "acquisition": "optional",
                "ceagent": "optional",
                "reconstruction": "optional"
            }
        },
        {
            "suffixes": [
                "defacemask"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "run": "optional",
                "acquisition": "optional",
                "ceagent": "optional",
                "reconstruction": "optional",
                "modality": "optional"
            }
        },
        {
            "suffixes": [
                "MESE",
                "MEGRE"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "run": "optional",
                "acquisition": "optional",
                "ceagent": "optional",
                "reconstruction": "optional",
                "echo": "required",
                "part": "optional"
            }
        },
        {
            "suffixes": [
                "VFA"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "run": "optional",
                "acquisition": "optional",
                "ceagent": "optional",
                "reconstruction": "optional",
                "flip": "required",
                "part": "optional"
            }
        },
        {
            "suffixes": [
                "IRT1"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "run": "optional",
                "acquisition": "optional",
                "ceagent": "optional",
                "reconstruction": "optional",
                "inversion": "required",
                "part": "optional"
            }
        },
        {
            "suffixes": [
                "MP2RAGE"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "run": "optional",
                "acquisition": "optional",
                "ceagent": "optional",
                "reconstruction": "optional",
                "echo": "optional",
                "flip": "optional",
                "inversion": "required",
                "part": "optional"
            }
        },
        {
            "suffixes": [
                "MPM",
                "MTS"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "run": "optional",
                "acquisition": "optional",
                "ceagent": "optional",
                "reconstruction": "optional",
                "echo": "optional",
                "flip": "required",
                "mtransfer": "required",
                "part": "optional"
            }
        },
        {
            "suffixes": [
                "MTR"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "run": "optional",
                "acquisition": "optional",
                "ceagent": "optional",
                "reconstruction": "optional",
                "mtransfer": "required",
                "part": "optional"
            }
        }
    ],
    "eeg": [
        {
            "suffixes": [
                "eeg"
            ],
            "extensions": [
                ".json",
                ".edf",
                ".vhdr",
                ".vmrk",
                ".eeg",
                ".set",
                ".fdt",
                ".bdf"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "task": "required",
                "acquisition": "optional",
                "run": "optional"
            }
        },
        {
            "suffixes": [
                "channels"
            ],
            "extensions": [
                ".json",
                ".tsv"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "task": "required",
                "acquisition": "optional",
                "run": "optional"
            }
        },
        {
            "suffixes": [
                "coordsystem"
            ],
            "extensions": [
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "acquisition": "optional",
                "space": "optional"
            }
        },
        {
            "suffixes": [
                "electrodes"
            ],
            "extensions": [
                ".json",
                ".tsv"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "acquisition": "optional",
                "space": "optional"
            }
        },
        {
            "suffixes": [
                "events"
            ],
            "extensions": [
                ".json",
                ".tsv"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "task": "required",
                "acquisition": "optional",
                "run": "optional"
            }
        },
        {
            "suffixes": [
                "photo"
            ],
            "extensions": [
                ".jpg"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "acquisition": "optional"
            }
        }
    ],
    "func": [
        {
            "suffixes": [
                "bold",
                "cbv",
                "sbref"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "task": "required",
                "acquisition": "optional",
                "ceagent": "optional",
                "reconstruction": "optional",
                "direction": "optional",
                "run": "optional",
                "echo": "optional",
                "part": "optional"
            }
        },
        {
            "suffixes": [
                "phase"
            ],
            "extensions": [
                ".nii.gz",
                ".nii",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "task": "required",
                "acquisition": "optional",
                "ceagent": "optional",
                "reconstruction": "optional",
                "direction": "optional",
                "run": "optional",
                "echo": "optional"
            }
        },
        {
            "suffixes": [
                "events"
            ],
            "extensions": [
                ".tsv",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "task": "required",
                "acquisition": "optional",
                "ceagent": "optional",
                "reconstruction": "optional",
                "direction": "optional",
                "run": "optional"
            }
        },
        {
            "suffixes": [
                "physio",
                "stim"
            ],
            "extensions": [
                ".tsv.gz",
                ".json"
            ],
            "entities": {
                "subject": "required",
                "session": "optional",
                "task": "required",
                "acquisition": "optional",
                "ceagent": "optional",
                "reconstruction": "optional",
                "direction": "optional",
                "run": "optional",
                "recording": "optional"
            }
        }
    ]
}

entities = {
    "subject": {
        "name": "Subject",
        "entity": "sub",
        "description": "A person or animal participating in the study.\n",
        "format": "label"
    },
    "session": {
        "name": "Session",
        "entity": "ses",
        "description": "A logical grouping of neuroimaging and behavioral data consistent across\nsubjects.\nSession can (but doesn't have to) be synonymous to a visit in a\nlongitudinal study.\nIn general, subjects will stay in the scanner during one session.\nHowever, for example, if a subject has to leave the scanner room and then\nbe re-positioned on the scanner bed, the set of MRI acquisitions will still\nbe considered as a session and match sessions acquired in other subjects.\nSimilarly, in situations where different data types are obtained over\nseveral visits (for example fMRI on one day followed by DWI the day after)\nthose can be grouped in one session.\nDefining multiple sessions is appropriate when several identical or similar\ndata acquisitions are planned and performed on all -or most- subjects,\noften in the case of some intervention between sessions\n(for example, training).\n",
        "format": "label"
    },
    "task": {
        "name": "Task",
        "entity": "task",
        "format": "label",
        "description": "Each task has a unique label that MUST only consist of letters and/or\nnumbers (other characters, including spaces and underscores, are not\nallowed).\nThose labels MUST be consistent across subjects and sessions.\n"
    },
    "acquisition": {
        "name": "Acquisition",
        "entity": "acq",
        "description": "The `acq-<label>` key/value pair corresponds to a custom label the\nuser MAY use to distinguish a different set of parameters used for\nacquiring the same modality.\nFor example this should be used when a study includes two T1w images - one\nfull brain low resolution and and one restricted field of view but high\nresolution.\nIn such case two files could have the following names:\n`sub-01_acq-highres_T1w.nii.gz` and `sub-01_acq-lowres_T1w.nii.gz`, however\nthe user is free to choose any other label than highres and lowres as long\nas they are consistent across subjects and sessions.\nIn case different sequences are used to record the same modality\n(for example, RARE and FLASH for T1w)\nthis field can also be used to make that distinction.\nAt what level of detail to make the distinction (for example,\njust between RARE and FLASH, or between RARE, FLASH, and FLASHsubsampled)\nremains at the discretion of the researcher.\n",
        "format": "label"
    },
    "ceagent": {
        "name": "Contrast Enhancing Agent",
        "entity": "ce",
        "description": "The `ce-<label>` key/value can be used to distinguish\nsequences using different contrast enhanced images.\nThe label is the name of the contrast agent.\nThe key `ContrastBolusIngredient` MAY also be added in the JSON file,\nwith the same label.\n",
        "format": "label"
    },
    "tracer": {
        "name": "Tracer",
        "entity": "trc",
        "description": "The `trc-<label>` key/value can be used to distinguish\nsequences using different tracers.\nThe key `TracerName` MUST also be included in the associated JSON file,\nalthough the label may be different.\n",
        "format": "label"
    },
    "reconstruction": {
        "name": "Reconstruction",
        "entity": "rec",
        "description": "The `rec-<label>` key/value can be used to distinguish\ndifferent reconstruction algorithms (for example ones using motion\ncorrection).\n",
        "format": "label"
    },
    "direction": {
        "name": "Phase-Encoding Direction",
        "entity": "dir",
        "description": "The `dir-<label>` key/value can be set to an arbitrary alphanumeric label\n(for example, `dir-LR` or `dir-AP`) to distinguish different phase-encoding\ndirections.\n",
        "format": "label"
    },
    "vp": {
        "name": "Recording",
        "entity": "vp",
        "description": "ax/cor/sag\n",
        "format": "label"
    },
    "run": {
        "name": "Run",
        "entity": "run",
        "description": "If several scans with the same acquisition parameters are acquired in the same session,\nthey MUST be indexed with the [`run-<index>`](../99-appendices/09-entities.md#run) entity:\n`_run-1`, `_run-2`, `_run-3`, and so on (only nonnegative integers are allowed as\nrun labels).\n\nIf different entities apply,\nsuch as a different session indicated by [`ses-<label>`](../99-appendices/09-entities.md#ses),\nor different acquisition parameters indicated by\n[`acq-<label>`](../99-appendices/09-entities.md#acq),\nthen `run` is not needed to distinguish the scans and MAY be omitted.\n",
        "format": "index"
    },
    "modality": {
        "name": "Corresponding Modality",
        "entity": "mod",
        "description": "The `mod-<label>` key/value pair corresponds to modality label for defacing\nmasks, for example, T1w, inplaneT1, referenced by a defacemask image.\nFor example, `sub-01_mod-T1w_defacemask.nii.gz`.\n",
        "format": "label"
    },
    "echo": {
        "name": "Echo",
        "entity": "echo",
        "description": "If files belonging to an entity-linked file collection are acquired at different\necho times, the `_echo-<index>` key/value pair MUST be used to distinguish\nindividual files.\nThis entity represents the `EchoTime` metadata field. Please note that the `<index>`\ndenotes the number/index (in the form of a nonnegative integer), not the\n`EchoTime` value which needs to be stored in the field `EchoTime` of the separate\nJSON file.\n",
        "format": "index"
    },
    "flip": {
        "name": "Flip Angle",
        "entity": "flip",
        "description": "If files belonging to an entity-linked file collection are acquired at different\nflip angles, the `_flip-<index>` key/value pair MUST be used to distinguish\nindividual files.\nThis entity represents the `FlipAngle` metadata field. Please note that the `<index>`\ndenotes the number/index (in the form of a nonnegative integer), not the `FlipAngle`\nvalue which needs to be stored in the field `FlipAngle` of the separate JSON file.\n",
        "format": "index"
    },
    "inversion": {
        "name": "Inversion Time",
        "entity": "inv",
        "description": "If files belonging to an entity-linked file collection are acquired at different\ninversion times, the `_inv-<index>` key/value pair MUST be used to distinguish\nindividual files.\nThis entity represents the `InversionTime` metadata field. Please note that the `<index>`\ndenotes the number/index (in the form of a nonnegative integer), not the `InversionTime`\nvalue which needs to be stored in the field `InversionTime` of the separate JSON file.\n",
        "format": "index"
    },
    "mtransfer": {
        "name": "Magnetization Transfer",
        "entity": "mt",
        "description": "If files belonging to an entity-linked file collection are acquired at different\nmagnetization transfer (MT) states, the `_mt-<label>` key/value pair MUST be used to\ndistinguish individual files.\nThis entity represents the `MTState` metadata field. Allowed label values for this\nentity are `on` and `off`, for images acquired in presence and absence of an MT pulse,\nrespectively.\n",
        "format": "label"
    },
    "part": {
        "name": "Part",
        "entity": "part",
        "description": "This entity is used to indicate which component of the complex\nrepresentation of the MRI signal is represented in voxel data.\nThe `part-<label>` key/value pair is associated with the DICOM tag\n`0008,9208`.\nAllowed label values for this entity are `phase`, `mag`, `real` and `imag`,\nwhich are typically used in `part-mag`/`part-phase` or\n`part-real`/`part-imag` pairs of files.\n\nPhase images MAY be in radians or in arbitrary units.\nThe sidecar JSON file MUST include the units of the `phase` image.\nThe possible options are `rad` or `arbitrary`.\n\nWhen there is only a magnitude image of a given type, the `part` key MAY be\nomitted.\n",
        "format": "label"
    },
    "recording": {
        "name": "Recording",
        "entity": "recording",
        "description": "More than one continuous recording file can be included (with different\nsampling frequencies).\nIn such case use different labels.\nFor example: `_recording-contrast`, `_recording-saturation`.\n",
        "format": "label"
    },
    "processing": {
        "name": "Processed (on device)",
        "entity": "proc",
        "description": "The proc label is analogous to rec for MR and denotes a variant of a file\nthat was a result of particular processing performed on the device.\nThis is useful for files produced in particular by Elekta\u2019s MaxFilter\n(for example, sss, tsss, trans, quat or mc),\nwhich some installations impose to be run on raw data because of active\nshielding software corrections before the MEG data can actually be\nexploited.\n",
        "format": "label"
    },
    "space": {
        "name": "Space",
        "entity": "space",
        "description": "The space entity can be used to indicate\nthe way in which electrode positions are interpreted\n(for EEG/MEG/iEEG data) or\nthe spatial reference to which a file has been aligned (for MRI data).\nThe space `<label>` MUST be taken from one of the modality specific lists in\n[Appendix VIII](../99-appendices/08-coordinate-systems.md).\nFor example for iEEG data, the restricted keywords listed under\n[iEEG Specific Coordinate Systems](../99-appendices/08-coordinate-systems.md#ieeg-specific-coordinate-systems)\nare acceptable for `<label>`.\n\nFor EEG/MEG/iEEG data, this entity can be applied to raw data, but\nfor other data types, it is restricted to derivative data.\n",
        "format": "label"
    },
    "split": {
        "name": "Split",
        "entity": "split",
        "description": "In the case of long data recordings that exceed a file size of 2Gb, the\n.fif files are conventionally split into multiple parts.\nEach of these files has an internal pointer to the next file.\nThis is important when renaming these split recordings to the BIDS\nconvention.\n\nInstead of a simple renaming, files should be read in and saved under their\nnew names with dedicated tools like [MNE-Python](https://mne.tools/),\nwhich will ensure that not only the file names, but also the internal file\npointers will be updated.\nIt is RECOMMENDED that .fif files with multiple parts use the\n`split-<index>` entity to indicate each part.\nIf there are multiple parts of a recording and the optional `scans.tsv` is provided,\nremember to list all files separately in `scans.tsv` and that the entries for the\n`acq_time` column in `scans.tsv` MUST all be identical, as described in\n[Scans file](../03-modality-agnostic-files.md#scans-file).\n",
        "format": "index"
    },
    "resolution": {
        "name": "Resolution",
        "entity": "res",
        "description": "Resolution of regularly sampled N-dimensional data.\nMUST have a corresponding `Resolution` metadata field to provide\ninterpretation.\n\nThis entity is only applicable to derivative data.\n",
        "format": "label"
    },
    "density": {
        "name": "Density",
        "entity": "den",
        "description": "Density of non-parametric surfaces.\nMUST have a corresponding `Density` metadata field to provide\ninterpretation.\n\nThis entity is only applicable to derivative data.\n",
        "format": "label"
    },
    "label": {
        "name": "Label",
        "entity": "label",
        "description": "Tissue-type label, following a prescribed vocabulary.\nApplies to binary masks and probabilistic/partial volume segmentations\nthat describe a single tissue type.\n\nThis entity is only applicable to derivative data.\n",
        "format": "label"
    },
    "description": {
        "name": "Description",
        "entity": "desc",
        "description": "When necessary to distinguish two files that do not otherwise have a\ndistinguishing entity, the `_desc-<label>` keyword-value SHOULD be used.\n\nThis entity is only applicable to derivative data.\n",
        "format": "label"
    }
}

def get_bidsname_custom(subid: str, sesid: str, run: dict, runtime: bool=False) -> str:
    """
    Composes a filename as it should be according to the BIDS standard using the BIDS keys in run. The bids values are
    dynamically updated and cleaned, and invalid bids keys are ignored

    :param subid:       The subject identifier, i.e. name of the subject folder (e.g. 'sub-001' or just '001')
    :param sesid:       The optional session identifier, i.e. name of the session folder (e.g. 'ses-01' or just '01'). Can be left empty
    :param run:         The run mapping with the BIDS key-value pairs
    :param runtime:     Replaces <<dynamic>> bidsvalues if True
    :return:            The composed BIDS file-name (without file-extension)
    """

    # Try to update the sub/ses-ids
    if not subid.startswith('sub-'):
        subid = f"sub-{cleanup_value(subid)}"
    if sesid and not sesid.startswith('ses-'):
        sesid = f"ses-{cleanup_value(sesid)}"

    # Compose a bidsname from valid BIDS entities only
    bidsname = f"{subid}{add_prefix('_', sesid)}"                               # Start with the subject/session identifier
    # Get default entity keys
    defaultkeys = [default_entities[entity]['entity'] for entity in default_entities]
    for entitykey in [entities[entity]['entity'] for entity in entities]:
        if entitykey not in defaultkeys:
            LOGGER.warning(f"{entitykey} not in bids standard, but added as part of new file name, please confirm this is inentional")
        bidsvalue = run['bids'].get(entitykey)                                  # Get the entity data from the run
        if isinstance(bidsvalue, list):
            bidsvalue = bidsvalue[bidsvalue[-1]]                                # Get the selected item
        else:
            bidsvalue = get_dynamicvalue(bidsvalue, Path(run['provenance']), runtime)
        if bidsvalue:
            bidsname = f"{bidsname}_{entitykey}-{cleanup_value(bidsvalue)}"     # Append the key-value data to the bidsname
    bidsname = f"{bidsname}{add_prefix('_', run['bids']['suffix'])}"            # And end with the suffix

    return bidsname

    
def load_bidsmap_custom(yamlfile: Path, folder: Path=Path(), report: Union[bool,None]=True) -> Tuple[dict, Path]:
    """
    Read the mapping heuristics from the bidsmap yaml-file. If yamlfile is not fullpath, then 'folder' is first searched before
    the default 'heuristics'. If yamfile is empty, then first 'bidsmap.yaml' is searched for, them 'bidsmap_template.yaml'. So fullpath
    has precendence over folder and bidsmap.yaml has precedence over bidsmap_template.yaml

    :param yamlfile:    The full pathname or basename of the bidsmap yaml-file. If None, the default bidsmap_template.yaml file in the heuristics folder is used
    :param folder:      Only used when yamlfile=basename or None: yamlfile is then first searched for in folder and then falls back to the ./heuristics folder (useful for centrally managed template yaml-files)
    :param report:      Report log.info when reading a file
    :return:            Tuple with (1) ruamel.yaml dict structure, with all options, BIDS mapping heuristics, labels and attributes, etc and (2) the fullpath yaml-file
    """
    LOGGER.info(f"Reloading bidsmap file")
    # Input checking
    if not folder.name:
        folder = heuristics_folder
    if not yamlfile.name:
        yamlfile = folder/'bidsmap.yaml'
        if not yamlfile.is_file():
            yamlfile = bidsmap_template

    # Add a standard file-extension if needed
    if not yamlfile.suffix:
        yamlfile = yamlfile.with_suffix('.yaml')

    # Get the full path to the bidsmap yaml-file
    if len(yamlfile.parents) == 1:
        if (folder/yamlfile).is_file():
            yamlfile = folder/yamlfile
        else:
            yamlfile = heuristics_folder/yamlfile

    if not yamlfile.is_file():
        if report:
            LOGGER.info(f"No existing bidsmap file found: {yamlfile}")
        return dict(), yamlfile
    elif report:
        LOGGER.info(f"Reading: {yamlfile}")

    # Read the heuristics from the bidsmap file
    with yamlfile.open('r') as stream:
        bidsmap = yaml.load(stream)

    # Issue a warning if the version in the bidsmap YAML-file is not the same as the bidscoin version
    if 'bidscoin' in bidsmap['Options'] and 'version' in bidsmap['Options']['bidscoin']:
        bidsmapversion = bidsmap['Options']['bidscoin']['version']
    elif 'version' in bidsmap['Options']:                       # Handle legacy bidsmaps
        bidsmapversion = bidsmap['Options']['version']
    else:
        bidsmapversion = 'Unknown'
    if bidsmapversion != bidscoin.version() and report:
        LOGGER.warning(f'BIDScoiner version conflict: {yamlfile} was created using version {bidsmapversion}, but this is version {bidscoin.version()}')

    # Add missing provenance info, run dictionaries and bids entities
    run_ = get_run_()
    for dataformat in bidsmap:
        if dataformat in ('Options','PlugIns'): continue        # Handle legacy bidsmaps (-> 'PlugIns')
        if not bidsmap[dataformat]:             continue
        for datatype in bidsmap[dataformat]:
            if not isinstance(bidsmap[dataformat][datatype], list): continue
            for index, run in enumerate(bidsmap[dataformat][datatype]):

                # Add missing provenance info
                if not run.get('provenance'):
                    run['provenance'] = f"sub-unknown/ses-unknown/{dataformat}_{datatype}_id{index+1:03}"

                # Add missing run dictionaries (e.g. "meta" or "filesystem")
                for key, val in run_.items():
                    if key not in run or not run[key]:
                        run[key] = val

                # Add missing bids entities
                for typegroup in bidsdatatypes.get(datatype,[]):
                    if run['bids']['suffix'] in typegroup['suffixes']:      # run_found = True
                        for entityname in typegroup['entities']:
                            entitykey = entities[entityname]['entity']
                            if entitykey not in run['bids'] and entitykey not in ('sub','ses'):
                                LOGGER.debug(f"Adding missing {dataformat}/{datatype} entity key: {entitykey}")
                                run['bids'][entitykey] = ''

    # Make sure we get a proper dictionary with plugins
    if not bidsmap['Options'].get('plugins'):
        bidsmap['Options']['plugins'] = {}
    for plugin, options in bidsmap['Options']['plugins'].items():
        if not bidsmap['Options']['plugins'].get(plugin):
            bidsmap['Options']['plugins'][plugin] = {}

    # Validate the bidsmap entries
    check_bidsmap(bidsmap, report)

    return bidsmap, yamlfile


def test(options) -> bool:
    """
    Performs shell tests dcm2niix

    :return:        True if the tool generated the expected result, False if there was a tool error, None if not tested
    """

    LOGGER.info('Testing the dcm2niix installation:')

    if 'path' not in options:
        LOGGER.error(f"The expected 'path' key is not defined in the custom_pancreas options")
    if 'args' not in options:
        LOGGER.error(f"The expected 'args' key is not defined in the custom_pancreas options")

    command = f"{options.get('path')}dcm2niix -u"

    return bidscoin.run_command(command)


def bidscoiner_plugin(session: Path, bidsmap: dict, bidsfolder: Path, personals: dict, subprefix: str, sesprefix: str) -> None:
    """
    The bidscoiner plugin to convert the session DICOM and PAR/REC source-files into BIDS-valid nifti-files in the
    corresponding bidsfolder and extract personals (e.g. Age, Sex) from the source header

    :param session:     The full-path name of the subject/session source file/folder
    :param bidsmap:     The full mapping heuristics from the bidsmap YAML-file
    :param bidsfolder:  The full-path name of the BIDS root-folder
    :param personals:   The dictionary with the personal information
    :param subprefix:   The prefix common for all source subject-folders
    :param sesprefix:   The prefix common for all source session-folders
    :return:            Nothing
    """

    # reload bidsmap
    # bidsmap, _ = load_bidsmap_custom(Path('bidsmap.yaml'), bidsfolder/'code'/'bidscoin')
    print(f"Running custom pancreas plugin")
    # See what dataformat we have
    dataformat = bids.get_dataformat(session)

    # Get valid BIDS subject/session identifiers from the (first) DICOM- or PAR/XML source file
    sourcefile   = Path()
    manufacturer = 'UNKNOWN'
    if dataformat == 'DICOM':
        sources = bidscoin.lsdirs(session)
        for source in sources:
            sourcefile = bids.get_dicomfile(source)
            if sourcefile.name:
                manufacturer = bids.get_dicomfield('Manufacturer', sourcefile)
                break
        if not sourcefile.name:
            LOGGER.info(f"No data found in: {session}")
            return

    elif dataformat == 'PAR':
        sources = bids.get_parfiles(session)
        if sources:
            sourcefile   = sources[0]
            manufacturer = 'Philips Medical Systems'
        else:
            LOGGER.info(f"No data found in: {session}")
            return

    else:
        LOGGER.info(f"Session {session} cannot be processed by {__name__}")
        return

    subid, sesid = bids.get_subid_sesid(sourcefile,
                                        bidsmap[dataformat]['subject'],
                                        bidsmap[dataformat]['session'],
                                        subprefix, sesprefix)

    if subid == subprefix:
        LOGGER.error(f"No valid subject identifier found for: {session}")
        return

    # Create the BIDS session-folder and a scans.tsv file
    bidsses = bidsfolder/subid/sesid
    if bidsses.is_dir():
        LOGGER.warning(f"Existing BIDS output-directory found, which may result in duplicate data (with increased run-index). Make sure {bidsses} was cleaned-up from old data before (re)running the bidscoiner")
    bidsses.mkdir(parents=True, exist_ok=True)
    scans_tsv = bidsses/f"{subid}{bids.add_prefix('_',sesid)}_scans.tsv"
    if scans_tsv.is_file():
        scans_table = pd.read_csv(scans_tsv, sep='\t', index_col='filename')
    else:
        scans_table = pd.DataFrame(columns=['acq_time'], dtype='str')
        scans_table.index.name = 'filename'

    # Process all the source files or run subfolders
    for source in sources:

        # Get a source-file
        if dataformat == 'DICOM':
            sourcefile = bids.get_dicomfile(source)
        elif dataformat == 'PAR':
            sourcefile = source
        if not sourcefile.name:
            continue

        # Get a matching run from the bidsmap
        run, datatype, index = bids.get_matching_run(sourcefile, bidsmap, dataformat)
        # Check if we should ignore this run
        if datatype == bids.ignoredatatype:
            LOGGER.info(f"Leaving out: {source}")
            continue

        # Check if we already know this run
        if index is None:
            # LOGGER.error(f"Skipping unknown '{datatype}' run: {sourcefile}\n-> Re-run the bidsmapper and delete {bidsses} to solve this warning")
            continue

        LOGGER.info(f"Processing: {source}")

        # Create the BIDS session/datatype output folder
        if run['bids']['suffix'] in bids.get_derivatives(datatype):
            outfolder = bidsfolder/'derivatives'/manufacturer.replace(' ','')/subid/sesid/datatype
        else:
            outfolder = bidsses/datatype
        outfolder.mkdir(parents=True, exist_ok=True)

        # Modify echo time before run
        bidsname  = get_bidsname_custom(subid, sesid, run, runtime=True)
        echoindex = run['bids'].get('echo', '')
        if echoindex.startswith('<<') and echoindex.endswith('>>'):
            while list(outfolder.glob(bidsname + '.*')):
                currentechoindex = bids.get_bidsvalue(bidsname, 'echo')
                if runindex:
                    bidsname = get_bidsvalue(bidsname, 'echo', str(int(currentechoindex) + 1))

        # Compose the BIDS filename using the matched run
        runindex  = run['bids'].get('run', '')
        if runindex.startswith('<<') and runindex.endswith('>>'):
            bidsname = bids.increment_runindex(outfolder, bidsname)
        jsonfiles = [(outfolder/bidsname).with_suffix('.json')]      # List -> Collect the associated json-files (for updating them later) -- possibly > 1

        # Check if file already exists (-> e.g. when a static runindex is used)
        if (outfolder/bidsname).with_suffix('.json').is_file():
            LOGGER.warning(f"{outfolder/bidsname}.* already exists and will be deleted -- check your results carefully!")
            for ext in ('.nii.gz', '.nii', '.json', '.bval', '.bvec', '.tsv.gz'):
                (outfolder/bidsname).with_suffix(ext).unlink(missing_ok=True)

        # Convert physiological log files (dcm2niix can't handle these)
        if run['bids']['suffix'] == 'physio':
            if bids.get_dicomfile(source, 2).name:
                LOGGER.warning(f"Found > 1 DICOM file in {source}, using: {sourcefile}")
            physiodata = physio.readphysio(sourcefile)
            physio.physio2tsv(physiodata, outfolder/bidsname)

        # Convert the source-files in the run folder to nifti's in the BIDS-folder
        else:
            command = '{path}dcm2niix {args} -f "{filename}" -o "{outfolder}" "{source}"'.format(
                path      = bidsmap['Options']['plugins']['custom_pancreas'].get('path',''),
                args      = bidsmap['Options']['plugins']['custom_pancreas'].get('args',''),
                filename  = bidsname,
                outfolder = outfolder,
                source    = source)
            if not bidscoin.run_command(command):
                continue

            # Replace uncropped output image with the cropped one
            if '-x y' in bidsmap['Options']['plugins']['custom_pancreas']['args']:
                for dcm2niixfile in sorted(outfolder.glob(bidsname + '*_Crop_*')):                              # e.g. *_Crop_1.nii.gz
                    ext         = ''.join(dcm2niixfile.suffixes)
                    newbidsfile = str(dcm2niixfile).rsplit(ext,1)[0].rsplit('_Crop_',1)[0] + ext
                    LOGGER.info(f"Found dcm2niix _Crop_ postfix, replacing original file\n{dcm2niixfile} ->\n{newbidsfile}")
                    dcm2niixfile.replace(newbidsfile)

            # Rename all files that got additional postfixes from dcm2niix. See: https://github.com/rordenlab/dcm2niix/blob/master/FILENAMING.md
            dcm2niixpostfixes = ('_c', '_i', '_Eq', '_real', '_imaginary', '_MoCo', '_t', '_Tilt', '_e', '_ph')
            dcm2niixfiles     = sorted(set([dcm2niixfile for dcm2niixpostfix in dcm2niixpostfixes for dcm2niixfile in outfolder.glob(f"{bidsname}*{dcm2niixpostfix}*.nii*")]))
            if not jsonfiles[0].is_file() and dcm2niixfiles:                                                    # Possibly renamed by dcm2niix, e.g. with multi-echo data (but not always for the first echo)
                jsonfiles.pop(0)
            for dcm2niixfile in dcm2niixfiles:
                ext         = ''.join(dcm2niixfile.suffixes)
                postfixes   = str(dcm2niixfile).split(bidsname)[1].rsplit(ext)[0].split('_')[1:]
                newbidsname = dcm2niixfile.name                                                                 # Strip the additional postfixes and assign them to bids entities in the for-loop below
                for postfix in postfixes:                                                                       # dcm2niix postfixes _c%d, _e%d and _ph (and any combination of these in that order) are for multi-coil data, multi-echo data and phase data
                    newbidsname = newbidsname.replace(postfix, '').replace('__.', '.').replace('_.', '.')
                    newbidsfile = outfolder/newbidsname
                LOGGER.info(f"Found dcm2niix {postfixes} postfixes, renaming\n{dcm2niixfile} ->\n{newbidsfile}")
                if newbidsfile.is_file():
                    LOGGER.warning(f"Overwriting existing {newbidsfile} file -- check your results carefully!")
                dcm2niixfile.replace(newbidsfile)

                # Rename all associated files (i.e. the json-, bval- and bvec-files)
                oldjsonfile = dcm2niixfile.with_suffix('').with_suffix('.json')
                newjsonfile = newbidsfile.with_suffix('').with_suffix('.json')
                if not oldjsonfile.is_file():
                    LOGGER.warning(f"Unexpected file conversion result: {oldjsonfile} not found")
                else:
                    if oldjsonfile in jsonfiles:
                        jsonfiles.remove(oldjsonfile)
                    if newjsonfile not in jsonfiles:
                        jsonfiles.append(newjsonfile)
                for oldfile in outfolder.glob(dcm2niixfile.with_suffix('').stem + '.*'):
                    oldfile.replace(newjsonfile.with_suffix(''.join(oldfile.suffixes)))

        # Loop over and adapt all the newly produced json files and write to the scans.tsv file (NB: assumes every nifti-file comes with a json-file)
        for jsonfile in sorted(set(jsonfiles)):

            # Add a dummy b0 bval- and bvec-file for any file without a bval/bvec file (e.g. sbref, b0 scans)
            if datatype == 'dwi':
                bvecfile = jsonfile.with_suffix('.bvec')
                bvalfile = jsonfile.with_suffix('.bval')
                if not bvecfile.is_file():
                    LOGGER.info(f"Adding dummy bvec file: {bvecfile}")
                    bvecfile.write_text('0\n0\n0\n')
                if not bvalfile.is_file():
                    LOGGER.info(f"Adding dummy bval file: {bvalfile}")
                    bvalfile.write_text('0\n')

            # Load the json meta-data
            with jsonfile.open('r') as json_fid:
                jsondata = json.load(json_fid)

            # Add the TaskName to the meta-data
            if datatype == 'func' and 'TaskName' not in jsondata:
                jsondata['TaskName'] = run['bids']['task']

            # Add the TracerName and TaskName to the meta-data
            elif datatype == 'pet' and 'TracerName' not in jsondata:
                jsondata['TracerName'] = run['bids']['trc']

            # Add all the meta data to the json-file except `IntendedFor`, which is handled separately later
            for metakey, metaval in run['meta'].items():
                if metakey != 'IntendedFor':
                    LOGGER.info(f"Adding '{metakey}: {metaval}' to: {jsonfile}")
                    metaval = bids.get_dynamicvalue(metaval, sourcefile, cleanup=False, runtime=True)
                jsondata[metakey] = metaval
            with jsonfile.open('w') as json_fid:
                json.dump(jsondata, json_fid, indent=4)

            # Parse the acquisition time from the json file or else from the source header (NB: assuming the source file represents the first acquisition)
            outputfile = [file for file in jsonfile.parent.glob(jsonfile.stem + '.*') if file.suffix in ('.nii','.gz')]     # Find the corresponding nifti/tsv.gz file (there should be only one, let's not make assumptions about the .gz extension)
            if not outputfile:
                LOGGER.exception(f"No data-file found with {jsonfile} when updating {scans_tsv}")
            elif datatype not in bidsmap['Options']['bidscoin']['bidsignore'] and not run['bids']['suffix'] in bids.get_derivatives(datatype):
                if 'AcquisitionTime' not in jsondata or not jsondata['AcquisitionTime']:
                    jsondata['AcquisitionTime'] = bids.get_sourcevalue('AcquisitionTime', sourcefile)       # DICOM
                if not jsondata['AcquisitionTime']:
                    jsondata['AcquisitionTime'] = bids.get_sourcevalue('exam_date', sourcefile)             # PAR/XML
                try:
                    acq_time = dateutil.parser.parse(jsondata['AcquisitionTime'])
                    acq_time = '1925-01-01T' + acq_time.strftime('%H:%M:%S')                                # Privacy protection (see BIDS specification)
                except Exception as jsonerror:
                    LOGGER.warning(f"Could not parse the acquisition time from: {sourcefile}\n{jsonerror}")
                    acq_time = 'n/a'
                scanpath = outputfile[0].relative_to(bidsses)
                scans_table.loc[scanpath.as_posix(), 'acq_time'] = acq_time

    # Write the scans_table to disk
    LOGGER.info(f"Writing acquisition time data to: {scans_tsv}")
    scans_table.sort_values(by=['acq_time','filename'], inplace=True)
    scans_table.to_csv(scans_tsv, sep='\t', encoding='utf-8')

    # Add IntendedFor search results and TE1+TE2 meta-data to the fieldmap json-files. This has been postponed until all datatypes have been processed (i.e. so that all target images are indeed on disk)
    if (bidsses/'fmap').is_dir():
        for jsonfile in sorted((bidsses/'fmap').glob('sub-*.json')):

            # Load the existing meta-data
            with jsonfile.open('r') as json_fid:
                jsondata = json.load(json_fid)

            # Search for the imaging files that match the IntendedFor search criteria
            niifiles    = []
            intendedfor = jsondata.get('IntendedFor')
            if intendedfor:
                # Search with multiple patterns in all runs and store the relative path to the subject folder
                if intendedfor.startswith('<') and intendedfor.endswith('>'):
                    intendedfor = intendedfor[2:-2].split('><')
                elif not isinstance(intendedfor, list):
                    intendedfor = [intendedfor]
                for selector in intendedfor:
                    niifiles.extend([Path(niifile).relative_to(bidsfolder/subid) for niifile in sorted(bidsses.rglob(f"*{selector}*.nii*")) if selector])

                # Add the IntendedFor data
                if niifiles:
                    LOGGER.info(f"Adding IntendedFor to: {jsonfile}")
                    jsondata['IntendedFor'] = [niifile.as_posix() for niifile in niifiles]  # The path needs to use forward slashes instead of backward slashes
                else:
                    LOGGER.warning(f"Empty 'IntendedFor' fieldmap value in {jsonfile}: the search for {intendedfor} gave no results")
                    jsondata['IntendedFor'] = ''
            else:
                LOGGER.warning(f"Empty 'IntendedFor' fieldmap value in {jsonfile}: the IntendedFor value of the bidsmap entry was empty")

            # Extract the echo times from magnitude1 and magnitude2 and add them to the phasediff json-file
            if jsonfile.name.endswith('phasediff.json'):
                json_magnitude = [None, None]
                echotime       = [None, None]
                for n in (0,1):
                    json_magnitude[n] = jsonfile.parent/jsonfile.name.replace('_phasediff', f"_magnitude{n+1}")
                    if not json_magnitude[n].is_file():
                        LOGGER.error(f"Could not find expected magnitude{n+1} image associated with: {jsonfile}")
                    else:
                        with json_magnitude[n].open('r') as json_fid:
                            data = json.load(json_fid)
                        echotime[n] = data['EchoTime']
                jsondata['EchoTime1'] = jsondata['EchoTime2'] = None
                if None in echotime:
                    LOGGER.error(f"Cannot find and add valid EchoTime1={echotime[0]} and EchoTime2={echotime[1]} data to: {jsonfile}")
                elif echotime[0] > echotime[1]:
                    LOGGER.error(f"Found invalid EchoTime1={echotime[0]} > EchoTime2={echotime[1]} for: {jsonfile}")
                else:
                    jsondata['EchoTime1'] = echotime[0]
                    jsondata['EchoTime2'] = echotime[1]
                    LOGGER.info(f"Adding EchoTime1: {echotime[0]} and EchoTime2: {echotime[1]} to {jsonfile}")

            # Save the collected meta-data to disk
            with jsonfile.open('w') as json_fid:
                json.dump(jsondata, json_fid, indent=4)

    # Collect personal data from a source header (PAR/XML does not contain personal info)
    if dataformat=='DICOM' and sourcefile.name:
        personals['participant_id'] = subid
        if sesid:
            if 'session_id' not in personals:
                personals['session_id'] = sesid
            else:
                return                                              # Only take data from the first session -> BIDS specification
        age = bids.get_dicomfield('PatientAge', sourcefile)         # A string of characters with one of the following formats: nnnD, nnnW, nnnM, nnnY
        if age.endswith('D'):
            personals['age'] = str(int(float(age.rstrip('D'))/365.2524))
        elif age.endswith('W'):
            personals['age'] = str(int(float(age.rstrip('W'))/52.1775))
        elif age.endswith('M'):
            personals['age'] = str(int(float(age.rstrip('M'))/12))
        elif age.endswith('Y'):
            personals['age'] = str(int(float(age.rstrip('Y'))))
        elif age:
            personals['age'] = age
        personals['sex']     = bids.get_dicomfield('PatientSex',    sourcefile)
        personals['size']    = bids.get_dicomfield('PatientSize',   sourcefile)
        personals['weight']  = bids.get_dicomfield('PatientWeight', sourcefile)
