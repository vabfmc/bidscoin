## ************************************************************************
##
## Provides the acquisition slice order for an FFE-EPI fMRI sequence
##
## Validated for R5.3.1, assumes 1 package
##
## Guillaume Gilbert
## Philips Healthcare Canada
## 2017,2018
##
##
## Bug fix 2018-12-05: MB-SENSE with default ordering
## Bug fix 2019-05-21: MB-SENSE with default ordering and odd number of slices per package
## ************************************************************************
import numpy as np
import numpy.matlib
import snoop
import logging

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S'
)
#clear('all')
## 1 7 13 19 25 31 37 2 8 14 20 26 32 38 3 28 34 40 5 11 17 23 29 35 41 6 12 18 24 30 36 42
MB_SENSE_factor = None
Slices_per_band = None
Slice_order = None

# @snoop
def test():
    global Slices_per_band, MB_SENSE_factor, Slice_order

    Slices = 69
    steps = np.round(np.sqrt(Slices))
    slice = 1
    sliceind = 1
    startind = 1
    sliceorder = np.zeros(Slices)
    while sliceind < Slices:

        sliceorder[sliceind] = slice
        sliceind = sliceind + 1
        slice = slice + steps
        if slice > Slices:
            startind = startind + 1
            slice = startind


    ## ************************************************************************
    ## Sequence parameters
    MB_SENSE = 'yes'
    MB_SENSE_factor = 3
    Slice_scan_order = 'FH'

    TR = 1.5
    ## ************************************************************************

    #We first perform a sanity check
    if (np.mod(Slices,MB_SENSE_factor) != 0 and str(MB_SENSE) == str('yes')):
        raise Exception('These parameters are impossible')

    Slice_order = []
    #FH
    if (str(Slice_scan_order) == str('FH')):
        if (str(MB_SENSE) == str('yes')):
            Slices_per_band = Slices / MB_SENSE_factor
            for k in np.arange(1,MB_SENSE_factor+1).reshape(-1):
                Slice_order[k:] = np.array([np.arange(Slices_per_band * (k - 1),Slices_per_band * (k - 1) + (Slices_per_band) - 1+1,1)])
        else:
            Slice_order = np.array([np.arange(0,Slices - 1+1,1)])
        #HF
    else:
        if (str(Slice_scan_order) == str('HF')):
            if (str(MB_SENSE) == str('yes')):
                Slices_per_band = Slices / MB_SENSE_factor
                for k in np.arange(1,MB_SENSE_factor+1).reshape(-1):
                    Slice_order[k:] = np.array([np.arange(Slices_per_band * (k - 1) + (Slices_per_band) - 1,Slices_per_band * (k - 1)+- 1,- 1)])
            else:
                Slice_order = np.array([np.arange(Slices - 1,0+- 1,- 1)])
            #rev. central
        else:
            if (str(Slice_scan_order) == str('rev. central')):
                if (str(MB_SENSE) == str('yes')):
                    Slices_per_band = Slices / MB_SENSE_factor
                    for k in np.arange(1,MB_SENSE_factor+1).reshape(-1):
                        up = np.array([np.arange(Slices_per_band * (k - 1),Slices_per_band * (k - 1) + int(np.floor((Slices_per_band - 1) / 2))+1,1)])
                        down = np.array([np.arange(Slices_per_band * (k - 1) + Slices_per_band - 1,Slices_per_band * (k - 1) + np.ceil((Slices_per_band) / 2)+- 1,- 1)])
                        Slice_order[k,np.arange[1,Slices_per_band+2,2]] = up
                        Slice_order[k,np.arange[2,Slices_per_band+2,2]] = down
                else:
                    up = np.array([np.arange(0,int(np.floor((Slices - 1) / 2))+1,1)])
                    down = np.array([np.arange(Slices - 1,np.ceil((Slices) / 2)+- 1,- 1)])
                    Slice_order[np.arange[1,Slices+2,2]] = up
                    Slice_order[np.arange[2,Slices+2,2]] = down
                #interleaved
            else:
                if (str(Slice_scan_order) == str('interleaved')):
                    if (str(MB_SENSE) == str('yes')):
                        Slices_per_band = Slices / MB_SENSE_factor
                        SliceGroup = 1
                        temp_order = 1
                        step = np.round(np.sqrt(double(Slices_per_band)))
                        for k in np.arange(2,Slices_per_band+1).reshape(-1):
                            current = temp_order(k - 1) + step
                            if (current > Slices_per_band):
                                SliceGroup = SliceGroup + 1
                                current = SliceGroup
                            temp_order = np.array([temp_order,current])
                        for k in np.arange(1,MB_SENSE_factor+1).reshape(-1):
                            Slice_order[k,:] = (temp_order - 1) + (k - 1) * Slices_per_band
                    else:
                        SliceGroup = 1
                        Slice_order = 1
                        step = np.round(np.sqrt(double(Slices)))
                        for k in np.arange(2,Slices+1).reshape(-1):
                            current = Slice_order(k - 1) + step
                            if (current > Slices):
                                SliceGroup = SliceGroup + 1
                                current = SliceGroup
                            Slice_order = np.array([Slice_order,current])
                        Slice_order = Slice_order - 1
                    #default
                else:
                    if (str(Slice_scan_order) == str('default')):
                        if (str(MB_SENSE) == str('yes')):
                            Slices_per_band = Slices / MB_SENSE_factor
                            # The are a few special cases to consider here
                            if (Slices_per_band <= 6):
                                step = 2
                                temp_order = np.zeros((1,Slices_per_band))
                                half_locs_per_package = (Slices_per_band + 1) / step
                                low_part_start_loc = 0
                                low_part_act_loc = low_part_start_loc
                                high_part_start_loc = Slices_per_band - 1
                                high_part_act_loc = high_part_start_loc
                                order = 1
                                temp_order = temp_order + 1
                            else:
                                if (Slices_per_band == 8):
                                    factor = 1
                                    SliceGroup = 1
                                    temp_order = 1
                                    step = np.round(np.sqrt(double(Slices_per_band)))
                                    for k in np.arange(2,Slices_per_band+1).reshape(-1):
                                        current = temp_order(k - 1) + step
                                        if (current > Slices_per_band):
                                            SliceGroup = SliceGroup + 1
                                            current = SliceGroup
                                        temp_order = np.array([temp_order,current])
                                else:
                                    if (np.mod(Slices_per_band,2) == 0):
                                        step = 2
                                        temp_order = np.zeros((1,Slices_per_band))
                                        half_locs_per_package = (Slices_per_band + 1) / step
                                        low_part_start_loc = 0
                                        low_part_act_loc = low_part_start_loc
                                        high_part_start_loc = Slices_per_band - 1
                                        high_part_act_loc = high_part_start_loc
                                        order = 1
                                        temp_order = temp_order + 1
                                    else:
                                        part1 = np.array([np.arange(0,(Slices_per_band - 1)+2,2)])
                                        part2 = np.array([np.arange(1,(Slices_per_band - 1)+2,2)])
                                        temp_order = cat(2,part1,part2)
                                        temp_order = temp_order + 1
                            for k in np.arange(1,MB_SENSE_factor+1).reshape(-1):
                                Slice_order[k,:] = (temp_order - 1) + (k - 1) * Slices_per_band
                        else:
                            part1 = np.array([np.arange(0,(Slices - 1)+2,2)])
                            part2 = np.array([np.arange(1,(Slices - 1)+2,2)])
                            Slice_order = cat(2,part1,part2)

    slice_timings = slicetime()
    for slice_timing in slice_timings:
        for element in slice_timings:
            for k in element:
                print(k)


    logging.info(slice_timings)

# @snoop
def slicetime(TRsec=1.5, Slices=69, isAscending=True, isSequential=False, DelayBetweenVolumesSec=0):

    #compute slice timing
    # TRsec : sampling rate (TR) in seconds
    # nSlices: number of slices in each 3D volume
    # isAscending: ascending (true) or descending (false) order
    # isSequential: interleaved (false) or sequential (true) order
    # DelayBetweenVolumesSec: pause between final slice of volume and start of next
    #Examples
    # sliceTime(2.0, 10, true, true); #ascending, sequential
    # sliceTime(2.0, 10, false, true); #descending, sequential
    # sliceTime(2.0, 10, true, false); #ascending, interleaved
    # sliceTime(2.0, 10, false, false); #descending, interleaved
    # sliceTime(3.0, 10, true, true, 1.0); #ascending, sequential, sparse


    TA = (TRsec - DelayBetweenVolumesSec) / Slices

    bidsSliceTiming = np.array([np.arange(0,TRsec - TA+TA,TA)])

    if not isAscending :
        bidsSliceTiming = flip(bidsSliceTiming)

    # if not isSequential :
    #     if not np.mod(Slices,2) :
    #         logging.info('Timings for Philips or GE. Siemens volume with even number of slices differs https://www.mccauslandcenter.sc.edu/crnl/tools/stc)\n' % ())
    #     order = np.array([np.arange(1,Slices+2,2),np.arange(2,Slices+2,2)])
    #     #order = Slice_order;#MM
    #     bidsSliceTiming[order] = bidsSliceTiming

    # #report results
    # logging.info('"SliceTiming": [\n' % ())
    # for i in np.arange(1,Slices+1).reshape(-1):
    #     logging.info('		%.3f' % (bidsSliceTiming(i)))
    #     if (i < Slices):
    #         logging.info(',\n' % ())
    #     else:
    #         logging.info('	],\n' % ())

    TA = (TRsec - DelayBetweenVolumesSec) / (Slices / MB_SENSE_factor)

    bidsSliceTiming = np.array([np.arange(0,TRsec - TA+TA,TA)])

    if not isAscending :
        bidsSliceTiming = flip(bidsSliceTiming)

    if not isSequential :

        # logging.info(f"{Slices} {np.mod(Slices, 2)}")

        if not np.mod(Slices,2) == 0 :
            logging.info('Timings for Philips or GE. Siemens volume with even number of slices differs https://www.mccauslandcenter.sc.edu/crnl/tools/stc)\n' % ())
        
        else:
            order = np.array([np.arange(1,Slices_per_band+1,1)])
            #order = Slice_order;#MM
            bidsSliceTiming[order] = bidsSliceTiming

    bidsSliceTiming = np.matlib.repmat(bidsSliceTiming,1,MB_SENSE_factor)
    for idx, array_ in enumerate(Slice_order):
        logging.info(f"row {idx}: {array_}")


    logging.info(bidsSliceTiming.shape)

    # #report results MB
    # logging.info('"SliceTimingMB": [\n' % ())
    # for i in np.arange(0,Slices+1).reshape(-1):
    #     logging.info('		%.3f' % (bidsSliceTiming[i]))
    #     if (i < Slices):
    #         logging.info(',\n' % ())
    #     else:
    #         logging.info('	],\n' % ())

    return bidsSliceTiming


if __name__ == "__main__":
    test()
