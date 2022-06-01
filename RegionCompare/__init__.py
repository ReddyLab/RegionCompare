import sys
from math import isnan

import h5py
import numpy as np
import pyBigWig

from RegionCompare.ChromoRegions import ChromoRegionSet


def setGlobalVariables(args):
    global REGIONS
    global BIGWIGS
    global HDFS
    global HDF_CHROMO
    global OFFSET
    global OUTPUT

    BIGWIGS = False
    HDFS = False

    if args.bigwigs is not None:
        BIGWIGS = args.bigwigs

    if args.hdfs is not None:
        HDFS = args.hdfs

        if args.chromo:
            HDF_CHROMO = args.chromo
        else:
            exit("--chromo argument required to use --hdfs option")

    regionSet = ChromoRegionSet.loadBed(args.regions)
    blacklistRegionSet = ChromoRegionSet.loadBed(args.blacklist) if args.blacklist else None
    REGIONS = setAnlaysisRegion(regionSet, blacklistRegionSet)

    OFFSET = args.offset

    OUTPUT = args.output
    if OUTPUT is None:
        OUTPUT = sys.stdout


def setAnlaysisRegion(regionSet, blacklistRegionSet):
    regionSet.mergeRegions()

    if blacklistRegionSet is not None:
        blacklistRegionSet.mergeRegions()
        regionSet = regionSet - blacklistRegionSet

    return regionSet

def compareBigWigs(bw1, bw2, regions: ChromoRegionSet, outputFile):
    with pyBigWig.open(bw1) as bigwig1, pyBigWig.open(bw2) as bigwig2:
        for region in regions:
            chromo = region.chromo
            start = region.start
            end = region.end
            bigwigValues1 = bigwig1.values(chromo, start, end)
            bigwigValues2 = bigwig2.values(chromo, start, end)

            for i in range(end - start):
                if not (isnan(bigwigValues1[i]) or isnan(bigwigValues2[i])) and bigwigValues1[i] != bigwigValues2[i]:
                    try:
                        outputFile.write(f"{chromo}:{start + i}\t{bigwigValues1[i]}\t{bigwigValues2[i]}\t{abs(bigwigValues1[i] - bigwigValues2[i])}\n")
                    except BrokenPipeError:
                        return

def compareHdfs(hdf1, hdf2, chromo, regions, offset, outputFile):
    with h5py.File(hdf1, "r") as hdffile1, h5py.File(hdf2, "r") as hdffile2:
        dataset1 = hdffile1["covari"]
        dataset2 = hdffile2["covari"]

        for region in regions:
            if region.chromo != chromo:
                continue

            start = region.start - offset
            end = region.end - offset
            values1 = dataset1[start:end]
            values2 = dataset2[start:end]

            for i in range(end - start):
                if not np.array_equal(values1[i], values2[i]):
                    try:
                        outputFile.write(f"{region.chromo}:{start + i}\t{values1[i]}\t{values2[i]}\t{abs(values1[i] - values2[i])}\n")
                    except BrokenPipeError:
                        return

def init(args):
    setGlobalVariables(args)

def run(args):
    init(args)

    if BIGWIGS:
        compareBigWigs(*BIGWIGS, REGIONS, OUTPUT)

    if HDFS:
        np.set_printoptions(linewidth=np.inf)
        compareHdfs(*HDFS, HDF_CHROMO, REGIONS, OFFSET, OUTPUT)
