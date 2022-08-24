# input:
# intersect bed file with names of bspq and dle files that intersect
# dle and bspq files that contain all insertions (output from line1nickdetection_{DLE/BSPQ}.py)
# output: 1 file per overlap unless there are no common samples found in both files
# author: Tanaya Jadhav

import pandas as pd
import glob
import os.path
import re

def get_orientation(querylist):
    if sorted(querylist) == querylist:
        return '+'
    else:
        return '-'


def determine_dystances(dlesample, bspqsample, dle_df, bspq_df, relativeposlist, locationdict):
    outputdict = {}
    outputdict['Sample'] = dlesample + ',' + bspqsample
    outputdict['Chrom'] = bspq_df.loc[bspqsample]['Chrom']
    outputdict['DLE_ContigID'] = dle_df.loc[dlesample]['ContigID']
    outputdict['BSPQ_ContigID'] = bspq_df.loc[bspqsample]['ContigID']
    outputdict['DLE_Region'] = str(dle_df.loc[dlesample]['RegionStart']) + '-' + str(dle_df.loc[dlesample]['RegionEnd'])
    outputdict['BSPQ_Region'] = str(bspq_df.loc[bspqsample]['RegionStart']) + '-' + str(bspq_df.loc[bspqsample]['RegionEnd'])

    # determine order of inserted nick sites
    relativeposlist = list(set(relativeposlist))
    relativeposlist.sort()

    # get distances betwen adjacent nick sites
    i = 0
    for currentpos in relativeposlist[:-1]:
        i = i+1
        nextpos = relativeposlist[i]
        distance = nextpos - currentpos
        firstnickstr = 'nick' + str(i) + 'pos'
        secondnickstr = 'nick' + str(i+1) + 'pos'
        diststr = 'D' + str(i)
        outputdict[firstnickstr] = (currentpos, locationdict[currentpos])
        outputdict[diststr] = distance
        outputdict[secondnickstr] = (nextpos, locationdict[nextpos])
    return outputdict


def get_bspq_locations(sample, bspq_df, bspq_orientation):
    # get relative locations for bspq nick sites
    if bspq_orientation == '+':
        bspq_startpos = bspq_df.loc[sample]['RefPositionStart']

        bspq_nick1pos = bspq_startpos + bspq_df.loc[sample]['D1']
        if bspq_df.loc[sample]['D2'] != 0:
            bspq_nick2pos = bspq_nick1pos + bspq_df.loc[sample]['D2']
        else:
            bspq_nick2pos = bspq_nick1pos
    else:
        bspq_startpos = bspq_df.loc[sample]['RefPositionEnd']
        bspq_nick1pos = bspq_startpos + bspq_df.loc[sample]['D3']
        if bspq_df.loc[sample]['D2'] != 0:
            bspq_nick2pos = bspq_nick1pos + bspq_df.loc[sample]['D2']
        else:
            bspq_nick2pos = bspq_nick1pos

    return bspq_startpos, bspq_nick1pos, bspq_nick2pos


def get_dle_locations(sample, dle_df, dle_orientation):
    # get relative locations for dle nick sites
    if dle_orientation == '+':
        dle_startpos = dle_df.loc[sample]['RefPositionStart']

        dle_nick1pos = dle_startpos + dle_df.loc[sample]['D1']
        dle_nick2pos = dle_nick1pos + dle_df.loc[sample]['D2']
        if dle_df.loc[sample]['D3'] != 0:
            dle_nick3pos = dle_nick2pos + dle_df.loc[sample]['D3']
        else:
            dle_nick3pos = dle_nick2pos
    else:
        dle_startpos = dle_df.loc[sample]['RefPositionEnd']
        dle_nick1pos = dle_startpos + dle_df.loc[sample]['D4']
        if dle_df.loc[sample]['D3'] != 0:
            dle_nick2pos = dle_nick1pos + dle_df.loc[sample]['D3']
        else:
            dle_nick2pos = dle_nick1pos
        dle_nick3pos = dle_nick2pos + dle_df.loc[sample]['D2']

    return dle_startpos, dle_nick1pos, dle_nick2pos, dle_nick3pos


def main():
    # bspq = '/Volumes/TOSHIBAEXT/line1/BSPQ/BSPQ_103921885_103943798_contig1_line1_onlyfiltered.tsv'
    # dle = '/Volumes/TOSHIBAEXT/line1/DLE/DLE_103760118_108453756_contig1_line1_onlyfiltered.tsv'
    # overlapbedfile = '565dle507bspq_overlap.bed'
    # bspqfiledir = '/Volumes/TOSHIBAEXT/line1/565BSPQ/'
    # dlefiledir = '/Volumes/TOSHIBAEXT/line1/565DLE/'
    # outdir = '/Volumes/TOSHIBAEXT/line1/565distances/'
    # overlapheader = ['chrDLE','startDLE', 'endDLE', 'DLEfile', 'chrBSPQ', 'startBSPQ', 'endBSPQ', 'BSPQfile', 'overlap']
    # overlapdf = pd.read_csv(overlapbedfile, sep='\t', names=overlapheader)
    # print(overlapdf.head(5))
    # dlefiles = list(set(list(overlapdf['DLEfile'])))
    # print(dlefiles, len(dlefiles))

    # loops through the overlap file to find pairs of dle and bspq files
    count = 0
    nocommonsamples = []
    # for index, row in overlapdf.iterrows():
    #     bspqfile = row['BSPQfile']
    #     # bspqfile = 'BSPQ_' + bspqfilename
    #     dlefile = row['DLEfile']
    #     bspq = bspqfiledir + bspqfile
    #     dle = dlefiledir + dlefile
    #     outfile =  outdir + dlefile.split('_')[0] + '-' + dlefile.split('_')[1] + '-' + dlefile.split('_')[2] + '-' + \
    #                dlefile.split('_')[3] + '_' + bspqfile.split('_')[0] + '-' + dlefile.split('_')[1] + '-' + \
    #                dlefile.split('_')[2] + '-' + dlefile.split('_')[3] + '_distancesNEW.tsv'
        # print(bspq, dle, outfile)
    for dlefile in glob.glob('NA12878_DLE/*_onlyfiltered.tsv', recursive=True):
        start = dlefile.split('_')[2]
        end = dlefile.split('_')[3]
        chrom = dlefile.split('_')[4]
        bspq = 'NA12878_BSPQ/BSPQ_' + start + '_' + end + '_' + chrom + '_line1_onlyfiltered.tsv'
        outfile = 'NA12878_distances/' + start + '_' + end + '_' + chrom + '_distances.tsv'
        outdir = 'NA12878_distances/'
        try:
            bspq_df = pd.read_csv(bspq, sep='\t', header=0, index_col='Sample')
        except FileNotFoundError:
            print(bspq, 'not found')
            continue
        try:
            dle_df = pd.read_csv(dlefile, sep='\t', header=0, index_col='Sample')
        except FileNotFoundError:
            print(dlefile, 'not found')
            continue
        bspq_samples = bspq_df.index.values.tolist()
        dle_samples = dle_df.index.values.tolist()
        bspq_samples = [re.sub("\D", "", sample) for sample in bspq_samples]
        dle_samples = [re.sub("\D", "", sample) for sample in dle_samples]
        common_samples = set(bspq_samples).intersection(set(dle_samples))
        nocommonsamplesdict = {}
        if len(common_samples) == 0:
            # print(bspqfile, 'bspq samples: ', bspq_samples, 'dle samples: ', dle_samples)
            count = count + 1
            nocommonsamplesdict['bspqfile'] = bspq
            nocommonsamplesdict['dlefile'] = dlefile
            nocommonsamplesdict['bspq_samples'] = bspq_samples
            nocommonsamplesdict['dle_samples'] = dle_samples
            nocommonsamples.append(nocommonsamplesdict)
        outputlist = []
        # loops through samples that are found in both files
        for commonsample in common_samples:
            locationdict = {}
            dle_sampledf = dle_df.query("Sample == 'NA'+@commonsample or Sample == 'HG'+@commonsample or Sample == 'GM'+@commonsample")
            bspq_sampledf = bspq_df.query("Sample == 'NA'+@commonsample or Sample == 'HG'+@commonsample or Sample == 'GM'+@commonsample")
            # print(commonsample, dle_sampledf)
            # loops through each insertion for the sample in dle
            dlesample = dle_sampledf.index.values.tolist()[0]
            bspqsample = bspq_sampledf.index.values.tolist()[0]
            for i, r in dle_sampledf.iterrows():
                dlecontig = r['ContigID']
                dlequerystart = r['QuerySiteStart']
                dleregionstart = r['RegionStart']
                dleregionend = r['RegionEnd']
                sub_dledf = dle_df.query("ContigID == @dlecontig and Sample ==@dlesample and QuerySiteStart == @dlequerystart and RegionStart == @dleregionstart and RegionEnd == @dleregionend")
                dlequerylist = sub_dledf.loc[dlesample]['QuerySites']
                dlequerylist = dlequerylist.strip('[|]').split(', ')
                dle_orientation = get_orientation(dlequerylist)
                dle_startpos, dle_nick1pos, dle_nick2pos, dle_nick3pos = get_dle_locations(dlesample, sub_dledf,
                                                                                           dle_orientation)
                locationdict[dle_startpos] = 'dle_refstart'
                locationdict[dle_nick1pos] = 'dle'
                locationdict[dle_nick2pos] = 'dle'
                locationdict[dle_nick3pos] = 'dle'

                # loops through each insertion for the sample in bspq
                for i, r in bspq_sampledf.iterrows():
                    bspqcontig = r['ContigID']
                    bspqquerystart = r['QuerySiteStart']
                    sub_bspqdf = bspq_df.query("ContigID == @bspqcontig and Sample == @bspqsample and QuerySiteStart == @bspqquerystart")
                    bspqquerylist = sub_bspqdf.loc[bspqsample]['QuerySites']
                    bspqquerylist = bspqquerylist.strip('[|]').split(', ')
                    bspq_orientation = get_orientation(bspqquerylist)
                    bspq_startpos, bspq_nick1pos, bspq_nick2pos = get_bspq_locations(bspqsample, sub_bspqdf,
                                                                                     bspq_orientation)
                    locationdict[bspq_startpos] = 'bspq_refstart'
                    locationdict[bspq_nick1pos] = 'bspq'
                    locationdict[bspq_nick2pos] = 'bspq'
                    relativeposlist = [dle_nick1pos, dle_nick2pos, dle_nick3pos, bspq_nick1pos, bspq_nick2pos]
                    outputdict = determine_dystances(dlesample, bspqsample, sub_dledf, sub_bspqdf, relativeposlist, locationdict)
                    outputlist.append(outputdict)

        if len(outputlist) != 0:
            output_df = pd.DataFrame(outputlist)
            # output_df.to_csv('testoutfile.tsv', sep='\t', index=None)
            output_df.to_csv(outfile, sep='\t', index=None)
    nocommonsamplesdf = pd.DataFrame(nocommonsamples)
    nocommonsamplesdf.to_csv(outdir+'nocommonsamplesfileNEW.tsv', sep='\t', index=None)


if __name__ == '__main__':
    main()
