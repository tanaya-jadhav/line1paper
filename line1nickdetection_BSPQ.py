# input:
# Bionano xmap, cmap and rmap files that have been filtered to only contain contigs that cross specified region
# output: 1 file per region containing details of insertions found in each sample and contig
# author: Tanaya Jadhav

import pandas as pd
import re
import datetime
import glob
import os


def closest(lst, K):
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - K))]


def all_non_consecutive(arr):
    ans = []
    start = arr[0]
    index = 0
    for number in arr:
        if start == number:
            start += 1
            index += 1
            continue

        ans.append({'i': index, 'n': number})
        start = number + 1
        index += 1
    return ans


def main():
    # xmaplist = 'xmapfilestorun.txt'
    # with open(xmaplist, 'r') as i:
    #     xmaps = i.readlines()
    # for file in xmaps:
    #     xmap = file.strip('\n')
    #     print(xmap)
    for xmap in glob.glob('*.xmap', recursive=True):
        print(xmap)
        timestart = datetime.datetime.now()
        # xmap = '73606580_73618557_OpenNew_contig1.xmap'
        # print(xmap)
        rmap = 'EXP_REFINEFINAL1_r.cmap'
        ourverybeginning = int(xmap.split('_')[0])
        ourveryend = int(xmap.split('_')[1])
        refcontig = xmap.split('_')[-1].split('.')[0]
        qmap = str(ourverybeginning) +'_'+ str(ourveryend) + '_OpenNew_q.cmap'
        filteredfile = 'BSPQ_' + str(ourverybeginning) + '_' + str(ourveryend) + '_' + refcontig + '_line1_onlyfiltered.tsv'
        fullsetfile = 'BSPQ_' + str(ourverybeginning) + '_' + str(ourveryend) + '_' + refcontig + '_line1_fullset.tsv'
        xmap_headercols = ['XmapEntryID', 'QryContigID', 'RefContigID', 'QryStartPos', 'QryEndPos', 'RefStartPos',
                           'RefEndPos', 'Orientation', 'Confidence', 'HitEnum', 'QryLen', 'RefLen', 'LabelChannel',
                           'Alignment', 'SampleID']
        xmap_df = pd.read_csv(xmap, sep='\t', comment='#', header=None, names=xmap_headercols)
        rmap_headercols = ['CMapId','ContigLength','NumSites','SiteID','LabelChannel','Position','StdDev','Coverage',
                           'Occurrence','ChimQuality','SegDupL','SegDupR','FragileL','FragileR','OutlierFrac','ChimNorm']
        rmap_df = pd.read_csv(rmap, sep='\t', comment='#', header=None, names=rmap_headercols)
        # qmap_headercols = ['CMapId','ContigLength','NumSites','SiteID','LabelChannel','Position','StdDev','Coverage',
        #                    'Occurrence','SampleID']
        qmap_headercols = ['CMapId', 'ContigLength', 'NumSites', 'SiteID', 'LabelChannel', 'Position', 'StdDev', 'Coverage',
                           'Occurrence', 'ChimQuality', 'SegDupL', 'SegDupR', 'FragileL', 'FragileR', 'OutlierFrac',
                           'ChimNorm']
                           # 'ChimNorm', 'Mask', 'SampleID']
        try:
            qmap_df = pd.read_csv(qmap, sep='\t', comment='#', header=None, names=qmap_headercols)
        except FileNotFoundError:
            print(xmap, 'does not have a qmap')
            # continue
        filteredlist = []
        fullsetlist = []
        for index, row in xmap_df.iterrows():
            querycontig = row['QryContigID']
            querylen = row['QryLen']
            alignment = row['Alignment']
            chrom = row['RefContigID']
            # contiglen = row['QryLen']
            # sub_qmap = qmap_df.query("CMapId == @querycontig and ContigLength == @contiglen")
            sample = row['SampleID']
            # print(row)
            # sub_qmap = qmap_df.query("CMapId == @querycontig and SampleID == @sample")
            sub_qmap = qmap_df.query("CMapId == @querycontig")
            chrom_ref = rmap_df.query("CMapId == @chrom")
            refdict = pd.Series(chrom_ref.SiteID.values, index=chrom_ref.Position).to_dict()
            querydict = pd.Series(sub_qmap.Position.values, index=sub_qmap.SiteID).to_dict()
            if sample == 'HG02623' and len(querydict) == 0:
                sub_qmap = qmap_df.query("CMapId == @querycontig and ChimQuality == @sample")
                querydict = pd.Series(sub_qmap.Position.values, index=sub_qmap.SiteID).to_dict()
            allquerysites = list(querydict.keys())
            keylist = list(refdict.keys())
            refstart_nicksite = closest(keylist, ourverybeginning)
            refend_nicksite = closest(keylist, ourveryend)

            alignmentlist = re.split('\(|\)', alignment)
            alignmentlist = list(filter(None, alignmentlist))
            alignment_df = pd.DataFrame(alignmentlist)
            # print(alignmentlist)
            alignment_df[['refsite', 'querysite']] = alignment_df[0].str.split(',',expand=True)
            alignment_df.drop(columns=0, inplace=True)
            alignment_df = alignment_df.set_index('refsite', drop=False)
            reflist = alignment_df.index.values.tolist()
            # print(reflist)
            reflist = [int(i) for i in reflist]
            refdictstart = refdict[refstart_nicksite] - 1
            refdictend = refdict[refend_nicksite] + 1
            # print(refdictstart, refdictend)
            # print(reflist, refdictstart, min(reflist))
            while refdictstart not in reflist and refdictstart > min(reflist):
                # print(refdictstart, 'in refstart while loop')
                refdictstart = refdictstart - 1
            if refdictstart not in reflist and refdictstart <= min(reflist):
                refdictstart = min(reflist)
            while refdictend not in reflist and refdictend < max(reflist):
                # print(refdictend, 'in refend while loop')
                refdictend = refdictend + 1
            if refdictend not in reflist and refdictend >= max(reflist):
                refdictend = max(reflist)
            # if refdictend < refdictstart:
            #     temp = refdictend
            #     refdictend = refdictstart
            #     refdictstart = temp
            alignment_df = alignment_df[str(refdictstart):str(refdictend)]
            # print(sample)
            # print(refdictstart, refdictend)
            # print(alignment_df)
            alignment_df = alignment_df.astype('int')
            alignment_df = alignment_df.drop_duplicates(subset='querysite', keep='first')
            alignment_df = alignment_df.set_index('querysite', drop=True)
            querylist = alignment_df.index.values.tolist()
            # print(alignment_df)
            sorted_querylist = querylist.copy()
            sorted_querylist.sort()

            skippedlist = all_non_consecutive(sorted_querylist)
            possible_line1 = []
            for skipped in skippedlist:
                prev = skipped['i'] - 1
                current = skipped['i']
                if abs(sorted_querylist[current]-sorted_querylist[prev]) == 2 or abs(sorted_querylist[current]-sorted_querylist[prev]) == 3:
                    region = (sorted_querylist[prev], sorted_querylist[current])
                    possible_line1.append(region)

            if len(possible_line1) >= 1:
                for region in possible_line1:
                    insertionstart = region[0]
                    insertionend = region[1]
                    # print(sample, region, querydict)
                    insertionstart_pos = querydict[insertionstart]
                    insertionend_pos = querydict[insertionend]
                    nick1 = insertionstart + 1
                    nick1_pos = querydict[nick1]
                    D1 = nick1_pos - insertionstart_pos
                    if abs(insertionend - insertionstart) == 2:
                        nick2 = ''
                        nick2_pos = ''
                        D2 = 0
                        D3 = insertionend_pos - nick1_pos
                    else:
                        nick2 = insertionstart + 2
                        nick2_pos = querydict[nick2]
                        D2 = nick2_pos - nick1_pos
                        D3 = insertionend_pos - nick2_pos
                    refsitestart = alignment_df.loc[insertionstart]['refsite']
                    refsiteend = alignment_df.loc[insertionend]['refsite']
                    refsitestart_pos = list(refdict.keys())[list(refdict.values()).index(refsitestart)]
                    refsiteend_pos = list(refdict.keys())[list(refdict.values()).index(refsiteend)]
                    confirmeddict = {}
                    # confirmeddict['Sample'] = 'Dummy'
                    confirmeddict['Sample'] = sample
                    confirmeddict['Chrom'] = chrom
                    confirmeddict['ContigID'] = querycontig
                    confirmeddict['ContigLength'] = querylen
                    confirmeddict['RegionStart'] = ourverybeginning
                    confirmeddict['RegionEnd']= ourveryend
                    confirmeddict['RefSites'] = list(alignment_df['refsite'])
                    confirmeddict['QuerySites'] = querylist
                    confirmeddict['RefSiteStart'] = refsitestart
                    confirmeddict['RefSiteEnd'] = refsiteend
                    confirmeddict['QuerySiteStart'] = insertionstart
                    confirmeddict['InsertionSite1'] = nick1
                    confirmeddict['InsertionSite2'] = nick2
                    confirmeddict['QuerySiteEnd'] = insertionend
                    confirmeddict['RefPositionStart'] = refsitestart_pos
                    confirmeddict['RefPositionEnd'] = refsiteend_pos
                    confirmeddict['RefDistance'] = refsiteend_pos - refsitestart_pos
                    confirmeddict['QueryPositionStart'] = insertionstart_pos
                    confirmeddict['InsertionPosition1'] = nick1_pos
                    confirmeddict['InsertionPosition2'] = nick2_pos
                    confirmeddict['QueryPositionEnd'] = insertionend_pos
                    confirmeddict['D1'] = D1
                    confirmeddict['D2'] = D2
                    confirmeddict['D3'] = D3
                    fullsetlist.append(confirmeddict)
                    filteredlist.append(confirmeddict)
                    # print(confirmeddict)
            else:
                confirmeddict = {}
                # confirmeddict['Sample'] = 'Dummy'
                confirmeddict['Sample'] = sample
                confirmeddict['Chrom'] = chrom
                confirmeddict['ContigID'] = querycontig
                confirmeddict['ContigLength'] = querylen
                confirmeddict['RegionStart'] = ourverybeginning
                confirmeddict['RegionEnd'] = ourveryend
                confirmeddict['RefSites'] = list(alignment_df['refsite'])
                confirmeddict['QuerySites'] = querylist
                fullsetlist.append(confirmeddict)
        fullsetdf = pd.DataFrame(fullsetlist)
        filtereddf = pd.DataFrame(filteredlist)
        fullsetdf.to_csv(fullsetfile, sep='\t', index=None)
        filtereddf.to_csv(filteredfile, sep='\t', index=None)
        timeend = datetime.datetime.now()
        print(xmap, "TimeTaken: ", timeend-timestart)

        filelist = [xmap, qmap, fullsetfile, filteredfile]
        for file in filelist:
            cmd = 'mv ' + file + ' completed/'
            os.system(cmd)


if __name__ == '__main__':
    main()
