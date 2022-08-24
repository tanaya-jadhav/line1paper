# input:
# Bionano xmap, cmap and rmap files that have been filtered to only contain contigs that cross specified region
# output: 1 file per region containing details of insertions found in each sample and contig filtered on D2 and D3
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
    for xmap in glob.glob('*.xmap', recursive=True):
        print(xmap)
    # xmap = '66559506_66604186_OpenNew_contig6.xmap'
        timestart = datetime.datetime.now()
        rmap = 'EXP_REFINEFINAL1_r.cmap'
        ourverybeginning = int(xmap.split('_')[0])
        ourveryend = int(xmap.split('_')[1])
        refcontig = xmap.split('_')[-1].split('.')[0]
        qmap = str(ourverybeginning) +'_'+ str(ourveryend) + '_OpenNew_q.cmap'
        filteredfile = 'DLE_' + str(ourverybeginning) + '_' + str(ourveryend) + '_' + refcontig + '_line1_onlyfiltered.tsv'
        # fullsetfile = str(ourverybeginning) + '_' + str(ourveryend) + '_' + refcontig + '_line1_fullset.tsv'
        xmap_headercols = ['XmapEntryID', 'QryContigID', 'RefContigID', 'QryStartPos', 'QryEndPos', 'RefStartPos',
                           'RefEndPos', 'Orientation', 'Confidence', 'HitEnum', 'QryLen', 'RefLen', 'LabelChannel',
                           'Alignment', 'SampleID']
        xmap_df = pd.read_csv(xmap, sep='\t', comment='#', header=None, names=xmap_headercols)
        rmap_headercols = ['CMapId','ContigLength','NumSites','SiteID','LabelChannel','Position','StdDev','Coverage',
                           'Occurrence','ChimQuality','SegDupL','SegDupR','FragileL','FragileR','OutlierFrac','ChimNorm']
        rmap_df = pd.read_csv(rmap, sep='\t', comment='#', header=None, names=rmap_headercols)
        qmap_headercols = ['CMapId','ContigLength','NumSites','SiteID','LabelChannel','Position','StdDev','Coverage',
                           'Occurrence','ChimQuality','SegDupL','SegDupR','FragileL','FragileR','OutlierFrac','ChimNorm',
                           'Mask']
                           # 'Mask', 'SampleID']
        # try:
        qmap_df = pd.read_csv(qmap, sep='\t', comment='#', header=None, names=qmap_headercols)
        # except FileNotFoundError:
        #     print(xmap, 'does not have a qmap')
        #     continue
        # print(xmap_df.head(5))
        # print(qmap_df.head(5))
        filteredlist = []
        # fullsetlist = []
        for index, row in xmap_df.iterrows():
            querycontig = row['QryContigID']
            querylen = row['QryLen']
            alignment = row['Alignment']
            chrom = row['RefContigID']
            # contiglen = row['QryLen']
            # sub_qmap = qmap_df.query("CMapId == @querycontig and ContigLength == @contiglen")
            sample = row['SampleID']
            # sub_qmap = qmap_df.query("CMapId == @querycontig and SampleID == @sample")
            sub_qmap = qmap_df.query("CMapId == @querycontig")
            chrom_ref = rmap_df.query("CMapId == @chrom")
            refdict = pd.Series(chrom_ref.SiteID.values, index=chrom_ref.Position).to_dict()
            querydict = pd.Series(sub_qmap.Position.values, index=sub_qmap.SiteID).to_dict()
            # print(querydict)
            allquerysites = list(querydict.keys())
            keylist = list(refdict.keys())
            refstart_nicksite = closest(keylist, ourverybeginning)
            refend_nicksite = closest(keylist, ourveryend)
            alignmentlist = re.split('\(|\)', alignment)
            alignmentlist = list(filter(None, alignmentlist))
            alignment_df = pd.DataFrame(alignmentlist)
            alignment_df[['refsite', 'querysite']] = alignment_df[0].str.split(',',expand=True)
            alignment_df.drop(columns=0, inplace=True)
            alignment_df = alignment_df.set_index('refsite', drop=False)
            reflist = alignment_df.index.values.tolist()
            refdictstart = refdict[refstart_nicksite]
            refdictend = refdict[refend_nicksite]
            while str(refdictstart) not in reflist:
                refdictstart = refdictstart - 1
            while str(refdictend) not in reflist:
                refdictend = refdictend + 1
            alignment_df = alignment_df[str(refdictstart):str(refdictend)]
            alignment_df = alignment_df.astype('int')
            alignment_df = alignment_df.drop_duplicates(subset='querysite', keep='first')
            alignment_df = alignment_df.set_index('querysite', drop=True)
            querylist = alignment_df.index.values.tolist()
            sorted_querylist = querylist.copy()
            sorted_querylist.sort()
            skippedlist = all_non_consecutive(sorted_querylist)
            possible_line1 = []
            for skipped in skippedlist:
                prev = skipped['i'] - 1
                current = skipped['i']
                if abs(sorted_querylist[current]-sorted_querylist[prev]) == 3 or abs(sorted_querylist[current]-sorted_querylist[prev]) == 4:
                    region = (sorted_querylist[prev], sorted_querylist[current])
                    possible_line1.append(region)
            if len(possible_line1) >= 1:
                for region in possible_line1:
                    insertionstart = region[0]
                    insertionend = region[1]
                    insertionstart_pos = querydict[insertionstart]
                    insertionend_pos = querydict[insertionend]
                    nick1 = insertionstart + 1
                    nick2 = insertionstart + 2
                    nick1_pos = querydict[nick1]
                    nick2_pos = querydict[nick2]
                    D1 = nick1_pos - insertionstart_pos
                    D2 = nick2_pos - nick1_pos
                    if insertionend - insertionstart == 3:
                        nick3 = ''
                        nick3_pos = ''
                        D3 = 0
                        D4 = insertionend_pos - nick2_pos
                    else:
                        nick3 = insertionstart + 3
                        nick3_pos = querydict[nick3]
                        D3 = nick3_pos - nick2_pos
                        D4 = insertionend_pos - nick3_pos
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
                    confirmeddict['InsertionSite3'] = nick3
                    confirmeddict['QuerySiteEnd'] = insertionend
                    confirmeddict['RefPositionStart'] = refsitestart_pos
                    confirmeddict['RefPositionEnd'] = refsiteend_pos
                    confirmeddict['RefDistance'] = refsiteend_pos - refsitestart_pos
                    confirmeddict['QueryPositionStart'] = insertionstart_pos
                    confirmeddict['InsertionPosition1'] = nick1_pos
                    confirmeddict['InsertionPosition2'] = nick2_pos
                    confirmeddict['InsertionPosition3'] = nick3_pos
                    confirmeddict['QueryPositionEnd'] = insertionend_pos
                    confirmeddict['D1'] = D1
                    confirmeddict['D2'] = D2
                    confirmeddict['D3'] = D3
                    confirmeddict['D4'] = D4
                    # fullsetlist.append(confirmeddict)
                    if 3000 <= D2 <= 4400 or 3000 <= D3 <= 4400:
                        filteredlist.append(confirmeddict)
                    # print(confirmeddict)
            # else:
            #     confirmeddict = {}
            #     # confirmeddict['Sample'] = 'Dummy'
            #     confirmeddict['Sample'] = sample
            #     confirmeddict['Chrom'] = chrom
            #     confirmeddict['ContigID'] = querycontig
            #     confirmeddict['ContigLength'] = querylen
            #     confirmeddict['RegionStart'] = ourverybeginning
            #     confirmeddict['RegionEnd'] = ourveryend
            #     confirmeddict['RefSites'] = list(alignment_df['refsite'])
            #     confirmeddict['QuerySites'] = querylist
            #     fullsetlist.append(confirmeddict)
            # fullsetdf = pd.DataFrame(fullsetlist)
            # fullsetdf.to_csv(fullsetfile, sep='\t', index=None)
            filtereddf = pd.DataFrame(filteredlist)
            filtereddf.to_csv(filteredfile, sep='\t', index=None)
        timeend = datetime.datetime.now()
        print(xmap, "TimeTaken: ", timeend-timestart)

        filelist = [xmap, qmap, filteredfile]
        for file in filelist:
            cmd = 'mv ' + file + ' completed/'
            os.system(cmd)


if __name__ == '__main__':
    main()
