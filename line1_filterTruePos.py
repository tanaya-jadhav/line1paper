"""
Script for Line1 Step 3 

Line1 pipeline:
1. Extract inserted nicking sites from xmaps and qmaps for each region
2. Overlap BSPQ and DLE results to calculate distances 
3. Filter true positives in tsvs from step 2 (this script)
4. Get list of unique regions (this)

True positives conditions: 
DLE - BSPQ - DLE pattern 
D1 = 2500 +- 500
D2 = 1000 +- 500 (or oppposite)
"""
import pandas as pd
import os
import sys

def filter_truepos(input_tsv, in_name):
    headercols = ['Sample', 'Chrom', 'DLE_ContigID', 'BSPQ_ContigID', 'DLE_Region', 'BSPQ_Region', 'nick1pos', 'D1', 'nick2pos', 'D2', 'nick3pos', 'D3', 'nick4pos']    
    df = pd.read_csv(input_tsv, sep='\t', header=None, names=headercols)
    
    truepos = []
    for index, row in df.iterrows():
        nick1label = row['nick1pos'].split('\'')[1]
        nick2label =row['nick2pos'].split('\'')[1]
        nick3label = row['nick3pos'].split('\'')[1]

        ## if nick4pos is not null, extract nick4label
        if pd.notna(row['nick4pos']):
            nick4label = row['nick4pos'].split('\'')[1]
        else: nick4label = 'None'

        ## check for dle, bspq, dle pattern in nick1pos-nick3pos
        if nick1label == 'dle' and nick2label == 'bspq' and nick3label == 'dle':
            if (row['D1'] >= 2000 and row['D1'] <= 3000) and (row['D2'] >= 500 and row['D2'] <= 1500):
                truepos.append(row)
            if (row['D2'] >= 2000 and row['D2'] <= 3000) and (row['D1'] >= 500 and row['D1'] <= 1500):
                truepos.append(row)
        
        ## check for dle, bspq, dle pattern in nick2pos-nick4pos
        if nick2label == 'dle' and nick3label == 'bspq' and nick4label == 'dle':
            if (row['D2'] >= 2000 and row['D2'] <= 3000) and (row['D3'] >= 500 and row['D3'] <= 1500):
                truepos.append(row)
            if (row['D3'] >= 2000 and row['D3'] <= 3000) and (row['D2'] >= 500 and row['D2'] <= 1500):
                truepos.append(row)

    
    truepos_df = pd.DataFrame(truepos,columns=['Sample', 'Chrom', 'DLE_ContigID', 'BSPQ_ContigID', 'DLE_Region', 'BSPQ_Region', 'nick1pos', 'D1', 'nick2pos', 'D2', 'nick3pos', 'D3', 'nick4pos'])
    truepos_df.sort_values(by=['Chrom', 'DLE_Region', 'Sample'], ascending=True)
    truepos_df.to_csv('{}_truepositives.tsv'.format(in_name), sep='\t', index=None)

    return truepos_df


def main(input_tsv, in_name):
    ## open true positives results file and get DLE regions
    indf = filter_truepos(input_tsv, in_name)

    ## true positives file
    regiondf = indf[['Chrom', 'DLE_Region']]
    regiondf[['RegionStart', 'RegionEnd']] = regiondf['DLE_Region'].str.split('-', 1, expand=True).astype(int)
    regiondf = regiondf[['Chrom', 'RegionStart', 'RegionEnd', 'DLE_Region']]

    regiondf['Chrom'] = 'chr' + regiondf['Chrom'].astype(str)
    regiondf.loc[(regiondf['Chrom'] == 'chr23'),'Chrom']='chrX'
    regiondf.loc[(regiondf['Chrom'] == 'chr24'),'Chrom']='chrY'
    regiondf.sort_values(['Chrom', 'RegionStart', 'RegionEnd'], inplace=True, ascending=True)
    print(regiondf.tail())

    # export regions to bed file
    regiondf.to_csv('temp.bed', sep='\t', index=None)

    ## use bedtools to merge regions for unique list of regions 
    uniq_regions_out = '{}_UniqRegions.bed'.format(in_name)
    merge_cmd = 'bedtools merge -i temp.bed > {}'.format(uniq_regions_out)
    os.system(merge_cmd)
    rm_cmd = 'rm temp.bed'
    os.system(rm_cmd)

    ## open the unique regions and convert to list for further processing
    headercols = ['Chrom', 'start', 'end']
    uniq_regiondf = pd.read_csv(uniq_regions_out, sep='\t', header=None, names=headercols)
    uniqregionslist = uniq_regiondf.values.tolist()
    print(uniqregionslist[:5])
    print(len(uniqregionslist))


if __name__ == "__main__":
    # input_tsv = '566alldistances_noheaders.tsv'
    input_tsv = sys.argv[1]
    in_name = input_tsv.split('_')[0]

    main(input_tsv, in_name)