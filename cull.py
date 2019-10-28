## Read input file.
## creates a dict with seqID as key and tRNAs and eVal, start, end as subdict(s)
## store lowest eVal for each seqID subdict under 'besteVal' key
import time
import sys
try:
    in1 = sys.argv[1] #in1=RNA_08-10_blast/180719Ded_D18-6932_1.CCAblast_paired
    in2 = sys.argv[2] #in2=D18-6932
except:
    print(__doc__)
    sys.exit(1)

def createdict(filename):
    start_time = time.time()
    fin=open(filename, 'rU').readlines()
    headerRow = fin.pop(0)
    headerRow = headerRow.split('\t')
    masterdict={}
    row=0
    while (row<len(fin)):
        fin[row]=fin[row].split('\t')
        seqID = fin[row][0]
        tRNAName = fin[row][1]
        start1 = min(int(fin[row][8]),int(fin[row][9]))
        end1 = max(int(fin[row][8]),int(fin[row][9]))
        start2 = min(int(fin[row][20]),int(fin[row][21]))
        end2 = max(int(fin[row][20]),int(fin[row][21]))
        eVal = float(fin[row][10])
        bScore = float(fin[row][11])+float(fin[row][23])

        #check if masterdict already has seqID as key
        if masterdict.has_key(seqID):
            #check if current eVal is smaller than besteVal
            if eVal < masterdict[seqID]['besteVal']: #eVal is smaller, discard existing tRNA
                masterdict[seqID]={}
                masterdict[seqID][tRNAName] = [eVal, min(start1,start2), max(end1,end2),bScore]
                masterdict[seqID]['besteVal'] = eVal #update 'besteVal' key
                masterdict[seqID]['bestbScore'] = bScore #update 'bestbScore' key
                masterdict[seqID]['longest'] = max(end1,end2)-min(start1,start2) # update longest length alignment
                row=row+1
            elif eVal == masterdict[seqID]['besteVal']: #eVal is equal,
                if bScore > masterdict[seqID]['bestbScore'] or max(end1,end2)-min(start1,start2) > masterdict[seqID]['longest']:
                #... but if current bScore is > stored OR current length is longer than longest stored, replace
                    masterdict[seqID]={}
                    masterdict[seqID][tRNAName] = [eVal, min(start1,start2), max(end1,end2),bScore]
                    masterdict[seqID]['besteVal'] = eVal #update 'besteVal' key
                    masterdict[seqID]['bestbScore'] = bScore #update 'bestbScore' key
                    masterdict[seqID]['longest'] = max(end1,end2)-min(start1,start2) # update longest length alignment
                    row=row+1
                elif bScore == masterdict[seqID]['bestbScore'] and max(end1,end2)-min(start1,start2)==masterdict[seqID]['longest']:
                #...but if bScore is equal AND max align length is the same, add current tRNA to existing
                    masterdict[seqID][tRNAName] = [eVal, min(start1,start2), max(end1,end2),bScore]
                    row=row+1
                else: # current bscore is < than stored bScore
                    row=row+1
            else: #current eVal is larger than existing, do nothing. go on to next row.
                row=row+1
        else:
            masterdict[seqID]={}
            masterdict[seqID][tRNAName] = [eVal, min(start1,start2), max(end1,end2),bScore]
            masterdict[seqID]['besteVal'] = eVal #create 'besteVal' key
            masterdict[seqID]['bestbScore'] = bScore #create 'bestbScore' key
            masterdict[seqID]['longest'] = max(end1,end2)-min(start1,start2) # update longest length alignment
            row=row+1

    masterdictkeys = masterdict.keys()
    for key in masterdictkeys:
        subdict = masterdict[key]
        if 'besteVal' in subdict:
            del subdict['besteVal']
        if 'bestbScore' in subdict:
            del subdict['bestbScore']
        if 'longest' in subdict:
            del subdict['longest']
    totaltime = time.time() - start_time
    print 'filedict script took', totaltime, 'seconds or', totaltime/float(60),'minutes to run'
    return [masterdict, len(masterdict)]

def dictitems(inputdict):
    start_time = time.time()
    inputdictkeys = inputdict.keys()
    itemcounter = 0
    for seqID in inputdictkeys:
        subdict = inputdict[seqID]
        subdictkeys = subdict.keys()
        for item in subdictkeys:
            if item == 'besteVal' or item in subdictkeys == 'bestbScore' or item in subdictkeys == 'longest':
                subdictkeys.remove(item)
        itemcounter = itemcounter + len(subdictkeys)
    totaltime = time.time() - start_time
    print 'dictitems script took', totaltime, 'seconds or', totaltime/float(60),'minutes to run'
    return itemcounter

suffix = ['.txt','_disc.txt','_nodupe.txt']
fileout=in2+'_culled' #D18-6932_culled
## begin writing files
print '----'
print 'for file:', in1
print '----'
filedict=createdict(in1)
print 'filedict has this many seqIDs:', filedict[1]
filedict_items = dictitems(filedict[0]) #excludes 'besteVal', 'longest','bestbScore' keys
print 'number of total items in filedict is:', filedict_items

# write all lines in filedict to folder 'RNA_11_culledfileout/.txt'
lines1=['SeqID(all)\ttRNA\teVal\tstart\tend\tbitScore']
# only write seqIDs that have been discarded to 'RNA_11_culledfileout/_disc.txt'
lines2=['SeqID(disc)\ttRNA\teVal\tstart\tend\tbitScore']
# write seqIDs with unique tRNA match to 'RNA_11_culledfileout/_nodupe.txt'
lines3=['SeqID(unique)\ttRNA\teVal\tstart\tend\tbitScore']

for seqID in filedict[0].keys():
    tRNAkeys = filedict[0][seqID].keys()
    for tRNAkey in tRNAkeys:
        lines1.append(seqID+'\t'+tRNAkey+'\t'+str(filedict[0][seqID][tRNAkey][0])+'\t'+str(filedict[0][seqID][tRNAkey][1])+'\t'+str(filedict[0][seqID][tRNAkey][2])+'\t'+str(filedict[0][seqID][tRNAkey][3]))
        if len(tRNAkeys) > 1:
            lines2.append(seqID+'\t'+tRNAkey+'\t'+str(filedict[0][seqID][tRNAkey][0])+'\t'+str(filedict[0][seqID][tRNAkey][1])+'\t'+str(filedict[0][seqID][tRNAkey][2])+'\t'+str(filedict[0][seqID][tRNAkey][3]))
        else:
            lines3.append(seqID+'\t'+tRNAkey+'\t'+str(filedict[0][seqID][tRNAkey][0])+'\t'+str(filedict[0][seqID][tRNAkey][1])+'\t'+str(filedict[0][seqID][tRNAkey][2])+'\t'+str(filedict[0][seqID][tRNAkey][3]))
print 'number of lines in total', len(lines1)-1 #minus header line
print 'number of lines in discarded', len(lines2)-1 #minus header line
print 'number of lines in non-duplicated', len(lines3)-1 #minus header line
with open('RNA_11_culledfileout/'+fileout+suffix[0], 'w') as fo:
    fo.write('\n'.join(lines1))
with open('RNA_11_culledfileout/'+fileout+suffix[1], 'w') as fo:
   fo.write('\n'.join(lines2))
with open('RNA_11_culledfileout/'+fileout+suffix[2], 'w') as fo:
   fo.write('\n'.join(lines3))

