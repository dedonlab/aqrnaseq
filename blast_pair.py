## Inputs are blast results from read1 and read2
## Read1 and read2 are merged by integrating their start and end positions to
## generate new start and end positions that reflect their combined coverage

import sys
try:
    in1 = sys.argv[1] #in1=RNA_08-10_blast/180719Ded_D18-6932_1.RNACCAblast
    in2 = sys.argv[2] #in2=RNA_08-10_blast/180719Ded_D18-6932_2.RNACCAblast
    out = sys.argv[3] #out=RNA_08-10_blast/180719Ded_D18-6932_1.CCAblast_paired
except:
    print(__doc__)
    sys.exit(1)

output=['queryId\tsubjectId\tpercIdentity\talnLength\tmismatchCount\tgapOpenCount\tqueryStart\tqueryEnd\tsubjectStart\tsubjectEnd\teVal\tbitScore\tqueryId\tsubjectId\tpercIdentity\talnLength\tmismatchCount\tgapOpenCount\tqueryStart\tqueryEnd\tsubjectStart\tsubjectEnd\teVal\tbitScore']
shendure={}
fin2=open(in2, 'rU').readlines()
row2=0
fin1=open(in1, 'rU').readlines()
row1=0

while (row2<len(fin2)):
	fin2[row2]=fin2[row2].rstrip()
	s2=fin2[row2].split('/')
	seqname=s2[0]
	t2=fin2[row2].split('\t')
	tRNA=t2[1]
	pos=seqname+'\t'+tRNA
	shendure[pos] = fin2[row2]
	row2=row2+1

while(row1<len(fin1)):
	fin1[row1]=fin1[row1].rstrip()
	s1=fin1[row1].split('/')
	seqname=s1[0]
	t1=fin1[row1].split('\t')
	tRNA=t1[1]
	index=seqname+'\t'+tRNA
	if shendure.has_key(index):
		s2=shendure[index].split('\t')
		if (t1[3]>=20 and s2[3]>=20):
			output.append(fin1[row1]+'\t'+shendure[index])
		elif t1[8]!=s2[9]:
			output.append(fin1[row1]+'\t'+shendure[index])
		elif t1[9]!=s2[8]:
			output.append(fin1[row1]+'\t'+shendure[index])
	row1=row1+1;

with open(out, 'w') as fo:
    fo.write('\n'.join(output))
