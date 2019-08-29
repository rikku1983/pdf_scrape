#! /Library/Frameworks/Python.framework/Versions/3.7/bin/python3

from pdfminer.pdfparser import PDFParser
from pdfminer.pdfdocument import PDFDocument
from pdfminer.layout import LAParams
from pdf2data.pdf import PageIterator
import re, sys, os, subprocess

# cmd: ./task.py NIHMS69173-supplement-Supplementary_Appendix.pdf outdir
# Take file as argument
inpdf = sys.argv[1]
outdir = sys.argv[2]

if not os.path.isfile(inpdf):
    print('PDF file not exist.')
    sys.exit()

subprocess.call(["mkdir", "-p", outdir])

# clustrows tabke sorted list of y positions of all text entities, 
# group them into groups based on differences in values
def clustrows(y0s):
	r = [[y0s[0]]]
	for i in range(1,len(y0s)):
		diffs = [abs(y0s[i] - sum(x)/len(x)) for x in r]
		minidx = diffs.index(min(diffs))
		if min(diffs) > 2:
			r = r + [[y0s[i]]]
		else:
			r[minidx] = r[minidx] + [y0s[i]]
	return(r)

# parsetable function take 1 page of pdf in the format of pdf2data.pdf.PageIterator and
# return a list of rows of the table contents, each row is text entities separated by tab
def parsetable(page, header1=''):
	tls = page.lines
	if header1 != '':
		hd1tl = [x for x in tls if x.text==header1][0]
		headerline = [x for x in tls if abs(x.y0 - hd1tl.y0) < 2]
		headerline = sorted(headerline, key=lambda h: h.x0)
		headers = '\t'.join([x.text for x in headerline])
		tbrows = [x for x in tls if x.y1 < hd1tl.y0]
	else:
		tbrows = [x for x in tls]
	ry0s = sorted([x.y0 for x in tbrows], reverse=True)
	rows = clustrows(ry0s)
	tbrows2 = [sorted([x for x in tbrows if x.y0 in ri], key=lambda x: x.x0) for ri in rows]
	tbrows3 = ['\t'.join([x.text for x in r]) for r in tbrows2]
	if header1 != '':
		return([headers] + tbrows3)
	else:
		return(tbrows3)


infile = open(inpdf, 'rb')
document = PDFDocument(PDFParser(infile))
page_it = PageIterator(document, LAParams(char_margin=0.2))

tbs5=[]
tbs6=[]
tbs7=[]
tbs8=[]
pg = 1
while pg < 264:
	if pg == 153:
		tbs5 = tbs5 + parsetable(page_it, header1='SAMPLE')
	elif pg > 153 and pg < 224:
		tbs5 = tbs5 + parsetable(page_it)
	elif pg == 224:
		tbs6 = tbs6 + parsetable(page_it, header1='PDID')
	elif pg > 224 and pg < 248:
		tbs6 = tbs6 + parsetable(page_it)
	elif pg == 248:
		tbs7 = tbs7 + parsetable(page_it, header1='Pairs')
	elif pg > 248 and pg < 257:
		tbs7 = tbs7 + parsetable(page_it)
	elif pg == 257:
		tbs8 = tbs8 + parsetable(page_it, header1='Triplet')
	elif pg > 257:
		tbs8 = tbs8 + parsetable(page_it)
	page_it.advance()
	pg = pg + 1

infile.close()

############### Following section Cleans fusion columns, this is based on manual checking
## table s7 and s8 are clean
## table s6
tbs6b = [x.split('\t') for x in tbs6]
# column fusion in headers
tbs6b[0] = tbs6b[0][0:7] + \
	[tbs6b[0][7][:7]] + [tbs6b[0][7][7:24]] + [tbs6b[0][7][24:30]] + [tbs6b[0][7][30:]] + \
	[tbs6b[0][8][:11]] + [tbs6b[0][8][11:]] + tbs6b[0][9:]
## padding last two summary columns
tbs6b[-2] = [''] + tbs6b[-2]
tbs6b[-1] = [''] + tbs6b[-1] + ['']*22
tbs6c = ['\t'.join(x) for x in tbs6b]

## table s5 has many fused columns
tbs5b = [x.split('\t') for x in tbs5]
# header missing one
tbs5b[0] = tbs5b[0][:8] + ['gene_variant'] + [tbs5b[0][8]]
# check
[x for x in tbs5b[1:] if not re.match('^PD[0-9]{4,5}a$',x[0])] # column 1 no fusion
[x for x in tbs5b[1:] if not x[1] in ['Sub','I','D','ID','PTD','ITD']] # column 2 no fusion
[x for x in tbs5b[1:] if not x[2] in ([str(c) for c in range(1,23)] + ['X', 'Y'])] # column 3 no fusion
[x for x in tbs5b[1:] if (not re.match('^[0-9,]+$',x[3])) and (x[3] != 'na')] # column 4 no fusion
[x[4] for x in tbs5b[1:] if re.search('[0-9]$',x[4])] # no case where col 5, 6, and 7 fuse
[x[4] for x in tbs5b[1:] if re.match('[0-9\.]+$',x[5])] # 8 cases col 5 and 6 fused, e.g, tcgctaW
# Correct fused column 5 and 6
print('Correcting fused column 5 and 6 for following rows...')
for i in range(1,len(tbs5b)):
	r = tbs5b[i]
	if re.match('[0-9\.]+$',r[5]):
		print(i)
		tbs5b[i] = r[:4] + [r[4][:-1]] + [r[4][-1]] + r[5:]
[x[5] for x in tbs5b[1:] if re.search('[0-9\.]$',x[5])] # 6 cases col 6 and 7 fused, e.g, GAAAACTG36.4
# Correct column 6
print('Correcting fused column 6 and 7 for following rows...')
for i in range(1,len(tbs5b)):
	r = tbs5b[i]
	if re.search('[0-9\.]$',r[5]):
		print(i)
		tbs5b[i] = r[:5] + [re.sub('^([ATCG]+)([0-9\.]+)$','\g<1>',r[5])] + [re.sub('^([ATCG]+)([0-9\.]+)$','\g<2>',r[5])] + r[6:]
[x for x in tbs5b[1:] if (not re.match('^[0-9\.]+$',x[6])) and (x[6] != 'na')] # column 7 no fusion now
[x for x in tbs5b[1:] if (not re.match('^[0-9\.]+$',x[7])) and (x[7] != 'na')] # column 8 no fusion
[x for x in tbs5b[1:] if '_' in x[-1]] # column 9 and column 10 fused
# Correct column 9
print('Correcting fused column 9 and 10 for following rows...')
for i in range(1,len(tbs5b)):
	r = tbs5b[i]
	if r[-1] not in ['ONCOGENIC', 'POSSIBLE']:
		print(i)
		tbs5b[i] = r[:-1] + [re.sub('(ONCOGENIC|POSSIBLE)','',r[-1])] + [re.sub('.*((ONCOGENIC|POSSIBLE))','\g<1>',r[-1])]
set([len(x) for x in tbs5b]) # correction done!
print('Correction is done!')
################################ Clean done

## separate column 9 to be two columns
## 93 cases have two '_' and it seems that only the first part is gene name
tbs5c = [x[:8] + [x[8].split('_')[0]] + ['_'.join(x[8].split('_')[1:])] + [x[9]] for x in tbs5b]
tbs5d = ['\t'.join(x) for x in tbs5c]

############ Calculation of frequencies
## combinecol is takeing idx as a list of index to subset list l and join using sep
def combinecol(l, idx, sep='::'):
	newl = sep.join([l[i] for i in idx])
	return(newl)
## getfreq function count unique combinations of vars1 for each unique combinations of vars2
## It takes a list of rows as data input, and it assumes the first row is col names
## vars1 and vars2 need to be type of list
def getfreq(ll, vars1, vars2):
	vars1idx = [ind for ind, x in enumerate(ll[0]) if x in vars1]
	vars2idx = [ind for ind, x in enumerate(ll[0]) if x in vars2]
	gs = list(set([combinecol(x, vars2idx) for x in ll[1:]]))
	gs = sorted(gs)
	gs = ['+'.join(vars2)] + gs
	freq = ['+'.join(vars1)]
	for g in gs[1:]:
		freq.append(len(set([combinecol(x, vars1idx) for x in ll if combinecol(x, vars2idx)==g])))
	gf = [gs[i] + '\t' + str(freq[i]) for i in range(len(gs))]
	return(gf)

# Indication level freq: number of unique gene + variant for each sample
freq1 = getfreq(tbs5c, ['gene', 'variant'], ['SAMPLE'])
# gene level freq: number of unique gene for each sample
freq2 = getfreq(tbs5c, ['gene'], ['SAMPLE'])
# variant level freq: number of unique variant for each sample + gene
freq3 = getfreq(tbs5c, ['variant'], ['SAMPLE','gene'])


####### Write all results to tab delimited txt file
## writefile function take list of rows, each row is '\t' separated, write to file f
def writefile(lr, f):
	with open(outdir + '/' + f, 'w') as of:
		for r in lr:
			of.write(r)
			of.write('\n')

## Table s5
writefile(tbs5d, 'table_s5.txt')
writefile(tbs6c, 'table_s6.txt')
writefile(tbs7, 'table_s7.txt')
writefile(tbs8, 'table_s8.txt')
writefile(freq1, 'indication_freq.txt')
writefile(freq2, 'gene_freq.txt')
writefile(freq3, 'variant_freq.txt')


# End





