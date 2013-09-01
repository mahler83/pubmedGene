# -*- coding: utf8 -*-
# Entrez pubmed search
from Bio import Entrez
from Bio import Medline
import re
from collections import Counter
import csv
import string

# Settings
yearsplit = 2

################ Step 1
# Get input
q = raw_input("Type your search term and press enter\n: ")
retmax = raw_input("How many abstracts? (Default = 200)\n: ")
retmax = int(retmax)
if retmax <= 0:
	retmax = 200

# Start search
Entrez.email = "mahler83@gmail.com"
handle = Entrez.esearch(db="pubmed", term=q, retmax=retmax)
record = Entrez.read(handle)
idlist = record["IdList"]
f = open('PMIDs - ' + q + '.txt', 'w')

for pmid in idlist:
	f.write(pmid + "\n")

f.close()


############### Step 2
# Get file
try:
	f = open('PMIDs - ' + q + '.txt', 'r')
except:
	print 'No such file. Make sure you did the search first!'
	raise SystemExit(0)
f2 = open('Abstracts - ' + q + '.txt', 'w')

i = 0
for pmid in f:
	handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
	records = Medline.parse(handle)
	records = list(records)
	records = records[0]
	try:
		records['PMID']
		records['DA']
		records['TI']
		records['TA']
		records['AB']
	except KeyError:
		continue
	f2.write('PMID: ' + records['PMID'] + "\n")
	f2.write('Date: ' + records['DA'] + "\n")
	f2.write('Title: ' + records['TI'] + "\n")
	f2.write('Journal: ' + records['TA'] + "\n")
	f2.write('Abstract: ' + records['AB'] + "\n")
	f2.write("\n")

	i += 1
	print '\rFetching from Pubmed - Abstract #' + str(i),

f.close()
f2.close()
print



############## Step 3
f = open('protein-coding_gene.txt', 'r')

genelist = []
matching = {}

next(f) #skip headings
reader = csv.reader(f, delimiter='\t')
i = 0
for hid,sym,nam,sta,prs,prn,syn,chr,sccessionnum,refseq in reader:
	#print sym + ': ' + prs + syn
	sym = sym.strip().upper()
	if sym == '':
		continue
	genelist.append(sym)
	prs = prs.strip().upper()
	if prs:
		genelist.append(prs)
		matching[prs] = sym
	syn = syn.strip()
	if syn:
		synList = syn.split(', ')
		for s in synList:
			s = s.strip().upper()
			genelist.append(s)
			matching[s] = sym

try:
	f = open('Abstracts - ' + q + '.txt', 'r')
except:
	print 'No such file. Make sure you did the abstract retrieval first!'
	raise SystemExit(0)
print 'Search term = ' + q
f2 = open('Wordcount - ' + q + '.txt', 'w')
f3 = open('Monthly wordcount - ' + q + '.txt', 'w')

def myrepl(matchobj):
	return matchobj.group(0).upper()

countPerMonth = {}
countTotal = {}
months = []

stopwords = ['THE', 'AND', 'OF', 'IN', 'TO', 'WITH', 'WERE', 'WAS', 'FOR', 'IN', 'BY', 'ON', 'AS', 'THAT', 'THIS', 'OR', 'BE', 'ARE', 'AN', 'WE', 'AT', 'NOT', 'THAN', 'HAVE', 'IS', 'FROM', 'HAD', 'NO', 'THERE', 'BUT', 'THESE', 'MORE', 'FROM', 'IS', 'MAY', 'ALL', 'HAS', 'ALSO', 'BOTH', 'BEEN', 'CAN', 'IT', 'MOST', 'USE', 'OTHER', 'ITS', 'THOSE', 'VS', 'OUR', 'INTO', 'COULD', 'WHEN', 'WILL', 'MIGHT', 'UP', 'DOWN', 'OFTEN', 'IF', 'HOW', 'THUS', 'NEW', 'AFTER']
stopwords2 = ['MICE']

switch = 0
for line in f:
	#print '{0}\r'.format(Searching),
	print '\rSearching Abstract #' + str(switch+1),
	if line.startswith('Date: '):
		pubdate = line[6:]
		#pubYM = pubdate[0:6]
		if yearsplit == 4:
			if int(pubdate[4:6]) > 9:
				pubYM = pubdate[0:4] + '.10-12'
			elif int(pubdate[4:6]) > 6:
				pubYM = pubdate[0:4] + '.07-09'
			elif int(pubdate[4:6]) > 3:
				pubYM = pubdate[0:4] + '.04-06'
			else:
				pubYM = pubdate[0:4] + '.01-03'
		elif yearsplit == 2:
			if int(pubdate[4:6]) > 6:
				pubYM = pubdate[0:4] + '.07-12'
			else:
				pubYM = pubdate[0:4] + '.01-06'
		else:
			pubYM = pubdate[0:4]
		if pubYM not in months:
			months.append(pubYM)
		continue
	elif line.startswith('Abstract: '):
		abstract = line[9:]
		txt = re.sub('[,():";\[\]=]', ' ', abstract)
		txt = re.sub('  ', ' ', txt)
		txt = re.sub('  ', ' ', txt)
		#txt = re.sub('([A-Z])([a-z]+)\W', myrepl, txt)
		txt = txt.upper()
		#txt = re.sub('(\w+)\.', r'\1', txt)
		#t = Counter(txt.split()).most_common()
		inThisAbstract = []
		for word in txt.split():
			if len(word) < 3:
				continue
			#word = re.sub('^([A-Z])([a-z]+)$', myrepl, word)
			if word in stopwords:
				continue
			if word in stopwords2:
				continue
			if word in inThisAbstract:
				continue
			else:
				inThisAbstract.append(word)
			try:
				word = matching[word]
			except:
				word
			word = re.sub('^(\w+)\.$', r'\1', word)
			try:
				countPerMonth[word][pubYM] += 1
			except KeyError:
				try:
					countPerMonth[word][pubYM] = 1
				except KeyError:
					countPerMonth[word] = {}
					countPerMonth[word][pubYM] = 1
			try:
				countTotal[word] += 1
			except KeyError:
				countTotal[word] = 1

		switch += 1
	else:
		continue

print

f2.write("Gene\tOccurance\n")
f3.write("Month")
for pubYM in months:
	f3.write("\t" + pubYM)
f3.write("\n")

i = 0
for w in sorted(countTotal, key=countTotal.get, reverse=True):
	print '\rVerifying word #' + str(i+1),
	#print w, countTotal[w]
	if w not in genelist:
		continue
	f2.write(w + "\t" + str(countTotal[w]) + "\n")
	f3.write(w)
	for pubYM in months:
		try:
			f3.write("\t" + str(countPerMonth[w][pubYM]))
		except:
			f3.write("\t0")
	f3.write("\n")

	i += 1
	if i >= 1000:
		break
print

raw_input('Press Enter to end.\n: ')
