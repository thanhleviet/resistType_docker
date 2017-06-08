#!/usr/bin/python
from Bio import SeqIO
import urllib
import re
import HTMLParser
from urlparse import urlparse
import StringIO
import sys

class MLSTFetch:

	def __init__(self):
		self.genes={}
		self.profiles=[]
		self.mlstpost=None

	def run(self):
		pass

	def extractFasta(self,txt):
	    return txt

	def extractTab(self,txt):
	    return txt

	def readAlleles(self,regex):
		for gene in self.geneSet:
			fasta=StringIO.StringIO(self.extractFasta(urllib.urlopen(self.alleleUrl.format(gene)).read()))
			seqs=[i for i in SeqIO.parse(fasta,"fasta")]
			for i in seqs:
				gname,n=re.findall(regex,i.name)[0]
				self.genes.setdefault(gname,{})[int(n)]=str(i.seq)
    
	def getOrder(self,firstline):
		pass

	def checkInfo(self):
		if not self.genes or not self.profiles: print self.code,"Fetch failed"
		else: print self.code,[(i,len(j)) for i,j in self.genes.items()],len(self.profiles)

	def readProfile(self,headersize=0):
		if not self.mlstpost:
			profile=StringIO.StringIO(self.extractTab(urllib.urlopen(self.mlstProfileUrl).read()))

		else:
			profile=StringIO.StringIO(self.extractTab(urllib.urlopen(self.mlstProfileUrl,self.mlstpost).read()))

		profile=[i.strip().split() for i in profile]
		order=self.getOrder(profile[0])
		res_profile=[]

		for i in profile[headersize:]:
			if len(i)<len(self.geneSet)+1: continue
			st=i[0]
			alleles=dict(zip(order,[int(j) for j in i[1:1+len(self.genes)]]))
			badst=False
			for gn,n in alleles.items():
				if gn not in self.genes or n not in self.genes[gn]:
					badst=True
					break
				
			if badst: 
				print "bad st",st,alleles
				continue
			self.profiles.append([st,alleles])
	def store(self):
		fn1 = self.code + "-alleles.fa"
		fn2 = self.code + "-stfile.tsv"
		f = open(fn1, "w")
		for gene,alleles in self.genes.items():
			for n,seq in alleles.items():
				alcode="{0}-{1}-{2}".format(self.code,gene,n)
				f.write(">{0}\n{1}\n".format(alcode, seq))
		f.close()
		f =open(fn2, "w")
		f.write("stid\ta1\ta2\ta3\ta4\ta5\ta6\ta7\n")
		for st,alleles in self.profiles:
			stcode="{0}-ST-{1}".format(self.code,st)
			alcodes=["{0}-{1}-{2}".format(self.code,i,j) for i,j in alleles.items()]
			f.write("{0}\t{1}\n".format(stcode, "\t".join(alcodes)))
		f.close()
				
		return
	

class PubMLST(MLSTFetch):
	def __init__(self):
		MLSTFetch.__init__(self)
		up=urlparse(self.mlstProfileUrl)
		self.mlstProfileUrl="{0}://{1}{2}".format(up.scheme,up.netloc,up.path)
		self.mlstpost=up.query

	def getOrder(self,firstline):
		return firstline[1:1+len(self.genes)]

	def run(self):
		self.readAlleles("^([^_]+)_([0-9]+)$")
		self.readProfile(1)
		self.checkInfo()

class MLSTNet(MLSTFetch):
	def __init__(self):
		MLSTFetch.__init__(self)

	def getOrder(self,firstline):
		return self.geneSet

	def extractFasta(self,txt):
		c=[i for i in re.finditer("</?textarea[^>]*>",txt)]
		txt= txt[c[0].end():c[1].start()]
		h = HTMLParser.HTMLParser()
		return h.unescape(txt)

	def extractTab(self,txt):
		txt=txt.replace("\n"," ")
		txt=txt.replace("\r"," ")
		txt=txt.replace("<br>","\n")
		txt=re.sub(" +"," ",txt)
		txt=re.sub("^ "," ",txt)
		return txt

	def run(self):
		self.readAlleles("^([^0-9]+)([0-9]+)$")
		self.readProfile()
		self.checkInfo()

class MLSTWarwick(MLSTFetch):
	def __init__(self):
		MLSTFetch.__init__(self)
		self.dropColumns.sort(reverse=True)

	def extractTab(self,txt):
		txt2=[]
		txt=txt.split("\n")
		for i in txt:
			if not i: continue
			i=i.split("\t")
			i=[re.sub("NUM$","",j) for j in i]
			[i.pop(j) for j in self.dropColumns]

			txt2.append(i)
		return "\n".join(["\t".join(i) for i in txt2])

	def run(self):
		self.readAlleles("^([^0-9]+)([0-9]+)$")
		self.readProfile(1)
		self.checkInfo()

	def getOrder(self,firstline):
		return firstline[1:1+len(self.genes)]

class CdifMLST(PubMLST):
	def __init__(self):
		self.code="Cdif"
		self.name="Clostridium difficile"
		self.geneSet=["adk","atpA","dxr","glyA","recA","sodA","tpi"]
		self.alleleUrl="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_cdifficile_seqdef&page=downloadAlleles&locus={0}"
		self.mlstProfileUrl="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_cdifficile_seqdef&page=downloadProfiles&scheme_id=1"
		self.mlstProfileUrl="https://pubmlst.org/data/profiles/cdifficile.txt"
		PubMLST.__init__(self)
		

class SsuiMLST(PubMLST):
	def __init__(self):
		self.code="Ssui"
		self.name="Streptococcus suis"
		self.geneSet=["aroA","cpn60","dpr","gki","mutS","recA","thrA"]
		self.alleleUrl="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_ssuis_seqdef&page=downloadAlleles&locus={0}"
		self.mlstProfileUrl="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_ssuis_seqdef&page=downloadProfiles&scheme_id=1"
		self.mlstProfileUrl="https://pubmlst.org/data/profiles/ssuis.txt"
		PubMLST.__init__(self)

		
class SaurMLST(MLSTNet):
	def __init__(self):
		self.code="Saur"
		self.name="Staphylococcus aureus"
		self.geneSet=["arcc","aroe","glpf","gmk_","pta_","tpi_","yqil"]
		self.alleleUrl="http://saureus.mlst.net/sql/fasta.asp?allele={0}"
		self.mlstProfileUrl="http://saureus.mlst.net/sql/st_tab.asp"
		MLSTNet.__init__(self)

class SepiMLST(MLSTNet):
	def __init__(self):
		self.code="Sepi"
		self.name="Staphylococcus epidermis"
		self.geneSet=["arcC","aroE","gtr","mutS","pyr","tpi","yqiL"]
		self.alleleUrl="http://sepidermidis.mlst.net/sql/fasta.asp?allele={0}"
		self.mlstProfileUrl="http://sepidermidis.mlst.net/sql/st_tab.asp"
		MLSTNet.__init__(self)

class SpneMLST(PubMLST):
	def __init__(self):
		self.code="Spne"
		self.name="Streptococcus pneumoniae"
		self.geneSet=["aroE","gdh","gki","recP","spi","xpt","ddl"]
		self.alleleUrl="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_spneumoniae_seqdef&page=downloadAlleles&locus={0}"
		self.mlstProfileUrl="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_spneumoniae_seqdef&page=downloadProfiles&scheme_id=1"
		self.mlstProfileUrl="https://pubmlst.org/data/profiles/spneumoniae.txt"
		PubMLST.__init__(self)
	
class KpneMLST(PubMLST):
	def __init__(self):
		self.code="Kpne"
		self.name="Klebsiella pneumoniae"
		self.geneSet=["gapA","infB","mdh","pgi","phoE","rpoB","tonB"]
		self.alleleUrl="http://bigsdb.web.pasteur.fr/perl/bigsdb/bigsdb.pl?db=pubmlst_klebsiella_seqdef_public&page=downloadAlleles&locus={0}"
		self.mlstProfileUrl="http://bigsdb.pasteur.fr/perl/bigsdb/bigsdb.pl?db=pubmlst_klebsiella_seqdef_public&page=downloadProfiles&scheme_id=1"

		
		PubMLST.__init__(self)

class KoxyMLST(PubMLST):
	def __init__(self):
		self.code="Koxy"
		self.name="Klebsiella oxytoca"
		self.geneSet=["gapA","infB","mdh","pgi","phoE","rpoB","tonB"]
		self.alleleUrl="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_koxytoca_seqdef&page=downloadAlleles&locus={0}"
		self.mlstProfileUrl="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_koxytoca_seqdef&page=downloadProfiles&scheme_id=1"
		self.mlstProfileUrl="https://pubmlst.org/data/profiles/koxytoca.txt"
		PubMLST.__init__(self)

class EcloMLST(PubMLST):
	def __init__(self):
		self.code="Eclo"
		self.name="Enterobacter cloacae"
		self.geneSet=["dnaA","fusA","gyrB","leuS","pyrG","rplB","rpoB"]
		self.alleleUrl="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_ecloacae_seqdef&page=downloadAlleles&locus={0}"
		self.mlstProfileUrl="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_ecloacae_seqdef&page=downloadProfiles&scheme_id=1"
		self.mlstProfileUrl="https://pubmlst.org/data/profiles/ecloacae.txt"
		PubMLST.__init__(self)

class PaerMLST(PubMLST):
	def __init__(self):
		self.code="Paer"
		self.name="Pseudomonas aeruginosa"
		self.geneSet=["acs","aro","gua","mut","nuo","pps","trp"]
		self.alleleUrl="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_paeruginosa_seqdef&page=downloadAlleles&locus={0}"
		self.mlstProfileUrl="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_paeruginosa_seqdef&page=downloadProfiles&scheme_id=1"

		self.mlstProfileUrl="https://pubmlst.org/data/profiles/paeruginosa.txt"
		PubMLST.__init__(self)

class SentMLST(MLSTWarwick):
	def __init__(self):
		self.code="Sent"
		self.name="Salmonella enterica"
		self.geneSet=["AROC","DNAN","HEMD","HISD","PURE","SUCA","THRA"]
		self.alleleUrl="http://mlst.warwick.ac.uk/mlst/dbs/Senterica/handlers/getFileData/home/cbailster/mlst/zope/Extensions/gadfly/Senterica/DB/{0}.fas"
		self.mlstProfileUrl="http://mlst.warwick.ac.uk/mlst/dbs/Senterica/handlers/getFileData/home/cbailster/mlst/zope/Extensions/gadfly/Senterica/DB/publicSTs.txt"
		self.dropColumns=[1]
		MLSTWarwick.__init__(self)

class EcolMLST(MLSTWarwick):
	def __init__(self):
		self.code="Ecol"
		self.name="Escherichia coli"
		self.geneSet=["ADK","FUMC","GYRB","ICD","MDH","PURA","RECA"]
		self.alleleUrl="http://mlst.warwick.ac.uk/mlst/dbs/Ecoli/handlers/getFileData/home/cbailster/mlst/zope/Extensions/gadfly/Ecoli/DB/{0}.fas"
		self.mlstProfileUrl="http://mlst.warwick.ac.uk/mlst/dbs/Ecoli/handlers/getFileData/home/cbailster/mlst/zope/Extensions/gadfly/Ecoli/DB/publicSTs.txt"
		self.dropColumns=[1,2]
		MLSTWarwick.__init__(self)

class NgonMLST(PubMLST):
	def __init__(self):
		self.code="Ngon"
		self.name="Neisseria gonorrhoeae"
		self.geneSet=["abcZ","adk","aroE","fumC","gdh","pdhC","pgm"]
		self.alleleUrl="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_neisseria_seqdef&page=downloadAlleles&locus={0}"
		self.mlstProfileUrl="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_neisseria_seqdef&page=downloadProfiles&scheme_id=1"
		self.mlstProfileUrl="https://pubmlst.org/data/profiles/neisseria.txt"
		PubMLST.__init__(self)
class CfreMLST(PubMLST):
	def __init__(self):
		self.code="Cfre"
		self.name="Citrobacter freundii"
		self.geneSet=["aspC","clpX","fadD","mdh","arcA","dnaG","lysP"]
		self.alleleUrl="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_cfreundii_seqdef&page=downloadAlleles&locus={0}"
		self.mlstProfileUrl="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_cfreundii_seqdef&page=downloadProfiles&scheme_id=1"
		self.mlstProfileUrl="https://pubmlst.org/data/profiles/cfreundii.txt"
		PubMLST.__init__(self)
class CjejMLST(PubMLST):
	def __init__(self):
		self.code="Cjej"
		self.name="Campylobacter jejuni"
		self.geneSet=["aspA","glnA","gltA","glyA","pgm","tkt","uncA"]

		self.alleleUrl="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_campylobacter_seqdef&page=downloadAlleles&locus={0}"
		self.mlstProfileUrl="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_campylobacter_seqdef&page=downloadProfiles&scheme_id=1"
		self.mlstProfileUrl="https://pubmlst.org/data/profiles/campylobacter.txt"
		PubMLST.__init__(self)
a=SsuiMLST()
a.run()
a.store()

a=CdifMLST()
a.run()
a.store()


a=KpneMLST()
a.run()
a.store()
a=CjejMLST()
a.run()
a.store()
a=KoxyMLST()
a.run()
a.store()
a=SepiMLST()
a.run()
a.store()
a=CfreMLST()
a.run()
a.store()
a=SaurMLST()
a.run()
a.store()
a=SpneMLST()
a.run()
a.store()
a=EcloMLST()
a.run()
a.store()
a=PaerMLST()
a.run()
a.store()
a=SentMLST()
a.run()
a.store()
a=NgonMLST()
a.run()
a.store()
a=EcolMLST()
a.run()
a.store()
