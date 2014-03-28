import OBOParser
import os

def who():
	print "__file__ = ", __file__
	print "path = ", os.path.split(__file__)

class SlimGOGroups:

    def __init__(self, fullOBOFile=None, slimOBOFile=None):
        self.path = os.path.split(__file__)[0]
        if fullOBOFile==None:
            fullOBOFile = self.path+'/data/gene_ontology.1_2.obo'
        if slimOBOFile==None:
            slimOBOFile = self.path+'/data/goslim_generic.obo'

        parser = OBOParser.OBOparser()
        ont_full = parser.createOntologyFromOBOFile(fullOBOFile)
        ont_slim = parser.createOntologyFromOBOFile(slimOBOFile)
        all_terms = [x.id for x in ont_full.terms.itervalues() if ((x.getIsObsolete()=='')or(not x.getIsObsolete().upper()[0]=='T'))]
        slim_terms = [x.id for x in ont_slim.terms.itervalues()]
    
        self.slim2ancs = {}
        for term in slim_terms:
            self.slim2ancs[term] = [x.id for x in ont_slim.getAncestors(term)]
    
        self.all_slim_ancs = set()
        for term in self.slim2ancs:
            self.all_slim_ancs.add(term)
            self.all_slim_ancs.union(self.slim2ancs[term])
        all_slim_ancs = list(self.all_slim_ancs)
    
        self.all2slim = {}
        for term in all_terms:
            if (len(self.all2slim)%100==0):
                print "done %i terms" %len(self.all2slim)
            to_look = [term]
            while len(to_look)>0:
                t = ont_full.getTermById(to_look.pop())
                if t.id in slim_terms:
                    break
                for tt in t.isA:
                    to_look.insert(0,tt)
            if t.id not in slim_terms:
                print "Error! %s has no ancestor in slim" %term
            self.all2slim[term] = t.getId()
        self.term2genes = dict(((go,set()) for go in self.all_slim_ancs))
        self.badGos = []

    def add_1_annot(self, gene, go):
        if (go in self.all2slim):
            slim_go = self.all2slim[go]
        elif (go.upper() in self.all2slim):
            slim_go = self.all2slim[go.upper()]
        else:
            #raise StandardError('Annotation term %s not found in Ontology' %go)
            self.badGos.append(go)
            raise ValueError("Bad GO %s" %go)
        self.term2genes[slim_go].add(gene)
        for anc in self.slim2ancs[slim_go]:
            self.term2genes[anc].add(gene)

    def read_gene_assoc(self,fname=None):
        if fname==None:
            fname = self.path + '/data/gene_association.sgd'

        gaf = open(fname,'Ur')
        nline = 0
        for line in gaf:
            nline +=1
            if (nline % 100)==0:
                print "done %i genes" %nline
    
            if (line.rstrip()=='')or(line[0]=='!'):
                continue
            w = line.rstrip().split('\t')
            gene = w[10].split('|')[0]
            go = w[4]
            if (go in self.all2slim):
                slim_go = self.all2slim[go]
            elif (go.upper() in self.all2slim):
                slim_go = self.all2slim[go.upper()]
            else:
                #raise StandardError('Annotation term %s not found in Ontology' %go)
                self.badGos.append(go)
            self.term2genes[slim_go].add(gene)
            for anc in self.slim2ancs[slim_go]:
                self.term2genes[anc].add(gene)

    def write_sets(self,outfile='gosets.txt'):
        fo = open(outfile,'w')
        for t,gs in self.term2genes.iteritems():
        	fo.write(t+'\t'+'\t'.join(gs)+'\n')
        fo.close()


if __name__=='__main__':
	
	sg = SlimGOGroups()
	sg.read_gene_assoc()
	sg.write_sets()



        
