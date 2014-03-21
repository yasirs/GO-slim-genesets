import OBOParser

def writeGeneSets(fullOBOFile='data/gene_ontology.1_2.obo',slimOBOFile='data/goslim_generic.obo',geneAssociationFile='data/gene_association.sgd',outfile='gosets.txt'):
    parser = OBOParser.OBOparser()
    ont_full = parser.createOntologyFromOBOFile(fullOBOFile)
    ont_slim = parser.createOntologyFromOBOFile(slimOBOFile)
    all_terms = [x.id for x in ont_full.terms if ((x.getIsObsolete()=='')or(not x.getIsObsolete().upper()[0]=='T'))]
    slim_terms = [x.id for x in ont_slim.terms]

    slim2ancs = {}
    for term in slim_terms:
        slim2ancs[term] = [x.id for x in ont_slim.getAncestors(term)]

    all_slim_ancs = set()
    for term in slim2ancs:
        all_slim_ancs.add(term)
        all_slim_ancs.union(slim2ancs[term])
    all_slim_ancs = list(all_slim_ancs)

    all2slim = {}
    for term in all_terms:
        if (len(all2slim)%100==0):
	    print "done %i terms" %len(all2slim)
        to_look = [term]
	while len(to_look)>0:
	    t = ont_full.getTermById(to_look.pop())
	    if t.id in slim_terms:
	    	break
	    for tt in t.isA:
	    	to_look.insert(0,tt)
        if t.id not in slim_terms:
	    print "Error! %s has no ancestor in slim" %term
	all2slim[term] = t.getId()

	

    term2genes = dict(((go,set()) for go in all_slim_ancs))
    gaf = open(geneAssociationFile,'r')
    badGos = []
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
        if (go in all2slim):
            slim_go = all2slim[go]
        elif (go.upper() in all2slim):
            slim_go = all2slim[go.upper()]
        else:
            #raise StandardError('Annotation term %s not found in Ontology' %go)
            badGos.append(go)
        term2genes[slim_go].add(gene)
        for anc in slim2ancs[slim_go]:
            term2genes[anc].add(gene)
    fo = open(outfile,'w')
    for t,gs in term2genes.iteritems():
    	fo.write(t+'\t'+str(gs))
    fo.close()
    #return term2genes, badGos


if __name__=='__main__':
	writeGeneSets()
        
