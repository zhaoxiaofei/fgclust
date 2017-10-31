from __future__ import division
import sys
from collections import defaultdict
import xml.etree.ElementTree as ET
from math import log

EXPCODES = set(['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP'])

def init_term_to_reachableterms(term, term_to_nxtterms, term_to_reachableterms):
    if term not in term_to_reachableterms:
        reachableterms = set([term])
        for subterm in term_to_nxtterms[term]:
            reachableterms.update(init_term_to_reachableterms(subterm, term_to_nxtterms, term_to_reachableterms))
        term_to_reachableterms[term] = reachableterms
    return term_to_reachableterms[term]

alt_to_id = {}
term_to_supterms = defaultdict(list)
term_to_subterms = defaultdict(list)
term_to_namespace = []

root = ET.parse(sys.argv[1]).getroot()
for term in root.findall('term'):
    child = term.find('id')
    alts = term.findall('alt_id')
    for alt in alts: alt_to_id[alt.text] = child.text
    alt_to_id[child.text] = child.text
    if term.find('is_obsolete') is not None and term.find('is_obsolete').text == '1': alt_to_id[child.text] = None
    namespace = term.find('namespace')
    term_to_namespace[term.text] = namespace.text
    for parent in term.findall('is_a'):
         term_to_supterms[child.text].append(parent.text)
         term_to_subterms[parent.text].append(child.text)

term_to_offspring = {}
term_to_ancestors = {}
for term in term_to_namespace:
    term_to_offspring[term] = init_term_to_reachableterms(term, term_to_subterms, term_to_offspring)
    term_to_ancestors[term] = init_term_to_reachableterms(term, term_to_supterms, term_to_ancestors)

def gos_to_ancestor_specificity(gos):
    for go in gos:
        if None == go: return (None, 0)
    gos = [alt_to_id[go] for go in gos]
    ancestors = term_to_ancestors[gos[0]]
    for go in gos[1:]:
        ancestors = ancestors & term_to_ancestors[go]
    most_specific_ancestor = None
    highest_specificity = 0
    for ancestor in ancestors:
        if term_to_specificity[ancestor] > highest_specificity:
            most_specific_ancestor = ancestor
            highest_specificity = term_to_specificity[ancestor]
    return (most_specific_ancestor, highest_specificity)

acc_to_gos = {}
with open(sys.argv[2]) as sprotfile:
    acc = ''
    gos = set([])
    i = 0
    for line in sprotfile:
        tokens = line.strip().split()
        if len(tokens) > 1 and tokens[0] == 'AC':
            if acc != '': 
                if (i & (i-1)) == 0 : sys.stderr.write('acc, gos = {}, {}\n'.format(acc, gos))
                acc_to_gos[acc] = gos
                i+=1
            acc = tokens[1].strip(';')
            gos = set([])
        if len(tokens) > 2 and tokens[0] == 'DR' and tokens[1] == 'GO;':
            if (
                #tokens[3].startswith('F:') and 
                alt_to_id[tokens[2].strip(';')] and tokens[-1].split(':')[0] 
                #!= ''
                #!= 'IEA'
                in EXPCODES
                ):
                gos.add(alt_to_id[tokens[2].strip(';')])
    if acc != '': acc_to_gos[acc] = gos
rep_to_goslist = defaultdict(list)
rep_to_clusize = defaultdict(int)
rep_to_gos = {}
nseqs = 0
clusters = set([])
with open(sys.argv[3]) as clufile:
    acc = ''
    for i, line in enumerate(clufile):
        tokens = line.strip().split()
        inner = tokens[0].lstrip('UniRef100_')
        outer = tokens[1].lstrip('UniRef100_')
        clusters.add(inner)
        rep_to_clusize[inner] += 1
        nseqs += 1
        if outer in acc_to_gos and acc_to_gos[outer]: 
            #sys.stderr.write('The following accession has go info {}, {}\n'.format(outer, acc_to_gos[outer]))
            rep_to_goslist[inner].append(acc_to_gos[outer])
            if inner == outer:
                rep_to_gos[inner] = acc_to_gos[inner]
        else:
            if (i & (i-1)) == 0: sys.stderr.write('acc {} has no go info\n'.format(outer)) 
            #sys.stderr.write('The following accession has no go info : {}\n'.format(outer))

go_to_freq, n_gos = (defaultdict(int), 0)
for rep, goslist in rep_to_goslist.items():
    for gos in goslist:
        for go in gos:
            for anc in term_to_ancestors[go]:
                go_to_freq[anc] += 1
            n_gos += 1
go_to_freq[None] = n_gos

term_to_specificity = {}
for term in term_to_namespace:
    n_offsprings = len(term_to_offspring[term])
    n_ancestors = len(term_to_ancestors[term])
    if n_offsprings + n_ancestors == 2:
        specificity = 0
    else:
        specificity = 1.0 - float(n_offsprings - 1) / float(n_offsprings + n_ancestors - 2)
    term_to_specificity[term] = specificity


def sim11(go1, go2):
    p1 = go_to_freq[go1] / n_gos 
    p2 = go_to_freq[go2] / n_gos
    go0, _ = gos_to_ancestor_specificity([go1] + [go2])
    if not go0: return 0
    p0 = go_to_freq[go0] / n_gos
    assert 0 < p1 and p1 <= 1
    assert 0 < p2 and p2 <= 1
    assert 0 < p0 and p0 <= 1
    ret = 2 * log(p0) / (log(p1) + log(p2)) if log(p1) + log(p2) != 0 else 1
    assert -1e-9 <= ret and ret <= 1+1e-9, 'p0,p1,p2,ret = {},{},{},{}'.format(p0, p1, p2, ret)
    return max((min((ret, 1)), 0))

def sim12(go1, gos): return max([sim11(go1, go2) for go2 in gos])
    
def sim22(gos1, gos2): return (sum([sim12(go1, gos2) for go1 in gos1]) + sum([sim12(go2, gos1) for go2 in gos2])) / (len(gos1) + len(gos2))

uniclsims = []
uniclsims_worst = []
for rep, goslist in rep_to_goslist.items():
    for gos1 in goslist:
        simsc = sum([sim22(gos1, gos2) for gos2 in goslist]) / len(goslist)
        assert 0 < simsc and simsc <= 1
        uniclsims.append(simsc)
        simsc = min([sim22(gos1, gos2) for gos2 in goslist])
        assert 0 <= simsc and simsc <= 1, '{} is not inclusively between 0 and 1'.format(simsc)
        uniclsims_worst.append(simsc)

nrecalls = 0
uniclsims2 = []
for rep, gosinner in rep_to_gos.items(): 
    uniclsims2.extend([sim22(gosinner, gosouter) for gosouter in rep_to_goslist[rep]])
    nrecalls += rep_to_clusize[rep] 
avginfdists = []
for rep, goslist in rep_to_goslist.items():
    for gos1 in goslist:
        infdists = []
        anc1, spec1 = gos_to_ancestor_specificity(gos1)
        for gos2 in goslist:
            anc2, spec2 = gos_to_ancestor_specificity(gos2)
            anc0, spec0 = gos_to_ancestor_specificity([anc1] + [anc2])
            p0 = go_to_freq[anc0] / n_gos
            p1 = go_to_freq[anc1] / n_gos
            p2 = go_to_freq[anc2] / n_gos
            assert p0 > 0
            assert p1 > 0
            assert p2 > 0
            infdists.append(2 * log(p0) / (log(p1) + log(p2)) if (log(p1) + log(p2) != 0) else 1)
        avginfdists.append(sum(infdists) / len(infdists))

specificities = []
ncat1 = 0
ncat2 = 0
ncat3 = 0
ncat4 = 0
for rep, goslist in rep_to_goslist.items():
    assert goslist
    #if len(goslist) == 1: continue
    theancestor, thespecificity = gos_to_ancestor_specificity([go for gos in goslist for go in gos])
    for gos in goslist:
        ances, spec = gos_to_ancestor_specificity([go for go in gos])
        specificities.append(spec - thespecificity)
    gos1 = goslist[0]
    if all([gos1 == gos for gos in goslist]):
        ncat1 += len(goslist)
    else:
        mcalist = [gos_to_ancestor_specificity([go for go in gos])[0] for gos in goslist]
        if len(set(mcalist)) == 1:
            ncat2 += len(goslist)
        else:
            ancestor, specificity = gos_to_ancestor_specificity(mcalist)
            #assert specificity == thespecificity, '{} != {}'.format(specificity, thespecificity)
            outliers = [mca for mca in mcalist if (mca not in term_to_subterms[ancestor] and mca != ancestor)] 
            if len(outliers) == 0:
                ncat3 += len(goslist)
            else:
                ncat4 += len(goslist)
repcnt = ncat1 + ncat2 + ncat3 + ncat4
if repcnt > 0:
    print('categories123 {:.7f}\t{:.7f}\t{:.7f} out-of \t{} curated-members \t{} curated-clusters \t{} clusters \t{} seqs \t{} avgspec \t{} curated-centroids \t{} avg-infdist \t{} uniclustsims \t{} uniclustsims2 \t{} ncurated-recalls \t{} nrecalls \t{} uniclustsims_worst'
          .format(ncat1/repcnt, (ncat1+ncat2)/repcnt, (ncat1+ncat2+ncat3)/repcnt, repcnt, len(rep_to_goslist), len(clusters), nseqs, float(sum(specificities)) / len(specificities), len(rep_to_gos), sum(avginfdists) / len(avginfdists), sum(uniclsims) / len(uniclsims), sum(uniclsims2) / len(uniclsims2), len(uniclsims2), nrecalls, sum(uniclsims_worst) / len(uniclsims_worst) ))
else:
    print('timeout')

exit(0)
specificities = []
for rep, goslist in rep_to_goslist.items():
    ancestor, specificity = gos_to_ancestor_specificity([go for gos in goslist for go in gos])
    assert(specificity >= 0), 'The following list of GOs does not share any ancestor{}'.format(goslist)
    for _ in goslist: specificities.append(float(specificity))
for specificity in sorted(specificities):
    print('{}'.format(specificity))
print('Average={} out of {} specificities'.format(float(sum(specificities))/float(len(specificities)+1e-9), len(specificities)))



exit(0)

jaccardindexes = []
for rep, goslist in rep_to_goslist.items():
    #sys.stderr.write('The rep {} has {} elemnents with GO annotation\n'.format(rep, len(goslist)))    
    for i in range(0, len(goslist), 1):
        for j in range(i+1, len(goslist), 1):
            gos1, gos2 = (goslist[i], goslist[j])
            jaccardindex = float(len(gos1 & gos2)) / float(len(gos1 | gos2))
            jaccardindexes.append(jaccardindex)
for jaccardindex in sorted(jaccardindexes):
    print('{}'.format(jaccardindex))
print('\n++++++++++\nAverage = {} out of {} jaccardindexes'.format(float(sum(jaccardindexes))/float(len(jaccardindexes)+1e-9), len(jaccardindexes)))

