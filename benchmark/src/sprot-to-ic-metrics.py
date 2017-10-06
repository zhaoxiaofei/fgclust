from __future__ import division
import itertools
import sys
import xml.etree.ElementTree as ET

from collections import defaultdict
from math import log

EXPCODES = set(['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP'])
CAFACODES = set(['EXP', 'IDA', 'IMP', 'IGI', 'IEP', 'TAS'])
def gos_to_ancestor_specificity(gos, term_to_ancestors, term_to_freq, tot_freq):
    gos = [go for go in gos]
    ancestors = term_to_ancestors[gos[0]]
    for go in gos[1:]:
        ancestors = ancestors & term_to_ancestors[go]
    most_specific_ancestor = None
    max_ic = 0
    for ancestor in ancestors:
        if -log(term_to_freq[ancestor] / tot_freq) > max_ic:
            most_specific_ancestor = ancestor
            max_ic = -log(term_to_freq[ancestor] / tot_freq)
    return (most_specific_ancestor, max_ic)

def sim11(go1, go2, term_to_ancestors, term_to_freq, tot_freq):
    p1 = go_to_freq[go1] / n_gos 
    p2 = go_to_freq[go2] / n_gos
    go0, _ = gos_to_ancestor_specificity([go1] + [go2], term_to_ancestors, term_to_freq, tot_freq)
    if not go0: return 0
    p0 = go_to_freq[go0] / n_gos
    assert 0 < p1 and p1 <= 1
    assert 0 < p2 and p2 <= 1
    assert 0 < p0 and p0 <= 1
    ret = 2 * log(p0) / (log(p1) + log(p2)) if log(p1) + log(p2) != 0 else 1
    assert -1e-9 <= ret and ret <= 1+1e-9, 'p0,p1,p2,ret = {},{},{},{}'.format(p0, p1, p2, ret)
    return max((min((ret, 1)), 0))

def sim12(go1, gos, term_to_ancestors, term_to_freq, tot_freq): 
    return max([sim11(go1, go2, term_to_ancestors, term_to_freq, tot_freq) for go2 in gos])
    
def sim22(gos1, gos2, term_to_ancestors, term_to_freq, tot_freq): 
    return ((sum([sim12(go1, gos2, term_to_ancestors, term_to_freq, tot_freq) for go1 in gos1]) 
           + sum([sim12(go2, gos1, term_to_ancestors, term_to_freq, tot_freq) for go2 in gos2])) 
           / (len(gos1) + len(gos2)))

def parse_go(fname):
    
    def init_term_to_reachableterms(term, term_to_nxtterms, term_to_reachableterms):
        if term not in term_to_reachableterms:
            reachableterms = set([term])
            for subterm in term_to_nxtterms[term]:
                reachableterms.update(init_term_to_reachableterms(subterm, term_to_nxtterms, term_to_reachableterms))
            term_to_reachableterms[term] = reachableterms
        return term_to_reachableterms[term]
    
    term_to_supterms = defaultdict(list)
    term_to_subterms = defaultdict(list)
    alt_to_id = {}
    term_to_namespace = {}
    root = ET.parse(fname).getroot()
    for term in root.findall('term'):
        child = term.find('id')
        alts = term.findall('alt_id')
        for alt in alts: alt_to_id[alt.text] = child.text
        alt_to_id[child.text] = child.text
        if term.find('is_obsolete') is not None and term.find('is_obsolete').text == '1': alt_to_id[child.text] = None
        namespace = term.find('namespace')
        term_to_namespace[child.text] = namespace.text
        for parent in term.findall('is_a'):
             term_to_supterms[child.text].append(parent.text)
             term_to_subterms[parent.text].append(child.text)
    term_to_offspring = {}
    term_to_ancestors = {}
    for term in term_to_namespace:
        term_to_offspring[term] = init_term_to_reachableterms(term, term_to_subterms, term_to_offspring)
        term_to_ancestors[term] = init_term_to_reachableterms(term, term_to_supterms, term_to_ancestors)
    return (alt_to_id, term_to_namespace, term_to_offspring, term_to_ancestors)

def parse_sprot(fname, alt_to_id, gotypes):
    acc_to_gos = {}
    with open(fname) as sprotfile:
        acc, gos = ('', set([]))
        for line in sprotfile:
            tokens = line.strip().split()
            if len(tokens) > 1 and tokens[0] == 'AC':
                if acc != '': acc_to_gos[acc] = gos
                acc, gos = (tokens[1].strip(';'), set([]))
            if (len(tokens) > 2 and tokens[0] == 'DR' and tokens[1] == 'GO;' 
                    and alt_to_id[tokens[2].strip(';')] and tokens[-1].split(':')[0] in CAFACODES and tokens[3].split(':')[0] in gotypes):
                gos.add(alt_to_id[tokens[2].strip(';')])
        if acc != '': acc_to_gos[acc] = gos
    return acc_to_gos

def parse_clu(fname, acc_to_gos):
    inner_to_outers = defaultdict(set)
    with open(fname) as clufile:
        for line in clufile:
            tokens = line.strip().split()
            inner = tokens[0].lstrip('UniRef100_')
            outer = tokens[1].lstrip('UniRef100_')
            if inner in acc_to_gos and outer in acc_to_gos and acc_to_gos[inner] and acc_to_gos[outer]:
            inner_to_outers[inner].add(outer)
    return inner_to_outers

def init_go_to_freq_n_gos(gos, term_to_ancestors):
    go_to_freq, n_gos = (defaultdict(int), 0)
    for go in gos:
        for anc in term_to_ancestors[go]:
            go_to_freq[anc] += 1
            n_gos += 1
    return (go_to_freq, n_gos)

def remove_unreliable_gos(acc_to_gos, term_to_freq, inner_to_outers, thres=10):
    ret = {}
    for acc, gos in acc_to_gos.items():
        ret[acc] = set([go for go in gos if term_to_freq[go] >= thres])
    ret_inner_to_outers = {}
    for inner, outers in inner_to_outers.items():
        ret_inner_to_outers[inner] = set([outer for outer in outers if outer in ret and inner in ret and ret[inner] and ret[outer]])
    return ret, ret_inner_to_outers

def compute_sim22s(inner_to_outers, acc_to_gos, term_to_ancestors, term_to_freq, tot_freq):
    perlnk_sim22s, permem_sim22s, perclu_sim22s = ([], [], [])
    norclu_sim22s = []
    for i, (inner, outers) in enumerate(inner_to_outers.items()):
        locclu_sim22s = []
        excclu_sim22s = []
        for acc1 in outers:
            if acc1 not in acc_to_gos or not acc_to_gos[acc1]: continue
            locmem_sim22s = []
            for acc2 in outers:
                if acc2 not in acc_to_gos or not acc_to_gos[acc2]: continue 
                gos1, gos2 = (acc_to_gos[acc1], acc_to_gos[acc2])
                sim22elem = sim22(gos1, gos2, term_to_ancestors, term_to_freq, tot_freq)
                locmem_sim22s.append(sim22elem)
                locclu_sim22s.append(sim22elem)
                perlnk_sim22s.append(sim22elem)
                if gos1 == gos2: assert 0.99 < sim22elem and sim22elem < 1.01
                if acc1 != acc2: excclu_sim22s.append(sim22elem)
            if locmem_sim22s: permem_sim22s.append(sum(locmem_sim22s) / len(locmem_sim22s))
        if locclu_sim22s: perclu_sim22s.append(sum(locclu_sim22s) /len(locclu_sim22s))
        if excclu_sim22s:
            for _ in range(0, len(outers)-1, 1): norclu_sim22s.append(sum(excclu_sim22s) / len(excclu_sim22s))
        if 0 == (i & (i-1)): sys.stderr.write('perlnk_sim22s {}th inner is {}\n'.format(i, inner))
    return (sorted(norclu_sim22s), sorted(perlnk_sim22s), sorted(permem_sim22s), sorted(perclu_sim22s))

def compute_sim22s_inn(inner_to_outers, acc_to_gos, term_to_ancestors, term_to_freq, tot_freq):
    perinn_sim22s = []
    perclu_sim22s = []
    for i, (inner, outers) in enumerate(inner_to_outers.items()):
        for outer in outers:
            if inner != outer and inner in acc_to_gos and outer in acc_to_gos and acc_to_gos[inner] and acc_to_gos[outer]:
                gos1, gos2 = (acc_to_gos[inner], acc_to_gos[outer])
                sim22elem = sim22(gos1, gos2, term_to_ancestors, term_to_freq, tot_freq)
                perinn_sim22s.append(sim22elem)
        locclu_sim = []
        for acc1 in outers:
            if acc1 not in acc_to_gos or not acc_to_gos[acc1]: continue
            for acc2 in outers:
                if acc2 not in acc_to_gos or not acc_to_gos[acc2] or acc1 == acc2: continue


    return sorted(perinn_sim22s)

def compute_infoloss(inner_to_outers, acc_to_gos, term_to_ancestors, term_to_freq, tot_freq, namespace, term_to_namespace):
    infolosses, infolosses_nclus, infolosses_nmems, perclu_infolosses = ([], 0, 0, [])
    for i, (inner, outers) in enumerate(inner_to_outers.items()):
        outers = [outer for outer in outers if outer in acc_to_gos and acc_to_gos[outer]]
        goslist = [[go for go in acc_to_gos[outer] if term_to_namespace[go] == namespace] for outer in outers]
        goslist = [gos for gos in goslist if gos]
        if not goslist: continue
        infolosses_nclus += 1
        infolosses_nmems += len(anc_spec_list)
        anc_spec_list = [gos_to_ancestor_specificity(gos, term_to_ancestors, term_to_freq, tot_freq) for gos in goslist]    
        rootanc, rootspec = gos_to_ancestor_specificity([anc for (anc, spec) in anc_spec_list], term_to_ancestors, term_to_freq, tot_freq)
        permem_infolosses = [(spec - rootspec) for anc, spec in anc_spec_list]
        perclu_infolosses.append(sum(permem_infolosses) / len(permem_infolosses))
        infolosses.extend(permem_infolosses)
        if 0 == (i & (i-1)): sys.stderr.write('infolosses {}th inner is {}\n'.format(i, inner))
    return (sorted(infolosses), infolosses_nclus, infolosses_nmems, sorted(perclu_infolosses))

def avg(vals): return sum(vals) / (1e-9 + len(vals))
def percentile(vals, i): return vals[int((len(vals)-1)*i/100)] if len(vals) > 0 else 0

alt_to_id, term_to_namespace, term_to_offspring, term_to_ancestors = parse_go(sys.argv[1])
acc_to_gos = parse_sprot(sys.argv[2], alt_to_id, sys.argv[4])
inner_to_outers = parse_clu(sys.argv[3], acc_to_gos)
go_to_freq, n_gos = init_go_to_freq_n_gos([go for gos in acc_to_gos.values() for go in gos], term_to_ancestors)
acc_to_gos, inner_to_outers = remove_unreliable_gos(acc_to_gos, go_to_freq, inner_to_outers)

#finfolosses, finfolosses_nclus, finfolosses_nmems, perclu_finfolosses = compute_infoloss(
#        inner_to_outers, acc_to_gos, term_to_ancestors, go_to_freq, n_gos, 'molecular_function', term_to_namespace)
perinn_sim22s = compute_sim22s_inn(inner_to_outers, acc_to_gos, term_to_ancestors, go_to_freq, n_gos)
norclu_sim22s, perlnk_sim22s, permem_sim22s, perclu_sim22s = compute_sim22s(inner_to_outers, acc_to_gos, term_to_ancestors, go_to_freq, n_gos)
#assert len(perclu_finfolosses) == len(permem_sim22s), '{} == {} failed'.format(len(perclu_finfolosses), len(permem_sim22s))
print('n--clus-mems {} {} sim22s--per-inn-lnk-mem-clu--len-avg-25perc-10perc {} {:.4f} {:.4f} {:.4f} {} {:.4f} {:.4f} {:.4f} {} {:.4f} {:.4f} {:.4f} {} {:.4f} {:.4f} {:.4f} {} {:.4f} {:.4f} {:.4f}'
      .format(len(inner_to_outers), len([outer for outers in inner_to_outers.values() for outer in outers]),
              len(norclu_sim22s), avg(norclu_sim22s), percentile(norclu_sim22s, 25), percentile(norclu_sim22s, 10),
              len(perinn_sim22s), avg(perinn_sim22s), percentile(perinn_sim22s, 25), percentile(perinn_sim22s, 10),
              len(perlnk_sim22s), avg(perlnk_sim22s), percentile(perlnk_sim22s, 25), percentile(perlnk_sim22s, 10),
              len(permem_sim22s), avg(permem_sim22s), percentile(permem_sim22s, 25), percentile(permem_sim22s, 10),
              len(perclu_sim22s), avg(perclu_sim22s), percentile(perclu_sim22s, 25), percentile(perclu_sim22s, 10)))
