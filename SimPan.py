import subprocess, sys, os, numpy as np, ete3, re, glob, shutil
from multiprocessing import Pool

complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
def rc(seq, missingValue='N') :
    return ''.join([complement.get(s, missingValue) for s in reversed(seq.upper())])

externals = {
    'simbac' : shutil.which('SimBac') if shutil.which('SimBac') else os.path.join(os.path.dirname(os.path.realpath(__file__)), 'dependencies', 'SimBac'),
    'indel-seq-gen' : shutil.which('indel-seq-gen') if shutil.which('indel-seq-gen') else os.path.join(os.path.dirname(os.path.realpath(__file__)), 'dependencies', 'indel-seq-gen'),
}

def getEvent(trees, tipAcc, num) :
    results = []
    w = np.array([t[1] for t in trees]).astype(float)
    eventSources = np.random.choice(len(trees), size=num, p=w/np.sum(w))
    for tId, cnt in zip(*np.unique(eventSources, return_counts=True)) :
        tree = ete3.Tree(trees[tId][2])
        annotateTree(tree)
        indelTree(tree, tipAcc)

        nodes, weights = [], []
        for n in tree.get_descendants() :
            nodes.append(n.descendants)
            weights.append(n.dist)
        weights = np.array(weights)
        picks = np.random.choice(nodes, size=cnt, p=weights/np.sum(weights))
        results.extend(picks)
    return results
def annotateTree(tree) :
    for node in tree.traverse('postorder') :
        if node.is_leaf() :
            node.descendants = [int(node.name)]
            node.height = 0.0
        else :
            node.descendants = [n for c in node.get_children() for n in c.descendants]
            node.height = np.mean([c.height + c.dist for c in node.get_children()])
    return tree

def indelTree(tree, tipAccelerate) :
    lmbda = np.log(tipAccelerate)#/tree.height    
    for node in tree.traverse('postorder') :
        node.dist = lmbda*abs(np.exp(-lmbda*node.height) - np.exp(-lmbda*(node.height+node.dist)))

def borderedGeometric(mean, low, high, size) :
    data = np.random.geometric(p=1./mean, size=size)
    data[data < low] = low
    data[data > high] = high
    return data


def borderedNB(mean, low, high, size, r=2) :
    data = np.random.negative_binomial(n=r, p=(1. - mean/(2.+mean)), size=size)
    data[data < low] = low
    data[data > high] = high
    return data

def setGeneContent(items, itrees, dtrees, tipAccelerate, genomeNum, aveSize, nCore, backboneBlockSize, mobileBlockSize) :
    genomes = np.random.permutation(np.repeat(np.where(items.T[0] == 0)[0], genomeNum).reshape(-1, genomeNum))
    
    events = getEvent(dtrees, tipAccelerate, genomes.shape[0] )
    backBoneSite = np.random.choice(genomes.shape[0], genomes.shape[0])
    backBoneSize = borderedGeometric(backboneBlockSize[0], backboneBlockSize[1], backboneBlockSize[2], genomes.shape[0])
    
    coreIdx = np.ones(genomes.shape[0], dtype=bool)
    for b, s, l in zip(events, backBoneSite, backBoneSize) :
        e = s + l
        c1 = np.sum(coreIdx)
        
        if e <= genomes.shape[0] :
            coreIdx[s:e] = False
            c2 = np.sum(coreIdx)
            if abs(c1-nCore) < abs(c2-nCore) :
                break
            genomes[s:e, b] = -1
        else :
            l1 = genomes.shape[0] - s
            l2 = l - l1
            coreIdx[s:e] = False
            c2 = np.sum(coreIdx)
            if abs(c1-nCore) < abs(c2-nCore) :
                break            
            genomes[s:, b], genomes[:l2, b] = -1, -1
    
    genePool = np.random.permutation(np.where(items.T[0] == 1)[0])
    toAdd = aveSize*genomeNum - np.sum(genomes >= 0)
    events = getEvent(itrees, tipAccelerate, toAdd)
    mobileBlock = borderedGeometric(mobileBlockSize[0], mobileBlockSize[1], mobileBlockSize[2], toAdd)
    gInEvent = np.cumsum([n*len(e) for n, e in zip(mobileBlock, events)])
    
    events = events[:np.argmin(np.abs(gInEvent - toAdd))+1]
    mobiles = []
    for b, s, l in zip(events, np.random.choice(len(genePool), len(events)), mobileBlock) :
        mobile = np.empty([l, genomeNum], dtype=int)
        mobile[:] = -1
        e = s + l
        if e <= len(genePool) :
            mobile[:, b] = genePool[s:e, np.newaxis]
        else :
            l1 = len(genePool) - s
            l2 = l - l1
            mobile[:l1,  b] = genePool[s:,  np.newaxis]
            mobile[-l2:, b] = genePool[:l2, np.newaxis]
        mobiles.append(mobile)
        
    blk = (np.random.dirichlet(np.ones(len(mobiles)+1, dtype=int))*genomes.shape[0]+0.5).astype(int)
    genomes = np.split(genomes, np.cumsum(blk)[:-1])
    genomes = np.vstack([nn for n in zip(genomes[:-1], mobiles) for nn in n] + genomes[-1:])
    
    return genomes

def treeGenerator(prefix, items, genomeNum, rec, recLen) :
    genomeLength = int(np.sum(items[:, 2:]))
    
    r0 = float(recLen)/genomeLength
    r = np.log(1-rec)/np.log(1-r0)/2/genomeLength
    if r <= 0 :
        r = 0
    cmd = '{simbac} -N {0} -T 0 -R {1} -B {2} -c {3}.phylogeny -l {3}.local'.format(\
        genomeNum, r, genomeLength, prefix, **externals
    )

    sys.stdout.write(cmd+'\n')
    p = subprocess.Popen(cmd.split(), universal_newlines=True, stdout=subprocess.PIPE).wait()
    
    globaltrees, localtrees = [], []
    c = 0
    with open('{0}.phylogeny'.format(prefix), 'r') as fin :
        for line in fin :
            globaltrees.append([genomeLength, genomeLength, line.strip()])
    with open('{0}.local'.format(prefix), 'r') as fin :
        for line in fin :
            if line.startswith('[') :
                w, t = line.strip()[1:].split(']')
                w = int(w)
                c += w
                localtrees.append([c, w, t])

    os.unlink('{0}.local'.format(prefix))
    return globaltrees, localtrees
    
def itemShaper(items, geneLen, igrLen) :
    nOrtho = np.max((items.T[1]/10).astype(int))+1
    orthos = borderedNB(geneLen[0], geneLen[1], geneLen[2], nOrtho, r=2.)
    orthos *= 3
    oigrs1 = borderedGeometric(igrLen[0], igrLen[1], igrLen[2], nOrtho)
    oigrs2 = borderedGeometric(igrLen[0], igrLen[1], igrLen[2], nOrtho)
    genes = orthos[(items.T[1]/10).astype(int)]
    igrs1 = oigrs1[(items.T[1]/10).astype(int)]
    igrs2 = oigrs2[(items.T[1]/10).astype(int)]
    return np.hstack([items, igrs1[:, np.newaxis], genes[:, np.newaxis], igrs2[:, np.newaxis]])

def itemGenerator(nBackbone, nMobile, pBackbone, pMobile) :
    types = np.concatenate([np.zeros(nBackbone, dtype=int), np.ones(nMobile, dtype=int)])

    genes = []
    for n, p in zip([nBackbone, nMobile], [pBackbone, pMobile]) :
        p2 = p - np.sqrt(p*(p+3)) + 1
        nOrtho = int(n / (1+2*p2) + 0.5)
        orthoGene = np.arange(nOrtho)
        gene = np.sort(np.concatenate([np.random.choice(np.repeat(orthoGene, 2), n - nOrtho, replace=False), orthoGene]))*10
        for i, g in enumerate(gene[1:]) :
            if int(g/10) == int(gene[i]/10) :
                gene[i+1] = gene[i]+1
        if len(genes) :
            gene += np.max(genes[-1])+10
        genes.append(gene)
    genes = np.concatenate(genes)
    ortho = np.unique(genes)
    return np.vstack([types, genes]).T


def alignTrees(trees) :
    result = [[], [], []]
    for idx, ts in enumerate(zip(*trees)) :
        res = result[idx]
        indices = [ 0 for t in ts ]
        while indices[0] < len(ts[0]) :
            s = np.min([ t[i][0] for t, i in zip(ts, indices) ])
            res.append([s] + [ t[i][1] for t, i in zip(ts, indices) ])
            for t, i in zip(ts, indices) :
                t[i][0] -= s
            for ti, t in enumerate(ts) :
                i = indices[ti]
                if t[i][0] <= 0 :
                    indices[ti] += 1
    return result
def duplicateOrthoTree(trees1, trees2, distDuplication=0.005) :
    result = alignTrees([trees1, trees2])
    tips = np.random.permutation(result[1][0][1].get_leaf_names())
    dists = np.random.exponential(distDuplication, size=len(tips))
    events = [[t.split('_', 1)[0], t, d] for t, d in zip(tips.tolist(), dists)]
    for res in result :
        for r in res :
            w, t1, t2 = r
            tr, ts = t1.copy(), t2.copy()
            n1 = {}
            for n in tr.get_leaves() :
                o = int(n.name.split('_', 1)[0])
                if o not in n1 :
                    n1[o] = [n]
                else :
                    n1[o].append(n)
            n2 = { n.split('_', 1)[0]:n for n in ts.get_leaf_names() }
            for o, cn1, dist in events :
                c2 = ts.get_leaves_by_name(n2[o])
                if not c2 :
                    continue
                c2 = c2[0]
                c1 = tr.get_leaves_by_name(cn1)[0]
                d1 = dist
                while d1 - c1.dist > 0 and c1.up :
                    d1, c1 = d1 - c1.dist, c1.up
                d2 = dist
                while d2 - c2.dist > 0 and c2.up :
                    d2, c2 = d2 - c2.dist, c2.up
                p, l = c1.up, c1.dist
                
                n = ete3.TreeNode()
                if p is not None :
                    p.remove_child(c1)
                n.add_child(c1)
                c1.up, c1.dist = n, d1
                
                if p is not None :
                    p.add_child(n)
                    n.dist = l - c1.dist
                else :
                    tr = n
                    for _, x, _ in events :
                        tr.get_leaves_by_name(x)[0]
                p2 = c2.up
                if p2 is not None :
                    p2.remove_child(c2)
                n.add_child(c2)
                c2.up, c2.dist = n, d2
                if p2 is not None :
                    gp = p2.up
                    if gp :
                        gp.remove_child(p2)
                        for c in p2.get_children() :
                            p2.remove_child(c)
                            gp.add_child(c)
                            c.up, c.dist = gp, c.dist + p2.dist
                    elif len(p2.get_children()) == 1 :
                        ts = p2.get_children()[0]
                        ts.up = None
                else :
                    break
            r[:] = [w, tr]
            
    return result
                

def getOrthoTree(trees, distOrtholog, postfix, strand, s, e, t, v, i0=0) :
    res = [[], [], []]
    while v > trees[-1][0] :
        x = min(s, trees[-1][0])
        s, e, t, v = s-x, e-x, t-x, v-x
    for r, i, j in ((res[0], s, e), (res[1], e, t), (res[2], t, v)) :
        while trees[i0][0] <= i :
            i0 += 1
        for m in range(i0, len(trees)) :
            tree = trees[m]
            i1, j1 = max(tree[0]-tree[1], i), min(tree[0], j)
            if j1 - i1 > 0 :
                tre = ete3.Tree(tree[2])
                for n in tre.get_descendants() :
                    n.dist *= distOrtholog
                for n in tre.get_leaves() :
                    n.name = '{1}_{0}'.format(postfix, n.name)            
                r.append([int(j1-i1), tre])
            else :
                break
    if strand < 0 :
        res = [ r[::-1] for r in res[::-1] ]
    return res, i0

def mergeOrthoTrees(prefix, orthoTrees, distParalog) :
    genes = list(orthoTrees.values())
    trees = alignTrees(genes)

    cmd = '{simbac} -N {0} -T 0 -R 0 -c {1}.gene'.format(\
        len(genes), prefix, **externals )
    subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, universal_newlines=True).communicate()
    paratree = ete3.Tree('{0}.gene'.format(prefix))
    os.unlink('{0}.gene'.format(prefix))
    for n in paratree.get_descendants() :
        n.dist *= distParalog

    results = [[], [], []]
    for res, ts in zip(results, trees) :
        for tt in ts :
            tree = paratree.copy()
            for n in tree.get_leaves() :
                t = tt[int(n.name)+1]
                p = n.up
                p.remove_child(n)
                p.add_child(t)
                t.up, t.dist = p, n.dist
            res.append([tt[0], tree])
    return results

i00 = 0
def getHomologTree(prefix, homologs, trees, distOrtholog, distParalog, distDuplication) :
    orthoTrees = {}
    global i00
    i0s = []
    for gene in homologs :
        if gene[3] not in orthoTrees :
            orthoTrees[gene[3]], i0 = getOrthoTree(trees, distOrtholog, gene[0], gene[1], gene[7], gene[8], gene[9], gene[10], i00)
        else :
            tree, i0 = getOrthoTree(trees, distOrtholog, gene[0], gene[1], gene[7], gene[8], gene[9], gene[10], i00)
            orthoTrees[gene[3]] = duplicateOrthoTree(orthoTrees[gene[3]], tree, distDuplication)
        i0s.append(i0)
    i00 = np.min(i0s)
    if len(orthoTrees) > 1 :
        homoTree = mergeOrthoTrees(prefix, orthoTrees, distParalog)
    else :
        homoTree = list(orthoTrees.values())[0]
    return homoTree

def getHomoSequence(prefix, homologs, homologTree, indelRate, indelMax, freqStart, freqStop) :
    tips = [int(t.split('_')[0])-1 for t in homologTree[1][0][1].get_leaf_names()]
    results = { homolog:{t:[] for t in tips} for homolog in homologs.T[0] }
    for region, treeBlock in zip(['igr', 'cds', 'igr'], homologTree) :
        if region == 'cds' :
            for i, (w, t) in reversed(list(enumerate(treeBlock))) :
                if w % 3 > 0 :
                    if i > 0 :
                        if w % 3 == 1 :
                            treeBlock[i-1][0] += 1
                        else :
                            treeBlock[i-1][0] -= 1
                    if w % 3 == 1 :
                        treeBlock[i][0] -= 1
                    else :
                        treeBlock[i][0] += 1
        for i, (w, t) in reversed(list(enumerate(treeBlock))) :
            if w <= 1 :
                if i > 0 :
                    treeBlock[i-1][0] += w
                elif i < len(treeBlock) - 1 :
                    treeBlock[i+1][0] += w
                treeBlock[i:i+1] = []
        if len(treeBlock) :
            codon = '' if region == 'igr' else '-c 2,1,8'
            x = subprocess.Popen('{indel-seq-gen} -m HKY -e {0} {2} -m HKY -o p -u sp'.format(\
                prefix, region, codon, **externals
                ).split(), universal_newlines=True, stderr=subprocess.PIPE, \
                             stdin=subprocess.PIPE).communicate(\
                                 input='\n'.join(['[{0}]{{{2},{3}/{3},seqgen.{4}.indel}}{1}'.format(tree[0], tree[1].write(format=1, dist_formatter='%0.12f'),\
                                                                                       indelMax, indelRate/2, region
                                                                                       ) for tree in treeBlock])+'\n'\
                             )
            with open('{0}.ma'.format(prefix, region)) as fin :
                fin.readline()
                for line in fin :
                    part = line.strip().split()
                    if len(part) > 1 :
                        n, s = part
                        go, ge = np.array(n.split('_'), dtype=int)
                        if go-1 not in results[ge] :
                            results[ge][go-1] = [s]
                        else :
                            results[ge][go-1].append(s)
        else :
            for ge in results :
                for go in results[ge] :
                    results[ge][go].append('')
    for fname in glob.glob('{0}.*'.format(prefix, region)) :
        os.unlink(fname)

    starts = np.unique([v[1].replace('-', '')[:3] for gene_seq in results.values() for v in gene_seq.values()], return_counts=True)
    starts = starts[0][np.argsort(-starts[1])]
    n = int(np.ceil(float(len(starts))/len(freqStart[0])))
    if n <= 1 :
        startConvs = np.random.choice(freqStart[0], p=freqStart[1], size = starts.shape[0], replace=False)
    else :
        fs0 = freqStart[0] * n
        fs1 = np.array(freqStart[1].tolist() * n)/n
        startConvs = np.random.choice(fs0, p=fs1, size = starts.shape[0], replace=False)
    startConvs = dict(zip(starts, startConvs))
    
    stops = np.unique([v[1].replace('-', '')[-3:] for gene_seq in results.values() for v in gene_seq.values()], return_counts=True)
    stops = stops[0][np.argsort(-stops[1])]
    n = int(np.ceil(float(len(stops))/len(freqStop[0])))
    if n <= 1 :
        stopConvs = np.random.choice(freqStop[0], p=freqStop[1], size = stops.shape[0], replace=False)
    else :
        fs0 = freqStop[0] * n
        fs1 = np.array(freqStop[1].tolist() * n)/n
        stopConvs = np.random.choice(fs0, p=fs1, size = stops.shape[0], replace=False)
    stopConvs = dict(zip(stops, stopConvs))
    
    neg_strand = set(homologs[homologs.T[1] == -1, 0])    
    for gene, genomes in results.items() :
        for genome, seq in genomes.items() :
            seq[1] = re.sub(r'(.{3})(-*)$', lambda g:g.group(2) + stopConvs[g.group(1)], re.sub(r'^(-*)(.{3})', lambda g:startConvs[g.group(2)] + g.group(1), seq[1]))
            if gene in neg_strand :
                seq = [ rc(s, '-') for s in seq[::-1] ]
            m, n, k, m0, n0, k0 = len(seq[0]), len(seq[1]), len(seq[2]), len(seq[0].replace('-', '')), len(seq[1].replace('-', '')), len(seq[2].replace('-', ''))
            genomes[genome] = [''.join(seq), m+1, m+n, m+n+k, m0+1, m0+n0, m0+n0+k0, -1 if gene in neg_strand else 1]
    return results

def getSequences(prefix, genomes, items, operonBlock, trees, idenOrtholog, idenParalog, idenDuplication, indelRate, indelMax, freqStart, freqStop) :
    distOrtholog = -3./4.*np.log(1.-4./3.*(1.-idenOrtholog))/2.
    distParalog = -3./4.*np.log(1.-4./3.*(1.-idenParalog))/2. - distOrtholog
    distDuplication = -3./4.*np.log(1.-4./3.*(1.-idenDuplication))/2.
    genes = np.max(genomes, 1)
    blocks = np.cumsum(borderedGeometric(operonBlock[0], operonBlock[1], operonBlock[2], genes.size))
    blocks = blocks[blocks < genes.size]
    strand = np.random.choice([-1, 1], 1)
    strands = np.zeros(genes.size, dtype=int)
    for s, e in zip([0] + blocks.tolist(), blocks.tolist()+[genes.size]) :
        strands[s:e] = strand
        strand *= -1
    genes = np.hstack([np.arange(genomes.shape[0])[:, np.newaxis], strands[:, np.newaxis], items[genes]])
    genes[genes.T[1] == -1, 4:] = genes[genes.T[1] == -1, 4:][:, ::-1]
    sites = np.cumsum(genes[:, 4:]).reshape(-1, 3)
    genes = np.hstack([genes, np.concatenate([[0], sites[:-1, -1]])[:, np.newaxis], sites])
    
    alignments, genome_seqs, annotations, archives = [ [] for g in np.arange(genomes.shape[1]) ], [ [] for g in np.arange(genomes.shape[1]) ], [ [] for g in np.arange(genomes.shape[1]) ], {}
    with open('seqgen.igr.indel', 'w') as fout :
        fout.write('2627.66769806,743.784419906488,355.480238119831,210.534712476797,140.238546142805,100.621803469817,75.9965040586087,59.5936994260522,48.0906317747209,39.6957521563191,33.3714549098946,28.4818852014648,24.6191452491735,21.5114779269409,18.9719696342231,16.8685200150737,15.1055503871874,13.6124756885757,12.3362203381659,11.2362312830263,10.2810775422532,9.44608340313232,8.71165041315475,8.06204775362389,7.48452699660906,6.96866528490998,6.50587182208479,6.08901275569724,5.71212299608529,5.37018262974892,5.05894185533674,4.77478274111892,4.51460918994168,4.27575870434353,4.05593113961235,3.85313079770075,3.66561907437647,3.49187551182199,3.33056558937959,3.1805139489103,3.0406820287779,2.91014929377263,2.78809741335811,2.67379686920619,2.56659557376873,2.46590916110141,2.37121267416572,2.28203342306742,2.19794482894246,2.11856110061303,2.04353261735985,1.97254191246708,1.90530016958954,1.84154415824645,1.78103354647507,1.72354853836471,1.66888779222258,1.61686658180352,1.56731516861245,1.52007735795841,1.47500921536305,1.4319779232322,1.39086076049377,1.35154419027362,1.31392304269358,1.27789978159125,1.24338384542653,1.21029105389477,1.17854307284416,1.14806693102194,1.11879458297388,1.09066251311289,1.06361137657111,1.03758567297005,1.01253344969503,0.988406031654495,0.965157774848733,0.942745841373371,0.921129993746859,0.900272406682654,0.880137494630415,0.860691753589846,0.841903615859037,0.823743316519059,0.806182770580176,0.78919545982484,0.772756328479934,0.756841686937273,0.741429122818457,0.726497418748777,0.712026476266241,0.697997245346518,0.684391659073637,0.671192573030209,0.658383709020344,0.645949602773805,0.633875555311786,0.622147587683329,0.610752398807241,0.599677326177617,0.58891030921216,0.578439855041472,0.568255006554736,0.558345312532765,0.548700799713571,0.539311946648411,0.530169659217982,0.521265247688958,0.512590405200812,0.504137187581589,0.495897994399363,0.487865551163388,0.480032892595684,0.472393346899888,0.464940520959775,0.457668286405011,0.450570766486371,0.443642323706956,0.436877548159893,0.430271246526647,0.423818431693374,0.417514312945853,0.411354286706303,0.405333927778046,0.39944898106635,0.393695353745956,0.388069107847919,0.382566453240197,0.377183740978193,0.37191745700308,0.366764216167198,0.361720756567208,0.356783934166972,0.351950717693281,0.347218183788694,0.342583512406741,0.338043982435683,0.333596967537936,0.329239932193047,0.324970427932885,0.320786089758426,0.316684632728148,0.312663848708688,0.308721603278961,0.304855832779504,0.301064541499273,0.297345798992609,0.293697737519508,0.290118549602749,0.286606485695794,0.283159851955747,0.279777008115973,0.276456365453306,0.273196384845048,0.269995574911242,0.266852490237955,0.263765729677552,0.260733934722148,0.257755787946661,0.254830011518064,0.25195536576763,0.249130647823143,0.24635469029819,0.243626360035843,0.240944556904132,0.238308212640893,0.235716289745684,0.233167780416573,0.230661705529731,0.228197113659875,0.225773080139685,0.223388706156424,0.221043117884104,0.21873546564958,0.216464923131075,0.21423068658769,0.212031974118538,0.209868024950204,0.207738098751289,0.205641474972868,0.203577452213756,0.201545347609496,0.199544496244073,0.19757425058339,0.195633979929586,0.193723069895309,0.19184092189714,0.189986952667338,0.188160593783181,0.186361291213162,0.18458850487935,0.182841708235274,0.181120387858689,0.179424043058626,0.177752185496167,0.176104338818378,0.174480038304898,0.172878830526678,0.171300273016381,0.169743933950018,0.168209391839349,0.166696235234656,0.165204062437485,0.16373248122297,0.162281108571374,0.160849570408511,0.159437501354692,0.158044544481898,0.156670351078856,0.155314580423733,0.153976899564159,0.152656983104326,0.15135451299887,0.150069178353336,0.14880067523094,0.147548706465428,0.146312981479816,0.145093216110773,0.143889132438476,0.142700458621727,0.141526928738149,0.140368282629282,0.139224265750408,0.138094629024942,0.136979128703231,0.135877526225601,0.134789588089517,0.133715085720704,0.132653795348108,0.131605497882552,0.130569978798971,0.1295470280221,0.128536439815511,0.127538012673867,0.126551549218306,0.125576856094842,0.124613743875682,0.123662026963371,0.122721523497662,0.121792055265032,0.120873447610757,0.119965529353452,0.119068132702017,0.118181093174894,0.117304249521565,0.116437443646235,0.115580520533605,0.114733328176693,0.113895717506622,0.11306754232433,0.112248659234118,0.111438927579011,0.110638209377844,0.109846369264042,0.109063274426035,0.108288794549258,0.107522801759683,0.106765170568846,0.106015777820322,0.105274502637592,0.104541226373283,0.103815832559718,0.10309820686075,0.102388237024838,0.101685812839328,0.100990826085908,0.100303170497194,0.0996227417144249,0.0989494372462204,0.0982831564283899,0.0976238003847422,0.0969712719888816,0.0963254758269547,0.095686318161324,0.095053706895142,0.0944275515377979,0.0938077631712138,0.0931942544169696,0.092586939404226,0.0919857337384313,0.0913905544707846,0.0908013200684377,0.0902179503854141,0.0896403666342237,0.0890684913581586,0.0885022484042469,0.0879415628968477,0.0873863612118733,0.0868365709516161,0.0862921209201679,0.0857529410994147,0.0852189626255917,0.0846901177663816,0.0841663398985463,0.083647563486072,0.0831337240588201,0.0826247581916649,0.0821206034841111,0.0816211985403725,0.081126482949906')
    with open('seqgen.cds.indel', 'w') as fout :
        fout.write('0,0,3726.93235608632,0,0,451.395062089419,0,0,183.680835259382,0,0,101.549092267679,0,0,65.1025928103375,0,0,45.5865460908368,0,0,33.8535291634455,0,0,26.219781569911,0,0,20.9590641036038,0,0,17.1713183815315,0,0,14.3483337863973,0,0,12.1848206416566,0,0,10.488060175578,0,0,9.13134527146083,0,0,8.02848985633303,0,0,7.11915525833455,0,0,6.36003854691535,0,0,5.71938624030307,0,0,5.17346987706236,0,0,4.70425910837438,0,0,4.29784789908902,0,0,3.94336701455846,0,0,3.63221797216546,0,0,3.35752402710871,0,0,3.1137304992362,0,0,2.8963096478766,0,0,2.70153989505993,0,0,2.52633868596794,0,0,2.36813455888495,0,0,2.22476822850451,0,0,2.0944153806864,0,0,1.97552588482436,0,0,1.86677554180236,0,0,1.76702749043125,0,0,1.67530111880107,0,0,1.59074685355535,0,0,1.51262558718176,0,0,1.44029179065896,0,0,1.37317957385116,0,0,1.3107911183935,0,0,1.25268703134553,0,0,1.19847826259035,0,0,1.14781930206631,0,0,1.10040242973749,0,0,1.05595283564895,0,0,1.01422446238036,0,0,0.974996449884357,0,0,0.938070084715798,0,0,0.903266173271386,0,0,0.870422772818051,0,0,0.839393225525025,0,0,0.810044449994245,0,0,0.782255452346362,0,0,0.755916025108837,0,0,0.730925607238165,0,0,0.70719228280315,0,0,0.684631899329291,0,0,0.663167289690108,0,0,0.642727583837302,0,0,0.623247598674361,0,0,0.604667296067325,0,0,0.586931300408285,0,0,0.569988468347658,0,0,0.553791504327786,0,0,0.538296616413482,0,0,0.523463207649954,0,0,0.509253598805748,0,0,0.495632778895111,0,0,0.482568180334576,0,0,0.470029475984486,0,0,0.457988395667355,0,0,0.446418560049704,0,0,0.435295330029065,0,0,0.424595669989159,0,0,0.414298023478581,0,0,0.404382200035822,0,0,0.394829272029631,0,0,0.385621480511478,0,0,0.376742149188829,0,0,0.368175605726065,0,0,0.359907109666226,0,0,0.351922786342694,0,0,0.34420956621692,0,0,0.336755129137459,0,0,0.329547853067921,0,0,0.322576766877788,0,0,0.315831506831196,0,0,0.309302276445306,0,0,0.30297980942243,0,0,0.296855335389035,0,0,0.290920548200579,0,0,0.285167576594264,0,0,0.279588956992409,0,0,0.274177608277653,0,0,0.268926808377796,0,0,0.263830172512968,0,0,0.258881632971199,0,0,0.25407542029052,0,0,0.249406045736557,0,0,0.24486828497439')
    for gId, gene in enumerate(genes) :
        print(gId, gene)
        if gId not in archives :
            homologs = genes[(genes.T[3]/10).astype(int) == (gene[3]/10).astype(int)]
            homologTree = getHomologTree(prefix, np.random.permutation(homologs), trees, distOrtholog, distParalog, distDuplication)
            homologSeq = getHomoSequence('{0}.gene.{1}'.format(prefix, gId), homologs, homologTree, indelRate, indelMax, freqStart, freqStop)
            archives.update(homologSeq)
        gene_seqs = archives.pop(gId)
        gene_stat = genomes[gId]
        for genomeId, gseq in gene_seqs.items() :
            if gene_stat[genomeId] >= 0 :
                alignments[genomeId].append(gseq[0])
                genome_seqs[genomeId].append(gseq[0].replace('-', ''))
                annotations[genomeId].append([gene[3]]+gseq[1:])
            else :
                alignments[genomeId].extend('-'*len(gseq[0]))
                annotations[genomeId].append([gene[3]] + gseq[1:4] + [ 0, 0, 0, gseq[7]])
    alignments = [ ''.join(aln) for aln in alignments ]
    genome_seqs = [ ''.join(genome) for genome in genome_seqs ]
    annotations = np.array(annotations)
    aSite = np.cumsum(annotations[:, :, 3], axis=1) - annotations[:, :, 3]
    annotations[:, :, 1:4] += aSite[:, :, np.newaxis]
    gSite = np.cumsum(annotations[:, :, 6], axis=1) - annotations[:, :, 6]
    annotations[:, :, 4:7] += gSite[:, :, np.newaxis]
    return alignments, genome_seqs, annotations

def writeOut(prefix, alignments, genome_seqs, annotations) :
    with open(prefix+'.aligned.fasta', 'w') as fout :
        for id, aln in enumerate(alignments) :
            fout.write('>{0}\n{1}\n'.format(id, aln))
    with open(prefix+'.aligned.tbl', 'w') as fout :
        for annotation in annotations[0] :
            fout.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(annotation[1], annotation[2], '+' if annotation[7]>0 else '-', int(annotation[0]/10), annotation[0]%10))
    for id, (genome, ann) in enumerate(zip(genome_seqs, annotations)) :
        with open(prefix + '_{0}.gff'.format(id), 'w') as fout :
            fout.write('##gff-version 3\n##sequence-region {0} 1 {1}\n'.format(id, len(genome)))
            for gid, annotation in enumerate(ann) :
                if annotation[5] > annotation[4] :
                    fout.write('{0}\tGenomeSim\tgene\t{1}\t{2}\t.\t{3}\t.\tID={4}_gene;gene={5};locus_tag={4}\n'.format(\
                        id, annotation[4], annotation[5], '+' if annotation[7]>0 else '-', \
                        'G_{0}_{1}'.format(id, gid+1), \
                        '{0}_{1}'.format(int(annotation[0]/10), annotation[0]%10)) )
                    fout.write('{0}\tGenomeSim\tCDS\t{1}\t{2}\t.\t{3}\t0\tID={4};Parent={4}_gene;gene={5};locus_tag={4};transl_table=11;product=simulated_gene\n'.format(\
                        id, annotation[4], annotation[5], '+' if annotation[7]>0 else '-', \
                        'G_{0}_{1}'.format(id, gid+1), \
                        '{0}_{1}'.format(int(annotation[0]/10), annotation[0]%10)) )
            fout.write('##FASTA\n')
            fout.write('>{0}\n{1}\n'.format(id, '\n'.join([genome[i:i+80] for i in range(0, len(genome), 80)]) ))
    for id, genome in enumerate(genome_seqs) :
        with open(prefix+'_{0}.fna'.format(id), 'w') as fout :
            fout.write('>{0} [organism=Simulated Genome] [gcode=11]\n{1}\n'.format(id, '\n'.join([genome[i:i+80] for i in range(0, len(genome), 80)]) ))
    for id, ann in enumerate(annotations) :
        with open(prefix+'_{0}.tbl'.format(id), 'w') as fout :
            fout.write('>Feature {0} Table1\n'.format(id))
            for gid, annotation in enumerate( ann ) :
                if annotation[5] > annotation[4] :
                    if annotation[7] > 0 :
                        fout.write('{0}\t{1}\tgene\n\t\t\tgene\t{2}_{3}\n\t\t\tlocus_tag\t{4}\n{0}\t{1}\tCDS\n\t\t\tproduct\tsimulated_gene\n\t\t\ttransl_table\t11\n'.format(annotation[4], annotation[5], \
                                                                                                                       int(annotation[0]/10), annotation[0]%10, \
                                                                                                                       'G_{0}_{1}'.format(id, gid+1)))
                    else :
                        fout.write('{1}\t{0}\tgene\n\t\t\tgene\t{2}_{3}\n\t\t\tlocus_tag\t{4}\n{1}\t{0}\tCDS\n\t\t\tproduct\tsimulated_gene\n\t\t\ttransl_table\t11\n'.format(annotation[4], annotation[5], \
                                                                                                                       int(annotation[0]/10), annotation[0]%10, \
                                                                                                                       'G_{0}_{1}'.format(id, gid+1)))
                        
    
def runGenomeSim(args) :
    # prepare gene lengths and intergenic region
    items = itemGenerator(args.nBackbone, args.nMobile, 1-args.pBackbone, 1-args.pMobile )
    items = itemShaper(items, args.geneLen, args.igrLen)
    
    # simulate a global tree
    globaltrees, localtrees = treeGenerator(args.prefix, items, args.genomeNum, args.rec, args.recLen)
    seqtrees = localtrees if args.insRec == 1 else globaltrees
    itrees   = localtrees[:10000] if args.insRec == 1 else globaltrees
    dtrees   = itrees     if args.delRec == 1 else globaltrees
    # impose gene content variation on the tree
    genomes = setGeneContent(items, itrees, dtrees, args.tipAccelerate, args.genomeNum, args.aveSize, args.nCore, args.backboneBlock, args.mobileBlock)
    
    # write parameter and gene content
    with open(args.prefix+'.params.log', 'w') as fout :
        fout.write(str(args).replace(', ', ',\n'))
    with open('{0}.gene.content.tsv'.format(args.prefix), 'w') as fout :
        fout.write('#ID\t{0}\n'.format( '\t'.join(np.arange(genomes.shape[1]).astype(str)) ))
        for gId, gene in enumerate(genomes) :
            fout.write('{0}\t{1}\n'.format(gId+1, '\t'.join([ str(items[g, 1])[:-1]+'_'+str(items[g, 1]%10) if g >= 0 else '-' for g in gene ])))
    
    if not args.noSeq :
        # simulate sequences
        alignments, genome_seqs, annotations = getSequences(args.prefix, genomes, items, args.operonBlock, seqtrees, args.idenOrtholog, args.idenParalog, args.idenDuplication, args.indelRate, args.indelMax, args.freqStart, args.freqStop)
        # output
        writeOut(args.prefix, alignments, genome_seqs, annotations)

def genomeSim(a) :
    import argparse
    parser = argparse.ArgumentParser(description='''SimPan is a simulator for bacterial pan-genome. 
Global phylogeny and tree distortions are derived from SimBac and the gene and intergenic sequences are simulated using INDELile.''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-p', '--prefix', help='prefix for all intermediate files and outputs. {DEFAULT: SimPan]', default='SimPan')
    parser.add_argument('--genomeNum', help='No of genome in population [DEFAULT: 20]', default=20, type=int)

    parser.add_argument('--geneLen', help='[negative bionomial with r=2] mean,min,max sizes of genes [DEFAULT: 900,150,6000]', default='900,150,6000')
    parser.add_argument('--igrLen', help='[negative bionomial] mean,min,max sizes of intergenic regions [DEFAULT: 50,0,300]', default='50,0,300')
    parser.add_argument('--backboneBlock', help='[geometric] mean,min,max number of backbone genes per block [DEFAULT: 3,0,30]', default='3,0,30')
    parser.add_argument('--mobileBlock', help='[geometric] mean,min,max number of mobile genes per block [DEFAULT: 10,0,100]', default='10,0,100')
    parser.add_argument('--operonBlock', help='[geometric] mean,min,max number of continuous genes that share the same coding strand [DEFAULT: 3,0,15]', default='3,0,15')

    parser.add_argument('--aveSize', help='average gene number per genome (greater than nBackbone). [DEFAULT: 4500]', default=4500, type=int)
    parser.add_argument('--nBackbone', help='number of backbone genes (present in common ancestor) per genome. [DEFAULT: 4000]', default=4000, type=int)
    parser.add_argument('--nCore', help='sizea of core gene (smaller than the size of backbone genes). [DEFAULT: 3500]', default=3500, type=int)
    parser.add_argument('--nMobile', help='size of mobile gene pool for accessory genome. [DEFAULT: 20000]', default=20000, type=int)

    parser.add_argument('--pBackbone', help='propotion of paralogs in backbone (core) genes. [DEFAULT: 0.05]', default=0.05, type=float)
    parser.add_argument('--pMobile', help='propotion of paralogs in mobile (accessory) genes. [DEFAULT: 0.4]', default=0.4, type=float)

    parser.add_argument('--tipAccelerate', help='grandient increasing of gene indels in recent times. [DEFAULT: 100]', default=100., type=float)
    parser.add_argument('--rec', help='expected coverage of homoplastic events in pairwise comparisons. [DEFAULT: 0.05]', default=0.05, type=float)
    parser.add_argument('--recLen', help='expected size of homoplastic events. [DEFAULT: 1000]', default=1000, type=float)
    parser.add_argument('--seqRec', help='Use homoplastic events to infer sequences. Use 0 to disable [DEFAULT: 1]', default=1, type=int)
    parser.add_argument('--insRec', help='Use homoplastic events to infer gene insertions. Use 0 to disable [DEFAULT: 1]', default=1, type=int)
    parser.add_argument('--delRec', help='Use homoplastic events to infer gene deletions. Use 0 to disable [DEFAULT: 1]', default=1, type=int)
    
    parser.add_argument('--noSeq', help='Do not infer sequence but only the gene presence/absence. [DEFAULT: False]', default=False, action='store_true')
    
    parser.add_argument('--idenOrtholog', help='average nucleotide identities for orthologous genes. [DEFAULT: 0.98]', default=0.98, type=float)
    parser.add_argument('--idenParalog', help='average nucleotide identities for paralogous genes. [DEFAULT: 0.6]', default=0.6, type=float)
    parser.add_argument('--idenDuplication', help='average nucleotide identities for recent gene duplications. [DEFAULT: 0.995]', default=0.995, type=float)
    
    parser.add_argument('--indelRate', help='summarised average frequency of indel events. [DEFAULT: 0.01]', default=0.01, type=float)
    parser.add_argument('--indelMax', help='maximum size of short indel events within each gene (<=300). [DEFAULT: 30]', default=30, type=float)
    parser.add_argument('--freqStart', help='frequencies of start codons of ATG,GTG,TTG. DEFAULT: 0.83,0.14,0.03', default='0.83,0.14,0.03')
    parser.add_argument('--freqStop', help='frequencies of stop codons of TAA,TAG,TGA. DEFAULT: 0.63,0.08,0.29', default='0.63,0.08,0.29')

    args = parser.parse_args(a)
    args.geneLen = np.round(np.array(args.geneLen.split(',')).astype(float)/3)
    args.igrLen = np.array(args.igrLen.split(',')).astype(float)
    args.backboneBlock = np.array(args.backboneBlock.split(',')).astype(float)
    args.mobileBlock = np.array(args.mobileBlock.split(',')).astype(float)
    args.operonBlock = np.array(args.operonBlock.split(',')).astype(float)
    args.freqStart = ('ATG,GTG,TTG'.split(','), np.array(args.freqStart.split(',')).astype(float))
    args.freqStop = ('TAA,TAG,TGA'.split(','), np.array(args.freqStop.split(',')).astype(float))
    runGenomeSim(args)
    
if __name__ == '__main__' :
    genomeSim(sys.argv[1:])
