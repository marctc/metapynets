#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re
import subprocess
import sys
import redis
import functools
import xml.parsers.expat

from Bio import pairwise2

from mwmatching import maxWeightMatching
from reactions import react2enz


class PNMLParser(object):
    def __init__(self, f):
        self.pnml_file = f
        self._i_matrix = {}
        self._o_matrix = {}

    def _start_element(self, name, attrs):
        if name == 'place':
            self._i_matrix[attrs.values()[0]] = []
        elif name == 'transition':
            self._o_matrix[attrs.values()[0]] = []
        elif name == 'arc':
            if attrs['source'] in self._i_matrix:
                self._i_matrix[attrs['source']].append(attrs['target'])
            elif attrs['source'] in self._o_matrix:
                self._o_matrix[attrs['source']].append(attrs['target'])

    def parse_data(self):
        parser = xml.parsers.expat.ParserCreate()
        parser.StartElementHandler = self._start_element
        parser.Parse(self.pnml_file[0])
    
    def get_matrix(self):
        return self._i_matrix, self._o_matrix
    
    def __repr__(self):
        print "Input"
        for i in self._i_matrix:
            print i, self._i_matrix[i]
        print
        print "Output"
        for i in self._o_matrix:
            print i, self._o_matrix[i]


class ReactionSet(object):
    def __init__(self, i_matrix, o_matrix):
        self._i_matrix = i_matrix
        self._o_matrix = o_matrix
        self._result = []
        self._visited = {}
        self._p_forbidden = []

        for i in self._i_matrix.keys():
            self._visited[i] = False
        for i in self._o_matrix.keys():
            self._visited[i] = False
        for i in self._o_matrix.values():
            self._p_forbidden += i

    def calc_result(self):
        for x in self._i_matrix.keys():
            if x not in self._p_forbidden:
                l = []
                self._calc(l, x)
                for i in self._i_matrix.keys():
                    self._visited[i] = False
                for i in self._o_matrix.keys():
                    self._visited[i] = False
                        
    def _calc(self, l, p):
        self._visited[p] = True
        for x in self._i_matrix[p]:
            if not self._visited[x]:
                self._visited[x] = True
                l.append(x)
                for y in self._o_matrix[x]:
                    if not self._visited[y]:
                        self._calc(l, y)
                self._check_len(l[:])
                l.remove(x)
            
    def _check_len(self, l):
        exists = False
        for i in self._result:
            if set(l).issubset(set(i)):
                exists = True
                break
            if set(l).issuperset(set(i)):
                self._result[self._result.index(i)] = l
                exists = True
                break
        if not exists:
            self._result.append(l)
        
    def get_result(self):
        d_result = {}
        for i in range(len(self._result)):
            d_result[i] = self._result[i]
        return d_result


class Align(object):
    def __init__(self, i_matrix1, i_matrix2, o_matrix1, o_matrix2, r_dict1, r_dict2, gap_pen=0):
        self.r_dict1 = r_dict1.copy()
        self.r_dict2 = r_dict2.copy()
        self.c_dict1 = r_dict1.copy()
        self.c_dict2 = r_dict2.copy()
        self.i_matrix1 = i_matrix1
        self.i_matrix2 = i_matrix2
        self.o_matrix1 = o_matrix1
        self.o_matrix2 = o_matrix2
        self.react1 = set(o_matrix1.keys())
        self.react2 = set(o_matrix2.keys())
        self.react = self.react1.union(self.react2)
        self.gap_pen = gap_pen
        self.align_matrix = {}
        self.score_matrix = {}
        self.r_to_c = {}
        self.c_to_r = {}
        self.sub_matrix = {}
        self.compscore = {}
        self._calc_char_set()
        self._calc_sub_matrix()

    def _calc_sub_matrix(self):
        def _enz_score(a, b):
            b_score = 0
            score = 0
            for i in a:
                s1 = re.split(":", a[0])[1]
                s1 = re.split(".", s1)
                for j in b:
                    s2 = re.split(":", b[0])[1]
                    s2 = re.split(".", s2)
                    if s1 == s2:
                        score = 1
                    elif s1[0:2] == s2[0:2]:
                        score = 0.75
                    elif s1[0:1] == s2[0:1]:
                        score = 0.5
                    elif s1[0] == s2[0]:
                        score = 0.25
                    if score >= b_score:
                        b_score = score
            return score

        def _comp_score(c1, c2):
            bi_graph = []
            score_comp = {}
            n1 = 0
            for i in c1:
                n2 = len(c1)
                for j in c2:
                   #simcomp = calc_score_cache(i, j)
                    simcomp = calc_score(i, j)
                    bi_graph.append((n1, n2, simcomp))
                    score_comp[n1, n2] = simcomp
                    n2 += 1
                n1 += 1
            match = maxWeightMatching(bi_graph)
            i = 0
            n = 0
            score = 0
            while i < len(match) and n < min(len(c1), len(c2)):
                if match[i] != -1:
                    try:
                        score += score_comp[i, match[i]]
                    except:
                        pass
                    n += 1
                i += 1
            return score/max(len(c1), len(c2))
                    
        def _get_ocomp(om, r):
            ocomp = []
            for c in om[r]:
                ocomp.append(c.split(":C")[1])
            return ocomp

        def _get_icomp(im, r):
            l_items = im.items()
            icomp = []
            for c in l_items:
                if r in c[1]:
                    icomp.append(c[0].split(":C")[1])
            return icomp
            
        # wsdl_url = 'http://soap.genome.jp/KEGG.wsdl'
        # wsdl = WSDL.Proxy(wsdl_url)
        react_enz = {}
        react_icomp1 = {}
        react_icomp2 = {}
        react_ocomp1 = {}
        react_ocomp2 = {}
        for r in self.react1:
            react_icomp1[r] = _get_icomp(self.i_matrix1, r)
            react_ocomp1[r] = _get_ocomp(self.o_matrix1, r)
        for r in self.react2:
            react_icomp2[r] = _get_icomp(self.i_matrix2, r)
            react_ocomp2[r] = _get_ocomp(self.o_matrix2, r)
        for r in self.react:
            #react_enz[r] = wsdl.get_enzymes_by_reaction(r.split('#')[0])
            react_enz[r] = react2enz[r.split(' ')[0].split('#')[0]]
        for i in self.react1:
            a = self.r_to_c[i]
            e1 = react_enz[i]
            ic1 = react_icomp1[i]
            oc1 = react_ocomp1[i]
            for j in self.react2:
                b = self.r_to_c[j] 
                if i == j:
                    score = 1
                else:
                    e2 = react_enz[j]
                    ic2 = react_icomp2[j]
                    oc2 = react_ocomp2[j]
                    e_score = _enz_score(e1, e2)
                    ic_score = _comp_score(ic1, ic2)
                    oc_score = _comp_score(oc1, oc2)
                    score = e_score * 0.4 + ic_score * 0.3 + oc_score * 0.3
                self.sub_matrix[a, b] = score
                self.sub_matrix[b, a] = score
                
    def _calc_char_set(self):
        n = 0
        f = lambda x : re.sub('^'+i+'$', c, x)
        for i in self.react: 
            c = chr(32+n)
            self.r_to_c[i] = c
            self.c_to_r[c] = i
            for j in self.c_dict1:
                self.c_dict1[j] = map(f, self.c_dict1[j])
            for j in self.c_dict2:
                self.c_dict2[j] = map(f, self.c_dict2[j])
            n += 1
            if c == 45:
                n += 1
        for i in self.c_dict1:
            self.c_dict1[i] = "".join(self.c_dict1[i])
        for i in self.c_dict2:
            self.c_dict2[i] = "".join(self.c_dict2[i])

    def align(self):
        desp = len(self.c_dict1) + 1
        for i in self.c_dict1:
            a = self.c_dict1[i]
            for j in self.c_dict2:
                r0 = ""
                r1 = ""
                b = self.c_dict2[j]
                align = pairwise2.align.localds(a, b, self.sub_matrix, self.gap_pen, self.gap_pen)
                if not align:
                    score = 0
                else:
                    align = align[0]
                    for c in align[0]:
                        if c != '-':
                            r0 += self.c_to_r[c]
                        else:
                            r0 += c
                    for c in align[1]:
                        if c != '-':
                            r1 += self.c_to_r[c]
                        else:
                            r1 += c
                    align = list(align)
                    align[0] = r0
                    align[1] = r1
                    score = align[2]/max(len(a), len(b))
                self.align_matrix[i, j+desp] = align
                self.score_matrix[i, j+desp] = score
                    
    def get_score(self):
        bi_graph = []
        score = 0
        # v0, v1, v2
        # [2, -1, 0]
        for i in self.score_matrix:
            l = list(i)
            l.append(self.score_matrix[i])
            bi_graph.append(tuple(l))
        r = maxWeightMatching(bi_graph)
        for b in r:
            a = r.index(b)
            if self.score_matrix.has_key((a, b)):
                s = self.score_matrix[a, b]
                score += s
        if max(len(self.r_dict1), len(self.r_dict2)) == 0:
            return 0.0
        else:
            return score/max(len(self.r_dict1), len(self.r_dict2))


def redis_memoize_compound(func):
    """Needs to be associative"""
    @functools.wraps(func)
    def memoize(x, y):
        key = "compound"
        key1 = "{2}:{0}:{1}".format(x, y, key)
        key2 = "{2}:{1}:{0}".format(x, y, key)
        if REDIS.exists(key1):
            val = REDIS.get(key1)
        elif REDIS.exists(key2):
            val = REDIS.get(key2)
        else:
            val = func(x, y)
            REDIS.set(key1, val)
            REDIS.set(key2, val)
        return float(val)
    return memoize


@redis_memoize_compound
def calc_score(i, j):
    if i == j:
        return 1.0
    simcomp_call = subprocess.check_output(["./simcomp", "compound", "-e", i, "-f", j])
    simcomp = float(simcomp_call.split()[5])
    return simcomp


def preprocess_data(filename):
    f = open(filename).readlines()
    p = PNMLParser(f)
    p.parse_data()
    i, o = p.get_matrix()
    s = ReactionSet(i, o)
    s.calc_result()
    r = s.get_result()
    return i, o, r


def redis_memoize_file(func):
    @functools.wraps(func)
    def memoize(x, y):
        key = "file:{0}:{1}".format(x, y)
        print key
        if REDIS.exists(key):
            val = REDIS.get(key)
        else:
            val = func(x, y)
            REDIS.set(key, val)
        return float(val)
    return memoize


@redis_memoize_file
def calc_score_files(file1, file2):
    i0, o0, r0 = preprocess_data(file1)
    i1, o1, r1 = preprocess_data(file2)

    a = Align(i0, i1, o0, o1, r0, r1)
    a.align()
    score = a.get_score()
    return score

if __name__ == '__main__':
#   global REDIS
    REDIS = redis.Redis()
    try:
        file1, file2 = sys.argv[1], sys.argv[2]
        score = calc_score_files(file1, file2)
        print "similitude score between: {0} and {1}".format(file1, file2)
        print(score)
    except:
        print "Usage: python {0} [FILE1] [FILE2]".format(sys.argv[0])

