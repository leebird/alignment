# -*- coding: utf-8 -*-
import operator


class Alignment(object):
    SCORE_UNIFORM = 1
    SCORE_PROPORTION = 2

    def __init__(self, seqa, seqb, global_optimal=False):
        self.seqa = seqa
        self.seqb = seqb
        self.lena = len(self.seqa)
        self.lenb = len(self.seqb)
        self.scoreNull = 5
        self.scoreSub = -100
        self.scoreDel = -2
        self.scoreIns = -1
        self.separator = u'|'
        self.mode = Alignment.SCORE_UNIFORM
        self.global_optimal = global_optimal

    def set_score(self, scoreNull=None, scoreSub=None, scoreDel=None, scoreIns=None):
        if scoreNull is not None:
            self.scoreNull = scoreNull
        if scoreSub is not None:
            self.scoreSub = scoreSub
        if scoreDel is not None:
            self.scoreDel = scoreDel
        if scoreIns is not None:
            self.scoreIns = scoreIns

    def match(self, a, b):
        if a == b and self.mode == Alignment.SCORE_UNIFORM:
            return self.scoreNull
        elif self.mode == Alignment.SCORE_UNIFORM:
            return self.scoreSub
        elif a == b:
            return self.scoreNull * len(a)
        else:
            return self.scoreSub * len(a)

    def delete(self, a):
        """
        deleted elements are on seqa
        """
        if self.mode == Alignment.SCORE_UNIFORM:
            return self.scoreDel
        return self.scoreDel * len(a)

    def insert(self, a):
        """
        inserted elements are on seqb
        """
        if self.mode == Alignment.SCORE_UNIFORM:
            return self.scoreIns
        return self.scoreIns * len(a)

    def score(self, Z, W):
        score = 0
        for a, b in zip(Z, W):
            if a == b:
                score += self.scoreNull
            else:
                if a == self.separator:
                    score += self.scoreIns
                elif b == self.separator:
                    score += self.scoreDel
                else:
                    score += self.scoreSub
        return score

    def map_alignment(self, Z, W):
        wToz = []
        idx = 0
        for x, y in zip(Z, W):
            if x == y:
                wToz.append(idx)
                idx += 1
            elif x == self.separator:
                wToz.append(idx)
            elif y == self.separator:
                idx += 1
                continue
        return wToz


class Needleman(Alignment):
    def __init__(self, *args):
        super(Needleman, self).__init__(*args)

    def init_matrix(self):
        rows = self.lena + 1
        cols = self.lenb + 1
        self.matrix = [[0] * cols for i in range(rows)]

    def compute_matrix(self):
        seqa = self.seqa
        seqb = self.seqb
        lena = self.lena
        lenb = self.lenb

        for i in range(1, lena + 1):
            self.matrix[i][0] = self.delete(seqa[i - 1]) + self.matrix[i - 1][0]

        for i in range(1, lenb + 1):
            self.matrix[0][i] = self.insert(seqb[i - 1]) + self.matrix[0][i - 1]

        for i in range(1, lena + 1):
            for j in range(1, lenb + 1):
                """
                Note that rows = lena+1, cols = lenb+1
                """
                scoreSub = self.matrix[i - 1][j - 1] + self.match(seqa[i - 1], seqb[j - 1])
                scoreDel = self.matrix[i - 1][j] + self.delete(seqa[i - 1])
                scoreIns = self.matrix[i][j - 1] + self.insert(seqb[j - 1])
                self.matrix[i][j] = max(scoreSub, scoreDel, scoreIns)

    def backtrack(self):
        alignSeqa, alignSeqb = [], []
        seqa, seqb = self.seqa, self.seqb
        i, j = self.lena, self.lenb
        mat = self.matrix

        while i > 0 or j > 0:
            if (i > 0 and j > 0 and
                        mat[i][j] == mat[i - 1][j - 1] + self.match(seqa[i - 1], seqb[j - 1])):
                alignSeqa.insert(0, seqa[i - 1])
                alignSeqb.insert(0, seqb[j - 1])
                i -= 1
                j -= 1

            elif j > 0 and mat[i][j] == mat[i][j - 1] + self.insert(seqb[j - 1]):
                alignSeqa.insert(0, self.separator * len(seqb[j - 1]))
                alignSeqb.insert(0, seqb[j - 1])
                j -= 1

            elif i > 0 and mat[i][j] == mat[i - 1][j] + self.delete(seqa[i - 1]):
                alignSeqa.insert(0, seqa[i - 1])
                alignSeqb.insert(0, self.separator * len(seqa[i - 1]))
                i -= 1
            else:
                print(seqa)
                print(seqb)
                print(alignSeqa)
                print(alignSeqb)
                print(mat)
                raise Exception('backtrack error', i, j, seqa[i - 2:i + 1], seqb[j - 2:j + 1])
                pass

        return alignSeqa, alignSeqb

    def align(self, mode=None):
        if mode is not None:
            self.mode = mode
        self.init_matrix()
        self.compute_matrix()
        return self.backtrack()


class Hirschberg(Alignment):
    def __init__(self, *args):
        super(Hirschberg, self).__init__(*args)

    def last_row(self, seqa, seqb):
        lena = len(seqa)
        lenb = len(seqb)
        pre_row = [0] * (lenb + 1)
        cur_row = [0] * (lenb + 1)

        for j in range(1, lenb + 1):
            pre_row[j] = pre_row[j - 1] + self.insert(seqb[j - 1])

        for i in range(1, lena + 1):
            cur_row[0] = self.delete(seqa[i - 1]) + pre_row[0]
            for j in range(1, lenb + 1):
                score_sub = pre_row[j - 1] + self.match(seqa[i - 1], seqb[j - 1])
                score_del = pre_row[j] + self.delete(seqa[i - 1])
                score_ins = cur_row[j - 1] + self.insert(seqb[j - 1])
                cur_row[j] = max(score_sub, score_del, score_ins)

            pre_row = cur_row
            cur_row = [0] * (lenb + 1)

        return pre_row

    def align_rec(self, seq_a, seq_b):
        aligned_a, aligned_b = [], []
        len_a, len_b = len(seq_a), len(seq_b)

        if len_a == 0:
            for i in range(len_b):
                aligned_a.append(self.separator * len(seq_b[i]))
                aligned_b.append(seq_b[i])
        elif len_b == 0:
            for i in range(len_a):
                aligned_a.append(seq_a[i])
                aligned_b.append(self.separator * len(seq_a[i]))

        elif len(seq_a) == 1:
            needleman = Needleman(seq_a, seq_b)
            aligned_a, aligned_b = needleman.align()

        else:
            mid_a = int(len_a / 2)
            mid_b = int(len_b / 2)
            diff = abs(len_a - len_b) + 1

            if self.global_optimal:
                rowleft = self.last_row(seq_a[:mid_a], seq_b)
                rowright = self.last_row(seq_a[mid_a:][::-1], seq_b[::-1])
            else:
                rowleft = self.last_row(seq_a[:mid_a], seq_b[:mid_b + diff])
                rowright = self.last_row(seq_a[mid_a:][::-1], seq_b[len_b - mid_b - diff::-1])

            rowright.reverse()

            row = [l + r for l, r in zip(rowleft, rowright)]
            maxidx, maxval = max(enumerate(row), key=operator.itemgetter(1))

            mid_b = maxidx

            aligned_a_left, aligned_b_left = self.align_rec(seq_a[:mid_a], seq_b[:mid_b])
            aligned_a_right, aligned_b_right = self.align_rec(seq_a[mid_a:], seq_b[mid_b:])
            aligned_a = aligned_a_left + aligned_a_right
            aligned_b = aligned_b_left + aligned_b_right

        return aligned_a, aligned_b

    def align(self, mode=None):
        if mode is not None:
            self.mode = mode
        return self.align_rec(self.seqa, self.seqb)


if __name__ == '__main__':
    seqa = list('12345678')
    seqb = list('123478901')
    seqa = list(
        'Antiangiogenic role of miR-361 in human umbilical vein endothelial cells: functional interaction with the peptide somatostatin . Somatostatin (SRIF) acts as antiangiogenic factor, but its role in the regulation of microRNAs (miRNAs) targeting proangiogenic factors is unknown. We used human umbilical vein endothelial cells (HUVEC) to investigate whether (1) miRNAs targeting proangiogenic factors are influenced by hypoxia,(2) their expression is regulated by SRIF, and (3) SRIF-regulated miRNAs affect HUVEC angiogenic phenotype. The involvement of signal transducer and activator of transcription ~@~(STAT) 3! and hypoxia inducible factor (HIF) -1 in miRNA effects was studied. Quantitative real-time PCR, Western blot, cell proliferation assays, and enzyme-linked immunosorbent assay (ELISA) were used. Using specific algorithms, three miRNAs (miR-17, miR-18b, and miR-361) were predicted to bind angiogenesis-associated factors including STAT3, HIF-1alpha, and vascular endothelial growth factor (VEGF). Hypoxia downregulates miR-17 and miR-361 without affecting miR-18b . SRIF restored decreased levels of miR-361 acting at the SRIF receptor sst (1). Downregulated miR-361 was also restored by HIF-1alpha inhibition with YC-1. Combined application of SRIF did not influence YC-1-induced miR-361 downregulation, suggesting that YC-1 and SRIF modulate miR-361 through a common mechanism involving HIF-1alpha . This possibility was confirmed by the result that HIF-1alpha activation in normoxia-downregulated miR-361 and that this downregulation was prevented by SRIF. miR-361 overexpression reduced hypoxia-induced cell proliferation and VEGF release indicating miR-361 involvement in the acquisition of an angiogenic phenotype by HUVEC. miR-361 effects on VEGF were enhanced by the coadministration of SRIF. Our results suggest that (1) SRIF regulates miR-361 expression through a control on HIF-1,(2) miR-361 affects HUVEC angiogenic phenotype, and (3) SRIF and miR-361 act cooperatively in limiting hypoxia-induced VEGF release.')
    seqb = list(
        'Antiangiogenic role of miR-361 in human umbilical vein endothelial cells: functional interaction with the peptide somatostatin. Somatostatin (SRIF) acts as antiangiogenic factor, but its role in the regulation of microRNAs (miRNAs) targeting proangiogenic factors is unknown. We used human umbilical vein endothelial cells (HUVEC) to investigate whether (1) miRNAs targeting proangiogenic factors are influenced by hypoxia, (2) their expression is regulated by SRIF, and (3) SRIF-regulated miRNAs affect HUVEC angiogenic phenotype. The involvement of signal transducer and activator of transcription (STAT) 3 and hypoxia inducible factor (HIF)-1 in miRNA effects was studied. Quantitative real-time PCR, Western blot, cell proliferation assays, and enzyme-linked immunosorbent assay (ELISA) were used. Using specific algorithms, three miRNAs (miR-17, miR-18b, and miR-361) were predicted to bind angiogenesis-associated factors including STAT3, HIF-1α, and vascular endothelial growth factor (VEGF). Hypoxia downregulates miR-17 and miR-361 without affecting miR-18b. SRIF restored decreased levels of miR-361 acting at the SRIF receptor sst(1). Downregulated miR-361 was also restored by HIF-1α inhibition with YC-1. Combined application of SRIF did not influence YC-1-induced miR-361 downregulation, suggesting that YC-1 and SRIF modulate miR-361 through a common mechanism involving HIF-1α. This possibility was confirmed by the result that HIF-1α activation in normoxia-downregulated miR-361 and that this downregulation was prevented by SRIF. miR-361 overexpression reduced hypoxia-induced cell proliferation and VEGF release indicating miR-361 involvement in the acquisition of an angiogenic phenotype by HUVEC. miR-361 effects on VEGF were enhanced by the coadministration of SRIF. Our results suggest that (1) SRIF regulates miR-361 expression through a control on HIF-1, (2) miR-361 affects HUVEC angiogenic phenotype, and (3) SRIF and miR-361 act cooperatively in limiting hypoxia-induced VEGF release.')
    # n = Needleman(seqa, seqb)
    # a, b = n.align()
    # print(a)
    # print(b)

    h = Hirschberg(seqa, seqb, False)
    a, b = h.align()
    h = Hirschberg(seqa, seqb, True)
    a1, b1 = h.align()

    print(a == a1, b == b1)
    print(''.join(a))
    print(''.join(a1))
    print(''.join(b))
    print(''.join(b1))

    # print(a)
    # print(b)