# -*- coding: utf-8 -*-
import operator


class Alignment(object):
    SCORE_UNIFORM = 1
    SCORE_PROPORTION = 2

    def __init__(self, *args):
        self.seq_a = None
        self.seq_b = None
        self.len_a = None
        self.len_b = None
        self.score_null = 5
        self.score_sub = -100
        self.score_del = -2
        self.score_ins = -1
        self.separator = u'|'
        self.mode = Alignment.SCORE_PROPORTION

    def set_score(self, score_null=None, score_sub=None, score_del=None, score_ins=None):
        if score_null is not None:
            self.score_null = score_null
        if score_sub is not None:
            self.score_sub = score_sub
        if score_del is not None:
            self.score_del = score_del
        if score_ins is not None:
            self.score_ins = score_ins

    def match(self, a, b):
        if a == b and self.mode == Alignment.SCORE_UNIFORM:
            return self.score_null
        elif self.mode == Alignment.SCORE_UNIFORM:
            return self.score_sub
        elif a == b:
            return self.score_null * len(a)
        else:
            return self.score_sub * len(a)

    def delete(self, a):
        """
        deleted elements are on seqa
        """
        if self.mode == Alignment.SCORE_UNIFORM:
            return self.score_del
        return self.score_del * len(a)

    def insert(self, a):
        """
        inserted elements are on seqb
        """
        if self.mode == Alignment.SCORE_UNIFORM:
            return self.score_ins
        return self.score_ins * len(a)

    def score(self, Z, W):
        score = 0
        for a, b in zip(Z, W):
            if a == b:
                score += self.score_null
            else:
                if a == self.separator:
                    score += self.score_ins
                elif b == self.separator:
                    score += self.score_del
                else:
                    score += self.score_sub
        return score

    def map_alignment(self, aligned_seq_a, aligned_seq_b):
        map_b2a = []
        idx = 0
        for x, y in zip(aligned_seq_a, aligned_seq_b):
            if x == y:
                map_b2a.append(idx)
                idx += 1
            elif x == self.separator:
                map_b2a.append(idx)
            elif y == self.separator:
                idx += 1
                continue
        return map_b2a


class Needleman(Alignment):
    def __init__(self, *args):
        super(Needleman, self).__init__(*args)

    def init_matrix(self):
        rows = self.len_a + 1
        cols = self.len_b + 1
        self.matrix = [[0] * cols for i in range(rows)]

    def compute_matrix(self):
        seq_a = self.seq_a
        seq_b = self.seq_b
        len_a = self.len_a
        len_b = self.len_b

        for i in range(1, len_a + 1):
            self.matrix[i][0] = self.delete(seq_a[i - 1]) + self.matrix[i - 1][0]

        for i in range(1, len_b + 1):
            self.matrix[0][i] = self.insert(seq_b[i - 1]) + self.matrix[0][i - 1]

        for i in range(1, len_a + 1):
            for j in range(1, len_b + 1):
                """
                Note that rows = len_a+1, cols = len_b+1
                """

                score_sub = self.matrix[i - 1][j - 1] + self.match(seq_a[i - 1], seq_b[j - 1])
                score_del = self.matrix[i - 1][j] + self.delete(seq_a[i - 1])
                score_ins = self.matrix[i][j - 1] + self.insert(seq_b[j - 1])
                self.matrix[i][j] = max(score_sub, score_del, score_ins)

    def backtrack(self):
        aligned_seq_a, aligned_seq_b = [], []
        seq_a, seq_b = self.seq_a, self.seq_b
        i, j = self.len_a, self.len_b
        mat = self.matrix
        while i > 0 or j > 0:
            # from end to start, choose insert/delete over match for a tie
            if j > 0 and mat[i][j] == mat[i][j - 1] + self.insert(seq_b[j - 1]):
                aligned_seq_a.insert(0, self.separator * len(seq_b[j - 1]))
                aligned_seq_b.insert(0, seq_b[j - 1])
                j -= 1

            elif i > 0 and mat[i][j] == mat[i - 1][j] + self.delete(seq_a[i - 1]):
                aligned_seq_a.insert(0, seq_a[i - 1])
                aligned_seq_b.insert(0, self.separator * len(seq_a[i - 1]))
                i -= 1

            elif (i > 0 and j > 0 and
                          mat[i][j] == mat[i - 1][j - 1] + self.match(seq_a[i - 1], seq_b[j - 1])):
                aligned_seq_a.insert(0, seq_a[i - 1])
                aligned_seq_b.insert(0, seq_b[j - 1])
                i -= 1
                j -= 1

            else:
                print(seq_a)
                print(seq_b)
                print(aligned_seq_a)
                print(aligned_seq_b)
                print(mat)
                raise Exception('backtrack error', i, j, seq_a[i - 2:i + 1], seq_b[j - 2:j + 1])
                pass

        return aligned_seq_a, aligned_seq_b

    def align(self, seq_a, seq_b, mode=None):
        self.seq_a = seq_a
        self.seq_b = seq_b
        self.len_a = len(self.seq_a)
        self.len_b = len(self.seq_b)
        if mode is not None:
            self.mode = mode
        self.init_matrix()
        self.compute_matrix()
        return self.backtrack()


class Hirschberg(Alignment):
    def __init__(self, *args):
        super(Hirschberg, self).__init__(*args)
        self.needleman = Needleman()

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
            aligned_a, aligned_b = self.needleman.align(seq_a, seq_b)

        else:
            mid_a = int(len_a / 2)

            rowleft = self.last_row(seq_a[:mid_a], seq_b)
            rowright = self.last_row(seq_a[mid_a:][::-1], seq_b[::-1])

            rowright.reverse()

            row = [l + r for l, r in zip(rowleft, rowright)]
            maxidx, maxval = max(enumerate(row), key=operator.itemgetter(1))

            mid_b = maxidx

            aligned_a_left, aligned_b_left = self.align_rec(seq_a[:mid_a], seq_b[:mid_b])
            aligned_a_right, aligned_b_right = self.align_rec(seq_a[mid_a:], seq_b[mid_b:])
            aligned_a = aligned_a_left + aligned_a_right
            aligned_b = aligned_b_left + aligned_b_right

        return aligned_a, aligned_b

    def align(self, seq_a, seq_b, mode=None):
        self.seq_a = seq_a
        self.seq_b = seq_b
        self.len_a = len(self.seq_a)
        self.len_b = len(self.seq_b)
        if mode is not None:
            self.mode = mode
        return self.align_rec(self.seq_a, self.seq_b)


class SegmentAlignment(object):
    step = 50

    def __init__(self):
        pass

    @classmethod
    def align(cls, seq_a, seq_b):
        len_a = len(seq_a)
        len_b = len(seq_b)
        diff = abs(len_a - len_b)

        curr_a = 0
        curr_b = 0

        h = Hirschberg()
        aligned_a = []
        aligned_b = []

        insert_length = 0
        while curr_a < len_a and curr_b < len_b:
            sub_seq_a = seq_a[curr_a:curr_a + cls.step]
            sub_seq_b = seq_b[curr_b:curr_b + cls.step + diff]


            # print(curr_a, curr_b, insert_length)
            # print(''.join(sub_seq_a), ''.join(sub_seq_b), sep='\n')
            # print()
            # print(seq_a, sub_seq_b)
            aligned_sub_a, aligned_sub_b = h.align(sub_seq_a, sub_seq_b)
            # print(''.join(aligned_sub_a), ''.join(aligned_sub_b), sep='\n')
            # print()


            insert_length = 0
            for char in aligned_sub_a[::-1]:
                if char == '|':
                    insert_length += 1
                else:
                    break

            aligned_a += aligned_sub_a[:len(aligned_sub_a) - insert_length]
            aligned_b += aligned_sub_b[:len(aligned_sub_b) - insert_length]

            curr_a += len(sub_seq_a)
            curr_b += len(sub_seq_b) - insert_length

        if curr_b < len_b:
            aligned_a += ['|'] * (len_b - curr_b)
            aligned_b += seq_b[curr_b:]
        else:
            aligned_a += seq_a[curr_a:]
            aligned_b += ['|'] * (len_a - curr_a)

        return aligned_a, aligned_b


def test():
    import os

    h = Hirschberg()
    s = SegmentAlignment()

    for root, _, files in os.walk('data/raw'):
        for f in files:
            raw_text_file = open(os.path.join(root, f), 'r')
            raw_text = raw_text_file.read()
            raw_text_file.close()

            altered_text_file = open(os.path.join('data/altered', f), 'r')
            altered_text = altered_text_file.read()
            altered_text_file.close()

            ha, hb = h.align(list(raw_text), list(altered_text))

            sa, sb = s.align(list(raw_text), list(altered_text))

            print(f, ha == sa, hb == sb)

            res = open('data/aligned/' + f, 'w')
            res.write(''.join(ha) + '\n' + ''.join(hb) + '\n')
            res.write(''.join(sa) + '\n' + ''.join(sb))
            res.close()

if __name__ == '__main__':
    # seqa = list('12345678')
    # seqb = list('123478908')
    # seqa = list(
    # 'Antiangiogenic role of miR-361 in human umbilical vein endothelial cells: functional interaction with the peptide somatostatin . Somatostatin (SRIF) acts as antiangiogenic factor, but its role in the regulation of microRNAs (miRNAs) targeting proangiogenic factors is unknown. We used human umbilical vein endothelial cells (HUVEC) to investigate whether (1) miRNAs targeting proangiogenic factors are influenced by hypoxia,(2) their expression is regulated by SRIF, and (3) SRIF-regulated miRNAs affect HUVEC angiogenic phenotype. The involvement of signal transducer and activator of transcription ~@~(STAT) 3! and hypoxia inducible factor (HIF) -1 in miRNA effects was studied. Quantitative real-time PCR, Western blot, cell proliferation assays, and enzyme-linked immunosorbent assay (ELISA) were used. Using specific algorithms, three miRNAs (miR-17, miR-18b, and miR-361) were predicted to bind angiogenesis-associated factors including STAT3, HIF-1alpha, and vascular endothelial growth factor (VEGF). Hypoxia downregulates miR-17 and miR-361 without affecting miR-18b . SRIF restored decreased levels of miR-361 acting at the SRIF receptor sst (1). Downregulated miR-361 was also restored by HIF-1alpha inhibition with YC-1. Combined application of SRIF did not influence YC-1-induced miR-361 downregulation, suggesting that YC-1 and SRIF modulate miR-361 through a common mechanism involving HIF-1alpha . This possibility was confirmed by the result that HIF-1alpha activation in normoxia-downregulated miR-361 and that this downregulation was prevented by SRIF. miR-361 overexpression reduced hypoxia-induced cell proliferation and VEGF release indicating miR-361 involvement in the acquisition of an angiogenic phenotype by HUVEC. miR-361 effects on VEGF were enhanced by the coadministration of SRIF. Our results suggest that (1) SRIF regulates miR-361 expression through a control on HIF-1,(2) miR-361 affects HUVEC angiogenic phenotype, and (3) SRIF and miR-361 act cooperatively in limiting hypoxia-induced VEGF release.')
    # seqb = list(
    # 'Antiangiogenic role of miR-361 in human umbilical vein endothelial cells: functional interaction with the peptide somatostatin. Somatostatin (SRIF) acts as antiangiogenic factor, but its role in the regulation of microRNAs (miRNAs) targeting proangiogenic factors is unknown. We used human umbilical vein endothelial cells (HUVEC) to investigate whether (1) miRNAs targeting proangiogenic factors are influenced by hypoxia, (2) their expression is regulated by SRIF, and (3) SRIF-regulated miRNAs affect HUVEC angiogenic phenotype. The involvement of signal transducer and activator of transcription (STAT) 3 and hypoxia inducible factor (HIF)-1 in miRNA effects was studied. Quantitative real-time PCR, Western blot, cell proliferation assays, and enzyme-linked immunosorbent assay (ELISA) were used. Using specific algorithms, three miRNAs (miR-17, miR-18b, and miR-361) were predicted to bind angiogenesis-associated factors including STAT3, HIF-1α, and vascular endothelial growth factor (VEGF). Hypoxia downregulates miR-17 and miR-361 without affecting miR-18b. SRIF restored decreased levels of miR-361 acting at the SRIF receptor sst(1). Downregulated miR-361 was also restored by HIF-1α inhibition with YC-1. Combined application of SRIF did not influence YC-1-induced miR-361 downregulation, suggesting that YC-1 and SRIF modulate miR-361 through a common mechanism involving HIF-1α. This possibility was confirmed by the result that HIF-1α activation in normoxia-downregulated miR-361 and that this downregulation was prevented by SRIF. miR-361 overexpression reduced hypoxia-induced cell proliferation and VEGF release indicating miR-361 involvement in the acquisition of an angiogenic phenotype by HUVEC. miR-361 effects on VEGF were enhanced by the coadministration of SRIF. Our results suggest that (1) SRIF regulates miR-361 expression through a control on HIF-1, (2) miR-361 affects HUVEC angiogenic phenotype, and (3) SRIF and miR-361 act cooperatively in limiting hypoxia-induced VEGF release.')
    # # n = Needleman(seqa, seqb)
    # # a, b = n.align()
    # # print(a)
    # # print(b)
    #
    # # h = Hirschberg()
    # # a, b = h.align(seqa, seqb)
    #
    # s = SegmentAlignment()
    # sa, sb = s.align(seqa, seqb)
    # print(''.join(a), ''.join(b), sep='\n')
    # print(''.join(sa), ''.join(sb), sep='\n')
    # print(a == sa, b == sb)

    # print(a)
    # print(b)

    test()