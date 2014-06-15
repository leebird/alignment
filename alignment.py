# -*- coding: utf-8 -*-
import operator

class Alignment(object):
    def __init__(self,seqa,seqb):
        self.seqa = seqa
        self.seqb = seqb
        self.lena = len(self.seqa)
        self.lenb = len(self.seqb)
        self.scoreNull = 5
        self.scoreSub = -100
        self.scoreDel = -2
        self.scoreIns = -1
        self.separator = '|'

    def set_score(self,scoreNull=None,scoreSub=None,scoreDel=None,scoreIns=None):
        if scoreNull is not None:
            self.scoreNull = scoreNull
        if scoreSub is not None:
            self.scoreSub = scoreSub
        if scoreDel is not None:
            self.scoreDel = scoreDel
        if scoreIns is not None:
            self.scoreIns = scoreIns

    def match(self,a,b):
        if a == b:
            return self.scoreNull
        else:
            return self.scoreSub
    
    def delete(self,a):
        return self.scoreDel

    def insert(self,a):
        return self.scoreIns

    def score(self,Z,W):
        score = 0
        for a,b in zip(Z,W):            
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
    
    def map_alignment(self,Z,W):
        wToz = []
        idx = 0
        for x,y in zip(Z,W):
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
    def __init__(self,*args):
        super(Needleman,self).__init__(*args)
        
    def init_matrix(self):
        rows = self.lena + 1
        cols = self.lenb + 1
        self.matrix = [[0]*cols for i in xrange(rows)]

    def compute_matrix(self):
        seqa = self.seqa
        seqb = self.seqb
        lena = self.lena
        lenb = self.lenb
        
        for i in xrange(1,lena+1):
            self.matrix[i][0] = self.scoreDel * i
        
        for i in xrange(1,lenb+1):
            self.matrix[0][i] = self.scoreIns * i

        for i in xrange(1,lena+1):
            for j in xrange(1,lenb+1):
                scoreSub = self.matrix[i-1][j-1] + self.match(seqa[i-1],seqb[j-1])
                scoreDel = self.matrix[i-1][j] + self.delete(seqa[i-1])
                scoreIns = self.matrix[i][j-1] + self.insert(seqb[j-1])                
                self.matrix[i][j] = max(scoreSub,scoreDel,scoreIns)

    def backtrack(self):
        Z = ''
        W = ''
        seqa = self.seqa
        seqb = self.seqb
        i = self.lena
        j = self.lenb
        mat = self.matrix
        
        while i > 0 or j > 0:
            if (i > 0 and j > 0 and
                mat[i][j] == mat[i-1][j-1] + self.match(seqa[i-1],seqb[j-1])):

                Z = seqa[i-1] + Z
                W = seqb[j-1] + W
                i -= 1
                j -= 1

            elif j > 0 and mat[i][j] == mat[i][j-1] + self.insert(seqb[j-1]):
                Z = self.separator + Z
                W = seqb[j-1] + W
                j -= 1                

            elif i > 0 and mat[i][j] == mat[i-1][j] + self.delete(seqa[i-1]):
                Z = seqa[i-1] + Z
                W = self.separator + W
                i -= 1            

        return Z,W
    
    def align(self):
        self.init_matrix()
        self.compute_matrix()
        return self.backtrack()



class Hirschberg(Alignment):
    def __init__(self,*args):
        super(Hirschberg,self).__init__(*args)

    def last_row(self,seqa,seqb):
        lena = len(seqa)
        lenb = len(seqb)
        preRow = [0] * (lenb+1)
        curRow = [0] * (lenb+1)

        for j in xrange(1,lenb+1):
            preRow[j] = preRow[j-1] + self.insert(seqb[j-1])
        
        for i in xrange(1,lena+1):
            curRow[0] = i * self.delete(seqa[i-1])
            for j in xrange(1,lenb+1):
                scoreSub = preRow[j-1] + self.match(seqa[i-1],seqb[j-1])
                scoreDel = preRow[j] + self.delete(seqa[i-1])
                scoreIns = curRow[j-1] + self.insert(seqb[j-1])                
                curRow[j] = max(scoreSub,scoreDel,scoreIns)

            preRow = curRow
            curRow = [0] * (lenb+1)

        return preRow

    def align_rec(self,seqa,seqb):
        Z = ''
        W = ''
        lena = len(seqa)
        lenb = len(seqb)
        
        if lena == 0:
            for i in xrange(lenb):
                Z += self.separator
                W += seqb[i]
        elif lenb == 0:
            for i in xrange(lena):
                Z += seqa[i]
                W += self.separator

        elif len(seqa) == 1:
            needleman = Needleman(seqa,seqb)
            Z,W = needleman.align()

        else:
            mida = lena / 2
            rowLeft = self.last_row(seqa[:mida],seqb)            
            rowRight = self.last_row(seqa[mida:][::-1],seqb[::-1])
            rowRight.reverse()
            
            row = [l+r for l,r in zip(rowLeft,rowRight)]
            maxIdx,maxVal = max(enumerate(row),key=operator.itemgetter(1))

            midb = maxIdx
            
            ZL,WL = self.align_rec(seqa[:mida],seqb[:midb])            
            ZR,WR = self.align_rec(seqa[mida:],seqb[midb:])
            Z = ZL + ZR
            W = WL + WR

        return Z,W
            
    def align(self):
        return self.align_rec(self.seqa,self.seqb)
                
