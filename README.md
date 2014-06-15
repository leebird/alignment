Alignment
=========

Implementation of Needleman-Wunsch and Hirschberg's Algorithm.

Usage
-----
	from alignment import Needleman, Hirschberg
	seqa = list('12345678')
	seqb = list('123478901')
	n = Needleman(seqa,seqb)
	a,b = n.align()
	print a
	print b
	print

	h = Hirschberg(seqa,seqb)
	a,b = h.align()
	print a
	print b
	print
^
	Output:

	12345678|||
	1234||78901

	12345678|||
	1234||78901

Reference
---------
* http://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
* http://en.wikipedia.org/wiki/Hirschberg%27s_algorithm
* http://www.akira.ruc.dk/~keld/teaching/algoritmedesign_f03/Artikler/05/Hirschberg75.pdf