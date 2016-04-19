Alignment
=========

Implementation of Needleman-Wunsch and Hirschberg's Algorithm.

Usage
-----
	from alignment import Needleman, Hirschberg
	seqa = list('12345678')
	seqb = list('123478901')
	
	# Align using Needleman-Wunsch algorithm.
	n = Needleman()
	a,b = n.align(seqa, seqb)
	print a
	print b
	print

	# Align using Hirschberg's algorithm.
	h = Hirschberg()
	a,b = h.align(seqa,seqb)
	print a
	print b
	print

	# Score the alignment, the higher the score is,
	# the better the two sequences align.
	score = h.score(a, b)
	print score
	
	Output:
	12345678|||
	1234||78901

	12345678|||
	1234||78901
	
	# Score.
	20

Memory Usage
------------
Use `/usr/bin/time -v` to test the two algorithms time and memory usage.  
Input: A string of length 1400 and another string of length 1426, both encoded in UTF-8.  
  
Needleman-Wunsch  

	User time (seconds): 5.43
	System time (seconds): 0.04
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:05.48
	Maximum resident set size (kbytes): 33832
Hirschberg's  

	User time (seconds): 8.84
	System time (seconds): 0.00
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:08.89
	Maximum resident set size (kbytes): 3660

A Trick to Speed up
-------------------
Identify the maximum units which you want to align between original and altered text. For example, `[a-zA-Z]+|[0-9]+|\s+|[.,;!\(\)]+`. The alignment might be a little different from the one aligning every character, but it suffices my need and takes much less time. Use Hirschberg's algorithm with this trick on the same two strings,

	 User time (seconds): 1.19
	 System time (seconds): 0.00
	 Percent of CPU this job got: 99%
	 Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.20
	 Maximum resident set size (kbytes): 3608

Reference
---------
* http://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
* http://en.wikipedia.org/wiki/Hirschberg%27s_algorithm
* http://www.akira.ruc.dk/~keld/teaching/algoritmedesign_f03/Artikler/05/Hirschberg75.pdf
