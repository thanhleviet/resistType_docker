ReformatReads readme by Brian Bushnell
Last updated August 29, 2014.
Please contact me at bbushnell@lbl.gov if you have any questions or encounter any errors.

This is currently a stub.

Change Log:

Moved processing phase from constructor into a new function.
Added pigz (parallel gzip) support, at suggestion of Rob Egan.
Added ftr/ftl (force trim to a certain base position), at request of Alicia Clum.
Added shellscript support for java -Xmx flag (Suggested by James Han).
Fixed stdout.fa.gz writing uncompressed via ReadStreamWriter.
Added scarf input support.
Major refactoring.
Output quality switch (qout) was being ignored and ASCII33 was always used; this has been fixed.
TrimRead.testOptimal() mode added, and made default when quality trimming is performed; old mode can be used with 'otf=f' flag.
Fixed bug in ByteBuilder for reads with no quality values.  Noted by Brian Foster.
All program message information now defaults to stderr.
Added 'tbs' (trim bad sequence) flag to fix broken files from NCBI (like NT).  Works in test cases, but causes problems with NT for unknown reasons. 
Added 'rbb' (requirebothbad) flag for tossing pairs shorter than minlen.  Default: false.
Added 'qfake' flag for quality level of fasta -> fastq reformatting.
Moved parsing to new Parser class.
Fixed random seed setting crash (found by Michael Barton).
Added support for breaking long fastq reads into shorter reads (maxlength and minlength flags).
Added "overrideinterleaved" flag to allow unpaired input when specifying out1 and out2.  Requested by Vasanth Singan.
Added "def" (deleteempty) flag which deletes output files that did not get any reads.  Requested by Vasanth Singan.  TODO: Consider adding to bbmap/bbsplit.
Added "uniquenames" flag.  Allows reads to be renamed to ensure all have unique names.  Important for fasta files that will be mapped to.
Added mhist, bhist, qhist, and rcomp flags.
Fixed by in which qual file was being written to same destination as fasta file. Found by Brian Foster.
Added 'append' flag.
Added 'outsingle' for reads that have a mate that gets discard (generally because they are too short after trimming).
Fixed bug of ignoring "out=null" flag.  Probably existed ever since Parser class was added.
Added "maxlen" flag: discard reads over a certain length.  Old "maxlen" flag became "breaklen" flag.
Added "mingc" and "maxgc".
Fixed a bug in which tabs were not allowed in fasta headers.
Added mappedonly flag.
Added ehist, idhist, indelhist, lhist, gchist.
Added obqhist.

