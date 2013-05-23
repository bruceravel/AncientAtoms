

all:	atoms.info atoms.ps atoms.html

atoms.info:	atoms.texi
	makeinfo atoms.texi

atoms.ps:	atoms.dvi
	dvips -o atoms.ps atoms.dvi

atoms.dvi:	atoms.texi
	texi2dvi atoms.texi

atoms.html:	atoms.texi
	texi2html atoms.texi
