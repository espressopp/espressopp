PYTHON  = python
RST2S5  = ${PYTHON} ./rst2s5.py
RST2S5_FLAGS =
RST2S5_FLAGS +=--cloak-email-addresses
RST2S5_FLAGS +=--current-slide
#RST2S5_FLAGS +=--visible-controls

OUTDIR = slides

.PHONY: slides
slides:
	${RM} -r ${OUTDIR} && mkdir -p ${OUTDIR}
	${RST2S5} ${RST2S5_FLAGS} mpi4py-part1.rst ${OUTDIR}/mpi4py-part1.html
	${RST2S5} ${RST2S5_FLAGS} mpi4py-part2.rst ${OUTDIR}/mpi4py-part2.html
	cp logo.png ${OUTDIR}
	cp clockspeeds.jpg ${OUTDIR}
	unzip -q sendrecv1.png.zip -d ${OUTDIR}
	unzip -q sendrecv2.png.zip -d ${OUTDIR}

.PHONY: clean
clean:
	${RM} -r *.py[co] ${OUTDIR}
