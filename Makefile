all: clean colormap

colormap:
	mkdir -p ./bin
	cd spCorrection; make -j2; mv spCorrection ../bin/
	cd oeaCorrection; make -j2; mv oeaCorrection ../bin/; cp runBWAIndex.sh runBWAMemPacbio.sh runMinia.sh ../bin/

clean:
	cd spCorrection; make clean
	cd oeaCorrection; make clean

deps:
	mkdir -p ./bin
	cd utils/utils; make clean; make -j2; mv fastUtils ../../bin/
	cd utils/bwa; make clean; make -j2; mv bwa-proovread ../../bin/
	cd utils/samtools; make clean-all; make -j2; mv samtools ../../bin/
	cd utils/minia; rm -rf build; mkdir build; cd build; cmake .. -DSKIP_DOC=1; make -j2; mv bin/minia ../../../bin/
