# makefile for version2_thread
# from fileDisposal.o , the command is automatically derived by make

Object =  coordinaryTransform.o    fileDisposal.o    initialization.o \
          main.o    matrixFunction.o    MOPSOAidFunction.o    \
          MOPSOFunction.o    parameterSettings.o \
	setTime.o     preDispose.o  multiplyBetterInput.o\
       	disposedCodes.o	 seiveRep.o  inputParticles.o  createRep.o\
	checkSimilarity.o	Debug.o	freeUpSpace.o	applyVariables.o

GL : $(Object)
	g++ -o GL $(Object) -lpthread

coordinaryTransform.o : MOPSO.h
	g++ -c coordinaryTransform.cpp

fileDisposal.o : MOPSO.h

initialization.o : MOPSO.h

main.o : MOPSO.h

matrixFunction.o : MOPSO.h

MOPSOAidFunction.o : MOPSO.h

MOPSOFunction.o : MOPSO.h

parameterSettings.o : MOPSO.h

setTime.o : MOPSO.h

preDispose.o : MOPSO.h

multiplyBetterInput.o : MOPSO.h

disposedCodes.o : MOPSO.h

seiveRep.o : MOPSO.h

inputParticles.o : MOPSO.h

updateRep.o : MOPSO.h

checkSimilarity.o : MOPSO.h

Debug.o : MOPSO.h

freeUpSpace.o : MOPSO.h

applyVariables.o : MOPSO.h

.PHONY : clear
clear :
	-rm $(Object) GL
.PHONY : delete
delete :
	-rm -r /home/gengling_new/mopso_v1/data/answer/newAnswer

