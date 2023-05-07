all: CompMath

CompMath: sofa.a main.o GravPot.o Converter.o Vnm.o
	g++ main.o GravPot.o Converter.o Vnm.o sofa/libsofa_c.a -o CompMath

main.o: main.cpp main.h Math/ComplexNums.h 
	g++ -c main.cpp

GravPot.o: Formule/GravPot.cpp Formule/GravPot.h Math/LegFunc.h Math/ComplexNums.h Formule/DerivativesV.h Formule/Converter.h
	g++ -c Formule/GravPot.cpp
 
Converter.o: Formule/Converter.cpp Formule/Converter.h
	g++ -c Formule/Converter.cpp

Vnm.o: Formule/DerivativesV.cpp Formule/DerivativesV.h Math/ComplexNums.h Math/LegFunc.h
	g++ -c Formule/Vnm.cpp

sofa.a:
	$(MAKE) -C ./sofa

clean:
	rm -rf *.o
	$(MAKE) clean -C ./sofa
