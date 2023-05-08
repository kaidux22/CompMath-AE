all: CompMath

CompMath: sofa.a main.o GravPot.o Converter.o DerivativesV.o DistanceOrbits.o SupportMath.o LeastSquares.o
	g++ main.o GravPot.o Converter.o DerivativesV.o DistanceOrbits.o SupportMath.o LeastSquares.o sofa/libsofa_c.a -o CompMath

main.o: main.cpp main.h Math/ComplexNums.h Formule/Converter.h Math/Matrix.h
	g++ -c main.cpp

GravPot.o: Formule/GravPot.cpp Formule/GravPot.h Math/LegFunc.h Math/ComplexNums.h Formule/DerivativesV.h Formule/Converter.h
	g++ -c Formule/GravPot.cpp
 
Converter.o: Formule/Converter.cpp Formule/Converter.h
	g++ -c Formule/Converter.cpp

DerivativesV.o: Formule/DerivativesV.cpp Formule/DerivativesV.h Math/ComplexNums.h Math/LegFunc.h Math/SupportMath.h
	g++ -c Formule/DerivativesV.cpp

DistanceOrbits.o: InverseProblem/DistanceOrbits.cpp InverseProblem/DistanceOrbits.h
	g++ -c InverseProblem/DistanceOrbits.cpp

SupportMath.o: Math/SupportMath.h Math/SupportMath.cpp
	g++ -c Math/SupportMath.cpp

LeastSquares.o: InverseProblem/LeastSquares.cpp InverseProblem/LeastSquares.h Math/Matrix.h
	g++ -c InverseProblem/LeastSquares.cpp

sofa.a:
	$(MAKE) -C ./sofa

clean:
	rm -rf *.o
	rm CompMath
	$(MAKE) clean -C ./sofa
