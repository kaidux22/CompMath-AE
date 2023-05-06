all: comp_math

comp_math: sofa.a main.o Observatories.o DormandPrince.o GravPot.o AdvancedFunction.o Derivatives.o Converter.o SLE.o DoubleDerivatives.o
	g++ main.o GravPot.o  Observatories.o DormandPrince.o Derivatives.o AdvancedFunction.o Converter.o DoubleDerivatives.o SLE.o sofa/libsofa_c.a -o comp_math

Converter.o: Converter.cpp Converter.h
	g++ -c Converter.cpp

Derivatives.o: Derivatives.h Derivatives.cpp
	g++ -c Derivatives.cpp

DoubleDerivatives.o: DoubleDerivatives.h DoubleDerivatives.cpp
	g++ -c DoubleDerivatives.cpp

GravPot.o: GravPot.h GravPot.cpp Converter.h Derivatives.h sofa/sofa.h
	g++ -c GravPot.cpp

AdvancedFunction.o: AdvancedFunction.h AdvancedFunction.cpp sofa/sofa.h Converter.h Derivatives.h
	g++ -c AdvancedFunction.cpp

SLE.o: SLE.h SLE.cpp Converter.h
	g++ -c SLE.cpp

DormandPrince.o: DormandPrince.cpp DormandPrince.h GravPot.h AdvancedFunction.h
	g++ -c DormandPrince.cpp

Observatories.o: Observatories.cpp Observatories.h
	g++ -c Observatories.cpp

main.o: main.cpp main.h Observatories.h sofa/sofa.h DormandPrince.h SLE.h
	g++ -c main.cpp

sofa.a: sofa/sofa.h
	$(MAKE) -C ./sofa

clean:
	rm -rf *.o
	$(MAKE) clean -C ./sofa