all: comp_math

comp_math: sofa.a main.o Observatories.o DormandPrince.o GravPot.o Derivatives.o Converter.o
	g++ main.o GravPot.o  Observatories.o DormandPrince.o Derivatives.o Converter.o sofa/libsofa_c.a -o comp_math

Converter.o: Converter.cpp Converter.h
	g++ -c Converter.cpp

Derivatives.o: Derivatives.h Derivatives.cpp LegFunc.h
	g++ -c Derivatives.cpp

GravPot.o: GravPot.h GravPot.cpp Converter.h Derivatives.h sofa/sofa.h
	g++ -c GravPot.cpp

DormandPrince.o: DormandPrince.cpp DormandPrince.h GravPot.h
	g++ -c DormandPrince.cpp

Observatories.o: Observatories.cpp Observatories.h
	g++ -c Observatories.cpp

main.o: main.cpp main.h Observatories.h sofa/sofa.h DormandPrince.h
	g++ -c main.cpp

sofa.a: sofa/sofa.h
	$(MAKE) -C ./sofa

clean:
	rm -rf *.o
	$(MAKE) clean -C ./sofa