all: comp_math

comp_math: sofa.a main.o Observatories.o DormandPrince.o GravPot.o Vnm.o Converter.o
	g++ main.o GravPot.o  Observatories.o DormandPrince.o Vnm.o Converter.o sofa/libsofa_c.a -o comp_math

Converter.o: Converter.cpp Converter.h
	g++ -c Converter.cpp

Vnm.o: Vnm.h Vnm.cpp ComplexNums.h LegFunc.h
	g++ -c Vnm.cpp

GravPot.o: GravPot.h GravPot.cpp Converter.h Vnm.h sofa/sofa.h
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