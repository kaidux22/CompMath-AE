all: CompMath

CompMath: sofa.a main.o
	g++ main.o sofa/libsofa_c.a -o CompMath

main.o: main.cpp main.h Math/ComplexNums.h
	g++ -c main.cpp

sofa.a:
	$(MAKE) -C ./sofa

clean:
	rm -rf *.o
	$(MAKE) clean -C ./sofa
