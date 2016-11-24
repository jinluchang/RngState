all:
	g++ -o example -Wall -O2 example.C
	./example >example.out

clean:
	-rm example
