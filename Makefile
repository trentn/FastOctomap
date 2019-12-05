baseline:
	gcc -O2 benchmarks.c FastOctree.c -o benchmarks -lm

converted:
	gcc -O2 benchmarks.c convertedFastOctree.c -o benchmarks -lm

clean:
	rm benchmarks