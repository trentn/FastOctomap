baseline:
	gcc benchmarks.c FastOctree.c -o benchmarks -lm

converted:
	gcc benchmarks.c convertedFastOctree.c -o benchmarks -lm

clean:
	rm benchmarks