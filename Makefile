baseline:
	gcc -std=c99 -O2 benchmarks.c FastOctree.c -o benchmarks -lm -g

converted:
	gcc -std=c99 -O2 benchmarks.c convertedFastOctree.c -o benchmarks -lm -g

simd:
	gcc -march=native -std=c99 -O2 benchmarks.c SIMDFastOctree.c -o benchmarks_simd -lm -g

clean:
	rm benchmarks