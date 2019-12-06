baseline:
	gcc -std=c99 -O0 benchmarks.c FastOctree.c -o benchmarks -lm

converted:
	gcc -std=c99 -O0 benchmarks.c convertedFastOctree.c -o benchmarks -lm

simd:
	gcc -march=native -std=c99 -O2 benchmarks.c SIMDFastOctree.c -o benchmarks_simd -lm

clean:
	rm benchmarks