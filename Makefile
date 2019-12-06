all: baseline converted simd

baseline:
	gcc -std=c99 -O0 benchmarks.c FastOctree.c -o benchmarks_baseline -lm

converted:
	gcc -std=c99 -O0 benchmarks.c convertedFastOctree.c -o benchmarks_scalar -lm

simd:
	gcc -march=native -std=c99 -O2 benchmarks.c SIMDFastOctree.c -o benchmarks_simd -lm

run:
	./benchmarks_baseline
	./benchmarks_scalar
	./benchmarks_simd

clean:
	rm benchmarks_*