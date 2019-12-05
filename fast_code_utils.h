/*
Returns the processors clock value.
The clock ticks at nominal frequency,
which may differ from the clock frequency at which your programs execute.
*/

unsigned long long rdtsc() {
  unsigned a, d;
  __asm__ volatile("rdtsc" : "=a" (a), "=d" (d));
  return ((unsigned long long)a) | (((unsigned long long)d) << 32);
}