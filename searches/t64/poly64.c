/*
 * poly64: random search for APN functions of type 64 in univariate form on GF(2^10) (modulus x^10 + x^3 + 1):
 * F(x) = sum c_e x^e over the 25 exponents e = 2^j (2^k + 1), k odd (exactly the quadratic F with F(w x) = F(x)).
 * Each trial: T terms (exponents distinct, random), coefficients random nonzero in GF(2^10).
 * APN test for quadratic F: for every a != 0 the vectors beta(a, e_i) = F(a + e_i) + F(a) + F(e_i) + F(0) span 9 dims.
 * APN hits: fingerprint = differential spectrum of the ortho-derivative (as sw2.c); first hit per fingerprint printed
 * as "K<key> T c_e:e ..." and its LUT on the next line ("L v0,v1,..").
 * usage: poly64 TRIALS TMIN TMAX SEED [Q]   (Q: coefficients in the subfield GF(Q), Q in {4, 32, 1024})
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 1024
static int lg[N], ex[2 * N];
static int gmul(int a, int b) { if (!a || !b) return 0; return ex[lg[a] + lg[b]]; }
static int E[25];
static uint64_t rs;
static uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static int powtab[25][N];
static inline int par(uint32_t x) { return __builtin_parity(x); }
static int ech(int *b, int nb, int v) { for (int t = 0; t < nb; t++) if ((v ^ b[t]) < v) v ^= b[t]; return v; }
static int add_ech(int *b, int nb, int v) {
    v = ech(b, nb, v); if (!v) return nb;
    int i = nb; b[nb++] = v;
    while (i > 0 && b[i] > b[i - 1]) { int t = b[i]; b[i] = b[i - 1]; b[i - 1] = t; i--; }
    return nb;
}
static int is_apn(const int *F) {
    for (int a = 1; a < N; a++) {
        int B[10], nb = 0;
        for (int i = 0; i < 10; i++) nb = add_ech(B, nb, F[a ^ (1 << i)] ^ F[a] ^ F[1 << i] ^ F[0]);
        if (nb != 9) return 0;
    }
    return 1;
}
static uint64_t fingerprint(const int *F) {
    static int pi[N], cnt[N]; static uint32_t hist[N + 1];
    pi[0] = 0;
    for (int a = 1; a < N; a++) {
        int B[10], nb = 0;
        for (int i = 0; i < 10; i++) nb = add_ech(B, nb, F[a ^ (1 << i)] ^ F[a] ^ F[1 << i] ^ F[0]);
        int c = 0; for (int l = 1; l < N; l++) { int ok = 1; for (int t = 0; t < nb; t++) if (par(l & B[t])) { ok = 0; break; } if (ok) { c = l; break; } }
        pi[a] = c;
    }
    memset(hist, 0, sizeof hist);
    for (int a = 1; a < N; a++) { memset(cnt, 0, sizeof cnt); for (int x = 0; x < N; x++) cnt[pi[x ^ a] ^ pi[x]]++; for (int b = 0; b < N; b++) hist[cnt[b]]++; }
    uint64_t h = 1469598103934665603ULL; for (int i = 0; i <= N; i++) if (hist[i]) { h ^= (uint64_t)i * 1000003u + hist[i]; h *= 1099511628211ULL; }
    return h;
}
#define MAXK 100000
static uint64_t seen[MAXK]; static int nseen;

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: poly64 TRIALS TMIN TMAX SEED\n"); return 1; }
    int Q = argc > 5 ? atoi(argv[5]) : 1024; int step = 1023 / (Q - 1);
    long trials = atol(argv[1]); int tmin = atoi(argv[2]), tmax = atoi(argv[3]); rs = 0x9E3779B97F4A7C15ULL ^ (strtoull(argv[4], 0, 10) * 2654435761ULL); if (!rs) rs = 1;
    int x = 1; for (int i = 0; i < N - 1; i++) { ex[i] = x; lg[x] = i; x <<= 1; if (x & N) x ^= 0x409; }
    for (int i = N - 1; i < 2 * N; i++) ex[i] = ex[i - (N - 1)];
    int ne = 0;
    for (int k = 1; k <= 9; k += 2) for (int j = 0; j < 10; j++) {
        int e = (int)(((1L << j) * ((1L << k) + 1)) % 1023); int dup = 0; for (int t = 0; t < ne; t++) if (E[t] == e) dup = 1;
        if (!dup) E[ne++] = e;
    }
    if (ne != 25) { fprintf(stderr, "exponents %d\n", ne); return 1; }
    for (int t = 0; t < 25; t++) for (int v = 0; v < N; v++) powtab[t][v] = v ? ex[(long)lg[v] * E[t] % 1023] : 0;
    long hits = 0;
    static int F[N];
    for (long tr = 0; tr < trials; tr++) {
        int T = tmin + (int)(rnd() % (tmax - tmin + 1));
        int idx[25], cf[25];
        for (int t = 0; t < T; t++) { int e; do { e = (int)(rnd() % 25); int d = 0; for (int s = 0; s < t; s++) if (idx[s] == e) d = 1; if (!d) break; } while (1); idx[t] = e; cf[t] = 0; }
        for (int t = 0; t < T; t++) cf[t] = t == 0 ? 1 : ex[(int)(rnd() % (Q - 1)) * step];
        for (int v = 0; v < N; v++) { int s = 0; for (int t = 0; t < T; t++) s ^= gmul(cf[t], powtab[idx[t]][v]); F[v] = s; }
        if (!is_apn(F)) continue;
        hits++;
        uint64_t key = fingerprint(F); int sn = 0; for (int t = 0; t < nseen; t++) if (seen[t] == key) { sn = 1; break; }
        if (sn || nseen >= MAXK) continue;
        seen[nseen++] = key;
        printf("K%016llx %d", (unsigned long long)key, T); for (int t = 0; t < T; t++) printf(" %d:%d", cf[t], E[idx[t]]); printf("\nL");
        for (int v = 0; v < N; v++) printf("%s%d", v ? "," : "", F[v]); printf("\n"); fflush(stdout);
    }
    fprintf(stderr, "trials %ld, APN %ld, fingerprints %d\n", trials, hits, nseen);
    return 0;
}
