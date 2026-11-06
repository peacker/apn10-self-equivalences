/*
 * tri10: exhaustive search for quadratic APN functions of GF(2^10) of the "canonical triplicate" form
 *     F(x) = x^b + sum_{k=1..t} c_k x^(e_k),
 * where b in {3, 9, 33} is the base exponent, the e_k are distinct quadratic exponents 2^i + 2^j divisible by 3
 * (j - i odd; 25 such exponents at n = 10, b excluded) and c_k ranges over the nonzero elements of a subfield K
 * (F_4: 3 values; F_32: 31 values). This extends the search of Budaghyan-Ivkovic-Kaleyski (Cryptogr. Commun. 2023,
 * Sec. 9.2: F_4 coefficients, 5-7 added terms) to more terms and to F_32.
 *
 * APN test for quadratic F: for every a != 0, x -> B(a, x) has rank 9. F is evaluated only at the points needed,
 * from precomputed tables c * x^e, and the test stops at the first failing a. Hits are printed as lookup tables.
 * Sharding: combinations of exponents with index % nshards == shard.
 * usage: tri10 -base 3|9|33 -terms t -field 4|32 [-shard i -nshards k]
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 1024
static int expt[N], logt[N];
static inline int mul(int a, int b) { if (!a || !b) return 0; return expt[(logt[a] + logt[b]) % 1023]; }
static inline int pw(int a, int e) { if (!a) return e ? 0 : 1; return expt[(long)logt[a] * e % 1023]; }

static int ex[32], nex;          /* candidate exponents (excluding the base) */
static int kel[64], nk;          /* nonzero subfield elements */
static uint16_t *T;              /* T[(k * nk + ci) * N + x] = kel[ci] * x^ex[k] */
static uint16_t Bt[N];           /* base term x^b */
static int t, sel[16], csel[16];
static long long candidates, hits;

static inline int Fval(int x) {
    int y = Bt[x];
    for (int q = 0; q < t; q++) y ^= T[((long)sel[q] * nk + csel[q]) * N + x];
    return y;
}

static int is_apn(void) {
    int F0 = Fval(0), Fe[10];
    for (int m = 0; m < 10; m++) Fe[m] = Fval(1 << m);
    for (int a = 1; a < N; a++) {
        int Fa = Fval(a), rows[10], r = 0;
        for (int m = 0; m < 10; m++) rows[m] = Fval(a ^ (1 << m)) ^ Fa ^ Fe[m] ^ F0;
        for (int bit = 9; bit >= 0; bit--) {
            int p = -1;
            for (int q = r; q < 10; q++) if ((rows[q] >> bit) & 1) { p = q; break; }
            if (p < 0) continue;
            int tmp = rows[p]; rows[p] = rows[r]; rows[r] = tmp;
            for (int q = 0; q < 10; q++) if (q != r && ((rows[q] >> bit) & 1)) rows[q] ^= rows[r];
            r++;
        }
        if (r != 9) return 0;
    }
    return 1;
}

/* cheap invariant: differential spectrum of the ortho-derivative pi (pi(a) = the nonzero u orthogonal to the image
 * of x -> B(a, x)); printed as "v:count" pairs. Only the first hit with a given spectrum is printed (per process). */
#define HSIZE (1 << 20)
static uint64_t seen_keys[HSIZE];
static uint64_t spectrum_key(const int *F, char *txt) {
    static int pi[N], cnt[N + 1];
    pi[0] = 0;
    for (int a = 1; a < N; a++) {
        int rows[10];
        for (int m = 0; m < 10; m++) rows[m] = F[a ^ (1 << m)] ^ F[a] ^ F[1 << m] ^ F[0];
        /* u with parity(u & rows[m]) = 0 for all m: null space of the 10x10 matrix M[m][k] = bit k of rows[m] */
        int R[10], piv[10], nr = 0;
        for (int m = 0; m < 10; m++) {
            int vv = rows[m];
            for (int q = 0; q < nr; q++) if ((vv >> piv[q]) & 1) vv ^= R[q];
            if (!vv) continue;
            int pb = 31 - __builtin_clz(vv);
            for (int q = 0; q < nr; q++) if ((R[q] >> pb) & 1) R[q] ^= vv;
            R[nr] = vv; piv[nr++] = pb;
        }
        int pivmask = 0;
        for (int q = 0; q < nr; q++) pivmask |= 1 << piv[q];
        int fcol = 0;
        while ((pivmask >> fcol) & 1) fcol++;   /* the free column (rank 9) */
        int u = 1 << fcol;
        for (int q = 0; q < nr; q++) if ((R[q] >> fcol) & 1) u |= 1 << piv[q];
        pi[a] = u;
    }
    memset(cnt, 0, sizeof cnt);
    static int d[N];
    for (int a = 1; a < N; a++) {
        memset(d, 0, sizeof d);
        for (int x = 0; x < N; x++) d[pi[x ^ a] ^ pi[x]]++;
        for (int b = 0; b < N; b++) cnt[d[b]]++;
    }
    uint64_t h = 1469598103934665603ULL;
    char *p = txt;
    for (int v = 0; v <= N; v++) if (cnt[v]) { h = (h ^ (uint64_t)(v * 1000003 + cnt[v])) * 1099511628211ULL; p += sprintf(p, "%d:%d,", v, cnt[v]); }
    return h ? h : 1;
}

static void coeffs(int q) {
    if (q == t) {
        candidates++;
        if (is_apn()) {
            hits++;
            static int F[N];
            for (int x = 0; x < N; x++) F[x] = Fval(x);
            static char txt[8192];
            uint64_t h = spectrum_key(F, txt);
            uint64_t slot = h & (HSIZE - 1);
            while (seen_keys[slot] && seen_keys[slot] != h) slot = (slot + 1) & (HSIZE - 1);
            if (!seen_keys[slot]) {
                seen_keys[slot] = h;
                printf("%s\t", txt);
                for (int x = 0; x < N; x++) printf("%d%c", F[x], x == N - 1 ? '\n' : ',');
                fflush(stdout);
            }
        }
        return;
    }
    for (int c = 0; c < nk; c++) { csel[q] = c; coeffs(q + 1); }
}

static long long comb_index;
static int shard, nshards;
static void combos(int q, int start) {
    if (q == t) {
        if (comb_index++ % nshards == shard) coeffs(0);
        return;
    }
    for (int k = start; k <= nex - (t - q); k++) { sel[q] = k; combos(q + 1, k + 1); }
}

int main(int argc, char **argv) {
    int base = 3, field = 4;
    t = 5; shard = 0; nshards = 1;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-base")) base = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-terms")) t = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-field")) field = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-shard")) shard = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-nshards")) nshards = atoi(argv[++i]);
    }
    int v = 1;
    for (int i = 0; i < 1023; i++) { expt[i] = v; logt[v] = i; v <<= 1; if (v & N) v ^= 0x409; }
    for (int i = 0; i < 10; i++) for (int j = i + 1; j < 10; j++)
        if ((j - i) % 2 == 1) { int e = (1 << i) + (1 << j); if (e != base) ex[nex++] = e; }
    nk = field - 1;
    for (int i = 0; i < nk; i++) kel[i] = expt[(long)i * (1023 / nk) % 1023];
    T = malloc(sizeof(uint16_t) * nex * nk * N);
    for (int k = 0; k < nex; k++) for (int c = 0; c < nk; c++) for (int x = 0; x < N; x++)
        T[((long)k * nk + c) * N + x] = (uint16_t)mul(kel[c], pw(x, ex[k]));
    for (int x = 0; x < N; x++) Bt[x] = (uint16_t)pw(x, base);
    fprintf(stderr, "tri10: base x^%d, %d candidate exponents, %d added terms, field F_%d, shard %d/%d\n",
            base, nex, t, field, shard, nshards);
    combos(0, 0);
    fprintf(stderr, "done: candidates %lld, APN %lld\n", candidates, hits);
    return 0;
}
