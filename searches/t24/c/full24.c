/* Full classification of type 24 (A = B = J5; quadratic functions on F_{2^(2m)} with coefficients in F_{2^m}).
   Normal form (../form24.py): F(p, q) = (G(p) + S(p, p) + mu(p, q), D(p) + S(q, p) + G(q)), G = x^3 (or x^e),
   mu its polar, S bilinear, D quadratic (polar Dp alternating).
   S is taken modulo the translations T = {mu(q, Y p) + Y' mu(q, p)} and satisfies the diag conditions
   h_b(S(b, b)) = 1 (APN in the directions (0, b)): S = S0 + span(C), C a complement of T in the diag solution space.
   For each S a DFS over Dp: level k chooses Dp(e_k, e_j), j > k; the active directions (a, b) with top bit of a = k
   keep an echelon basis of beta((a, b), X) on X = span(e_1..e_k, chosen e_j) x F^m, whose rank must stay dim X - 1.
   usage: full24 mode [args]   mode "all": every S (small m);  mode "rand n seed": n random S;  env GE = Gold exponent.
   Layout: x = p | q << M (p low), value = top | bot << M (top low), as in form24.py. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#ifndef M
#define M 5
#endif
#define Q (1 << M)
#define NV (M * M * M)
typedef unsigned __int128 u128;
static int POLY;
static int G[Q], mu[Q][Q], hb[Q];
static int St[Q][Q];             /* S(q, p) table */
static uint8_t dpv[M][M];         /* Dp(e_i, e_j) */
static long long nodes[M + 2], leaves, apn, nS, reach[M + 1]; static int reached[M + 1];
static int fm(int a, int b) { int r = 0; while (b) { if (b & 1) r ^= a; b >>= 1; a <<= 1; if (a & Q) a ^= POLY; } return r; }
static int fpow(int a, int e) { int r = 1; while (e) { if (e & 1) r = fm(r, a); a = fm(a, a); e >>= 1; } return r; }
static int par(unsigned x) { return __builtin_parity(x); }
static u128 B1(int i) { return (u128)1 << i; }
static int bit(u128 v, int i) { return (int)((v >> i) & 1); }
static int vidx(int i, int j, int k) { return (i * M + j) * M + k; }
/* ---------- linear algebra on 128-bit vectors ---------- */
static int low(u128 v) { return (uint64_t)v ? __builtin_ctzll((uint64_t)v) : 64 + __builtin_ctzll((uint64_t)(v >> 64)); }
static int par128(u128 v) { return __builtin_parityll((uint64_t)v) ^ __builtin_parityll((uint64_t)(v >> 64)); }
static int ech_ins(u128 *bs, int *piv, int n, u128 v) {   /* reduce v by bs (pivot = lowest set bit); insert if nonzero */
    for (int i = 0; i < n; i++) if (bit(v, piv[i])) v ^= bs[i];
    if (!v) return n;
    int p = low(v);
    for (int i = 0; i < n; i++) if (bit(bs[i], p)) bs[i] ^= v;
    bs[n] = v; piv[n] = p; return n + 1;
}
static u128 S0, Cb[NV]; static int nC;
static void setup_space(void) {
    /* diag rows: for b != 0, linear functional on S bits; rhs 1 */
    u128 rows[Q]; int nr = 0;
    for (int b = 1; b < Q; b++) {
        u128 r = 0;
        for (int i = 0; i < M; i++) if (b >> i & 1) for (int j = 0; j < M; j++) if (b >> j & 1)
            for (int k = 0; k < M; k++) if (hb[b] >> k & 1) r ^= B1(vidx(i, j, k));
        rows[nr++] = r;
    }
    /* solve rows . S = 1: RREF with rhs bit at position NV */
    u128 R[Q]; int piv[Q], rk = 0;
    for (int i = 0; i < nr; i++) {
        u128 v = rows[i] | B1(NV);
        for (int t = 0; t < rk; t++) if (bit(v, piv[t])) v ^= R[t];
        u128 w = v & (B1(NV) - 1);
        if (!w) { if (v) { fprintf(stderr, "diag system inconsistent\n"); exit(1); } continue; }
        int p = 0; while (!bit(w, p)) p++;
        for (int t = 0; t < rk; t++) if (bit(R[t], p)) R[t] ^= v;
        R[rk] = v; piv[rk++] = p;
    }
    /* particular solution: free vars 0 */
    S0 = 0; for (int t = 0; t < rk; t++) if (bit(R[t], NV)) S0 |= B1(piv[t]);
    /* nullspace basis W */
    static u128 W[NV]; int nw = 0; int isp[NV] = {0}; for (int t = 0; t < rk; t++) isp[piv[t]] = 1;
    for (int f = 0; f < NV; f++) if (!isp[f]) {
        u128 v = B1(f); for (int t = 0; t < rk; t++) if (bit(R[t], f)) v |= B1(piv[t]);
        W[nw++] = v;
    }
    /* T */
    static u128 T[2 * M * M]; int nt = 0;
    for (int r = 0; r < M; r++) for (int s = 0; s < M; s++) {
        u128 v = 0, w = 0;
        for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) {
            int y = (j == s) ? mu[1 << i][1 << r] : 0;               /* mu(e_i, Y e_j), Y = E_rs */
            int z = (mu[1 << i][1 << j] >> s & 1) ? (1 << r) : 0;     /* Y' mu(e_i, e_j) */
            for (int k = 0; k < M; k++) { if (y >> k & 1) v |= B1(vidx(i, j, k)); if (z >> k & 1) w |= B1(vidx(i, j, k)); }
        }
        T[nt++] = v; T[nt++] = w;
    }
    u128 eb[NV]; int ep[NV], ne = 0;
    for (int i = 0; i < nt; i++) {
        for (int t = 0; t < nr; t++) if (par128(rows[t] & T[i])) {
            fprintf(stderr, "T not in diag space\n"); exit(1); }
        ne = ech_ins(eb, ep, ne, T[i]);
    }
    int dimT = ne;
    nC = 0; for (int i = 0; i < nw; i++) { int n2 = ech_ins(eb, ep, ne, W[i]); if (n2 > ne) { ne = n2; Cb[nC++] = W[i]; } }
    fprintf(stderr, "m=%d: S bits %d, diag rank %d, dim T %d, dim C %d\n", M, NV, rk, dimT, nC);
}
static void build_S(u128 s) {
    int sv[M][M];
    for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) { int v = 0; for (int k = 0; k < M; k++) v |= bit(s, vidx(i, j, k)) << k; sv[i][j] = v; }
    for (int q = 0; q < Q; q++) for (int p = 0; p < Q; p++) {
        int v = 0; for (int i = 0; i < M; i++) if (q >> i & 1) for (int j = 0; j < M; j++) if (p >> j & 1) v ^= sv[i][j];
        St[q][p] = v;
    }
}
/* image of basis input vectors under beta((a, b), .), with Dp(a, e_j) given separately */
static int img_p(int a, int b, int j, int dz) {  /* input (e_j, 0) */
    int e = 1 << j;
    int top = mu[a][e] ^ St[a][e] ^ St[e][a] ^ mu[b][e];
    int bot = dz ^ St[b][e];
    return top | bot << M;
}
static int img_q(int a, int b, int j) {          /* input (0, e_j) */
    int e = 1 << j;
    return mu[a][e] | (St[e][a] ^ mu[b][e]) << M;
}
static int dpa(int a, int j) { int r = 0; for (int i = 0; i < M; i++) if (a >> i & 1) r ^= dpv[i][j]; return r; }
/* echelon basis of 2M-bit vectors, pivot = highest bit; returns 0 if v dependent */
typedef struct { uint16_t v[2 * M]; uint8_t n; } Ech;
/* reduced echelon basis sorted by decreasing leading bit */
static int ered(const Ech *e, int v) { for (int i = 0; i < e->n; i++) if ((v ^ e->v[i]) < v) v ^= e->v[i]; return v; }
static int emem(const Ech *e, int v) { return ered(e, v) == 0; }
static int eins(Ech *e, int v) {
    v = ered(e, v);
    if (!v) return 0;
    int hv = 31 - __builtin_clz(v);
    for (int i = 0; i < e->n; i++) if (e->v[i] >> hv & 1) e->v[i] ^= v;
    int k = e->n++; while (k > 0 && e->v[k - 1] < v) { e->v[k] = e->v[k - 1]; k--; } e->v[k] = v;
    return 1;
}
static FILE *fout;
static int Dq(int p) { int r = 0; for (int i = 0; i < M; i++) if (p >> i & 1) for (int j = i + 1; j < M; j++) if (p >> j & 1) r ^= dpv[i][j]; return r; }
static void leaf(void) {
    static int F[Q * Q]; leaves++;
    for (int q = 0; q < Q; q++) for (int p = 0; p < Q; p++)
        F[p | q << M] = (G[p] ^ St[p][p] ^ mu[p][q]) | (Dq(p) ^ St[q][p] ^ G[q]) << M;
    int ok = 1; static int cnt[Q * Q];
    for (int d = 1; d < Q * Q && ok; d++) {
        memset(cnt, 0, sizeof cnt);
        for (int x = 0; x < Q * Q; x++) if (++cnt[F[x] ^ F[x ^ d]] > 2) { ok = 0; break; }
    }
    if (ok) apn++;
    if (fout) { fprintf(fout, "%d", ok); for (int x = 0; x < Q * Q; x++) fprintf(fout, " %d", F[x]); fprintf(fout, "\n"); }
}
/* active directions at level k: a = (1 << k) ^ t, t < (1 << k); index t * Q + b */
#define MAXACT ((Q / 2) * Q)
static Ech stk[M + 1][M + 1][MAXACT];
static void dfs(int k, int j, Ech *E) {
    /* level k (0-based), next column j (> k) to choose; E holds echelons for active directions on current X */
    int na = 1 << k;
    if (j == M) {
        if (k == M - 1) { leaf(); return; }
        /* go to level k+1: build initial echelons on X = span(e_0..e_{k+1}) x F^m */
        int k1 = k + 1, na1 = 1 << k1; Ech *E1 = stk[k1][k1 + 1];
        nodes[k1]++;
        for (int t = 0; t < na1; t++) {
            int a = (1 << k1) ^ t;
            for (int b = 0; b < Q; b++) {
                Ech *e = &E1[t * Q + b]; e->n = 0;
                for (int jj = 0; jj < M; jj++) eins(e, img_q(a, b, jj));
                for (int i = 0; i <= k1; i++) eins(e, img_p(a, b, i, dpa(a, i)));
                if (e->n != M + k1) return;
            }
        }
        if (!reached[k1]) { reached[k1] = 1; reach[k1]++; }
        dfs(k1, k1 + 1, E1);
        return;
    }
    /* choose z = Dp(e_k, e_j) */
    unsigned mask = (Q == 32) ? 0xffffffffu : ((1u << Q) - 1);
    int base[MAXACT];
    for (int t = 0; t < na && mask; t++) {
        int a = (1 << k) ^ t; int known = dpa(t, j);
        for (int b = 0; b < Q && mask; b++) {
            Ech *e = &E[t * Q + b];
            int v0 = img_p(a, b, j, known);
            base[t * Q + b] = v0;
            /* forbidden z: red(v0) = red(z << M) (gray code) */
            int r0 = ered(e, v0); int rz[M]; for (int t = 0; t < M; t++) rz[t] = ered(e, 1 << (M + t));
            unsigned forb = 0; int acc = 0;
            for (int z = 0; z < Q; z++) { if (z) acc ^= rz[__builtin_ctz(z)]; if (acc == r0) forb |= 1u << (z ^ (z >> 1)); }
            mask &= ~forb;
        }
    }
    Ech *E2 = stk[k][j + 1];
    while (mask) {
        int z = __builtin_ctz(mask); mask &= mask - 1;
        nodes[k]++;
        dpv[k][j] = dpv[j][k] = z;
        for (int t = 0; t < na; t++) for (int b = 0; b < Q; b++) {
            E2[t * Q + b] = E[t * Q + b];
            eins(&E2[t * Q + b], base[t * Q + b] ^ (z << M));
        }
        dfs(k, j + 1, E2);
    }
    dpv[k][j] = dpv[j][k] = 0;
}
static void run_S(u128 s) {
    nS++;
    build_S(s);
    /* directions (0, b): full check */
    for (int b = 1; b < Q; b++) {
        Ech e; e.n = 0;
        for (int j = 0; j < M; j++) eins(&e, img_q(0, b, j));
        for (int j = 0; j < M; j++) eins(&e, img_p(0, b, j, 0));
        if (e.n != 2 * M - 1) { fprintf(stderr, "diag check failed b=%d\n", b); exit(1); }
    }
    /* level 0: a = e_0, X = span(e_0) x F^m */
    memset(dpv, 0, sizeof dpv);
    Ech *E = stk[0][1];
    nodes[0]++;
    for (int b = 0; b < Q; b++) {
        Ech *e = &E[b]; e->n = 0;
        for (int jj = 0; jj < M; jj++) eins(e, img_q(1, b, jj));
        eins(e, img_p(1, b, 0, 0));
        if (e->n != M) return;
    }
    nodes[M + 1]++;   /* S passing level-0 */
    memset(reached, 0, sizeof reached);
    dfs(0, 1, E);
}
static uint64_t rs = 88172645463325252ull;
static uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
int main(int argc, char **argv) {
    POLY = (M == 3) ? 0xB : (M == 5) ? 0x25 : (M == 7) ? 0x83 : 0;
    int ge = getenv("GE") ? atoi(getenv("GE")) : 3;
    for (int x = 0; x < Q; x++) G[x] = fpow(x, ge);
    for (int a = 0; a < Q; a++) for (int b = 0; b < Q; b++) mu[a][b] = G[a ^ b] ^ G[a] ^ G[b] ^ G[0];
    for (int b = 1; b < Q; b++) {
        Ech e; e.n = 0; for (int y = 0; y < Q; y++) eins(&e, mu[b][y]);
        if (e.n != M - 1) { fprintf(stderr, "G not APN\n"); return 1; }
        hb[b] = 0;
        for (int h = 1; h < Q && !hb[b]; h++) { int ok = 1; for (int y = 0; y < Q; y++) if (par(h & mu[b][y])) { ok = 0; break; } if (ok) hb[b] = h; }
    }
    setup_space();
    if (getenv("OUT")) fout = fopen(getenv("OUT"), "w");
    clock_t c0 = clock();
    if (!strcmp(argv[1], "all")) {
        for (uint64_t c = 0; c < (1ull << nC); c++) {
            u128 s = S0; for (int i = 0; i < nC; i++) if (c >> i & 1) s ^= Cb[i];
            run_S(s);
        }
    } else {
        long n = atol(argv[2]); if (argc > 3) rs ^= atoll(argv[3]) * 0x9E3779B97F4A7C15ull;
        for (long i = 0; i < n; i++) {
            u128 s = S0; for (int t = 0; t < nC; t++) if (rnd() & 1) s ^= Cb[t];
            run_S(s);
        }
    }
    double sec = (double)(clock() - c0) / CLOCKS_PER_SEC;
    printf("S %lld (pass level0 %lld), nodes per level:", nS, nodes[M + 1]);
    for (int k = 0; k < M; k++) printf(" %lld", nodes[k]);
    printf(", S reaching level:"); for (int k = 1; k < M; k++) printf(" %lld", reach[k]);
    printf(", leaves %lld, APN %lld, %.2f s (%.2f us per S)\n", leaves, apn, sec, 1e6 * sec / nS);
    if (fout) fclose(fout);
    return 0;
}
