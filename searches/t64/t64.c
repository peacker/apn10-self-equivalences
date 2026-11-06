/*
 * t64: classification search for type 64 (F(omega x) = F(x) + linear on F_2^10 = F_4^5).
 * F is given by a 10-dim F_2-space V of Hermitian 5x5 matrices over F_4 (components x^* H x). F is APN iff every
 * nonzero vector a of F_4^5 is killed (H a = 0) by exactly one nonzero H in V (herm.py checks this against the
 * dual Hermitian rank-metric code). Kernels are F_4-subspaces, so this is a condition on the 341 points of PG(4, 4).
 *
 * Search: V_1 = <H_1>, H_1 = diag(0, 1, 1, 1, 1) (every rank-4 form is congruent to it; the classes with no element
 * of rank 4 are excluded by Schmidt's forced rank distribution, and checked separately). At level k, p = the first
 * point not killed by V_k; the new basis element is the unique H in V with H p = 0, so it runs over K_p = {H : Hp = 0}
 * (16-dim) with the condition H q not in V_k q for every point q already killed (then no point is killed twice;
 * uncovered points are automatically killed at most once). Every V containing H_1 is produced exactly once.
 * At dimension 10 every point is killed (dim V + dim K_q = 26 > 25), so every leaf is APN.
 *
 * usage: t64 -knuth PROBES [-seed S]     estimate the tree size
 *        t64 -dfs [-limit NODES] [-shard I K -sdepth D]   list leaves (10 matrices as 25-bit codes per line)
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NP 341
static const int F4M[4][4] = {{0,0,0,0},{0,1,2,3},{0,2,3,1},{0,3,1,2}};
static const int CONJ[4] = {0, 1, 3, 2};
static int pts[NP];                         /* normalised points: 5 F_4 coords, 2 bits each */
static int ptidx[1024];                     /* vector -> point index */
static uint16_t img[NP][25];                /* image of basis matrix b at point q: (E_b q) as 10-bit vector */
static uint16_t tab[NP][4][128];            /* H q by 7-bit chunks of the 25-bit code (25 = 7+7+7+4) */

static int coord(int v, int i) { return (v >> (2 * i)) & 3; }

/* basis matrix b (0..24) -> 5x5 F_4 matrix */
static void basis_mat(int b, int M[5][5]) {
    memset(M, 0, sizeof(int) * 25);
    if (b < 5) { M[b][b] = 1; return; }
    int k = 5;
    for (int i = 0; i < 5; i++) for (int j = i + 1; j < 5; j++) {
        if (b == k) { M[i][j] = 1; M[j][i] = 1; return; }          /* h = 1 */
        if (b == k + 1) { M[i][j] = 2; M[j][i] = 3; return; }      /* h = omega */
        k += 2;
    }
}
static int matvec(int M[5][5], int v) {
    int r = 0;
    for (int i = 0; i < 5; i++) { int s = 0; for (int j = 0; j < 5; j++) s ^= F4M[M[i][j]][coord(v, j)]; r |= s << (2 * i); }
    return r;
}
static inline int hq(uint32_t H, int q) {
    return tab[q][0][H & 127] ^ tab[q][1][(H >> 7) & 127] ^ tab[q][2][(H >> 14) & 127] ^ tab[q][3][(H >> 21) & 15];
}

static void setup(void) {
    int np = 0;
    memset(ptidx, -1, sizeof ptidx);
    for (int v = 1; v < 1024; v++) {
        int lead = 0;
        for (int i = 0; i < 5; i++) if (coord(v, i)) { lead = coord(v, i); break; }
        if (lead != 1) continue;
        pts[np] = v;
        for (int c = 1; c < 4; c++) { int w = 0; for (int i = 0; i < 5; i++) w |= F4M[c][coord(v, i)] << (2 * i); ptidx[w] = np; }
        np++;
    }
    if (np != NP) { fprintf(stderr, "points %d\n", np); exit(1); }
    for (int b = 0; b < 25; b++) { int M[5][5]; basis_mat(b, M); for (int q = 0; q < NP; q++) img[q][b] = (uint16_t)matvec(M, pts[q]); }
    for (int q = 0; q < NP; q++) for (int c = 0; c < 4; c++) for (int x = 0; x < 128; x++) {
        int r = 0; for (int t = 0; t < 7; t++) if ((x >> t) & 1 && 7 * c + t < 25) r ^= img[q][7 * c + t];
        tab[q][c][x] = (uint16_t)r;
    }
}

/* state per level */
static uint32_t basisV[10];
static int dimV;
static uint8_t covered[NP];                 /* number of nonzero elements of V killing the point (0 or 1) */
static uint8_t Sq[11][NP][128];             /* bitmap (1024 bits) of V_k q for each point, per level */
static int order_cov[11][NP], ncov[11];     /* covered points per level */

static void kernel_basis(int p, uint32_t *K, int *nk) {   /* basis of {H : H p = 0} */
    int pimg[10]; uint32_t ppre[10]; int has[10] = {0};   /* pivot by leading bit */
    *nk = 0;
    for (int b = 0; b < 25; b++) {
        int x = img[p][b]; uint32_t pr = 1u << b;
        for (int t = 9; t >= 0; t--) if ((x >> t) & 1) {
            if (has[t]) { x ^= pimg[t]; pr ^= ppre[t]; }
            else { has[t] = 1; pimg[t] = x; ppre[t] = pr; x = -1; break; }
        }
        if (x == 0) K[(*nk)++] = pr;
    }
}

static long long nodes, leaves, node_limit = -1;
static int shard = 0, nshards = 1, sdepth = 3; static long long scount = 0;
static int knuth_mode = 0; static double est[12];
static uint64_t rs = 88172645463325252ULL;
static uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

static void add_element(int k, uint32_t H) {   /* V_{k+1} = V_k + <H>; update covered and Sq */
    basisV[k] = H; dimV = k + 1;
    memcpy(Sq[k + 1], Sq[k], sizeof Sq[k]);
    ncov[k + 1] = ncov[k]; memcpy(order_cov[k + 1], order_cov[k], sizeof(int) * ncov[k]);
    /* new elements H + X: their kernels; and Sq update: V_{k+1} q = V_k q + {Hq} */
    for (int q = 0; q < NP; q++) {
        int h = hq(H, q);
        uint8_t *S = Sq[k + 1][q], *S0 = Sq[k][q];
        for (int x = 0; x < 1024; x++) if (S0[x >> 3] >> (x & 7) & 1) { int y = x ^ h; S[y >> 3] |= 1 << (y & 7); }
        if (!covered[q] && (S0[h >> 3] >> (h & 7) & 1)) { covered[q] = (uint8_t)(k + 1); order_cov[k + 1][ncov[k + 1]++] = q; }
    }
}
static void remove_element(int k) {
    for (int t = ncov[k]; t < ncov[k + 1]; t++) covered[order_cov[k + 1][t]] = 0;
    dimV = k;
}

static int mrv = 0;
/* candidates at point p; stop after cap (returns cap) */
static int cands_at(int k, int p, uint32_t *out, int cap) {
    uint32_t K[25]; int nk; kernel_basis(p, K, &nk);
    int nc = 0; int *cov = order_cov[k]; int ncv = ncov[k];
    uint32_t H = 0;
    for (uint32_t g = 0; g < (1u << nk); g++) {
        if (g) H ^= K[__builtin_ctz(g)];       /* Gray code */
        if (!H) continue;
        int ok = 1;
        for (int t = 0; t < ncv; t++) { int q = cov[t]; int h = hq(H, q); if (Sq[k][q][h >> 3] >> (h & 7) & 1) { ok = 0; break; } }
        if (ok) { if (out) out[nc] = H; if (++nc >= cap) return nc; }
    }
    return nc;
}
static int candidates(int k, uint32_t *out) {  /* all H in K_p (p = chosen uncovered point) with Hq not in V_k q for covered q */
    int p = 0; while (p < NP && covered[p]) p++;
    if (p == NP) return -1;
    if (mrv) {                                 /* most constrained uncovered point; 0 candidates anywhere: dead */
        int best = 1 << 30, bp = p;
        for (int q = p; q < NP; q++) if (!covered[q]) {
            int c = cands_at(k, q, 0, best);
            if (c < best) { best = c; bp = q; if (!c) return 0; }
        }
        p = bp;
    }
    return cands_at(k, p, out, 1 << 30);
}

static uint32_t *candbuf[11];

static void print_leaf(void) {
    for (int i = 0; i < 10; i++) printf("%u%c", basisV[i], i == 9 ? '\n' : ' ');
    fflush(stdout);
}

static void dfs(int k) {
    if (node_limit >= 0 && nodes >= node_limit) return;
    nodes++;
    if (k == 10) { leaves++; print_leaf(); return; }
    int nc = candidates(k, candbuf[k]);
    for (int c = 0; c < nc; c++) {
        if (k + 1 == sdepth && nshards > 1 && scount++ % nshards != shard) continue;
        add_element(k, candbuf[k][c]);
        dfs(k + 1);
        remove_element(k);
    }
}

int main(int argc, char **argv) {
    long long probes = 0;
    for (int a = 1; a < argc; a++) {
        if (!strcmp(argv[a], "-knuth")) { knuth_mode = 1; probes = atoll(argv[++a]); }
        else if (!strcmp(argv[a], "-seed")) rs ^= strtoull(argv[++a], 0, 10) * 0x9E3779B97F4A7C15ULL;
        else if (!strcmp(argv[a], "-dfs")) knuth_mode = 0;
        else if (!strcmp(argv[a], "-mrv")) mrv = 1;
        else if (!strcmp(argv[a], "-limit")) node_limit = atoll(argv[++a]);
        else if (!strcmp(argv[a], "-shard")) { shard = atoi(argv[++a]); nshards = atoi(argv[++a]); }
        else if (!strcmp(argv[a], "-sdepth")) sdepth = atoi(argv[++a]);
    }
    setup();
    for (int k = 0; k < 11; k++) candbuf[k] = malloc(sizeof(uint32_t) << 16);
    memset(Sq[0], 0, sizeof Sq[0]);
    for (int q = 0; q < NP; q++) Sq[0][q][0] = 1;      /* V_0 q = {0} */
    /* H_1 = diag(0, 1, 1, 1, 1): kills exactly point 0 = (1, 0, 0, 0, 0) */
    uint32_t H1 = (1u << 1) | (1u << 2) | (1u << 3) | (1u << 4);
    if (ptidx[1] != 0) { fprintf(stderr, "point order\n"); return 1; }
    add_element(0, H1);
    if (ncov[1] != 1) { fprintf(stderr, "H1 kills %d points\n", ncov[1]); return 1; }
    if (knuth_mode) {
        double leafest = 0;
        for (long long r = 0; r < probes; r++) {
            double w = 1; int k = 1;
            while (k < 10) {
                int nc = candidates(k, candbuf[k]);
                if (nc <= 0) break;
                w *= nc; est[k + 1] += w;
                add_element(k, candbuf[k][rnd() % nc]);
                k++;
            }
            if (k == 10) { leafest += w; print_leaf(); }
            for (int j = k - 1; j >= 1; j--) remove_element(j);
        }
        double tot = 1;
        for (int k = 2; k <= 10; k++) { fprintf(stderr, "  level %2d: %.3g\n", k, est[k] / probes); tot += est[k] / probes; }
        fprintf(stderr, "total nodes %.3g, leaves %.3g\n", tot, leafest / probes);
        return 0;
    }
    dfs(1);
    fprintf(stderr, "done: nodes %lld, leaves %lld\n", nodes, leaves);
    return 0;
}
