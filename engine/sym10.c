/*
 * sym10: exhaustive search for quadratic APN functions of F_2^n (n <= 10) in a linear space W of quadratic forms,
 * typically the forms with a prescribed self-equivalence (input from sym_setup.py).
 *
 * F is given by beta(e_i, e_j) = sum_k t_k w_k(i, j), t in F_2^d (w_k: reduced row echelon basis of W in search
 * order). Search order: input variables e_1, e_2, .. (level s adds e_s), inside a level the pairs (0, s), (1, s), ..
 * At pair (i, s) the bits t_k whose pivot lies in that pair are chosen (all 2^(#) combinations); the value
 * y = beta(e_i, e_s) is then final. APN on V' = <e_0..e_s> (variable approach): L(a) = beta(a, e_s) must avoid
 * Im_V beta(a, .) for every a != 0 in V = <e_0..e_{s-1}>; L(a) is checked as soon as y for the highest bit i of a is
 * known. A leaf (all levels done) is a quadratic APN function; it is printed as a lookup table
 * F(x) = sum_{i<j} x_i x_j beta(e_i, e_j) (F agrees with the symmetric function up to linear terms).
 * Output normalisation: beta(e_0, e_1) is restricted to the orbit representatives listed after the basis (see
 * sym_setup.py); -nonorm disables it.
 * usage: sym10 -knuth PROBES [-seed S] < W.txt
 *        sym10 -dfs [-print MAX | -dedupe] [-limit NODES] [-shard I K [-sdepth D]] [-pshard R K0 -pdepth D0] [-symdepth P] < W.txt
 * -symdepth P: output-side symmetry reduction (stabiliser chain in C(B)) on the first P pairs; replaces the
 * beta(e_0, e_1) normalisation. Needs the Z(B) basis at the end of W.txt (sym_setup.py).
 * -cp K1 K2: input-side checkpoint pruning when beta is complete on V_k, K1 <= k <= K2 (needs the Z(A) basis);
 * -nsin M: number of sampled input maps per k (default 24).
 * -limit NODES: stop after NODES nodes; the run is then NOT exhaustive: stderr gets "INCOMPLETE: node limit",
 * the -progress file ends with "INCOMPLETE" instead of "DONE", and the exit status is 3.
 * -progress FILE: every 60 s write elapsed time, nodes, estimated explored fraction and time left to FILE.
 * -offset FILE: lines "pair value": affine search space, beta = offset + span(W) (from intertwine.py).
 * -force FILE: lines "pair value": the value of beta on that pair is fixed (e.g. beta on Fix(A) x Fix(A) fixed to a
 * class representative). Combine with -cp only together with -cpfix U, when every forced pair lies in
 * <e_0..e_{U-1}> x <e_0..e_{U-1}>: the input maps are then the identity there, forced pairs are compared for
 * equality (not minimised) in the output canonical form, so only maps that keep the forced prefix are used.
 * -dedupe: print every leaf whose ortho-derivative spectrum fingerprint was not seen before in this process
 * ("K<fingerprint> <lookup table>"); the fingerprint is a cheap invariant, so classify the printed tables further.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXN 10
#define MAXP 45
#define MAXD 460
static int n, N, d, P, pairs_i[MAXP], pairs_j[MAXP], pidx[MAXN][MAXN];
static uint16_t wv[MAXD][MAXP];                    /* basis vectors: value at each pair */
static int pivpair[MAXD];
static uint8_t allowed01[1 << MAXN];               /* allowed values of beta(e_0, e_1) (output normalisation) */
static int forced[MAXP];                            /* -force: forced[p] = value + 1 (0 = free) */
static uint16_t offset[MAXP];                       /* -offset: beta = offset + sum t_k w_k (affine search space) */
static int kfirst[MAXP], kcount[MAXP];             /* basis vectors with pivot in pair p: kfirst[p] .. +kcount[p] */
static uint16_t val[MAXP];                         /* current pair values (sum of chosen t_k w_k) */
static uint16_t hb[1 << MAXN][MAXN];               /* hb[a][j] = beta(a, e_j) for a in the current span */
static uint16_t fb[MAXN][1 << MAXN][MAXN];         /* per level s and a in V: reduced echelon basis of */
static uint8_t fbn[MAXN][1 << MAXN];               /* Im_V beta(a, .) (pivot = highest bit, decreasing) */
static uint16_t Lv[MAXN][1 << MAXN];
static long long nodes, leaves, node_limit = -1;
static int print_max = 0, printed = 0, nonorm = 0, dedupe = 0;
static int shard = 0, nshards = 1, shard_depth = 6;
static long long shard_count = 0;
/* -pshard R K0 -pdepth D0: a first shard level (keep the depth-D0 nodes whose counter is R mod K0), so that one heavy
 * shard of an earlier run (same R, K0, D0) can be split again at a deeper -sdepth */
static int pshard = 0, npshards = 1, pshard_depth = 0;
static long long pshard_count = 0;
static uint64_t rs = 88172645463325252ULL;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static uint64_t rs2 = 0x2545F4914F6CDD1DULL;          /* separate generator for the Knuth probes (the pruning one is */
static inline uint64_t rnd2(void) { rs2 ^= rs2 << 13; rs2 ^= rs2 >> 7; rs2 ^= rs2 << 17; return rs2; }  /* reseeded per node) */

/* forbidden values of L(a) at level s: the subspace Im_V beta(a, .) = span{beta(a, e_j) : j < s}; stored as a
 * reduced echelon basis (O(2^s s^2) per level instead of a 2^s x 2^s table) */
static void level_start(int s) {
    int V = 1 << s;
    for (int a = 1; a < V; a++) {
        uint16_t *B = fb[s][a]; int nb = 0;
        for (int j = 0; j < s; j++) {
            uint16_t v = hb[a][j];
            for (int q = 0; q < nb; q++) if ((v >> (31 - __builtin_clz(B[q]))) & 1) v ^= B[q];
            if (!v) continue;
            int pv = 31 - __builtin_clz(v);
            for (int q = 0; q < nb; q++) if ((B[q] >> pv) & 1) B[q] ^= v;
            int q = nb++;                                  /* insert keeping pivots decreasing */
            while (q > 0 && B[q - 1] < v) { B[q] = B[q - 1]; q--; }
            B[q] = v;
        }
        fbn[s][a] = nb;
    }
    Lv[s][0] = 0;
}

static inline int forbidden(int s, int a, uint16_t L) {
    const uint16_t *B = fb[s][a]; int nb = fbn[s][a];
    for (int q = 0; q < nb; q++) if ((L >> (31 - __builtin_clz(B[q]))) & 1) L ^= B[q];
    return L == 0;
}

static void level_end(int s) {
    int V = 1 << s;
    for (int a = 0; a < V; a++) {
        hb[a][s] = Lv[s][a];
        for (int j = 0; j < s; j++) hb[V | a][j] = hb[a][j] ^ Lv[s][1 << j];
        hb[V | a][s] = Lv[s][a];
    }
}

/* check y = beta(e_i, e_s) at level s: L(a) for a in [2^i, 2^(i+1)) */
static int check(int s, int i, uint16_t y) {
    int lo = 1 << i, hi = 2 << i;
    for (int a = lo; a < hi; a++) {
        uint16_t L = Lv[s][a ^ lo] ^ y;
        if (forbidden(s, a, L)) return 0;
    }
    return 1;
}

static void set_y(int s, int i, uint16_t y) {
    int lo = 1 << i, hi = 2 << i;
    for (int a = lo; a < hi; a++) Lv[s][a] = Lv[s][a ^ lo] ^ y;
}

static void apply_t(int k) {                        /* toggle t_k: val += w_k */
    for (int p = pivpair[k]; p < P; p++) val[p] ^= wv[k][p];
}

static void print_leaf(void) {
    for (int x = 0; x < N; x++) {
        int y = 0;
        for (int p = 0; p < P; p++) if (((x >> pairs_i[p]) & (x >> pairs_j[p])) & 1) y ^= val[p];
        printf("%d%c", y, x == N - 1 ? '\n' : ',');
    }
}


/* cheap class fingerprint of a leaf: differential spectrum of the ortho-derivative pi (pi(a) = the nonzero u
 * orthogonal to the image of x -> beta(a, x)); only the first leaf with a given fingerprint is printed */
#define HSIZE (1 << 20)
static uint64_t seen_keys[HSIZE];
static long long distinct = 0;
static uint64_t leaf_key(void) {
    static int pi[1 << MAXN], cnt[(1 << MAXN) + 1], dd[1 << MAXN];
    pi[0] = 0;
    for (int a = 1; a < N; a++) {
        int R[MAXN], piv[MAXN], nr = 0;
        for (int m = 0; m < n; m++) {
            int vv = hb[a][m];
            for (int q = 0; q < nr; q++) if ((vv >> piv[q]) & 1) vv ^= R[q];
            if (!vv) continue;
            int pb = 31 - __builtin_clz(vv);
            for (int q = 0; q < nr; q++) if ((R[q] >> pb) & 1) R[q] ^= vv;
            R[nr] = vv; piv[nr++] = pb;
        }
        int pivmask = 0;
        for (int q = 0; q < nr; q++) pivmask |= 1 << piv[q];
        int fcol = 0;
        while ((pivmask >> fcol) & 1) fcol++;
        int u = 1 << fcol;
        for (int q = 0; q < nr; q++) if ((R[q] >> fcol) & 1) u |= 1 << piv[q];
        pi[a] = u;
    }
    memset(cnt, 0, sizeof(int) * (N + 1));
    for (int a = 1; a < N; a++) {
        memset(dd, 0, sizeof(int) * N);
        for (int x = 0; x < N; x++) dd[pi[x ^ a] ^ pi[x]]++;
        for (int b = 0; b < N; b++) cnt[dd[b]]++;
    }
    uint64_t h = 1469598103934665603ULL;
    for (int v = 0; v <= N; v++) if (cnt[v]) h = (h ^ (uint64_t)(v * 1000003 + cnt[v])) * 1099511628211ULL;
    return h ? h : 1;
}

/* ---- output-side symmetry: stabiliser chain in the centraliser C(B) ----
 * Z(B) = {X : XB = BX} has basis Mz[0..z-1]. The stabiliser of the prefix values y_0..y_{p-1} in C(B) is the set of
 * invertible X = sum c_j Mz[j] with X y_i = y_i: an affine system in c, kept in reduced echelon form per pair index.
 * At pair p only one value per orbit (under sampled stabiliser elements; union-find, so orbits can only be too fine:
 * always sound) is tried. Used for pairs p < symdepth while the stabiliser is not trivial. */
#define MAXZ 128
typedef unsigned __int128 u128;
static int z = 0, symdepth = 0, nsamp = 12;
static uint16_t Mz[MAXZ][MAXN];
typedef struct { int rank; u128 row[MAXZ]; uint8_t rhs[MAXZ]; int piv[MAXZ]; } Sys;
static Sys sys[MAXP + 1];
static int uf[1 << MAXN];
static long long sym_skipped = 0;

static int ctz128(u128 v) { uint64_t lo = (uint64_t)v; return lo ? __builtin_ctzll(lo) : 64 + __builtin_ctzll((uint64_t)(v >> 64)); }
static int par128(u128 v) { return __builtin_parityll((uint64_t)v) ^ __builtin_parityll((uint64_t)(v >> 64)); }

static void sys_add(Sys *S, u128 coef, int rhs) {
    for (int q = 0; q < S->rank; q++) if ((coef >> S->piv[q]) & 1) { coef ^= S->row[q]; rhs ^= S->rhs[q]; }
    if (!coef) return;                                  /* dependent (consistent: X = I is always a solution) */
    int pv = ctz128(coef);
    for (int q = 0; q < S->rank; q++) if ((S->row[q] >> pv) & 1) { S->row[q] ^= coef; S->rhs[q] ^= rhs; }
    S->row[S->rank] = coef; S->rhs[S->rank] = rhs; S->piv[S->rank] = pv; S->rank++;
}

static void sys_fix_vector(Sys *S, uint16_t y) {      /* add X y = y */
    for (int r = 0; r < n; r++) {
        u128 coef = 0;
        for (int j = 0; j < z; j++) {
            uint16_t im = 0;
            for (int i = 0; i < n; i++) if ((y >> i) & 1) im ^= Mz[j][i];
            if ((im >> r) & 1) coef |= (u128)1 << j;
        }
        sys_add(S, coef, (y >> r) & 1);
    }
}

static int rank10(const uint16_t *cols) {
    uint16_t r[MAXN]; memcpy(r, cols, sizeof(uint16_t) * n);
    int k = 0;
    for (int b = 0; b < n; b++) {
        int pv = -1;
        for (int i = k; i < n; i++) if ((r[i] >> b) & 1) { pv = i; break; }
        if (pv < 0) continue;
        uint16_t t = r[pv]; r[pv] = r[k]; r[k] = t;
        for (int i = 0; i < n; i++) if (i != k && ((r[i] >> b) & 1)) r[i] ^= r[k];
        k++;
    }
    return k;
}

/* deterministic randomness: the samples used for pruning at a node depend only on the prefix (pairs 0..p-1) and the
 * base seed, so that all shards make identical pruning decisions above the shard depth (otherwise the shards would
 * not partition the same tree). */
static uint64_t base_seed = 0x9E3779B97F4A7C15ULL;
static void reseed(int p) {
    uint64_t h = base_seed ^ (uint64_t)(p + 1) * 0xBF58476D1CE4E5B9ULL;
    for (int q = 0; q < p; q++) { h ^= val[q]; h *= 0x100000001B3ULL; h ^= h >> 29; }
    rs = h ? h : 1;
    for (int q = 0; q < 4; q++) rnd();
}

static int uf_find(int x) { while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; } return x; }

/* orbits of F_2^n under sampled invertible elements of the stabiliser described by S; returns 0 if trivial */
static int orbits(const Sys *S) {
    if (S->rank >= z) return 0;
    for (int x = 0; x < N; x++) uf[x] = x;
    int got = 0;
    for (int tries = 0; tries < 8 * nsamp && got < nsamp; tries++) {
        u128 c = 0;
        uint8_t isp[MAXZ] = {0};
        for (int q = 0; q < S->rank; q++) isp[S->piv[q]] = 1;
        for (int j = 0; j < z; j++) if (!isp[j] && (rnd() & 1)) c |= (u128)1 << j;
        for (int q = 0; q < S->rank; q++) {
            int v = S->rhs[q] ^ par128(S->row[q] & c & ~((u128)1 << S->piv[q]));
            if (v) c |= (u128)1 << S->piv[q];
        }
        uint16_t X[MAXN];
        for (int i = 0; i < n; i++) { X[i] = 0; for (int j = 0; j < z; j++) if ((c >> j) & 1) X[i] ^= Mz[j][i]; }
        if (rank10(X) < n) continue;
        got++;
        for (int x = 0; x < N; x++) {
            int y = 0;
            for (int i = 0; i < n; i++) if ((x >> i) & 1) y ^= X[i];
            int a = uf_find(x), b = uf_find(y);
            if (a != b) { if (a < b) uf[b] = a; else uf[a] = b; }
        }
    }
    return got > 0;
}

/* ---- input side: checkpoint pruning at the end of a level ----
 * When beta is known on V_k = <e_0..e_{k-1}>, take sampled invertible S in Z(A) = {X : XA = AX} with S V_k = V_k
 * (precomputed per k). Q = beta(S., S.) on V_k is again a partial solution of the same search. Bring Q to the greedy
 * lexicographic form under C(B) (orbits from sampled stabiliser elements, the output map built explicitly as a
 * product of generators). If the result is lexicographically smaller than the current prefix, the branch is pruned.
 * Sound: every applied map is a genuine symmetry, and the lexicographically least member of each orbit under
 * C(A) x C(B) is never pruned (by this rule or by the output chain), so each orbit keeps a representative. */
#define MAXS 4096
static int za = 0, cp_min = 0, cp_max = -1, ns_in = 24, cp_fix = 0;     /* cp_fix: bit mask of fixed e_i */
static Sys sys_forced;                               /* output maps fixing every forced value */
static uint16_t Ma[MAXZ][MAXN];
static uint16_t Sk[MAXN + 1][MAXS][MAXN];
static int nSk[MAXN + 1];
static long long cp_calls = 0, cp_pruned = 0;

static void precompute_inputs(void) {
    for (int k = 2; k <= n; k++) {
        /* homogeneous conditions on c in F_2^za: for i < k, bits r >= k of (sum c_j Ma[j][i]) vanish */
        Sys S; memset(&S, 0, sizeof S);
        for (int i = 0; i < k; i++) for (int r = k; r < n; r++) {
            u128 coef = 0;
            for (int j = 0; j < za; j++) if ((Ma[j][i] >> r) & 1) coef |= (u128)1 << j;
            if (coef) sys_add(&S, coef, 0);
        }
        /* -cpfix U / -cpmask M: the input maps are the identity on these e_i (they then keep forced U x U pairs) */
        for (int i = 0; i < n; i++) if ((cp_fix >> i) & 1) for (int r = 0; r < n; r++) {
            u128 coef = 0;
            for (int j = 0; j < za; j++) if ((Ma[j][i] >> r) & 1) coef |= (u128)1 << j;
            sys_add(&S, coef, r == i);
        }
        nSk[k] = 0;
        uint8_t isp[MAXZ] = {0};
        for (int q = 0; q < S.rank; q++) isp[S.piv[q]] = 1;
        for (int tries = 0; tries < 200 * ns_in && nSk[k] < ns_in; tries++) {
            u128 c = 0;
            for (int j = 0; j < za; j++) if (!isp[j] && (rnd() & 1)) c |= (u128)1 << j;
            for (int q = 0; q < S.rank; q++)
                if (S.rhs[q] ^ par128(S.row[q] & c & ~((u128)1 << S.piv[q]))) c |= (u128)1 << S.piv[q];
            uint16_t X[MAXN];
            for (int i = 0; i < n; i++) { X[i] = 0; for (int j = 0; j < za; j++) if ((c >> j) & 1) X[i] ^= Ma[j][i]; }
            int okfix = 1;
            for (int i = 0; i < n; i++) if (((cp_fix >> i) & 1) && X[i] != (1 << i)) okfix = 0;
            if (!okfix) continue;
            if (rank10(X) < n) continue;
            memcpy(Sk[k][nSk[k]++], X, sizeof(uint16_t) * n);
        }
    }
}

static inline uint16_t matvec(const uint16_t *X, uint16_t v) {
    uint16_t r = 0;
    for (int i = 0; i < n; i++) if ((v >> i) & 1) r ^= X[i];
    return r;
}

/* greedy form of L (len values) under C(B), compared with P: 1 if a genuinely smaller image is found */
static int canon_smaller(uint16_t *L, int len, const uint16_t *P) {
    static Sys S; S = sys_forced;
    static uint16_t G[64][MAXN];
    static int prevv[1 << MAXN], prevg[1 << MAXN], qu[1 << MAXN];
    static uint8_t vis[1 << MAXN];
    for (int p = 0; p < len; p++) {
        if (forced[p]) {                                /* forced pair: no minimisation; images must agree */
            if (L[p] != P[p]) return 0;
            sys_fix_vector(&S, L[p]);
            continue;
        }
        if (S.rank < z) {
            int ng = 0;
            uint8_t isp[MAXZ] = {0};
            for (int q = 0; q < S.rank; q++) isp[S.piv[q]] = 1;
            for (int tries = 0; tries < 8 * nsamp && ng < nsamp; tries++) {
                u128 c = 0;
                for (int j = 0; j < z; j++) if (!isp[j] && (rnd() & 1)) c |= (u128)1 << j;
                for (int q = 0; q < S.rank; q++) {
                    int v = S.rhs[q] ^ par128(S.row[q] & c & ~((u128)1 << S.piv[q]));
                    if (v) c |= (u128)1 << S.piv[q];
                }
                uint16_t X[MAXN];
                for (int i = 0; i < n; i++) { X[i] = 0; for (int j = 0; j < z; j++) if ((c >> j) & 1) X[i] ^= Mz[j][i]; }
                if (rank10(X) < n) continue;
                memcpy(G[ng++], X, sizeof(uint16_t) * n);
            }
            if (ng) {
                /* BFS orbit of L[p] under the generators, keep the path to the minimum */
                memset(vis, 0, N);
                int h = 0, t = 0, start = L[p], mn = start;
                vis[start] = 1; prevv[start] = -1; qu[t++] = start;
                while (h < t) {
                    int v = qu[h++];
                    for (int g = 0; g < ng; g++) {
                        int w = matvec(G[g], (uint16_t)v);
                        if (!vis[w]) { vis[w] = 1; prevv[w] = v; prevg[w] = g; qu[t++] = w; if (w < mn) mn = w; }
                    }
                }
                if (mn < start) {
                    int path[1 << MAXN], np = 0;
                    for (int v = mn; v != start; v = prevv[v]) path[np++] = prevg[v];
                    for (int q = p; q < len; q++) {
                        uint16_t v = L[q];
                        for (int e = np - 1; e >= 0; e--) v = matvec(G[path[e]], v);
                        L[q] = v;
                    }
                }
            }
        }
        if (L[p] < P[p]) return 1;
        if (L[p] > P[p]) return 0;
        sys_fix_vector(&S, L[p]);
    }
    return 0;
}

/* checkpoint after level k - 1 is complete (beta known on V_k): 1 = prune */
static int checkpoint(int k) {
    if (k < cp_min || k > cp_max || !za || !nSk[k]) return 0;
    cp_calls++;
    reseed(k * (k - 1) / 2);
    int Pk = k * (k - 1) / 2;
    static uint16_t L[MAXP];
    for (int q = 0; q < nSk[k]; q++) {
        const uint16_t *X = Sk[k][q];
        for (int p = 0; p < Pk; p++) {
            uint16_t x = X[pairs_i[p]], y = X[pairs_j[p]], b = 0;
            for (int j = 0; j < k; j++) if ((y >> j) & 1) b ^= hb[x][j];
            L[p] = b;
        }
        if (canon_smaller(L, Pk, val)) { cp_pruned++; return 1; }
    }
    return 0;
}

static void dfs(int p, int depth);

/* -gbasis F C0 C1 ..: canonical basis of the intertwiner space under GL(<e_F .. e_{n-1}>). The intertwiner
 * T_s = beta(e_s, .) (s >= F) is identified with the tuple (beta(e_C0, e_s), beta(e_C1, e_s), ..) (it is determined
 * by it, e.g. C = 0 2 4 for F_4-linear T on F_4^3), compared lexicographically. Kept only if T_F < T_F+1 < .. and
 * each T_s is the least element of T_s + span(T_F .. T_{s-1}): the greedy (lexicographically least) basis of the
 * span, one representative per GL orbit of injective T. In the search order this is exactly "lexicographically
 * least in its GL(<e_F..>) orbit" (T_F decides first, then T_F+1, ..), so like -symdepth and -cp it keeps the least
 * member of every orbit: sound together with them and with -force on U x U pairs (U = indices < F). */
static int gb_F = -1, gb_nc = 0, gb_c[MAXN];
static long long gb_pruned = 0;
static int gb_ok(int p, uint16_t y) {
    int s = pairs_j[p], i = pairs_i[p];
    if (gb_F < 0 || s < gb_F) return 1;
    int cc = -1;                                   /* column index of this pair, if any */
    for (int c = 0; c < gb_nc; c++) if (gb_c[c] == i) cc = c;
    if (cc < 0) return 1;
    int last = (cc == gb_nc - 1);                  /* all components known: ties also decide */
    uint16_t t[MAXN][MAXN];                        /* t[k][c] = component c of T_{F+k}, c <= cc */
    int m = s - gb_F + 1;
    for (int k = 0; k < m; k++) for (int c = 0; c <= cc; c++)
        t[k][c] = (k == m - 1 && c == cc) ? y : val[pidx[gb_c[c]][gb_F + k]];
    if (m >= 2) {                                  /* T_{s-1} < T_s on the known prefix */
        int c = 0;
        while (c <= cc && t[m - 2][c] == t[m - 1][c]) c++;
        if ((c <= cc && t[m - 2][c] > t[m - 1][c]) || (c > cc && last)) { gb_pruned++; return 0; }
    }
    for (int v = 1; v < (1 << (m - 1)); v++) {    /* T_s <= T_s + v on the known prefix */
        int c = 0;
        for (; c <= cc; c++) {
            uint16_t w = t[m - 1][c];
            for (int k = 0; k < m - 1; k++) if ((v >> k) & 1) w ^= t[k][c];
            if (t[m - 1][c] < w) break;
            if (t[m - 1][c] > w) { gb_pruned++; return 0; }
        }
    }
    return 1;
}

static void after_pair(int p, int depth) {
    int s = pairs_j[p], i = pairs_i[p];
    if (i + 1 < s) { dfs(p + 1, depth); return; }
    level_end(s);                                   /* level s complete */
    if (checkpoint(s + 1)) return;
    if (s + 1 == n) {
        leaves++;
        if (dedupe) {
            uint64_t h = leaf_key(), slot = h & (HSIZE - 1);
            while (seen_keys[slot] && seen_keys[slot] != h) slot = (slot + 1) & (HSIZE - 1);
            if (!seen_keys[slot]) { seen_keys[slot] = h; distinct++; printf("K%016llx ", (unsigned long long)h); print_leaf(); fflush(stdout); }
        } else if (printed < print_max) { print_leaf(); printed++; fflush(stdout); }
        return;
    }
    level_start(s + 1);
    dfs(p + 1, depth);
}

/* progress: for each level of the current path, the index of the child being explored and the number of children;
 * the explored fraction is the mixed-radix number sum_d idx_d / cnt_d * prod_{e<d} 1 / cnt_e (assumes subtrees of
 * similar size: an estimate). Written every 60 s to the -progress file: elapsed s, nodes, fraction, estimated s left. */
static int pidx_path[MAXP + 1], pcnt_path[MAXP + 1];
static const char *progress_file = 0;
static time_t t_start, t_last;
static void progress_write(int final) {
    double f = 0, mult = 1;
    for (int d = 0; d <= MAXP && pcnt_path[d] > 0; d++) { f += mult * pidx_path[d] / pcnt_path[d]; mult /= pcnt_path[d]; }
    if (final == 1) f = 1;
    time_t now = time(0);
    double el = difftime(now, t_start), left = f > 0 ? el * (1 - f) / f : -1;
    FILE *fp = fopen(progress_file, "w");
    if (fp) {
        fprintf(fp, "elapsed_s %.0f nodes %lld leaves %lld fraction %.6f est_left_s %.0f%s\n", el, nodes, leaves, f, left,
                final == 1 ? " DONE" : final == 2 ? " INCOMPLETE (node limit reached)" : "");
        fclose(fp);
    }
    t_last = now;
}

static void dfs(int p, int depth) {
    if (node_limit >= 0 && nodes >= node_limit) return;
    nodes++;
    if (progress_file && (nodes & 0x3FFF) == 0 && difftime(time(0), t_last) >= 60) progress_write(0);
    int s = pairs_j[p], i = pairs_i[p], c = kcount[p], k0 = kfirst[p];
    if (p < symdepth) reseed(p);
    int use_sym = (p < symdepth) && orbits(&sys[p]);
    /* pass 1: the valid choices of the new bits, in Gray-code order (all filtering happens here, before recursing,
     * so the orbit table is not overwritten by deeper levels while in use) */
    static long long codes_buf[MAXP][1 << MAXN];
    long long *codes = codes_buf[p];
    int nv = 0;
    long long cur = 0;
    for (long long q = 0; q < (1LL << c); q++) {
        if (q) { int b = __builtin_ctzll(q); apply_t(k0 + b); cur ^= 1LL << b; }
        uint16_t y = val[p];
        if (forced[p]) { if (y != forced[p] - 1) continue; }
        else {
            if (p == 0 && !allowed01[y]) continue;     /* both filters keep the orbit minimum: sound together */
            if (use_sym && uf_find(y) != y) { sym_skipped++; continue; }
        }
        if (!gb_ok(p, y) || !check(s, i, y)) continue;
        if (depth + 1 == pshard_depth && npshards > 1 && pshard_count++ % npshards != pshard) continue;
        if (depth + 1 == shard_depth && nshards > 1 && shard_count++ % nshards != shard) continue;
        codes[nv++] = cur;
    }
    /* pass 2: explore them in the same order */
    pcnt_path[p] = nv;
    for (int k = 0; k < nv; k++) {
        if (node_limit >= 0 && nodes >= node_limit) break;
        long long diff = cur ^ codes[k];
        for (int b = 0; b < c; b++) if ((diff >> b) & 1) apply_t(k0 + b);
        cur = codes[k];
        pidx_path[p] = k;
        uint16_t y = val[p];
        if (p + 1 <= symdepth) { sys[p + 1] = sys[p]; sys_fix_vector(&sys[p + 1], y); }
        set_y(s, i, y);
        after_pair(p, depth + 1);
    }
    for (int b = 0; b < c; b++) if ((cur >> b) & 1) apply_t(k0 + b);    /* back to code 0 */
    pcnt_path[p] = 0;
}

int main(int argc, char **argv) {
    long long probes = 0;
    int mode_dfs = 0;
    for (int a = 1; a < argc; a++) {
        if (!strcmp(argv[a], "-knuth")) probes = atoll(argv[++a]);
        else if (!strcmp(argv[a], "-dfs")) mode_dfs = 1;
        else if (!strcmp(argv[a], "-print")) print_max = atoi(argv[++a]);
        else if (!strcmp(argv[a], "-limit")) node_limit = atoll(argv[++a]);
        else if (!strcmp(argv[a], "-shard")) { shard = atoi(argv[++a]); nshards = atoi(argv[++a]); }
        else if (!strcmp(argv[a], "-sdepth")) shard_depth = atoi(argv[++a]);
        else if (!strcmp(argv[a], "-pshard")) { pshard = atoi(argv[++a]); npshards = atoi(argv[++a]); }
        else if (!strcmp(argv[a], "-pdepth")) pshard_depth = atoi(argv[++a]);
        else if (!strcmp(argv[a], "-nonorm")) nonorm = 1;
        else if (!strcmp(argv[a], "-dedupe")) dedupe = 1;
        else if (!strcmp(argv[a], "-symdepth")) symdepth = atoi(argv[++a]);
        else if (!strcmp(argv[a], "-nsamp")) nsamp = atoi(argv[++a]);
        else if (!strcmp(argv[a], "-cp")) { cp_min = atoi(argv[++a]); cp_max = atoi(argv[++a]); }
        else if (!strcmp(argv[a], "-nsin")) { ns_in = atoi(argv[++a]); if (ns_in > MAXS) ns_in = MAXS; }
        else if (!strcmp(argv[a], "-cpfix")) cp_fix = (1 << atoi(argv[++a])) - 1;
        else if (!strcmp(argv[a], "-cpmask")) cp_fix = atoi(argv[++a]);
        else if (!strcmp(argv[a], "-gbasis")) {
            gb_F = atoi(argv[++a]);
            while (a + 1 < argc && argv[a + 1][0] != '-') gb_c[gb_nc++] = atoi(argv[++a]);
        }
        else if (!strcmp(argv[a], "-progress")) progress_file = argv[++a];
        else if (!strcmp(argv[a], "-offset")) {
            FILE *ff = fopen(argv[++a], "r"); int pp, vv;
            while (ff && fscanf(ff, "%d %d", &pp, &vv) == 2) if (pp >= 0 && pp < MAXP) offset[pp] = (uint16_t)vv;
            if (ff) fclose(ff);
        }
        else if (!strcmp(argv[a], "-force")) {
            FILE *ff = fopen(argv[++a], "r"); int pp, vv;
            while (ff && fscanf(ff, "%d %d", &pp, &vv) == 2) if (pp >= 0 && pp < MAXP) forced[pp] = vv + 1;
            if (ff) fclose(ff);
        }
        else if (!strcmp(argv[a], "-seed")) rs ^= (uint64_t)atoll(argv[++a]) * 0x9E3779B97F4A7C15ULL;
    }
    for (int q = 0; q < 20; q++) rnd();
    rs2 ^= rs; for (int q = 0; q < 20; q++) rnd2();
    if (scanf("%d %d", &n, &d) != 2) return 1;
    N = 1 << n; P = 0;
    for (int j = 0; j < n; j++) for (int i = 0; i < j; i++) { pairs_i[P] = i; pairs_j[P] = j; pidx[i][j] = P++; }
    for (int p = 0; p < P; p++) { kfirst[p] = -1; kcount[p] = 0; }
    for (int k = 0; k < d; k++) {
        int pp, pb;
        if (scanf("%d %d", &pp, &pb) != 2) return 1;
        pivpair[k] = pp;
        for (int p = 0; p < P; p++) { int x; if (scanf("%d", &x) != 1) return 1; wv[k][p] = (uint16_t)x; }
        if (kfirst[pp] < 0) kfirst[pp] = k;
        kcount[pp]++;
    }
    for (int p = 0; p < P; p++) if (kfirst[p] < 0) kfirst[p] = 0;
    int nrep = 0;                  /* the representative block is always consumed; it is used only without -nonorm */
    if (scanf("%d", &nrep) == 1 && nrep >= 0) {
        for (int q = 0; q < nrep; q++) { int x; if (scanf("%d", &x) != 1 || x < 0 || x >= N) return 1; if (!nonorm) allowed01[x] = 1; }
        if (nonorm) for (int x = 0; x < N; x++) allowed01[x] = 1;
    } else for (int x = 0; x < N; x++) allowed01[x] = 1;
    if (scanf("%d", &z) == 1 && z > 0 && z <= MAXZ) {
        for (int j = 0; j < z; j++) for (int i = 0; i < n; i++) { int x; if (scanf("%d", &x) != 1) return 1; Mz[j][i] = (uint16_t)x; }
    } else z = 0;
    if (!z) symdepth = 0;
    if (scanf("%d", &za) == 1 && za > 0 && za <= MAXZ) {
        for (int j = 0; j < za; j++) for (int i = 0; i < n; i++) { int x; if (scanf("%d", &x) != 1) return 1; Ma[j][i] = (uint16_t)x; }
    } else za = 0;
    if (za && cp_max >= cp_min) precompute_inputs();
    memset(&sys_forced, 0, sizeof(Sys));
    if (z) for (int p = 0; p < P; p++) if (forced[p]) sys_fix_vector(&sys_forced, (uint16_t)(forced[p] - 1));
    sys[0] = sys_forced;
    fprintf(stderr, "sym10: n=%d, d=%d free bits; per level:", n, d);
    for (int s = 1; s < n; s++) { int c = 0; for (int p = 0; p < P; p++) if (pairs_j[p] == s) c += kcount[p]; fprintf(stderr, " %d", c); }
    fprintf(stderr, "\n");
    memset(hb, 0, sizeof hb);
    memcpy(val, offset, sizeof val);

    if (probes > 0) {
        double est[MAXP + 1] = {0}, leafest = 0;
        for (long long pr = 0; pr < probes; pr++) {
            memcpy(val, offset, sizeof val); memset(hb, 0, sizeof hb);
            level_start(1);
            double w = 1;
            int p = 0, alive = 1;
            while (alive) {
                int s = pairs_j[p], i = pairs_i[p], c = kcount[p], k0 = kfirst[p];
                /* collect the valid choices */
                static uint32_t ok[1 << 12]; int nok = 0;
                uint16_t base = val[p];
                int use_sym = (p < symdepth) && orbits(&sys[p]);
                for (long long q = 0; q < (1LL << c); q++) {
                    uint16_t y = base;
                    for (int b = 0; b < c; b++) if ((q >> b) & 1) y ^= wv[k0 + b][p];
                    if (forced[p]) { if (y != forced[p] - 1) continue; }
                    else {
                        if (p == 0 && !allowed01[y]) continue;     /* both filters keep the orbit minimum: sound together */
                        if (use_sym && uf_find(y) != y) continue;
                    }
                    if (gb_ok(p, y) && check(s, i, y)) ok[nok++] = (uint32_t)q;
                }
                if (!nok) break;
                w *= nok;
                uint32_t q = ok[rnd2() % nok];
                for (int b = 0; b < c; b++) if ((q >> b) & 1) apply_t(k0 + b);
                if (p + 1 <= symdepth) { sys[p + 1] = sys[p]; sys_fix_vector(&sys[p + 1], val[p]); }
                set_y(s, i, val[p]);
                est[p + 1] += w;
                if (i + 1 == s) {
                    level_end(s);
                    if (checkpoint(s + 1)) break;
                    if (s + 1 == n) { leafest += w; break; }
                    level_start(s + 1);
                }
                p++;
            }
        }
        double total = 1;
        fprintf(stderr, "Knuth estimate (%lld probes):\n", probes);
        for (int p = 1; p <= P; p++) {
            total += est[p] / probes;
            fprintf(stderr, "  pair %2d (%d,%d): %.3g\n", p, pairs_i[p - 1], pairs_j[p - 1], est[p] / probes);
        }
        fprintf(stderr, "total nodes %.3g, leaves (APN functions in W) %.3g; checkpoints %lld calls, %lld pruned\n", total, leafest / probes, cp_calls, cp_pruned);
        return 0;
    }
    if (mode_dfs) {
        level_start(1);
        t_start = t_last = time(0);
        dfs(0, 0);
        int limited = node_limit >= 0 && nodes >= node_limit;
        if (progress_file) progress_write(limited ? 2 : 1);
        if (gb_F >= 0) fprintf(stderr, "greedy basis: %lld pruned\n", gb_pruned);
        fprintf(stderr, "symmetry: skipped %lld values (symdepth %d, z = %d); checkpoints %lld calls, %lld pruned (za = %d)\n",
                sym_skipped, symdepth, z, cp_calls, cp_pruned, za);
        fprintf(stderr, "shard counters: depth %d: %lld, depth %d: %lld\n", pshard_depth, pshard_count, shard_depth, shard_count);
        fprintf(stderr, "done: nodes %lld, leaves %lld, distinct fingerprints %lld%s\n", nodes, leaves, distinct,
                limited ? " (node limit reached)" : "");
        if (limited) { fprintf(stderr, "INCOMPLETE: node limit %lld reached, the search is not exhaustive\n", node_limit); return 3; }
    }
    return 0;
}
