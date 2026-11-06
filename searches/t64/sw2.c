/*
 * sw2: codimension-2 Hermitian switching for type 64 (see switch64.py for codimension 1).
 * Input: V (10 Hermitian 5x5 matrices over F_4, 25-bit codes as in herm.py / switch64.py) on the command line.
 * For every 8-dim subspace U of V (174,251) find all planes <H', H''> mod U such that V' = U + <H', H''> kills every
 * point of PG(4, 4) exactly once:
 *   point q killed by U (U q is 7-dim in q^perp): pi_q(H'), pi_q(H'') independent in q^perp / U q (2-dim);
 *   point q not killed by U (U q is 8-dim):       (l_q(H' q), l_q(H'' q)) != (0, 0) in q^perp / U q (1-dim).
 * Prints each V' != V as a line of 10 codes (U basis, then H', H''). Exact for each U.
 * usage: sw2 C1 .. C10 [-umax N] [-ustart S -ustep T]
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const int F4M[4][4] = {{0,0,0,0},{0,1,2,3},{0,2,3,1},{0,3,1,2}};
static int pts[341]; static uint16_t img[341][25];
static int coord(int v, int i) { return (v >> (2 * i)) & 3; }
static void basis_mat(int b, int M[5][5]) {
    memset(M, 0, sizeof(int) * 25);
    if (b < 5) { M[b][b] = 1; return; }
    int k = 5;
    for (int i = 0; i < 5; i++) for (int j = i + 1; j < 5; j++) {
        if (b == k) { M[i][j] = 1; M[j][i] = 1; return; }
        if (b == k + 1) { M[i][j] = 2; M[j][i] = 3; return; }
        k += 2;
    }
}
static int matvec(int M[5][5], int v) {
    int r = 0;
    for (int i = 0; i < 5; i++) { int s = 0; for (int j = 0; j < 5; j++) s ^= F4M[M[i][j]][coord(v, j)]; r |= s << (2 * i); }
    return r;
}
static inline int hq(uint32_t H, int q) { int r = 0; while (H) { int b = __builtin_ctz(H); r ^= img[q][b]; H &= H - 1; } return r; }
static inline int par(uint32_t x) { return __builtin_parity(x); }

/* trace pairing on F_4^5: <q, y> = Tr(sum conj(q_i) y_i), Tr(z) = 1 iff z in {w, w^2} */
static int trpair(int q, int y) { int s = 0; for (int i = 0; i < 5; i++) { int c = coord(q, i); int cq = c == 2 ? 3 : c == 3 ? 2 : c; s ^= F4M[cq][coord(y, i)]; } return s >= 2; }

/* echelon basis helpers over 10-bit vectors */
static int ech(int *b, int nb, int v) { for (int t = 0; t < nb; t++) if ((v ^ b[t]) < v) v ^= b[t]; return v; }
static int add_ech(int *b, int nb, int v) {
    v = ech(b, nb, v); if (!v) return nb;
    int i = nb; b[nb++] = v;
    while (i > 0 && b[i] > b[i - 1]) { int t = b[i]; b[i] = b[i - 1]; b[i - 1] = t; i--; }
    return nb;
}
/* functionals l (10-bit) vanishing on span(B): solve */
static int annihilator(int *B, int nb, int *out) {
    int no = 0;
    /* brute force small: 1024 candidates */
    int basis[10]; int nbas = 0;
    for (int l = 1; l < 1024; l++) {
        int ok = 1; for (int t = 0; t < nb; t++) if (par(l & B[t])) { ok = 0; break; }
        if (!ok) continue;
        int r = ech(basis, nbas, l);
        if (r) { nbas = add_ech(basis, nbas, l); out[no++] = l; }
    }
    return no;
}


/* fingerprint: differential spectrum of the ortho-derivative of F(x) = (x^* H_k x)_k */
static uint8_t qf[25][1024];
static void qf_setup(void) {
    for (int b = 0; b < 25; b++) { int M[5][5]; basis_mat(b, M);
        for (int x = 0; x < 1024; x++) { /* x^* M x as Tr of off-diagonal + diagonal norms */
            int s = 0;
            for (int i = 0; i < 5; i++) { int ci = coord(x, i); int cc = ci == 2 ? 3 : ci == 3 ? 2 : ci; int r = 0; for (int j = 0; j < 5; j++) r ^= F4M[M[i][j]][coord(x, j)]; s ^= F4M[cc][r]; }
            if (s > 1) { fprintf(stderr, "x*Mx not in F2\n"); exit(1); }
            qf[b][x] = (uint8_t)s; } }
}
static uint64_t fingerprint(const uint32_t *W) {
    static int F[1024], pi[1024]; static int cnt[1024];
    for (int x = 0; x < 1024; x++) { int v = 0; for (int k = 0; k < 10; k++) { uint32_t H = W[k]; int s = 0; while (H) { int b = __builtin_ctz(H); s ^= qf[b][x]; H &= H - 1; } v |= s << k; } F[x] = v; }
    pi[0] = 0;
    for (int a = 1; a < 1024; a++) {
        int B[10]; int nb = 0;
        for (int i = 0; i < 10; i++) { int y = F[a ^ (1 << i)] ^ F[a] ^ F[1 << i] ^ F[0]; nb = add_ech(B, nb, y); }
        int c = 0; for (int l = 1; l < 1024; l++) { int ok = 1; for (int t = 0; t < nb; t++) if (par(l & B[t])) { ok = 0; break; } if (ok) { c = l; break; } }
        pi[a] = c;
    }
    static uint32_t hist[1025]; memset(hist, 0, sizeof hist);
    for (int a = 1; a < 1024; a++) { memset(cnt, 0, sizeof cnt); for (int x = 0; x < 1024; x++) cnt[pi[x ^ a] ^ pi[x]]++; for (int b = 0; b < 1024; b++) hist[cnt[b]]++; }
    uint64_t h = 1469598103934665603ULL; for (int i = 0; i <= 1024; i++) if (hist[i]) { h ^= (uint64_t)i * 1000003u + hist[i]; h *= 1099511628211ULL; }
    return h;
}
#define MAXK 100000
static uint64_t seenk[MAXK]; static int nseen = 0;

static int cmpu(const void *x, const void *y) { uint32_t a = *(const uint32_t *)x, b = *(const uint32_t *)y; return a < b ? -1 : a > b; }

int main(int argc, char **argv) {
    uint32_t V[10]; int nv = 0; long umax = -1, ustart = 0, ustep = 1;
    for (int a = 1; a < argc; a++) {
        if (!strcmp(argv[a], "-umax")) umax = atol(argv[++a]);
        else if (!strcmp(argv[a], "-ustart")) ustart = atol(argv[++a]);
        else if (!strcmp(argv[a], "-ustep")) ustep = atol(argv[++a]);
        else V[nv++] = (uint32_t)strtoul(argv[a], 0, 10);
    }
    if (nv != 10) { fprintf(stderr, "need 10 codes\n"); return 1; }
    qf_setup();
    int np = 0;
    for (int v = 1; v < 1024; v++) { int lead = 0; for (int i = 0; i < 5; i++) if (coord(v, i)) { lead = coord(v, i); break; } if (lead == 1) pts[np++] = v; }
    for (int b = 0; b < 25; b++) { int M[5][5]; basis_mat(b, M); for (int q = 0; q < 341; q++) img[q][b] = (uint16_t)matvec(M, pts[q]); }
    /* kill vectors t_q */
    int tq[341];
    for (int q = 0; q < 341; q++) {
        int im[10]; for (int k = 0; k < 10; k++) im[k] = hq(V[k], q);
        tq[q] = 0;
        for (int t = 1; t < 1024; t++) { int s = 0; for (int k = 0; k < 10; k++) if ((t >> k) & 1) s ^= im[k]; if (!s) { tq[q] = t; break; } }
        if (!tq[q]) { fprintf(stderr, "input is not a type-64 APN space (point %d)\n", q); return 1; }
    }
    long ucount = 0, found = 0;
    /* 2-dim subspaces of the dual: canonical pair f1 < f2 with f1 < (f1 ^ f2), f2 < (f1 ^ f2) */
    for (int f1 = 1; f1 < 1024; f1++) for (int f2 = f1 + 1; f2 < 1024; f2++) {
        int f3 = f1 ^ f2; if (f3 < f2) continue;
        if (ucount++ % ustep != ustart % ustep) continue;
        if (umax >= 0 && ucount > umax * ustep) goto done;
        /* U = {t : f1.t = f2.t = 0}: basis in V-coordinates */
        int ub[8]; int nu = 0;
        {
            int cand[10]; int nc = 0;
            for (int t = 1; t < 1024 && nc < 8; t++) if (!par(f1 & t) && !par(f2 & t)) { int r = ech(cand, nc, t); if (r) nc = add_ech(cand, nc, t); }
            for (int i = 0; i < 8; i++) ub[i] = cand[i]; nu = nc;
        }
        if (nu != 8) { fprintf(stderr, "U dim %d\n", nu); return 1; }
        uint32_t U[8];
        for (int i = 0; i < 8; i++) { U[i] = 0; for (int k = 0; k < 10; k++) if ((ub[i] >> k) & 1) U[i] ^= V[k]; }
        /* complement of U in Herm (25-bit): echelon of U, free bits */
        uint32_t ue[8]; int ne = 0;
        for (int i = 0; i < 8; i++) {
            uint32_t v = U[i];
            for (int t = 0; t < ne; t++) if ((v ^ ue[t]) < v) v ^= ue[t];
            if (v) { int j = ne++; ue[j] = v; while (j > 0 && ue[j] > ue[j - 1]) { uint32_t tmp = ue[j]; ue[j] = ue[j - 1]; ue[j - 1] = tmp; j--; } }
        }
        uint32_t lead = 0; for (int t = 0; t < ne; t++) lead |= 1u << (31 - __builtin_clz(ue[t]));
        int cb[17]; int ncb = 0; for (int b = 0; b < 25; b++) if (!((lead >> b) & 1)) cb[ncb++] = b;
        /* per point: functionals */
        int kq[341], l1[341], l2[341], nkill = 0;
        int order[341];
        for (int q = 0; q < 341; q++) {
            int Uq[10]; int nb = 0;
            for (int i = 0; i < 8; i++) nb = add_ech(Uq, nb, hq(U[i], q));
            int killed = !par(f1 & tq[q]) && !par(f2 & tq[q]);
            kq[q] = killed;
            /* functionals vanishing on Uq, restricted to be independent on q^perp: annihilator of Uq + <y0> */
            int y0 = 1; while (!trpair(pts[q], y0)) y0++;
            int B[11]; int nb2 = nb; memcpy(B, Uq, sizeof(int) * nb); nb2 = add_ech(B, nb2, y0);
            int ann[10]; int na = annihilator(B, nb2, ann);
            if (killed) { if (na != 2) { fprintf(stderr, "ann %d\n", na); return 1; } l1[q] = ann[0]; l2[q] = ann[1]; order[nkill++] = q; }
            else { if (na != 1) { fprintf(stderr, "ann1 %d\n", na); return 1; } l1[q] = ann[0]; l2[q] = 0; }
        }
        int nall = nkill; for (int q = 0; q < 341; q++) if (!kq[q]) order[nall++] = q;
        /* candidates x (17 bits over cb): pi_q(x) != 0 for all killed q */
        static uint32_t cx[1 << 17];
        int ncand = 0;
        for (uint32_t m = 1; m < (1u << 17); m++) {
            uint32_t H = 0; for (int k = 0; k < 17; k++) if ((m >> k) & 1) H |= 1u << cb[k];
            int ok = 1;
            for (int t = 0; t < nkill; t++) { int q = order[t]; int y = hq(H, q); if (!par(l1[q] & y) && !par(l2[q] & y)) { ok = 0; break; } }
            if (ok) cx[ncand++] = H;
        }
        qsort(cx, ncand, sizeof(uint32_t), cmpu);
        /* pairs */
        for (int i = 0; i < ncand; i++) for (int j = i + 1; j < ncand; j++) {
            uint32_t a = cx[i], b = cx[j], c = a ^ b;
            /* plane <a, b>: count once (a < b < a^b in candidate terms) */
            if (c < b) continue;
            int ok = 1;
            for (int t = 0; t < nall && ok; t++) {
                int q = order[t]; int ya = hq(a, q), yb = hq(b, q);
                if (kq[q]) {
                    int pa = par(l1[q] & ya) | (par(l2[q] & ya) << 1), pb = par(l1[q] & yb) | (par(l2[q] & yb) << 1);
                    if (!pb || pa == pb) ok = 0;
                } else {
                    if (!par(l1[q] & ya) && !par(l1[q] & yb)) ok = 0;
                }
            }
            if (!ok) continue;
            /* skip V itself */
            int inV = 1; uint32_t vb[10]; int nvb = 0;
            for (int k = 0; k < 10; k++) { uint32_t v = V[k]; for (int t = 0; t < nvb; t++) if ((v ^ vb[t]) < v) v ^= vb[t]; if (v) { int jj = nvb++; vb[jj] = v; while (jj > 0 && vb[jj] > vb[jj - 1]) { uint32_t tmp = vb[jj]; vb[jj] = vb[jj - 1]; vb[jj - 1] = tmp; jj--; } } }
            for (int s = 0; s < 2; s++) { uint32_t v = s ? b : a; for (int t = 0; t < nvb; t++) if ((v ^ vb[t]) < v) v ^= vb[t]; if (v) inV = 0; }
            if (inV) continue;
            found++;
            { uint32_t W[10]; for (int k = 0; k < 8; k++) W[k] = U[k]; W[8] = a; W[9] = b;
              uint64_t key = fingerprint(W); int seen = 0; for (int t = 0; t < nseen; t++) if (seenk[t] == key) { seen = 1; break; }
              if (!seen && nseen < MAXK) { seenk[nseen++] = key;
                printf("K%016llx", (unsigned long long)key); for (int k = 0; k < 10; k++) printf(" %u", W[k]); printf("\n"); fflush(stdout); } }
        }
    }
done:
    fflush(stdout);
    fprintf(stderr, "subspaces %ld, new spaces %ld\n", ucount, found);
    return 0;
}
