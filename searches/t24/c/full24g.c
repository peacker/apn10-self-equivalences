/* Type 24 full classification: plane formulation + canonical plane (see ../../REDESIGN.md, "Type 24: plane
   formulation"). Normal form (form24.py): F(p, q) = (G(p) + S(p, p) + mu(p, q), D(p) + S(q, p) + G(q)).
   F is APN iff the diag conditions hold and, for every plane U = <a, p> of F^m, Dp(a, p) avoids Forb(U)
   (Forb depends only on the columns S(., u), u in U).
   Orbit reduction: Aut(G) = {x -> l x^(2^k)} acts regularly on the planes of F_32 (planes.py); Forb transforms by an
   affine map, so key(U) = (|Forb|, affine dim, sorted sizes of the three parts) is invariant. We keep only
   solutions where U0 = <e_0, e_1> has the maximal key among all planes (ties kept). Sound for any group action that
   is transitive on planes.
   Order: S columns 0..3 (S only, canonical check on every new plane), then column-4 options (S only), then Dp on
   span(e_0..e_3) with stored plane masks, then Dp(., e_4) for each column-4 option.
   usage: full24d   env GE = Gold exponent, OUT = LUT file, PROG = progress seconds,
          KNUTH = probes with KDEPTH (random descent through columns < KDEPTH, exhaustive below),
          PART = "r/n": only column-1 options with index = r mod n (for splitting a run).
   Layout: x = p | q << M, value = top | bot << M. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>
#ifndef M
#define M 5
#endif
#define Q (1 << M)
#define NV (M * M * M)
typedef unsigned __int128 u128;
typedef unsigned long long u64;
static int POLY;
static int G[Q], mu[Q][Q], hb[Q], solv[Q][Q];
static int idx(int i, int j, int k) { return (j * M + i) * M + k; }
static int fm(int a, int b) { int r = 0; while (b) { if (b & 1) r ^= a; b >>= 1; a <<= 1; if (a & Q) a ^= POLY; } return r; }
static int fpow(int a, int e) { int r = 1; while (e) { if (e & 1) r = fm(r, a); a = fm(a, a); e >>= 1; } return r; }
static int bit(u128 v, int i) { return (int)((v >> i) & 1); }
static u128 B1(int i) { return (u128)1 << i; }
static int ispiv[NV]; static u128 Teb[2 * M * M]; static int Tep[2 * M * M], Tne;
static void setup_T(void) {
    u128 T[2 * M * M]; int nt = 0;
    for (int r = 0; r < M; r++) for (int s = 0; s < M; s++) {
        u128 v = 0, w = 0;
        for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) {
            int y = (j == s) ? mu[1 << i][1 << r] : 0;
            int z = (mu[1 << i][1 << j] >> s & 1) ? (1 << r) : 0;
            for (int k = 0; k < M; k++) { if (y >> k & 1) v |= B1(idx(i, j, k)); if (z >> k & 1) w |= B1(idx(i, j, k)); }
        }
        T[nt++] = v; T[nt++] = w;
    }
    u128 eb[2 * M * M]; int ep[2 * M * M], ne = 0;
    for (int t = 0; t < nt; t++) {
        u128 v = T[t];
        for (int i = 0; i < ne; i++) if (bit(v, ep[i])) v ^= eb[i];
        if (!v) continue;
        int p = 0; while (!bit(v, p)) p++;
        for (int i = 0; i < ne; i++) if (bit(eb[i], p)) eb[i] ^= v;
        eb[ne] = v; ep[ne++] = p;
    }
    for (int i = 0; i < ne; i++) { ispiv[ep[i]] = 1; Teb[i] = eb[i]; Tep[i] = ep[i]; } Tne = ne;
}
/* ---------- state ---------- */
static int sv[M][M];
static int St[Q][Q];
static int dpv[M][M];
static unsigned pm[M][Q][Q];          /* pm[j][w][up]: Forb mask of plane <w, e_j + up> */
static u64 key0; static int pid[Q][Q]; static unsigned pmask[200]; static int useMin;
static long long cnt_col[M + 1], cnt_canon[M + 1], dnodes[M + 1], leaves, apn;
static FILE *fout;
static volatile sig_atomic_t tick; static int progsec;
static void on_alarm(int s) { (void)s; tick = 1; }
static int knuth, kdepth = 99; static double kw = 1, kest[M + 1]; static long long subn;
static u64 rs = 88172645463325252ull;
static u64 rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static int part_r = 0, part_n = 1; static int forced, fsv[M][M];
static int col_ok(int j, const int *o) { if (!forced) return 1; for (int i = 0; i < M; i++) if (o[i] != fsv[i][j]) return 0; return 1; }
/* ---------- planes ---------- */
static unsigned forb1(int a, int p) {
    unsigned f = 0;
    int c0 = mu[a][p] ^ St[a][p] ^ St[p][a];
    for (int q = 0; q < Q; q++) {
        int b = solv[p][c0 ^ mu[a][q]]; if (b < 0) continue;
        f |= 1u << (St[b][p] ^ St[q][a] ^ mu[b][q]);
        b ^= p;
        f |= 1u << (St[b][p] ^ St[q][a] ^ mu[b][q]);
    }
    return f;
}
static int adim(unsigned m) {
    if (!m) return -1;
    int x0 = __builtin_ctz(m); int bs[M], nb = 0;
    for (unsigned t = m & (m - 1); t; t &= t - 1) {
        int v = __builtin_ctz(t) ^ x0;
        for (int i = 0; i < nb; i++) if ((v ^ bs[i]) < v) v ^= bs[i];
        if (v) { bs[nb++] = v; for (int i = nb - 1; i > 0 && bs[i] > bs[i - 1]; i--) { int x = bs[i]; bs[i] = bs[i - 1]; bs[i - 1] = x; } }
    }
    return nb;
}
static unsigned forb1c(int a, int p, uint8_t *cnt) {
    unsigned f = 0;
    int c0 = mu[a][p] ^ St[a][p] ^ St[p][a];
    for (int q = 0; q < Q; q++) {
        int b = solv[p][c0 ^ mu[a][q]]; if (b < 0) continue;
        int v = St[b][p] ^ St[q][a] ^ mu[b][q]; f |= 1u << v; cnt[v]++;
        b ^= p;
        v = St[b][p] ^ St[q][a] ^ mu[b][q]; f |= 1u << v; cnt[v]++;
    }
    return f;
}
static u64 hmix(u64 h, u64 x) { h ^= x + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2); return h; }
/* invariant key of a plane (finer): |F|, affine dim, sorted pairwise intersections, triple intersection,
   hash of the multiplicity histograms (total and per part, parts sorted) */
static u64 pkey(int a, int p, unsigned *mask) {
    uint8_t c[3][Q]; memset(c, 0, sizeof c);
    unsigned f1 = forb1c(a, p, c[0]), f2 = forb1c(a, a ^ p, c[1]), f3 = forb1c(a ^ p, p, c[2]), f = f1 | f2 | f3;
    *mask = f;
    int in[3] = {__builtin_popcount(f1 & f2), __builtin_popcount(f1 & f3), __builtin_popcount(f2 & f3)};
    for (int x = 0; x < 2; x++) for (int y = 0; y < 2 - x; y++) if (in[y] < in[y + 1]) { int t = in[y]; in[y] = in[y + 1]; in[y + 1] = t; }
    int ht[3 * 2 * Q + 1] = {0}; u64 hp[3];
    for (int k = 0; k < 3; k++) { int h[2 * Q + 1] = {0}; for (int v = 0; v < Q; v++) h[c[k][v]]++; u64 x = 0; for (int t = 1; t <= 2 * Q; t++) x = hmix(x, (u64)h[t] << 8 | t); hp[k] = x; }
    for (int v = 0; v < Q; v++) ht[c[0][v] + c[1][v] + c[2][v]]++;
    for (int x = 0; x < 2; x++) for (int y = 0; y < 2 - x; y++) if (hp[y] < hp[y + 1]) { u64 t = hp[y]; hp[y] = hp[y + 1]; hp[y + 1] = t; }
    u64 hh = 0; for (int t = 1; t <= 6 * Q; t++) if (ht[t]) hh = hmix(hh, (u64)ht[t] << 10 | t);
    hh = hmix(hmix(hmix(hh, hp[0]), hp[1]), hp[2]);
    return (u64)__builtin_popcount(f) << 58 | (u64)(adim(f) + 1) << 54 | (u64)in[0] << 48 | (u64)in[1] << 42 | (u64)in[2] << 36
         | (u64)__builtin_popcount(f1 & f2 & f3) << 30 | (hh & 0x3fffffffu);
}
/* new planes of column j: <w, e_j + up>, w in [1, 2^j), up < 2^j with bit top(w) clear; returns 0 if non-canonical */
static int planes_col(int j, unsigned pmj[Q][Q]) {
    for (int w = 1; w < (1 << j); w++) {
        int tb = 31 - __builtin_clz(w);
        for (int up = 0; up < (1 << j); up++) {
            if (up >> tb & 1) continue;
            unsigned m; u64 k = pkey(w, (1 << j) | up, &m);
            if (useMin) k = ~k;
            if (j == 1) key0 = k; else if (k > key0) return 0;
            pmj[w][up] = m; pmask[pid[w][(1 << j) | up]] = m;
        }
    }
    return 1;
}
/* ---------- column options ---------- */
typedef struct { int n; int val[1 << 16][M]; } ColList;   /* up to 2^16 options, each M row values */
static void col_options(int j, ColList *L) {
    int fv[M * M], nf = 0;
    for (int i = 0; i < M; i++) for (int k = 0; k < M; k++) if (!ispiv[idx(i, j, k)]) fv[nf++] = i * M + k;
    uint32_t rows[Q]; int nr = 0;
    for (int b = 1 << j; b < 2 << j; b++) {
        int rest = 0;
        for (int i = 0; i < M; i++) if (b >> i & 1) for (int l = 0; l < j; l++) if (b >> l & 1) rest ^= sv[i][l];
        uint32_t r = 0;
        for (int t = 0; t < nf; t++) { int i = fv[t] / M, k = fv[t] % M; if ((b >> i & 1) && (hb[b] >> k & 1)) r |= 1u << t; }
        rows[nr++] = r | (uint32_t)(1 ^ __builtin_parity(hb[b] & rest)) << 31;
    }
    int piv[Q], rk = 0;
    for (int t = 0; t < nf; t++) {
        int q = -1; for (int s = rk; s < nr; s++) if (rows[s] >> t & 1) { q = s; break; }
        if (q < 0) continue;
        uint32_t x = rows[q]; rows[q] = rows[rk]; rows[rk] = x;
        for (int s = 0; s < nr; s++) if (s != rk && (rows[s] >> t & 1)) rows[s] ^= rows[rk];
        piv[rk++] = t;
    }
    L->n = 0;
    for (int s = rk; s < nr; s++) if (rows[s]) return;
    int fr[M * M], nfr = 0, isp[M * M] = {0};
    for (int s = 0; s < rk; s++) isp[piv[s]] = 1;
    for (int t = 0; t < nf; t++) if (!isp[t]) fr[nfr++] = t;
    if (nfr > 16) { fprintf(stderr, "too many column options\n"); exit(1); }
    for (uint32_t c = 0; c < (1u << nfr); c++) {
        uint32_t x = 0;
        for (int t = 0; t < nfr; t++) if (c >> t & 1) x |= 1u << fr[t];
        for (int s = 0; s < rk; s++) { int v = rows[s] >> 31 & 1; v ^= __builtin_parity(rows[s] & x & 0x7fffffffu); if (v) x |= 1u << piv[s]; }
        int *o = L->val[L->n++];
        for (int i = 0; i < M; i++) o[i] = 0;
        for (int t = 0; t < nf; t++) if (x >> t & 1) o[fv[t] / M] |= 1 << (fv[t] % M);
    }
}
static void set_col(int j, const int *o) {
    for (int i = 0; i < M; i++) sv[i][j] = o[i];
    for (int q = 0; q < Q; q++) {
        int cj = 0; for (int i = 0; i < M; i++) if (q >> i & 1) cj ^= o[i];
        for (int p = 1 << j; p < 2 << j; p++) St[q][p] = St[q][p ^ (1 << j)] ^ cj;
    }
}
/* ---------- Dp ---------- */
static unsigned xshift(unsigned m, int z) {
    if (z & 1) m = ((m & 0x55555555u) << 1) | ((m >> 1) & 0x55555555u);
    if (z & 2) m = ((m & 0x33333333u) << 2) | ((m >> 2) & 0x33333333u);
    if (z & 4) m = ((m & 0x0f0f0f0fu) << 4) | ((m >> 4) & 0x0f0f0f0fu);
    if (z & 8) m = ((m & 0x00ff00ffu) << 8) | ((m >> 8) & 0x00ff00ffu);
    if (z & 16) m = (m << 16) | (m >> 16);
    return m;
}
static int Dpv(int a, int p) { int r = 0; for (int i = 0; i < M; i++) if (a >> i & 1) for (int j = 0; j < M; j++) if (p >> j & 1) r ^= dpv[i][j]; return r; }
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
    if (fout) { fprintf(fout, "%d", ok); for (int x = 0; x < Q * Q; x++) fprintf(fout, " %d", F[x]); fprintf(fout, "\n"); fflush(fout); }
    if (!ok) { fprintf(stderr, "BUG: leaf not APN\n"); exit(2); }
}
static ColList LL[M];
static int cheapkey;
/* cheap invariant key: |Forb|, sorted sizes of the three parts */
static u64 ckey(int a, int p, unsigned *mask) {
    unsigned f1 = forb1(a, p), f2 = forb1(a, a ^ p), f3 = forb1(a ^ p, p), f = f1 | f2 | f3; *mask = f;
    int s0 = __builtin_popcount(f1), s1 = __builtin_popcount(f2), s2 = __builtin_popcount(f3), t;
    if (s0 < s1) { t = s0; s0 = s1; s1 = t; } if (s1 < s2) { t = s1; s1 = s2; s2 = t; } if (s0 < s1) { t = s0; s0 = s1; s1 = t; }
    return (u64)__builtin_popcount(f) << 18 | (u64)s0 << 12 | (u64)s1 << 6 | (u64)s2;
}
static u64 key0f;
static int canon_col(int j, unsigned pmj[Q][Q]) {   /* masks of the new planes of column j; 0 if non-canonical */
    for (int w = 1; w < (1 << j); w++) {
        int tb = 31 - __builtin_clz(w);
        for (int up = 0; up < (1 << j); up++) {
            if (up >> tb & 1) continue;
            unsigned m; u64 k = ckey(w, (1 << j) | up, &m);
            if (useMin) k = ~k;
            if (j == 1) { key0 = k; unsigned mm; key0f = pkey(w, (1 << j) | up, &mm); if (useMin) key0f = ~key0f; }
            else if (k > key0) return 0;
            else if (k == key0 && cheapkey < 2) { unsigned mm; u64 f = pkey(w, (1 << j) | up, &mm); if (useMin) f = ~f; if (f > key0f) return 0; }
            pmj[w][up] = m;
        }
    }
    return 1;
}
/* bilinear alternating Dp on span(e_0..e_{k-1}) from coordinates d[i][j] */
static int dbil(int d[M][M], int a, int b) { int r = 0; for (int i = 0; i < M; i++) if (a >> i & 1) for (int j = 0; j < M; j++) if ((b >> j & 1) && i != j) r ^= d[i][j]; return r; }
/* ---- W-list: Dp on W = span(e_0, e_1, e_2) ---- */
typedef struct { uint8_t d01, d02, d12; uint8_t DW[8][8]; } WSol;
static WSol *wl; static int nwl;
static void build_wlist(void) {
    nwl = 0;
    int d[M][M]; memset(d, 0, sizeof d);
    unsigned full = 0xffffffffu;
    unsigned m01 = full & ~pm[1][1][0];
    for (unsigned a = m01; a; a &= a - 1) {
        int d01 = __builtin_ctz(a); d[0][1] = d[1][0] = d01;
        /* planes of column 2 with w = e0: <e0, e2 + up>, up in {0, e1}: Dp = d02 + Dp(e0, up) */
        unsigned m02 = full & ~xshift(pm[2][1][0], 0) & ~xshift(pm[2][1][2], d01);
        for (unsigned b2 = m02; b2; b2 &= b2 - 1) {
            int d02 = __builtin_ctz(b2); d[0][2] = d[2][0] = d02;
            /* w in {e1, e0+e1}, up in {0, e0} (bit 1 clear): Dp(w, e2 + up) = Dp(w, e2) + Dp(w, up), Dp(w, e2) = [w0] d02 + d12 */
            unsigned m12 = full;
            for (int w = 2; w < 4; w++) for (int up = 0; up < 4; up++) {
                if (up >> 1 & 1) continue;
                int known = ((w & 1) ? d02 : 0) ^ dbil(d, w, up);
                m12 &= ~xshift(pm[2][w][up], known);
            }
            for (unsigned c = m12; c; c &= c - 1) {
                int d12 = __builtin_ctz(c); d[1][2] = d[2][1] = d12;
                if (nwl >= 40000) { fprintf(stderr, "W-list overflow\n"); exit(3); }
                WSol *x = &wl[nwl++];
                x->d01 = d01; x->d02 = d02; x->d12 = d12;
                for (int u = 0; u < 8; u++) for (int v = 0; v < 8; v++) x->DW[u][v] = dbil(d, u, v);
            }
            d[1][2] = d[2][1] = 0;
        }
    }
}
/* ---- H0 = span(e_0..e_3) solutions ---- */
typedef struct { int d[M][M]; uint8_t D4[16][16]; } HSol;
static HSol *hl; static int nhl;
static long long cnt_h, cnt_w, cnt_c4, cnt_z;
static void extend_h0(void) {
    nhl = 0;
    static unsigned SH[8][8][Q];   /* forbidden y_w for plane <w, e3 + up>, shifted by Dp(w, up) */
    for (int w = 1; w < 8; w++) { int tb = 31 - __builtin_clz(w); for (int up = 0; up < 8; up++) { if (up >> tb & 1) continue; unsigned *T = SH[w][up]; T[0] = pm[3][w][up];
        for (int s2 = 1; s2 < Q; s2++) { int lb = s2 & -s2; unsigned m = T[s2 ^ lb];
            switch (lb) { case 1: m = ((m & 0x55555555u) << 1) | ((m >> 1) & 0x55555555u); break; case 2: m = ((m & 0x33333333u) << 2) | ((m >> 2) & 0x33333333u); break;
                case 4: m = ((m & 0x0f0f0f0fu) << 4) | ((m >> 4) & 0x0f0f0f0fu); break; case 8: m = ((m & 0x00ff00ffu) << 8) | ((m >> 8) & 0x00ff00ffu); break; default: m = (m << 16) | (m >> 16); }
            T[s2] = m; } } }
    for (int t = 0; t < nwl; t++) {
        WSol *x = &wl[t]; unsigned A[8];
        for (int w = 1; w < 8; w++) {
            int tb = 31 - __builtin_clz(w); unsigned f = 0;
            for (int up = 0; up < 8; up++) { if (up >> tb & 1) continue; f |= SH[w][up][x->DW[w][up]]; }
            A[w] = ~f;
            if (!A[w]) break;
            if (w == 3) { int okp = 0; for (unsigned a0 = A[1]; a0 && !okp; a0 &= a0 - 1) { int y0 = __builtin_ctz(a0); for (unsigned a1 = A[2]; a1; a1 &= a1 - 1) if (A[3] >> (y0 ^ __builtin_ctz(a1)) & 1) { okp = 1; break; } }
                if (!okp) { A[w] = 0; break; } }
        }
        int ok = 1; for (int w = 1; w < 8; w++) if (!A[w]) { ok = 0; break; }
        if (!ok) continue;
        for (unsigned a0 = A[1]; a0; a0 &= a0 - 1) { int y0 = __builtin_ctz(a0);
            for (unsigned a1 = A[2]; a1; a1 &= a1 - 1) { int y1 = __builtin_ctz(a1); if (!(A[3] >> (y0 ^ y1) & 1)) continue;
                for (unsigned a2 = A[4]; a2; a2 &= a2 - 1) { int y2 = __builtin_ctz(a2);
                    if (!(A[5] >> (y0 ^ y2) & 1) || !(A[6] >> (y1 ^ y2) & 1) || !(A[7] >> (y0 ^ y1 ^ y2) & 1)) continue;
                    if (nhl >= 200000) { fprintf(stderr, "H0-list overflow\n"); exit(3); }
                    HSol *h = &hl[nhl++]; memset(h->d, 0, sizeof h->d);
                    h->d[0][1] = h->d[1][0] = x->d01; h->d[0][2] = h->d[2][0] = x->d02; h->d[1][2] = h->d[2][1] = x->d12;
                    h->d[0][3] = h->d[3][0] = y0; h->d[1][3] = h->d[3][1] = y1; h->d[2][3] = h->d[3][2] = y2;
                    for (int u = 0; u < 16; u++) for (int v = 0; v < 16; v++) h->D4[u][v] = dbil(h->d, u, v);
                } } }
    }
    cnt_h += nhl;
}
/* ---- column 4: lazy plane masks, Dp(., e_4) ---- */
static unsigned m4[16][16]; static uint8_t have4[16][16];
static unsigned mask4(int w, int up) { if (!have4[w][up]) { m4[w][up] = forb1(w, 16 | up) | forb1(w, w ^ 16 ^ up) | forb1(w ^ 16 ^ up, 16 | up); have4[w][up] = 1; } return m4[w][up]; }
static unsigned A4(HSol *h, int w) {
    int tb = 31 - __builtin_clz(w); unsigned f = 0;
    for (int up = 0; up < 16; up++) { if (up >> tb & 1) continue; f |= xshift(mask4(w, up), h->D4[w][up]); if (f == 0xffffffffu) break; }
    return ~f;
}
static void finish(HSol *h, int z[4]) {
    memcpy(dpv, h->d, sizeof dpv);
    for (int i = 0; i < 4; i++) dpv[i][4] = dpv[4][i] = z[i];
    leaf();
}
static void do_col4(void) {
    int j = M - 1; col_options(j, &LL[j]);
    for (int t = 0; t < LL[j].n; t++) {
        if (!col_ok(j, LL[j].val[t])) continue;
        cnt_col[j]++; cnt_c4++;
        set_col(j, LL[j].val[t]);
        memset(have4, 0, sizeof have4);
        for (int k = 0; k < nhl; k++) {
            HSol *h = &hl[k]; unsigned A[16]; int ok = 1;
            for (int w = 1; w < 16; w++) { A[w] = A4(h, w); if (!A[w]) { ok = 0; break; } }
            if (!ok) continue;
            cnt_z++;
            int z[4];
            for (unsigned a0 = A[1]; a0; a0 &= a0 - 1) { z[0] = __builtin_ctz(a0);
              for (unsigned a1 = A[2]; a1; a1 &= a1 - 1) { z[1] = __builtin_ctz(a1); if (!(A[3] >> (z[0] ^ z[1]) & 1)) continue;
                for (unsigned a2 = A[4]; a2; a2 &= a2 - 1) { z[2] = __builtin_ctz(a2);
                  if (!(A[5] >> (z[0] ^ z[2]) & 1) || !(A[6] >> (z[1] ^ z[2]) & 1) || !(A[7] >> (z[0] ^ z[1] ^ z[2]) & 1)) continue;
                  for (unsigned a3 = A[8]; a3; a3 &= a3 - 1) { z[3] = __builtin_ctz(a3); int okz = 1;
                    for (int w = 9; w < 16 && okz; w++) { int s = 0; for (int i = 0; i < 4; i++) if (w >> i & 1) s ^= z[i]; if (!(A[w] >> s & 1)) okz = 0; }
                    if (okz) finish(h, z);
                  } } } }
        }
    }
}
static double tcan, th0, tc4;
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static void colS(int j) {
    if (tick) { tick = 0; fprintf(stderr, "progress: col/canon"); for (int k = 0; k < M; k++) fprintf(stderr, " %lld/%lld", cnt_col[k], cnt_canon[k]);
        fprintf(stderr, " | W %lld H0 %lld c4 %lld z %lld | leaves %lld\n", cnt_w, cnt_h, cnt_c4, cnt_z, leaves); alarm(progsec); }
    col_options(j, &LL[j]);
    int n = LL[j].n; if (!n) return;
    int from = 0, to = n;
    if (knuth && j < kdepth) { from = rnd() % n; to = from + 1; kw *= n; }
    if (knuth && j == kdepth) subn++;
    for (int t = from; t < to; t++) {
        if (j == 1 && part_n > 1 && t % part_n != part_r) continue;
        if (!col_ok(j, LL[j].val[t])) continue;
        cnt_col[j]++;
        set_col(j, LL[j].val[t]);
        double t0 = now(); int cc = (j >= 1) ? canon_col(j, pm[j]) : 1; if (j == 3) tcan += now() - t0;
        if (!cc) continue;
        cnt_canon[j]++;
        if (knuth && j < kdepth) kest[j] += kw;
        if (j < 2) colS(j + 1);
        else if (j == 2) { build_wlist(); cnt_w += nwl; if (nwl) colS(3); }
        else { double t1 = now(); extend_h0(); double t2 = now(); th0 += t2 - t1; if (nhl) do_col4(); tc4 += now() - t2; }
    }
}
/* ---------- planted test: F in type 24 given on F_1024 = K(w), w^2 = w + 1, x = q + w p ---------- */
typedef struct { int q, p; } E2;
static E2 mul2(E2 a, E2 b) { E2 r; int pp = fm(a.p, b.p); r.q = fm(a.q, b.q) ^ pp; r.p = fm(a.q, b.p) ^ fm(a.p, b.q) ^ pp; return r; }
static E2 add2(E2 a, E2 b) { E2 r = {a.q ^ b.q, a.p ^ b.p}; return r; }
static E2 pw2(E2 a, int e) { E2 r = {1, 0}; while (e) { if (e & 1) r = mul2(r, a); a = mul2(a, a); e >>= 1; } return r; }
static E2 plantF(int which, E2 x) {
    E2 y = pw2(x, 3);
    if (which == 2) { E2 z = pw2(x, 9), t = {0, 0}; for (int i = 0; i < 2 * M; i++) { t = add2(t, z); z = mul2(z, z); } y = add2(y, t); }
    if (which == 3) { E2 z = pw2(x, 9); y = add2(y, z); }   /* not APN in general; sanity */
    if (which == 4) y = pw2(x, 129);                           /* Gold k = 7: x^5 on F_32 (run with GE=5) */
    if (which == 5) { E2 z = pw2(x, 129 * 3 % 1023), t = {0, 0}; y = pw2(x, 129); (void)z; (void)t; }
    return y;
}
static int plant_test(int which) {
    static int top[Q][Q], bot[Q][Q];
    for (int p = 0; p < Q; p++) for (int q = 0; q < Q; q++) { E2 x = {q, p}; E2 y = plantF(which, x); top[p][q] = y.p; bot[p][q] = y.q; }
    int S0v[M][M];
    for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) { int qq = 1 << i, pp = 1 << j; S0v[i][j] = bot[pp][qq] ^ bot[pp][0] ^ bot[0][qq] ^ bot[0][0]; }
    /* check normal form */
    for (int p = 0; p < Q; p++) for (int q = 0; q < Q; q++) {
        int Spp = 0; for (int i = 0; i < M; i++) if (p >> i & 1) for (int j = 0; j < M; j++) if (p >> j & 1) Spp ^= S0v[i][j];
        if (top[p][q] != (G[p] ^ Spp ^ mu[p][q])) { fprintf(stderr, "plant %d: top not in normal form\n", which); return 1; }
        if (bot[0][q] != G[q]) { fprintf(stderr, "plant %d: bot(0, q) != G(q)\n", which); return 1; }
    }
    { /* Forb statistics for the planted S */
        for (int j = 0; j < M; j++) { int o[M]; for (int i = 0; i < M; i++) o[i] = S0v[i][j]; set_col(j, o); }
        int hs[40] = {0}, ha[8] = {0};
        for (int a = 1; a < Q; a++) for (int b = a + 1; b < Q; b++) { if ((a ^ b) < b) continue; unsigned m; pkey(a, b, &m); hs[__builtin_popcount(m)]++; ha[adim(m) + 1]++; }
        printf("plant %d |Forb| histogram:", which); for (int i = 0; i < 40; i++) if (hs[i]) printf(" %d:%d", i, hs[i]);
        printf("  affine dim:"); for (int i = 0; i < 8; i++) if (ha[i]) printf(" %d:%d", i - 1, ha[i]); printf("\n");
    }
    int npass = 0, nfound = 0;
    for (int l = 1; l < Q; l++) for (int k = 0; k < M; k++) {
        int L[Q], Ninv[Q], N[Q];
        for (int x = 0; x < Q; x++) { int y = x; for (int t = 0; t < k; t++) y = fm(y, y); L[x] = fm(l, y); }
        int l3 = fpow(l, getenv("GE") ? atoi(getenv("GE")) : 3); for (int y = 0; y < Q; y++) { int z = y; for (int t = 0; t < k; t++) z = fm(z, z); Ninv[y] = fm(l3, z); }
        for (int y = 0; y < Q; y++) N[Ninv[y]] = y;
        for (int x = 0; x < Q; x++) if (N[G[L[x]]] != G[x]) { fprintf(stderr, "Aut error\n"); return 1; }
        /* S_g(e_i, e_j) = N S(L e_i, L e_j) */
        u128 v = 0;
        for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) {
            int a = L[1 << i], b = L[1 << j], s = 0;
            for (int ii = 0; ii < M; ii++) if (a >> ii & 1) for (int jj = 0; jj < M; jj++) if (b >> jj & 1) s ^= S0v[ii][jj];
            s = N[s];
            for (int t = 0; t < M; t++) if (s >> t & 1) v |= B1(idx(i, j, t));
        }
        for (int i = 0; i < Tne; i++) if (bit(v, Tep[i])) v ^= Teb[i];
        for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) { int s = 0; for (int t = 0; t < M; t++) s |= bit(v, idx(i, j, t)) << t; fsv[i][j] = s; }
        forced = 1; long long lv = leaves; long long c3 = cnt_canon[M - 2];
        memset(dpv, 0, sizeof dpv);
        colS(0);
        int pass = cnt_canon[M - 2] > c3, found = leaves > lv;
        npass += pass; nfound += found;
    }
    printf("plant %d: %d of %d group images pass the canonical test, %d give leaves (APN total %lld)\n", which, npass, (Q - 1) * M, nfound, apn);
    return 0;
}

static int stat_mode(long n) {
    /* random complete S (no canonical pruning): key statistics over all planes */
    long top1 = 0, topmin = 0, tot = 0; double sumd = 0; long hist[200] = {0};
    for (long it = 0; it < n; it++) {
        for (int j = 0; j < M; j++) { col_options(j, &LL[j]); if (!LL[j].n) break; set_col(j, LL[j].val[rnd() % LL[j].n]); }
        u64 ks[200]; int nk = 0; unsigned m;
        for (int a = 1; a < Q; a++) for (int b = a + 1; b < Q; b++) { if ((a ^ b) < b) continue; /* plane {a, b, a^b} once: a < b < a^b */ ks[nk++] = pkey(a, b, &m); }
        u64 mx = 0, mn = ~0ull; for (int i = 0; i < nk; i++) { if (ks[i] > mx) mx = ks[i]; if (ks[i] < mn) mn = ks[i]; }
        int cmx = 0, cmn = 0, d = 0; for (int i = 0; i < nk; i++) { cmx += ks[i] == mx; cmn += ks[i] == mn; int seen = 0; for (int k = 0; k < i; k++) if (ks[k] == ks[i]) { seen = 1; break; } d += !seen; }
        top1 += cmx; topmin += cmn; tot += nk; sumd += d; hist[cmx]++;
    }
    printf("planes per S %ld, distinct keys %.1f, mean size of max class %.2f, of min class %.2f\n", tot / n, sumd / n, (double)top1 / n, (double)topmin / n);
    printf("max-class size histogram:"); for (int i = 0; i < 200; i++) if (hist[i]) printf(" %d:%ld", i, hist[i]); printf("\n");
    return 0;
}
int main(int argc, char **argv) {
    (void)argc; (void)argv;
    POLY = (M == 3) ? 0xB : (M == 5) ? 0x25 : 0;
    int ge = getenv("GE") ? atoi(getenv("GE")) : 3;
    for (int x = 0; x < Q; x++) G[x] = fpow(x, ge);
    for (int a = 0; a < Q; a++) for (int b = 0; b < Q; b++) mu[a][b] = G[a ^ b] ^ G[a] ^ G[b] ^ G[0];
    for (int p = 0; p < Q; p++) { for (int t = 0; t < Q; t++) solv[p][t] = -1; for (int b = 0; b < Q; b++) if (solv[p][mu[p][b]] < 0) solv[p][mu[p][b]] = b; }
    for (int b = 1; b < Q; b++) {
        int ci = 0; for (int t = 0; t < Q; t++) ci += solv[b][t] >= 0;
        if (ci != Q / 2) { fprintf(stderr, "G not APN\n"); return 1; }
        for (int h = 1; h < Q && !hb[b]; h++) { int ok = 1; for (int y = 0; y < Q; y++) if (__builtin_parity(h & mu[b][y])) { ok = 0; break; } if (ok) hb[b] = h; }
    }
    setup_T();
    { int n = 0; for (int a = 1; a < Q; a++) for (int b = a + 1; b < Q; b++) { if ((a ^ b) < b) continue; pid[a][b] = pid[b][a] = pid[a][a ^ b] = pid[a ^ b][a] = pid[b][a ^ b] = pid[a ^ b][b] = n++; } }
    useMin = getenv("MINKEY") != NULL;
    wl = malloc(sizeof(WSol) * 40000); hl = malloc(sizeof(HSol) * 200000); cheapkey = getenv("CHEAPONLY") ? 2 : 1;
    if (getenv("OUT")) fout = fopen(getenv("OUT"), "w");
    if (getenv("PART")) sscanf(getenv("PART"), "%d/%d", &part_r, &part_n);
    if (getenv("PROG")) { progsec = atoi(getenv("PROG")); signal(SIGALRM, on_alarm); alarm(progsec); }
    if (getenv("STAT")) return stat_mode(atol(getenv("STAT")));
    if (getenv("PLANT")) return plant_test(atoi(getenv("PLANT")));
    clock_t c0 = clock();
    if (getenv("KNUTH")) {
        long ns = atol(getenv("KNUTH")); knuth = 1; if (getenv("KDEPTH")) kdepth = atoi(getenv("KDEPTH"));
        if (getenv("SEED")) rs ^= atoll(getenv("SEED")) * 0x9E3779B97F4A7C15ull;
        for (long t = 0; t < ns; t++) { kw = 1; colS(0); }
        double sec = (double)(clock() - c0) / CLOCKS_PER_SEC;
        printf("Knuth %ld probes, KDEPTH %d: est. canonical nodes per column", ns, kdepth);
        for (int k = 0; k < kdepth && k < M; k++) printf(" %.3g", kest[k] / ns);
        printf(" | exhaustive subtrees %lld, %.4g s each | col/canon below:", subn, sec / (subn ? subn : 1));
        for (int k = kdepth; k < M; k++) printf(" %lld/%lld", cnt_col[k], cnt_canon[k]);
        printf(" | W %lld H0 %lld c4 %lld z %lld | leaves %lld | time canon3 %.1f h0 %.1f col4 %.1f\n", cnt_w, cnt_h, cnt_c4, cnt_z, leaves, tcan, th0, tc4);
        if (kdepth <= M) printf("estimated total: %.3g core-hours\n", kest[kdepth - 1] / ns * sec / (subn ? subn : 1) / 3600 * (subn ? 1 : 0));
        return 0;
    }
    colS(0);
    double sec = (double)(clock() - c0) / CLOCKS_PER_SEC;
    printf("GE %d PART %d/%d: col/canon", ge, part_r, part_n); for (int k = 0; k < M; k++) printf(" %lld/%lld", cnt_col[k], cnt_canon[k]);
    printf(" | W %lld H0 %lld c4 %lld z %lld | leaves %lld APN %lld | %.2f s\n", cnt_w, cnt_h, cnt_c4, cnt_z, leaves, apn, sec);
    if (fout) fclose(fout);
    return 0;
}
