/* stage2.c generalized: U = F_4^KU (bits 0..2KU-1, coordinate k = bits 2k, 2k+1 as 1, omega), Z = F_2^MZ.
   Same algorithm as stage2.c (see there); KU = 3, MZ = 4 is type 47 at n = 10, KU = 2, MZ = 2 the n = 6 analogue
   (A = B = Frobenius x -> x^4 on F_64), used to validate the solver on x^3 (gold6.py writes its U x U part).
   build: cc -O3 -DKU=3 -DMZ=4 -o stage2i stage2i.c        usage: stage2i cls_file */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>

#ifndef KU
#define KU 3
#endif
#ifndef MZ
#define MZ 4
#endif
#define NU (2 * KU)              /* bits of U, U' */
#define NB (NU + MZ)             /* n */
#define SU (1 << NU)
#define SZ (1 << MZ)
#define NPP ((SU - 1) / 3)       /* projective points of U */
#define NPAIR (MZ * (MZ - 1) / 2)
#define NCAND (1 << (NU * (KU - 1)))

static int M4[4][4];
static int sm[4][SU];
static int bUU[SU][SU];
static long long lvl[MZ + 1], nodes, leaves, found, napn;

static int coord(int u, int k) { return (u >> (2 * k)) & 3; }
static int ins(int *b, int r, int x) {   /* insert into descending reduced basis; returns new rank or -1 if dependent */
    for (int j = 0; j < r; j++) if ((x ^ b[j]) < x) x ^= b[j];
    if (!x) return -1;
    b[r++] = x; for (int j = r - 1; j > 0 && b[j] > b[j - 1]; j--) { int t = b[j]; b[j] = b[j - 1]; b[j - 1] = t; }
    return r;
}
static int indep(const int *v, int n) { int b[32], r = 0; for (int i = 0; i < n; i++) if ((r = ins(b, r, v[i])) < 0) return 0; return 1; }
static uint32_t qcur;
static int qval(uint32_t q, int a, int b) {
    int r = 0, p = 0;
    for (int j = 0; j < MZ; j++) for (int i = 0; i < j; i++, p++) {
        int c = ((a >> i & 1) & (b >> j & 1)) ^ ((a >> j & 1) & (b >> i & 1));
        if (c) r ^= (q >> (MZ * p)) & (SZ - 1);
    }
    return r;
}
static int tc[MZ][SU];
static int tau(int u, int z) { int r = 0; for (int j = 0; j < MZ; j++) if (z >> j & 1) r ^= tc[j][u]; return r; }
static int full_beta(int a, int b, uint32_t q) {
    int x = a & (SU - 1), z = a >> NU, y = b & (SU - 1), w = b >> NU;
    return bUU[x][y] ^ tau(x, w) ^ tau(y, z) ^ (qval(q, z, w) << NU);
}
static int dir_ok(int a, uint32_t q) {
    int b[NB], r = 0;
    for (int k = 0; k < NB; k++) { int t = ins(b, r, full_beta(a, 1 << k, q)); if (t > 0) r = t; }
    return r == NB - 1;
}
static int Wb[SU][NU], Wd[SU];
static int PP[NPP], UU[NPP + 1];
static int E0[MZ];
static void fill_col(int j, int c0, int rest) {
    for (int u = 0; u < SU; u++) {
        int r = sm[coord(u, 0)][c0];
        for (int k = 1; k < KU; k++) r ^= sm[coord(u, k)][(rest >> (NU * (k - 1))) & (SU - 1)];
        tc[j][u] = r;
    }
}
static int okB(int u, int nz) {
    int v[NU + MZ], n = 0;
    for (int i = 0; i < Wd[u]; i++) v[n++] = Wb[u][i];
    for (int j = 0; j < nz; j++) v[n++] = tc[j][u];
    return indep(v, n);
}
/* stage2i: AN[t][z0] holds MZ functionals (10-bit masks) whose common kernel is the image beta((u, z0), U), u = UU[t]
   (it has dimension NU by the good table); dir_part reduces with MZ parities. Tz[z][k] = T_z e_k, and
   Qv[z][k] = Q(z, z_k) are kept up to date incrementally. */
#define NBIT (NU + MZ)
static int Tz[SZ][NU], Qv[SZ][MZ], zcode[SZ];
static uint16_t (*ANT)[NPP + 1][MZ];        /* ANT[code][t]: the functionals for T = code, u = UU[t] (good codes only) */
static int setAB(int t, const int *tz, uint16_t *out) {
    int u = UU[t], row[NU], piv[NU];
    for (int k = 0; k < NU; k++) row[k] = bUU[u][1 << k] ^ tz[k];
    for (int k = 0; k < NU; k++) {              /* reduced row echelon form */
        int x = row[k];
        if (!x) return 0;
        int p = 31 - __builtin_clz(x);
        for (int i = 0; i < NU; i++) if (i != k && (row[i] >> p & 1)) row[i] ^= x;
        piv[k] = p;
        for (int i = k + 1; i < NU; i++) if (!row[i]) return 0;
    }
    int pm = 0; for (int k = 0; k < NU; k++) pm |= 1 << piv[k];
    for (int p = 0, n = 0; p < NBIT; p++) if (!(pm >> p & 1)) {
        int f = 1 << p;
        for (int k = 0; k < NU; k++) if (row[k] >> p & 1) f |= 1 << piv[k];
        out[n++] = (uint16_t)f;
    }
    return 1;
}
static int proj(const uint16_t *an, int x) { int y = 0; for (int i = 0; i < MZ; i++) y |= (__builtin_parity(x & an[i])) << i; return y; }
static int dir_part(int t, int z0, int nz) {
    int u = UU[t], b[MZ], r = 0; const uint16_t *an = ANT[zcode[z0]][t];
    for (int k = 0; k < nz; k++) {
        int y = proj(an, tc[k][u] ^ (Qv[z0][k] << NU));
        for (int j = 0; j < r; j++) if ((y ^ b[j]) < y) y ^= b[j];
        if (y) { b[r++] = y; for (int j = r - 1; j > 0 && b[j] > b[j - 1]; j--) { int tt = b[j]; b[j] = b[j - 1]; b[j - 1] = tt; } }
    }
    return r == nz - 1;
}
static void report(void) {
    uint32_t q = qcur;
    for (int a = 1; a < (1 << NB); a++) if (!dir_ok(a, q)) return;
    found++;
    static int F[1 << NB]; static uint8_t seen[1 << NB];
    for (int v = 0; v < (1 << NB); v++) {
        int r = 0;
        for (int j = 0; j < NB; j++) if (v >> j & 1) for (int i = 0; i < j; i++) if (v >> i & 1) r ^= full_beta(1 << i, 1 << j, q);
        F[v] = r;
    }
    int apn = 1;
    for (int a = 1; a < (1 << NB) && apn; a++) {
        memset(seen, 0, sizeof seen);
        for (int x = 0; x < (1 << NB); x++) { int d = F[x] ^ F[x ^ a]; if (seen[d] > 1) { apn = 0; break; } seen[d]++; }
    }
    napn += apn;
    printf("FOUND apn=%d LUT", apn);
    for (int v = 0; v < (1 << NB); v++) printf("%c%d", v ? ',' : ' ', F[v]);
    printf("\n"); fflush(stdout);
}
static int cand[MZ][NCAND], ncand[MZ], *candp[MZ];
/* PREFIX mode: the X representatives with the same first MZ - 1 echelon columns are searched together; they differ
   only in E0[MZ - 1], which enters at the last level: nalt alternatives with their own candidate lists. */
static int nalt, alt[SU * SU], (*candL)[NCAND], *ncandL;
static long long lastq;
/* stage2i: every single-column test (okB for one column, setAB) depends only on one intertwiner T_z, coded as
   (T e_0, T e_1, ..., T e_{KU-1}) in NU * KU bits (linear in T). good[code] = 1 iff T u is outside W_u for every
   point u and beta((u, z), .) restricted to U has full rank for u = 0 and every point u. The level-j test of a
   candidate then is good[code(T_zr) ^ code(T_c)] for every zr in span(z_0..z_{j-1}). */
#define NCODE (1 << (NU * KU))
static uint8_t *good;
static int tmp[SU];
static int gtest(int code) {
    for (int u = 0; u < SU; u++) {
        int r = sm[coord(u, 0)][code & (SU - 1)];
        for (int k = 1; k < KU; k++) r ^= sm[coord(u, k)][(code >> (NU * k)) & (SU - 1)];
        tmp[u] = r;
    }
    for (int i = 0; i < NPP; i++) {
        int u = PP[i], v[NU + 1], n = 0;
        for (int a = 0; a < Wd[u]; a++) v[n++] = Wb[u][a];
        v[n++] = tmp[u];
        if (!indep(v, n)) return 0;
    }
    for (int t = 0; t <= NPP; t++) {
        int u = UU[t], b[NU], r = 0;
        for (int k = 0; k < NU; k++) if ((r = ins(b, r, bUU[u][1 << k] ^ tmp[1 << k])) < 0) return 0;
    }
    return 1;
}
static int setAB(int t, const int *tz, uint16_t *out);
static void make_good(void) {
    good = malloc(NCODE); ANT = calloc(NCODE, sizeof *ANT);
    for (int c = 0; c < NCODE; c++) {
        good[c] = gtest(c);
        if (!good[c]) continue;
        int tz[NU]; for (int k = 0; k < NU; k++) tz[k] = tmp[1 << k];
        for (int t = 0; t <= NPP; t++) if (!setAB(t, tz, ANT[c][t])) { fprintf(stderr, "good table mismatch\n"); exit(4); }
    }
}
static int ccode[MZ];
static int pidx(int i, int j) { return j * (j - 1) / 2 + i; }
static void dfs(int j);
static void dfs_level(int j, const int *zc);
/* the images for z0 = zr + z_j are needed only in the last Q step of level j (i = j - 1): built lazily, once per candidate */

/* rank4[w]: rank of the MZ nibbles of w (MZ <= 4 vectors of F_2^MZ, MZ <= 4) */
static uint8_t rank4[1 << 16];
static void make_rank4(void) {
    for (int w = 0; w < (1 << 16); w++) {
        int b[4], r = 0;
        for (int k = 0; k < 4; k++) {
            int y = (w >> (4 * k)) & 15;
            for (int a = 0; a < r; a++) if ((y ^ b[a]) < y) y ^= b[a];
            if (y) { b[r++] = y; for (int a = r - 1; a > 0 && b[a] > b[a - 1]; a--) { int tt = b[a]; b[a] = b[a - 1]; b[a - 1] = tt; } }
        }
        rank4[w] = r;
    }
}
/* allowed values of Q(z_i, z_j) for direction (UU[t], z0) on U x span(z_0..z_j), given everything else: bit c set
   iff dir_part would pass with Q(z_i, z_j) = c (cur is the value currently stored in Qv) */
static int allowed(int t, int z0, int i, int j, int cur) {
    int u = UU[t], nz = j + 1, w = 0; const uint16_t *an = ANT[zcode[z0]][t];
    for (int k = 0; k < nz; k++) w |= proj(an, tc[k][u] ^ (Qv[z0][k] << NU)) << (4 * k);
    int dm = ((z0 >> i & 1) ? 15 << (4 * j) : 0) | ((z0 >> j & 1) ? 15 << (4 * i) : 0);
    int pd[SZ], m = 0; pd[0] = 0;
    for (int a = 0; a < MZ; a++) { int p = proj(an, 1 << (NU + a)); for (int d = 0; d < (1 << a); d++) pd[d | (1 << a)] = pd[d] ^ p; }
    for (int d = 0; d < SZ; d++) if (rank4[w ^ ((pd[d] * 0x1111) & dm)] == nz - 1) m |= 1 << (d ^ cur);
    return m;
}
static void pstep(int i, int j) {
    if (i == j) { dfs(j + 1); return; }
    uint32_t save = qcur;
    int cur = (qcur >> (MZ * pidx(i, j))) & (SZ - 1), c0 = cur, ok = (1 << SZ) - 1;
    for (int z0 = 1 << i; z0 < (2 << i) && ok; z0++)
        for (int t = 0; t <= NPP && ok; t++) ok &= allowed(t, z0, i, j, cur);
    if (ok && i == j - 1)
        for (int z0 = 1 << j; z0 < (2 << j) && ok; z0++)
            for (int t = 0; t <= NPP && ok; t++) ok &= allowed(t, z0, i, j, cur);
    for (uint32_t c = 0; c < SZ; c++) {
        if (!(ok >> c & 1)) continue;
        qcur = (save & ~((uint32_t)(SZ - 1) << (MZ * pidx(i, j)))) | (c << (MZ * pidx(i, j)));
        int d = (int)c ^ cur; cur = c;          /* Q(z_i, z_j) changed by d: rows of z containing z_i / z_j */
        for (int z = 0; z < SZ; z++) { if (z >> i & 1) Qv[z][j] ^= d; if (z >> j & 1) Qv[z][i] ^= d; }
        if (j == MZ - 1) lastq++;
        pstep(i + 1, j);
    }
    qcur = save;
    { int d = cur ^ c0; for (int z = 0; z < SZ; z++) { if (z >> i & 1) Qv[z][j] ^= d; if (z >> j & 1) Qv[z][i] ^= d; } }
}
static void dfs(int j) {
    nodes++; lvl[j]++;
    if (j == MZ) { leaves++; report(); return; }
    int zc[1 << MZ]; zc[0] = 0;
    for (int zr = 1; zr < (1 << j); zr++) zc[zr] = zc[zr & (zr - 1)] ^ ccode[__builtin_ctz(zr)];
    if (j == MZ - 1 && nalt) {
        for (int a = 0; a < nalt; a++) { E0[j] = alt[a]; candp[j] = candL[a]; ncand[j] = ncandL[a]; dfs_level(j, zc); }
        return;
    }
    dfs_level(j, zc);
}
static void dfs_level(int j, const int *zc) {
    for (int ci = 0; ci < ncand[j]; ci++) {
        int code = E0[j] | (candp[j][ci] << NU), ok = 1;
        for (int zr = 1; zr < (1 << j); zr++) if (!good[zc[zr] ^ code]) { ok = 0; break; }
        if (!ok) continue;
        ccode[j] = code;
        fill_col(j, E0[j], candp[j][ci]);
        for (int zr = 0; zr < (1 << j); zr++) zcode[zr | (1 << j)] = zc[zr] ^ code;
        if (j) pstep(0, j); else dfs(1);
    }
}
static void make_cand(void) {
    for (int j = 0; j < MZ; j++) {
        ncand[j] = 0; candp[j] = cand[j];
        for (int c = 0; c < NCAND; c++) {
            if (good[E0[j] | (c << NU)]) cand[j][ncand[j]++] = c;
        }
    }
}
static void read_class(const char *fn) {
    int b[NU][NU]; memset(b, 0, sizeof b);
    FILE *f = fopen(fn, "r"); if (!f) { perror(fn); exit(1); }
    int p, v;
    while (fscanf(f, "%d %d", &p, &v) == 2) {
        int j = 1; while (j * (j + 1) / 2 <= p) j++;
        int i = p - j * (j - 1) / 2;
        b[i][j] = b[j][i] = v;
    }
    fclose(f);
    for (int x = 0; x < SU; x++) for (int y = 0; y < SU; y++) {
        int r = 0;
        for (int i = 0; i < NU; i++) if (x >> i & 1) for (int j = 0; j < NU; j++) if (y >> j & 1) r ^= b[i][j];
        bUU[x][y] = r;
    }
}
static int setup_W(void) {
    for (int u = 1; u < SU; u++) {
        int ker = 0, b[NU], r = 0;
        for (int x = 0; x < SU; x++) {
            if (bUU[u][x] >> NU) continue;
            int c = bUU[u][x] & (SU - 1);
            if (!c) ker++;
            int t = ins(b, r, c); if (t > 0) r = t;
        }
        if (ker != 2) return 0;
        Wd[u] = r; memcpy(Wb[u], b, sizeof b);
        if (r + MZ > NU) return 0;
    }
    return 1;
}
static int GO[100000][KU], nGO;
static int applyP(const int *P, int v) { int r = 0; for (int k = 0; k < KU; k++) r ^= sm[coord(v, k)][P[k]]; return r; }
static void out_stab(void) {
    nGO = 0; int P[KU];
    for (long long code = 0; code < (1LL << (NU * KU)); code++) {
        for (int k = 0; k < KU; k++) P[k] = (code >> (NU * k)) & (SU - 1);
        static uint8_t hit[SU]; memset(hit, 0, sizeof hit); int ok = 1;
        for (int v = 0; v < SU; v++) { int w = applyP(P, v); if (hit[w]) { ok = 0; break; } hit[w] = 1; }
        if (!ok) continue;
        for (int x = 1; x < SU && ok; x <<= 1) for (int y = 1; y < SU; y <<= 1)
            if (applyP(P, bUU[x][y] & (SU - 1)) != (bUU[x][y] & (SU - 1))) { ok = 0; break; }
        if (!ok) continue;
        memcpy(GO[nGO++], P, sizeof P);
    }
}
static void spanset(const int *c, uint64_t *m) {
    memset(m, 0, sizeof(uint64_t) * ((SU + 63) / 64));
    for (int s = 0; s < SZ; s++) { int v = 0; for (int i = 0; i < MZ; i++) if (s >> i & 1) v ^= c[i]; m[v >> 6] |= 1ull << (v & 63); }
}
static int less(const uint64_t *a, const uint64_t *b) {
    for (int i = (SU + 63) / 64 - 1; i >= 0; i--) if (a[i] != b[i]) return a[i] < b[i];
    return 0;
}
static int xrep(const int *c) {
    uint64_t m[(SU + 63) / 64], m2[(SU + 63) / 64];
    spanset(c, m);
    for (int g = 0; g < nGO; g++) { int d[MZ]; for (int i = 0; i < MZ; i++) d[i] = applyP(GO[g], c[i]); spanset(d, m2); if (less(m2, m)) return 0; }
    return 1;
}
static long long nX, nXall;
static const char *xonly;
static int (*reps)[MZ], nreps;
static void runX(const int *c) {
    nXall++;
    if (!xrep(c)) return;
    nX++;
    if (getenv("PREFIX") || getenv("LISTP")) { memcpy(reps[nreps++], c, sizeof reps[0]); return; }
    if (xonly && atoi(xonly) != nX) return;
    memcpy(E0, c, sizeof E0);
    make_cand();
    if (getenv("SHOWX")) { fprintf(stderr, "X %lld:", nX); for (int i = 0; i < MZ; i++) fprintf(stderr, " %d", c[i]); fprintf(stderr, " cand"); for (int i = 0; i < MZ; i++) fprintf(stderr, " %d", ncand[i]); fprintf(stderr, "\n"); }
    dfs(0);
}
/* enumerate reduced echelon bases (c_0 < ... < c_{MZ-1}) of MZ-dim X with X + W_{f_0} direct */
static void enumX(int *c, int i, int lo) {
    if (i == MZ) {
        for (int a = 0; a < MZ; a++) for (int k = 0; k < MZ; k++)
            if (k != a && (c[k] >> (31 - __builtin_clz(c[a])) & 1)) return;
        int v[NU + MZ], n = 0;
        for (int a = 0; a < Wd[1]; a++) v[n++] = Wb[1][a];
        for (int a = 0; a < MZ; a++) v[n++] = c[a];
        if (indep(v, n)) runX(c);
        return;
    }
    for (int x = lo; x < SU; x++) {
        if (i && 31 - __builtin_clz(x) == 31 - __builtin_clz(c[i - 1])) continue;
        c[i] = x; enumX(c, i + 1, x + 1);
    }
}
int main(int argc, char **argv) {
    (void)argc;
    for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) {
        int r = 0; for (int i = 0; i < 2; i++) if (b >> i & 1) r ^= a << i;
        if (r & 4) r ^= 7;
        M4[a][b] = r;
    }
    for (int c = 0; c < 4; c++) for (int v = 0; v < SU; v++) {
        int r = 0; for (int k = 0; k < KU; k++) r |= M4[c][coord(v, k)] << (2 * k);
        sm[c][v] = r;
    }
    { int n = 0; for (int u = 1; u < SU; u++) { int k = 0; while (!coord(u, k)) k++; if (coord(u, k) == 1) PP[n++] = u; } }
    UU[0] = 0; for (int t = 1; t <= NPP; t++) UU[t] = PP[t - 1];
    xonly = getenv("XONLY");
    read_class(argv[1]);
    for (int x = 0; x < SU; x++) for (int y = 0; y < SU; y++) {
        if (bUU[x][x]) { fprintf(stderr, "not alternating\n"); return 1; }
        int l = bUU[sm[2][x]][sm[2][y]], r = bUU[x][y];
        if ((l >> NU) != (r >> NU) || (l & (SU - 1)) != sm[2][r & (SU - 1)]) { fprintf(stderr, "not equivariant\n"); return 1; }
    }
    if (!setup_W()) { printf("%s: dead on U x U\n", argv[1]); return 0; }
    out_stab();
    make_good(); make_rank4();
    reps = malloc(sizeof *reps * 200000);
    int c[MZ]; enumX(c, 0, 1);
    if (getenv("PREFIX") || getenv("LISTP")) {
        /* groups in order of first appearance (enumX is lexicographic, so a group is a run of consecutive reps) */
        int g = 0, want = getenv("PREFIX") ? atoi(getenv("PREFIX")) : 0, first = -1;
        for (int r = 0; r < nreps; r++) {
            if (r == 0 || memcmp(reps[r], reps[r - 1], sizeof(int) * (MZ - 1))) { g++; if (getenv("LISTP")) printf("%sgroup %d:", g > 1 ? "\n" : "", g); }
            if (getenv("LISTP")) printf(" %d", r + 1);
            if (g == want) { if (first < 0) first = r; alt[nalt++] = reps[r][MZ - 1]; }
        }
        if (getenv("LISTP")) { printf("\n%s: groups %d, X %d\n", argv[1], g, nreps); return 0; }
        if (first < 0) { fprintf(stderr, "no prefix %d\n", want); return 1; }
        memcpy(E0, reps[first], sizeof E0);
        make_cand();
        candL = malloc(sizeof *candL * nalt); ncandL = malloc(sizeof *ncandL * nalt);
        for (int a = 0; a < nalt; a++) {
            ncandL[a] = 0;
            for (int cc = 0; cc < NCAND; cc++) if (good[alt[a] | (cc << NU)]) candL[a][ncandL[a]++] = cc;
        }
        dfs(0);
        printf("%s: stab %d, prefix %d (%d X), nodes %lld, last-level Q passes %lld, tau leaves %lld, found %lld (APN %lld)\n", argv[1], nGO, want, nalt, nodes, lastq, leaves, found, napn);
        return 0;
    }
    printf("%s: stab %d, X %lld of %lld, nodes %lld, last-level Q passes %lld, tau leaves %lld, found %lld (APN %lld)\n", argv[1], nGO, nX, nXall, nodes, lastq, leaves, found, napn);
    return 0;
}
