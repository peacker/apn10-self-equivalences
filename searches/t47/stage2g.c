/* stage2.c generalized: U = F_4^KU (bits 0..2KU-1, coordinate k = bits 2k, 2k+1 as 1, omega), Z = F_2^MZ.
   Same algorithm as stage2.c (see there); KU = 3, MZ = 4 is type 47 at n = 10, KU = 2, MZ = 2 the n = 6 analogue
   (A = B = Frobenius x -> x^4 on F_64), used to validate the solver on x^3 (gold6.py writes its U x U part).
   build: cc -O3 -DKU=3 -DMZ=4 -o stage2g stage2g.c        usage: stage2g cls_file */
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
static int AB[NPP + 1][SZ][NU];
static int setAB(int t, int z0) {
    int u = UU[t], *b = AB[t][z0], r = 0;
    for (int k = 0; k < NU; k++) if ((r = ins(b, r, bUU[u][1 << k] ^ tau(1 << k, z0))) < 0) return 0;
    return 1;
}
static int dir_part(int t, int z0, int nz) {
    int u = UU[t], *ab = AB[t][z0], b[MZ], r = 0;
    for (int k = 0; k < nz; k++) {
        int x = tc[k][u] ^ (qval(qcur, z0, 1 << k) << NU);
        for (int j = 0; j < NU; j++) if ((x ^ ab[j]) < x) x ^= ab[j];
        int t2 = ins(b, r, x); if (t2 > 0) r = t2;
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
static int cand[MZ][NCAND], ncand[MZ];
static int pidx(int i, int j) { return j * (j - 1) / 2 + i; }
static void dfs(int j);
static void pstep(int i, int j) {
    if (i == j) { dfs(j + 1); return; }
    uint32_t save = qcur;
    for (uint32_t c = 0; c < SZ; c++) {
        qcur = (save & ~((uint32_t)(SZ - 1) << (MZ * pidx(i, j)))) | (c << (MZ * pidx(i, j)));
        int ok = 1;
        for (int z0 = 1 << i; z0 < (2 << i) && ok; z0++)
            for (int t = 0; t <= NPP; t++) if (!dir_part(t, z0, j + 1)) { ok = 0; break; }
        if (ok && i == j - 1)
            for (int z0 = 1 << j; z0 < (2 << j) && ok; z0++)
                for (int t = 0; t <= NPP; t++) if (!dir_part(t, z0, j + 1)) { ok = 0; break; }
        if (ok) pstep(i + 1, j);
    }
    qcur = save;
}
static void dfs(int j) {
    nodes++; lvl[j]++;
    if (j == MZ) { leaves++; report(); return; }
    for (int ci = 0; ci < ncand[j]; ci++) {
        fill_col(j, E0[j], cand[j][ci]);
        int ok = 1;
        if (j) for (int i = 0; i < NPP && ok; i++) if (!okB(PP[i], j + 1)) ok = 0;
        for (int zr = 0; zr < (1 << j) && ok; zr++)
            for (int t = 0; t <= NPP; t++) if (!setAB(t, zr | (1 << j))) { ok = 0; break; }
        if (ok) { if (j) pstep(0, j); else dfs(1); }
    }
}
static void make_cand(void) {
    for (int j = 0; j < MZ; j++) {
        ncand[j] = 0;
        for (int c = 0; c < NCAND; c++) {
            fill_col(0, E0[j], c);
            int ok = 1;
            for (int i = 0; i < NPP && ok; i++) if (!okB(PP[i], 1)) ok = 0;
            for (int t = 0; t <= NPP && ok; t++) if (!setAB(t, 1)) ok = 0;
            if (ok) cand[j][ncand[j]++] = c;
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
static void runX(const int *c) {
    nXall++;
    if (!xrep(c)) return;
    nX++;
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
    int c[MZ]; enumX(c, 0, 1);
    printf("%s: stab %d, X %lld of %lld, nodes %lld, tau leaves %lld, found %lld (APN %lld)\n", argv[1], nGO, nX, nXall, nodes, leaves, found, napn);
    return 0;
}
