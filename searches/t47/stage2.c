/* Type 47 stage 2, dedicated solver (replaces the generic sym10f run for one forced U x U class).
   Input U = F_4^3 (bits 0..5, coordinate k = bits 2k, 2k+1 as 1, omega), Z = Fix(A) = F_2^4; output U' (bits 0..5),
   Z' (bits 6..9).  Equivariance leaves exactly: beta on U x U (forced class (C, h), cls/c*.txt), the intertwiner
   tau(u, z) = T_z u (F_4-linear in u, U x Z -> U'), and Q: Z x Z -> Z' (U x Z -> Z' and Z x Z -> U' vanish).
   Necessary conditions without Q:
     (B)  direction (u, 0): {C(u, x) : h(u, x) = 0} has dim dim K_u - 1 and meets tau(u, Z) (injective) trivially;
     (A') direction (u, z0), z0 != 0: x -> (C(u, x) + T_z0 x, h(u, x)) injective on U (u = 0: T_z0 nonsingular).
   GL(Z) is used completely: E_0 = tau(f_0, .) (injective by (B)) is put in column echelon form, i.e. it is the
   echelon basis of its column space X (4-dim, X meets W_{f_0} trivially).  Then the columns j = 0..3 of
   E_1 = tau(f_1, .), E_2 = tau(f_2, .) are chosen (12 bits per column) with (A'), (B) checked as soon as defined.
   At a leaf every alternating Q with Q(a, .) of rank 3 (Gold type on F_16) is tried with the full APN test on all
   directions (u, z0), z0 != 0; survivors are checked again on the LUT and printed.
   usage: stage2 cls/cXX.txt          stage2 -test N   (random consistency test of (A'), (B) against the full test) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>

static int M4[4][4];
static int sm[4][64];            /* F_4 scalar times a vector of U' (or U) */
static int bUU[64][64];          /* beta on U x U, 10 bits */
static int nQ; static uint32_t *Qs;   /* Q packed: Q(z_i, z_j) for the 6 pairs (i < j), 4 bits each */
static long long lvl[5], nodes, leaves, qtests, found;

static int coord(int u, int k) { return (u >> (2 * k)) & 3; }
static int indep(const int *v, int n) {   /* are v[0..n-1] linearly independent */
    int b[16], r = 0;
    for (int i = 0; i < n; i++) {
        int x = v[i];
        for (int j = 0; j < r; j++) if ((x ^ b[j]) < x) x ^= b[j];
        if (!x) return 0;
        b[r++] = x; for (int j = r - 1; j > 0 && b[j] > b[j - 1]; j--) { int t = b[j]; b[j] = b[j - 1]; b[j - 1] = t; }
    }
    return 1;
}
static int qval(uint32_t q, int a, int b) {   /* Q(a, b) for a, b in F_2^4 */
    int r = 0, p = 0;
    for (int j = 0; j < 4; j++) for (int i = 0; i < j; i++, p++) {
        int c = ((a >> i & 1) & (b >> j & 1)) ^ ((a >> j & 1) & (b >> i & 1));
        if (c) r ^= (q >> (4 * p)) & 15;
    }
    return r;
}
/* tau as tables: tc[j][u] = tau(u, z_j) */
static int tc[4][64];
static int tau(int u, int z) { int r = 0; for (int j = 0; j < 4; j++) if (z >> j & 1) r ^= tc[j][u]; return r; }
static int full_beta(int a, int b) {   /* 10-bit input: bits 0..5 U, 6..9 Z */
    int x = a & 63, z = a >> 6, y = b & 63, w = b >> 6;
    return bUU[x][y] ^ tau(x, w) ^ tau(y, z);
}
static int dir_ok(int a, uint32_t q) {  /* ker beta(a, .) = {0, a} with Q = q */
    int v[10];
    for (int k = 0; k < 10; k++) v[k] = full_beta(a, 1 << k) ^ (qval(q, a >> 6, (1 << k) >> 6) << 6);
    /* kernel dim 1 <=> rank 9 */
    int b[10], r = 0;
    for (int i = 0; i < 10; i++) {
        int x = v[i];
        for (int j = 0; j < r; j++) if ((x ^ b[j]) < x) x ^= b[j];
        if (x) { b[r++] = x; for (int j = r - 1; j > 0 && b[j] > b[j - 1]; j--) { int t = b[j]; b[j] = b[j - 1]; b[j - 1] = t; } }
    }
    return r == 9;
}
/* per-u data for (B) */
static int Wb[64][6], Wd[64];
static int okB(int u, int nz) {    /* W_u basis + tau(u, z_0..z_{nz-1}) independent */
    int v[12], n = 0;
    for (int i = 0; i < Wd[u]; i++) v[n++] = Wb[u][i];
    for (int j = 0; j < nz; j++) v[n++] = tc[j][u];
    return indep(v, n);
}
static int okA(int u, int z0) {    /* x -> (C(u, x) + T_z0 x, h(u, x)) injective */
    int v[6];
    for (int k = 0; k < 6; k++) v[k] = bUU[u][1 << k] ^ tau(1 << k, z0);
    return indep(v, 6);
}
static int E0[4];   /* columns of E_0 (6-bit) */
static void fill_col(int j, int c0, int c1, int c2) {
    for (int u = 0; u < 64; u++) tc[j][u] = sm[coord(u, 0)][c0] ^ sm[coord(u, 1)][c1] ^ sm[coord(u, 2)][c2];
}
static int PP[21];   /* projective points of U (first nonzero F_4-coordinate 1): directions up to omega */
static uint32_t qcur;
static int checkQ(int z0) { for (int i = 0; i < 21; i++) if (!dir_ok(PP[i] | z0 << 6, qcur)) return 0; return 1; }
static void report(void) {
    uint32_t q = qcur;
    for (int a = 1; a < 1024; a++) if (!dir_ok(a, q)) return;
    found++;
    /* LUT: F(v) = sum_{i<j} v_i v_j beta(e_i, e_j); verify APN directly */
    static int F[1024]; static uint8_t seen[1024];
    for (int v = 0; v < 1024; v++) {
        int r = 0;
        for (int j = 0; j < 10; j++) if (v >> j & 1) for (int i = 0; i < j; i++) if (v >> i & 1)
            r ^= full_beta(1 << i, 1 << j) ^ (qval(q, (1 << i) >> 6, (1 << j) >> 6) << 6);
        F[v] = r;
    }
    int apn = 1;
    for (int a = 1; a < 1024 && apn; a++) {
        memset(seen, 0, sizeof seen);
        for (int x = 0; x < 1024; x++) { int d = F[x] ^ F[x ^ a]; if (seen[d] > 1) { apn = 0; break; } seen[d]++; }
    }
    printf("FOUND apn=%d q=%u E:", apn, q);
    for (int j = 0; j < 4; j++) printf(" %d,%d,%d", tc[j][1], tc[j][4], tc[j][16]);
    printf("\nLUT");
    for (int v = 0; v < 1024; v++) printf("%c%d", v ? ',' : ' ', F[v]);
    printf("\n"); fflush(stdout);
}
/* Q packed as in qval: pair p = (i, j), i < j, index j(j-1)/2 + i: (0,1)=0 (0,2)=1 (1,2)=2 (0,3)=3 (1,3)=4 (2,3)=5 */
static void leaf(void) {
    leaves++;
    if (getenv("PROG") && (leaves & (leaves - 1)) == 0) fprintf(stderr, "leaf %lld (nodes %lld %lld %lld %lld)\n", leaves, lvl[0], lvl[1], lvl[2], lvl[3]);
    if (getenv("NOQ")) return;
    for (uint32_t a = 0; a < 4096; a++) {          /* Q(z0, z1), Q(z0, z2), Q(z0, z3) */
        qcur = (a & 255) | ((a >> 8) << 12); qtests++;
        if (!checkQ(1)) continue;
        for (uint32_t b = 0; b < 256; b++) {       /* Q(z1, z2), Q(z1, z3) */
            qcur = (a & 255) | ((b & 15) << 8) | ((a >> 8) << 12) | ((b >> 4) << 16);
            if (!checkQ(2) || !checkQ(3)) continue;
            for (uint32_t c = 0; c < 16; c++) {     /* Q(z2, z3) */
                qcur = (qcur & 0xFFFFF) | (c << 20);
                int ok = 1;
                for (int z0 = 4; z0 < 16 && ok; z0++) ok = checkQ(z0);
                if (ok) report();
            }
        }
    }
}
static int cand[4][4096], ncand[4];
static int AB[22][16][6];   /* reduced basis (descending) of {C(u, x) + T_z0 x | h(u, x)} for t = 0 (u = 0), t = 1..21 (PP[t - 1]) */
static int UU[22];
static int setAB(int t, int z0) {
    int u = UU[t], *b = AB[t][z0], r = 0;
    for (int k = 0; k < 6; k++) {
        int x = bUU[u][1 << k] ^ tau(1 << k, z0);
        for (int j = 0; j < r; j++) if ((x ^ b[j]) < x) x ^= b[j];
        if (!x) return 0;
        b[r++] = x; for (int j = r - 1; j > 0 && b[j] > b[j - 1]; j--) { int tt = b[j]; b[j] = b[j - 1]; b[j - 1] = tt; }
    }
    return 1;
}
/* direction (u, z0) restricted to the domain U x span(z_0..z_{nz-1}) with the current tau, qcur: kernel = {0, a}? */
static int dir_part(int t, int z0, int nz) {   /* uses AB[t][z0] (x-part injective) */
    int u = UU[t], *ab = AB[t][z0], b[4], r = 0;
    for (int k = 0; k < nz; k++) {
        int x = tc[k][u] ^ (qval(qcur, z0, 1 << k) << 6);
        for (int j = 0; j < 6; j++) if ((x ^ ab[j]) < x) x ^= ab[j];
        for (int j = 0; j < r; j++) if ((x ^ b[j]) < x) x ^= b[j];
        if (x) { b[r++] = x; for (int j = r - 1; j > 0 && b[j] > b[j - 1]; j--) { int tt = b[j]; b[j] = b[j - 1]; b[j - 1] = tt; } }
    }
    return r == nz - 1;
}
static long long pst[4][4];
static int pidx(int i, int j) { return j * (j - 1) / 2 + i; }
static void dfs(int j);
/* pair step: choose Q(z_i, z_j) */
static void pstep(int i, int j) {
    if (i == j) { dfs(j + 1); return; }
    uint32_t save = qcur; pst[i][j]++;
    for (uint32_t c = 0; c < 16; c++) {
        qcur = (save & ~(15u << (4 * pidx(i, j)))) | (c << (4 * pidx(i, j)));
        int ok = 1;
        for (int z0 = 1 << i; z0 < (2 << i) && ok; z0++) {       /* top bit of z0 = i */
            for (int t = 0; t < 22; t++) if (!dir_part(t, z0, j + 1)) { ok = 0; break; }
        }
        if (ok && i == j - 1)
            for (int z0 = 1 << j; z0 < (2 << j) && ok; z0++) {   /* top bit j: Q(z_j, .) now known */
                for (int t = 0; t < 22; t++) if (!dir_part(t, z0, j + 1)) { ok = 0; break; }
            }
        if (ok) pstep(i + 1, j);
    }
    qcur = save;
}
static void dfs(int j) {
    nodes++; lvl[j]++;
    if (j == 4) { leaves++; report(); return; }
    for (int ci = 0; ci < ncand[j]; ci++) {
        int c12 = cand[j][ci];
        fill_col(j, E0[j], c12 & 63, c12 >> 6);
        int ok = 1;
        if (j) for (int i = 0; i < 21 && ok; i++) if (!okB(PP[i], j + 1)) ok = 0;
        for (int zr = 0; zr < (1 << j) && ok; zr++) {    /* (A') and cache; zr = 0 passed already in make_cand */
            int z0 = zr | (1 << j);
            for (int t = 0; t < 22; t++) if (!setAB(t, z0)) { ok = 0; break; }
        }
        if (ok) { if (j) pstep(0, j); else dfs(1); }
    }
}
/* single-column candidates for column j: B with tau(u, z_j) alone, A' at z0 = z_j */
static void make_cand(void) {
    for (int j = 0; j < 4; j++) {
        ncand[j] = 0;
        for (int c12 = 0; c12 < 4096; c12++) {
            fill_col(0, E0[j], c12 & 63, c12 >> 6);
            int ok = 1, sv = 0;
            for (int i = 0; i < 21 && ok; i++) if (!okB(PP[i], 1)) ok = 0;
            (void)sv;
            if (ok && !okA(0, 1)) ok = 0;
            for (int i = 0; i < 21 && ok; i++) if (!okA(PP[i], 1)) ok = 0;
            if (ok) cand[j][ncand[j]++] = c12;
        }
    }
}
static void setup_field(void) {
    for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) {
        int r = 0; for (int i = 0; i < 2; i++) if (b >> i & 1) r ^= a << i;
        if (r & 4) r ^= 7;
        M4[a][b] = r;
    }
    for (int c = 0; c < 4; c++) for (int v = 0; v < 64; v++) {
        int r = 0; for (int k = 0; k < 3; k++) r |= M4[c][coord(v, k)] << (2 * k);
        sm[c][v] = r;
    }
}
static void gen_Q(void) {
    Qs = malloc(sizeof(uint32_t) * (1u << 24)); nQ = 0;
    for (uint32_t q = 0; q < (1u << 24); q++) {
        int ok = 1;
        for (int a = 1; a < 16 && ok; a++) {
            int v[4], n = 0;
            for (int k = 0; k < 4; k++) if (k != __builtin_ctz(a)) v[n++] = qval(q, a, 1 << k);
            if (!indep(v, 3)) ok = 0;
        }
        if (ok) Qs[nQ++] = q;
    }
}
static void read_class(const char *fn) {
    int b[6][6]; memset(b, 0, sizeof b);
    FILE *f = fopen(fn, "r"); if (!f) { perror(fn); exit(1); }
    int p, v;
    while (fscanf(f, "%d %d", &p, &v) == 2) {
        int j = 1; while (j * (j + 1) / 2 <= p) j++;
        int i = p - j * (j - 1) / 2;
        b[i][j] = b[j][i] = v;
    }
    fclose(f);
    for (int x = 0; x < 64; x++) for (int y = 0; y < 64; y++) {
        int r = 0;
        for (int i = 0; i < 6; i++) if (x >> i & 1) for (int j = 0; j < 6; j++) if (y >> j & 1) r ^= b[i][j];
        bUU[x][y] = r;
    }
}
/* returns 0 if the class is dead already on U x U */
static int setup_W(void) {
    for (int u = 1; u < 64; u++) {
        int K[64], nk = 0;
        for (int x = 0; x < 64; x++) if (!(bUU[u][x] >> 6)) K[nk++] = x;
        int ker = 0, b[6], r = 0;
        for (int i = 0; i < nk; i++) {
            int c = bUU[u][K[i]] & 63;
            if (!c) ker++;
            for (int j = 0; j < r; j++) if ((c ^ b[j]) < c) c ^= b[j];
            if (c) { b[r++] = c; for (int j = r - 1; j > 0 && b[j] > b[j - 1]; j--) { int t = b[j]; b[j] = b[j - 1]; b[j - 1] = t; } }
        }
        if (getenv("DBG")) fprintf(stderr, "u %d |K| %d ker %d r %d\n", u, nk, ker, r);
        if (ker != 2) return 0;      /* {x in K_u : C(u, x) = 0} must be {0, u} */
        Wd[u] = r; memcpy(Wb[u], b, sizeof b);
        if (r + 4 > 6) return 0;     /* rank h(u, .) >= 3 needed */
    }
    return 1;
}
static uint64_t rs = 88172645463325252ull;
static uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static int curX;
static int GO[48 * 64][3], nGO;   /* P' in GL_3(F_4) with P' C = C (columns as 6-bit vectors) */
static int applyP(const int *P, int v) { return sm[coord(v, 0)][P[0]] ^ sm[coord(v, 1)][P[1]] ^ sm[coord(v, 2)][P[2]]; }
static void out_stab(void) {
    nGO = 0;
    for (int a = 1; a < 64; a++) for (int b = 1; b < 64; b++) for (int c = 1; c < 64; c++) {
        int P[3] = {a, b, c}, img[64], ok = 1; uint64_t seenv = 0;
        for (int v = 0; v < 64; v++) { img[v] = applyP(P, v); if (seenv >> img[v] & 1) { ok = 0; break; } seenv |= 1ull << img[v]; }
        if (!ok) continue;
        for (int x = 1; x < 64 && ok; x <<= 1) for (int y = 1; y < 64; y <<= 1) if (img[bUU[x][y] & 63] != (bUU[x][y] & 63)) { ok = 0; break; }
        if (!ok) continue;
        if (nGO >= 48 * 64) { fprintf(stderr, "GO overflow\n"); exit(1); }
        memcpy(GO[nGO++], P, sizeof P);
    }
}
static uint64_t spanset(const int *c) { uint64_t m = 0; for (int s = 0; s < 16; s++) { int v = 0; for (int i = 0; i < 4; i++) if (s >> i & 1) v ^= c[i]; m |= 1ull << v; } return m; }
static int xrep(const int *c) {   /* is X = span(c) minimal (as a 64-bit set) in its orbit under GO? */
    uint64_t m = spanset(c);
    for (int g = 0; g < nGO; g++) { int d[4]; for (int i = 0; i < 4; i++) d[i] = applyP(GO[g], c[i]); if (spanset(d) < m) return 0; }
    return 1;
}
static void onalrm(int sig) { (void)sig; fprintf(stderr, "[X %d] col %lld %lld %lld %lld leaves %lld | p01 %lld p02 %lld p12 %lld p03 %lld p13 %lld p23 %lld\n", curX, lvl[0], lvl[1], lvl[2], lvl[3], lvl[4], pst[0][1], pst[0][2], pst[1][2], pst[0][3], pst[1][3], pst[2][3]); alarm(20); }
int main(int argc, char **argv) {
    if (getenv("PROG")) { signal(SIGALRM, onalrm); alarm(20); }
    setup_field();
    { int n = 0; for (int u = 1; u < 64; u++) { int k = 0; while (!coord(u, k)) k++; if (coord(u, k) == 1) PP[n++] = u; } }
    UU[0] = 0; for (int t = 1; t < 22; t++) UU[t] = PP[t - 1];
    if (!strcmp(argv[1], "-test")) {   /* random consistency of (A'), (B) with the full kernel test */
        gen_Q(); int n = atoi(argv[2]); long long bad = 0, ca = 0, cb = 0;
        read_class(argv[3]);
        int wok = setup_W();
        for (int t = 0; t < n; t++) {
            for (int j = 0; j < 4; j++) fill_col(j, rnd() & 63, rnd() & 63, rnd() & 63);
            uint32_t q = Qs[rnd() % nQ];
            for (int u = 0; u < 64; u++) for (int z0 = 0; z0 < 16; z0++) {
                if (!u && !z0) continue;
                int full = dir_ok(u | z0 << 6, q);
                if (z0 && !okA(u, z0) && full) bad++;
                if (z0 && okA(u, z0)) ca++;
                if (!z0 && wok) { int pb = okB(u, 4); cb += pb; if (pb != full) bad++; }
                if (dir_ok(sm[2][u] | z0 << 6, q) != full) bad++;    /* omega symmetry */
            }
        }
        printf("test: bad %lld (A' holds %lld, B holds %lld, U x U alive %d)\n", bad, ca, cb, wok);
        return 0;
    }
    read_class(argv[1]);
    for (int x = 0; x < 64; x++) for (int y = 0; y < 64; y++) {   /* sanity: alternating, equivariant */
        if (bUU[x][x]) { fprintf(stderr, "not alternating\n"); return 1; }
        int l = bUU[sm[2][x]][sm[2][y]], r = bUU[x][y];
        if ((l >> 6) != (r >> 6) || (l & 63) != sm[2][r & 63]) { fprintf(stderr, "not equivariant\n"); return 1; }
    }
    if (!setup_W()) { printf("%s: dead on U x U\n", argv[1]); return 0; }
    out_stab(); fprintf(stderr, "output stabilizer of C: %d\n", nGO);
    long long nXall = 0;
    /* E_0 columns: echelon bases of 4-dim X in U' with X + W_{f_0} direct (f_0 = bit 0) */
    long long nX = 0;
    for (int c0 = 1; c0 < 64; c0++) for (int c1 = c0 + 1; c1 < 64; c1++) for (int c2 = c1 + 1; c2 < 64; c2++)
    for (int c3 = c2 + 1; c3 < 64; c3++) {
        int c[4] = {c0, c1, c2, c3};
        /* reduced echelon: pivots = top bits distinct, each pivot bit appears in no other vector */
        int piv = 0, ok = 1;
        for (int i = 0; i < 4 && ok; i++) { int t = 31 - __builtin_clz(c[i]); if (piv >> t & 1) ok = 0; piv |= 1 << t; }
        for (int i = 0; i < 4 && ok; i++) for (int k = 0; k < 4; k++) if (k != i && (c[k] >> (31 - __builtin_clz(c[i])) & 1)) ok = 0;
        if (!ok) continue;
        int v[10], n = 0;
        for (int i = 0; i < Wd[1]; i++) v[n++] = Wb[1][i];
        for (int i = 0; i < 4; i++) v[n++] = c[i];
        if (!indep(v, n)) continue;
        nXall++;
        if (!xrep(c)) continue;
        nX++; curX = nX;
        if (getenv("XONLY") && atoi(getenv("XONLY")) != nX) continue;
        if (getenv("XMOD") && nX % atoi(getenv("XMOD")) != atoi(getenv("XRES"))) continue;
        memcpy(E0, c, sizeof E0);
        make_cand();
        if (getenv("CANDONLY")) { printf("X %lld cand %d %d %d %d\n", nX, ncand[0], ncand[1], ncand[2], ncand[3]); continue; }
        dfs(0);
        if (getenv("PROG")) fprintf(stderr, "X %lld: cand %d %d %d %d levels %lld %lld %lld %lld %lld\n", nX, ncand[0], ncand[1], ncand[2], ncand[3], lvl[0], lvl[1], lvl[2], lvl[3], lvl[4]);
    }
    printf("%s: X %lld of %lld, nodes %lld, tau leaves %lld, Q tests %lld, found %lld\n", argv[1], nX, nXall, nodes, leaves, qtests, found);
    return 0;
}
