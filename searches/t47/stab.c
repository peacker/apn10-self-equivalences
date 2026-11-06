/* stab.c: full stabilizer of a type-47 U x U class (cls file) under input P in GL_3(F_4) (F_4-linear on U)
   with output (P' F_4-linear on U', N in GL(Z')). Prints |H| and, for each projective point f of U, the size of
   H_f = {P : P f in F_4^* f}.   build: cc -O3 -o stab stab.c   usage: stab cls_file */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define NU 6
#define SU 64
static int M4[4][4], sm[4][SU], bUU[SU][SU];
static int coord(int u, int k) { return (u >> (2 * k)) & 3; }
static int applyP(const int *P, int v) { int r = 0; for (int k = 0; k < 3; k++) r ^= sm[coord(v, k)][P[k]]; return r; }
/* try to extend the partial linear map src->dst (10-bit vectors); returns 0 on inconsistency */
static int addpair(int *bs, int *bd, int *r, int s, int d) {
    for (int j = 0; j < *r; j++) if ((s ^ bs[j]) < s) { s ^= bs[j]; d ^= bd[j]; }
    if (!s) return d == 0;
    bs[*r] = s; bd[*r] = d; (*r)++;
    for (int j = *r - 1; j > 0 && bs[j] > bs[j - 1]; j--) { int t = bs[j]; bs[j] = bs[j - 1]; bs[j - 1] = t; t = bd[j]; bd[j] = bd[j - 1]; bd[j - 1] = t; }
    return 1;
}
int main(int argc, char **argv) {
    for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) { int r = 0; for (int i = 0; i < 2; i++) if (b >> i & 1) r ^= a << i; if (r & 4) r ^= 7; M4[a][b] = r; }
    for (int c = 0; c < 4; c++) for (int v = 0; v < SU; v++) { int r = 0; for (int k = 0; k < 3; k++) r |= M4[c][coord(v, k)] << (2 * k); sm[c][v] = r; }
    int b[NU][NU]; memset(b, 0, sizeof b);
    FILE *f = fopen(argv[1], "r"); int p, v;
    while (fscanf(f, "%d %d", &p, &v) == 2) { int j = 1; while (j * (j + 1) / 2 <= p) j++; int i = p - j * (j - 1) / 2; b[i][j] = b[j][i] = v; }
    fclose(f);
    for (int x = 0; x < SU; x++) for (int y = 0; y < SU; y++) { int r = 0; for (int i = 0; i < NU; i++) if (x >> i & 1) for (int j = 0; j < NU; j++) if (y >> j & 1) r ^= b[i][j]; bUU[x][y] = r; }
    int PP[21], n = 0; for (int u = 1; u < SU; u++) { int k = 0; while (!coord(u, k)) k++; if (coord(u, k) == 1) PP[n++] = u; }
    long long H = 0, Hf[21] = {0};
    int P[3];
    for (P[0] = 1; P[0] < SU; P[0]++) for (P[1] = 1; P[1] < SU; P[1]++) for (P[2] = 1; P[2] < SU; P[2]++) {
        static int hit[SU]; memset(hit, 0, sizeof hit); int ok = 1;
        for (int u = 0; u < SU && ok; u++) { int w = applyP(P, u); if (hit[w]) ok = 0; hit[w] = 1; }
        if (!ok) continue;
        /* output map T on 10 bits with T beta(x,y) = beta(Px,Py); U' part must be F_4-linear: add omega-pairs */
        int bs[10], bd[10], r = 0;
        for (int x = 1; x < SU && ok; x <<= 1) for (int y = 1; y < SU && ok; y <<= 1) {
            int s = bUU[x][y], d = bUU[applyP(P, x)][applyP(P, y)];
            if ((s >> NU) == 0 && (d >> NU) != 0) { } /* allowed only if consistent; checked below via block test */
            ok = addpair(bs, bd, &r, s, d);
            int ws = sm[2][s & (SU - 1)] | (s & ~(SU - 1)), wd = sm[2][d & (SU - 1)] | (d & ~(SU - 1));
            if (ok) ok = addpair(bs, bd, &r, ws, wd);   /* omega acts on U', trivially on Z' (h(wx,wy) = h(x,y)) */
        }
        if (!ok) continue;
        /* block test: U' -> U' and Z' -> Z' */
        for (int j = 0; j < r && ok; j++) { /* reduce unit vectors */ }
        for (int e = 0; e < 10 && ok; e++) {
            int s = 1 << e, d = 0;
            for (int j = 0; j < r; j++) if ((s ^ bs[j]) < s) { s ^= bs[j]; d ^= bd[j]; }
            if (s) continue; /* outside span: free */
            if (e < NU && (d >> NU)) ok = 0;
            if (e >= NU && (d & (SU - 1))) ok = 0;
        }
        if (!ok) continue;
        H++;
        for (int t = 0; t < 21; t++) { int w = applyP(P, PP[t]); if (w == PP[t] || w == sm[2][PP[t]] || w == sm[3][PP[t]]) Hf[t]++; }
    }
    printf("%s |H| %lld  H_f:", argv[1], H);
    for (int t = 0; t < 21; t++) printf(" %lld", Hf[t]);
    printf("\n");
    return 0;
}
