#!/usr/bin/env python3
"""Variance in mutation spectra explained by Baltimore type (+ optional covariates).

Sequential (Type I) distance-based PERMANOVA: partitions multivariate variance
across terms in the order given, so each term's R2 is the fraction it explains
*beyond* the terms before it. Distance defaults to Aitchison (Euclidean on CLR).
Also runs PERMDISP on the primary group to flag dispersion (spread) differences
that could otherwise inflate a location result.

Inputs: a merged per-virus spectra TSV (12 mutation cols + group col), and an
optional annotation TSV joined on the group column to add covariates.
"""
import argparse, sys
import numpy as np, pandas as pd

BASES = "ACGT"


def find_mutation_columns(cols):
    return [c for c in cols
            if len(str(c).strip()) == 2
            and str(c).strip()[0] in BASES and str(c).strip()[1] in BASES
            and str(c).strip()[0] != str(c).strip()[1]]


def clr(X, pseudo):
    X = X / X.sum(axis=1, keepdims=True)
    if pseudo is None:
        pos = X[X > 0]
        pseudo = pos.min() * 0.5 if pos.size else 1e-6
    X = X + pseudo
    X = X / X.sum(axis=1, keepdims=True)
    L = np.log(X)
    return L - L.mean(axis=1, keepdims=True)


def dummies(series):
    cats = pd.Categorical(series)
    M = np.zeros((len(series), len(cats.categories) - 1))
    for j, c in enumerate(cats.categories[1:]):
        M[:, j] = (cats == c).astype(float)
    return M


def gower(D2):
    N = D2.shape[0]
    J = np.eye(N) - np.ones((N, N)) / N
    return -0.5 * J @ D2 @ J


def hat_trace(cols, G):
    N = G.shape[0]
    X = np.hstack([np.ones((N, 1))] + cols) if cols else np.ones((N, 1))
    H = X @ np.linalg.pinv(X.T @ X) @ X.T
    return np.einsum("ij,ji->", H, G), X.shape[1]


def seq_permanova(G, blocks, perms, rng):
    Gtot = np.trace(G)
    N = G.shape[0]

    def stats(Gm):
        cum, ss_cum, k_cum = [], [], []
        prev_ss, prev_k = 0.0, 1
        for _, arr in blocks:
            cum.append(arr)
            ss, k = hat_trace(cum, Gm)
            ss_cum.append(ss - prev_ss); k_cum.append(k - prev_k)
            prev_ss, prev_k = ss, k
        ss_res = np.trace(Gm) - prev_ss
        df_res = N - prev_k
        Fs = [(ss / df) / (ss_res / df_res) if df > 0 and ss_res > 0 else np.nan
              for ss, df in zip(ss_cum, k_cum)]
        return ss_cum, k_cum, Fs, ss_res, df_res

    ss_cum, df_cum, F_obs, ss_res, df_res = stats(G)
    ge = [1] * len(blocks)
    for _ in range(perms):
        p = rng.permutation(N)
        _, _, Fp, _, _ = stats(G[np.ix_(p, p)])
        for i in range(len(blocks)):
            if np.isfinite(Fp[i]) and Fp[i] >= F_obs[i]:
                ge[i] += 1
    out = []
    for i, (name, _) in enumerate(blocks):
        out.append((name, ss_cum[i] / Gtot, df_cum[i], F_obs[i], ge[i] / (perms + 1)))
    return out, ss_res / Gtot, df_res


def permdisp(Xclr, labels, perms, rng):
    groups = np.unique(labels)
    disp = np.zeros(len(labels))
    for g in groups:
        ix = np.where(labels == g)[0]
        disp[ix] = np.linalg.norm(Xclr[ix] - Xclr[ix].mean(axis=0), axis=1)

    def F(lab):
        grand = disp.mean(); ssb = ssw = 0.0
        for g in groups:
            d = disp[lab == g]
            ssb += len(d) * (d.mean() - grand) ** 2
            ssw += ((d - d.mean()) ** 2).sum()
        a, N = len(groups), len(lab)
        return (ssb / (a - 1)) / (ssw / (N - a)) if ssw > 0 else np.inf

    Fo = F(labels); ge = 1
    for _ in range(perms):
        if F(rng.permutation(labels)) >= Fo:
            ge += 1
    return {g: disp[labels == g].mean() for g in groups}, Fo, ge / (perms + 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True,
                    help="Merged per-virus spectra TSV (12 mutation cols + group col).")
    ap.add_argument("-a", "--annotations", default=None,
                    help="Annotation TSV joined on --group-col to add covariates.")
    ap.add_argument("--group-col", default="Baltimore_group")
    ap.add_argument("--terms", default=None,
                    help="Comma-separated predictors in sequential order "
                         "(default: group-col, then all annotation columns).")
    ap.add_argument("--annot-key", default="virus_group",
                    help="Join key in the annotation file (default: virus_group).")
    ap.add_argument("--distance", choices=["aitchison", "euclidean"], default="aitchison")
    ap.add_argument("--permutations", type=int, default=9999)
    ap.add_argument("--min-group-size", type=int, default=2,
                    help="Drop rows whose primary group is smaller than this.")
    ap.add_argument("--pseudo", type=float, default=None)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    df = pd.read_csv(args.input, sep="\t", keep_default_na=False, na_values=[""])
    if args.group_col not in df.columns:
        sys.exit(f"Group column '{args.group_col}' not found.")
    mut = find_mutation_columns(df.columns)
    if len(mut) != 12:
        sys.stderr.write(f"Warning: found {len(mut)} mutation columns\n")

    if args.annotations:
        ann = pd.read_csv(args.annotations, sep="\t", keep_default_na=False, na_values=[""])
        left = args.group_col if args.annot_key not in df.columns else args.annot_key
        df = df.merge(ann, left_on=left, right_on=args.annot_key, how="left",
                      suffixes=("", "_ann"))

    terms = ([t.strip() for t in args.terms.split(",")] if args.terms else
             [args.group_col] + ([c for c in ann.columns if c != args.annot_key]
                                 if args.annotations else []))
    for t in terms:
        if t not in df.columns:
            sys.exit(f"Term '{t}' not found. Columns: {list(df.columns)}")

    df = df.dropna(subset=mut + terms).copy()
    sizes = df[args.group_col].value_counts()
    keep = sizes[sizes >= args.min_group_size].index
    dropped = sizes[sizes < args.min_group_size]
    if len(dropped):
        sys.stderr.write("Dropped small primary groups: "
                         + ", ".join(f"{g}(n={n})" for g, n in dropped.items()) + "\n")
    df = df[df[args.group_col].isin(keep)].reset_index(drop=True)

    Xr = df[mut].to_numpy(float)
    Xc = clr(Xr, args.pseudo) if args.distance == "aitchison" else Xr / Xr.sum(1, keepdims=True)
    G_ = Xc @ Xc.T
    sq = np.diag(G_)
    D2 = np.maximum(sq[:, None] + sq[None, :] - 2 * G_, 0.0)
    G = gower(D2)

    blocks = [(t, dummies(df[t])) for t in terms]
    rng = np.random.default_rng(args.seed)
    res, r2_res, df_res = seq_permanova(G, blocks, args.permutations, rng)
    means, Fd, pd_ = permdisp(Xc, df[args.group_col].to_numpy(),
                              args.permutations, np.random.default_rng(args.seed + 1))

    print(f"N = {len(df)} viruses, distance = {args.distance}, "
          f"{args.permutations} permutations")
    print(f"sequential (Type I) PERMANOVA: spectrum ~ {' + '.join(terms)}\n")
    print(f"{'term':<18}{'df':>4}{'R2':>10}{'pseudo-F':>11}{'p':>10}")
    for name, r2, dfk, F, p in res:
        print(f"{name:<18}{dfk:>4}{r2:>10.4f}{F:>11.3f}{p:>10.4g}")
    print(f"{'Residual':<18}{df_res:>4}{r2_res:>10.4f}")
    print(f"\nPERMDISP on {args.group_col}: F = {Fd:.3f}, p = {pd_:.4g}   "
          f"({'unequal spread; interpret with care' if pd_ < 0.05 else 'homogeneous dispersion'})")


if __name__ == "__main__":
    main()