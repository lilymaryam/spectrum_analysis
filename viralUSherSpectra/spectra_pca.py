#!/usr/bin/env python3
"""Merge segmented-virus spectra, attach Baltimore type, and plot a PCA.

Steps:
  1. Read the opportunity-normalized spectrum file (one row per sample/segment).
  2. Collapse segments of the same virus into a single spectrum by weighted
     averaging (weight = number of mutations, so segments with more signal
     count more). A weighted average of rows that each sum to 1 also sums to 1.
  3. Attach the Baltimore group from the classification file.
  4. PCA on the 12 mutation-type features; scatter PC1 vs PC2 coloured by
     Baltimore group.

Segment groupings live in SEGMENT_PREFIXES below (plus the Influenza rule) and
are easy to edit. An optional --group-map TSV (sample<TAB>group) overrides them.
"""
import argparse
import re
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

BASES = "ACGT"

# Prefix -> shared virus group. Every sample whose name starts with the prefix
# is treated as one segment of that virus. These are single-isolate segment sets
# (or, for the numbered reo/orbi/rota families, the numbers are segments 1..N
# with some dropped by dedup, not serotypes).
# Ambiguous / worth a sanity check: Bluetongue_virus_* and Orbivirus_alphaequi_*
# (numbers = segments here, but these genera are also numbered by serotype);
# SFTSV appears under two names below and is left as two groups.
SEGMENT_PREFIXES = [
    ("Avian_orthoreovirus_AVS_B_", "Avian_orthoreovirus_AVS_B"),
    ("Mammalian_orthoreovirus_3_", "Mammalian_orthoreovirus_3"),
    ("Piscine_orthoreovirus_", "Piscine_orthoreovirus"),
    ("Rotavirus_A_", "Rotavirus_A"),
    ("Rotavirus_C_", "Rotavirus_C"),
    ("Rice_black_streaked_dwarf_virus_", "Rice_black_streaked_dwarf_virus"),
    ("Bluetongue_virus_", "Bluetongue_virus"),
    ("Orbivirus_alphaequi_", "Orbivirus_alphaequi"),
    ("Epizootic_hemorrhagic_disease_virus_serotype_1_strain_New_Jersey_",
     "Epizootic_hemorrhagic_disease_virus_serotype_1_strain_New_Jersey"),
    ("Cardamom_bushy_dwarf_virus_", "Cardamom_bushy_dwarf_virus"),
    ("Oropouche_virus_", "Oropouche_virus"),
    ("Rift_Valley_fever_virus_", "Rift_Valley_fever_virus"),
    ("Orthohantavirus_hantanense_", "Orthohantavirus_hantanense"),
    ("Orthohantavirus_seoulense_", "Orthohantavirus_seoulense"),
    ("Orthonairovirus_haemorrhagiae_", "Orthonairovirus_haemorrhagiae"),
    ("Orthotospovirus_tomatomaculae_", "Orthotospovirus_tomatomaculae"),
    ("Mammarenavirus_lassaense_", "Mammarenavirus_lassaense"),
    ("SFTS_virus_HB29_", "SFTS_virus_HB29"),
    ("Severe_fever_with_thrombocytopenia_syndrome_virus_",
     "Severe_fever_with_thrombocytopenia_syndrome_virus"),
    ("Isavirus_salaris_", "Isavirus_salaris"),
    ("Mogiana_tick_virus_", "Mogiana_tick_virus"),
    ("Infectious_pancreatic_necrosis_virus_", "Infectious_pancreatic_necrosis_virus"),
    ("Infectious_bursal_disease_virus_1977_from_Schobries_et_al_1977_",
     "Infectious_bursal_disease_virus_1977_from_Schobries_et_al_1977"),
]

# Baltimore group ordering + a colour for each.
BALT_ORDER = ["I", "II", "III", "IV", "V", "VI", "VII", "NA"]
BALT_COLORS = {
    "I": "#1f77b4", "II": "#ff7f0e", "III": "#2ca02c", "IV": "#d62728",
    "V": "#9467bd", "VI": "#8c564b", "VII": "#e377c2", "NA": "#7f7f7f",
}


def virus_group(name, override):
    if name in override:
        return override[name]
    if name.startswith("Influenza_"):          # segments are a trailing _<int>; keep strain
        return re.sub(r"_\d+$", "", name)
    for prefix, group in SEGMENT_PREFIXES:
        if name.startswith(prefix):
            return group
    return name


def convex_hull(pts):
    """Andrew's monotone chain; returns hull vertices. Needs >=3 points."""
    pts = sorted(map(tuple, pts))
    if len(pts) < 3:
        return None
    def cross(o, a, b):
        return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])
    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    hull = lower[:-1] + upper[:-1]
    return np.array(hull) if len(hull) >= 3 else None


def find_mutation_columns(cols):
    return [c for c in cols
            if len(str(c).strip()) == 2
            and str(c).strip()[0] in BASES and str(c).strip()[1] in BASES
            and str(c).strip()[0] != str(c).strip()[1]]


def write_html(path, merged, mut_cols, scores, evr, present, colors, labels):
    """Self-contained interactive HTML: PCA + 4x4 strips, hover tooltips. No deps."""
    import json
    meta_cols = [c for c in ("Baltimore_group", "genome_type", "n_segments",
                             "total_mutations", "number_tips", "segments")
                 if c in merged.columns]
    recs = []
    for i in range(len(merged)):
        r = {"name": str(merged["virus_group"].iat[i]),
             "pc1": float(scores[i, 0]), "pc2": float(scores[i, 1]),
             "g": str(merged["Baltimore_group"].iat[i])}
        for c in meta_cols:
            v = merged[c].iat[i]
            r[c] = None if pd.isna(v) else (float(v) if isinstance(v, (int, float, np.floating, np.integer)) else str(v))
        for m in mut_cols:
            r[m] = float(merged[m].iat[i])
        recs.append(r)

    payload = {
        "recs": recs, "muts": list(mut_cols), "groups": list(present),
        "colors": {g: colors.get(g, "#333") for g in present},
        "labels": {g: labels.get(g, g) for g in present},
        "evr": [float(x) for x in evr[:3]], "meta_cols": meta_cols,
    }

    html = """<!doctype html><meta charset="utf-8">
<title>Mutation spectra</title>
<style>
 body{font:13px/1.4 system-ui,sans-serif;margin:20px;color:#222}
 h2{font-size:15px;font-weight:600;margin:18px 0 6px}
 #tip{position:fixed;pointer-events:none;background:#fff;border:1px solid #bbb;
      border-radius:5px;padding:7px 9px;font-size:12px;box-shadow:0 2px 8px #0003;
      opacity:0;transition:opacity .08s;max-width:340px;z-index:9}
 #tip b{font-size:12.5px}
 #tip div{color:#555;margin-top:2px}
 circle{cursor:pointer}
 circle:hover{stroke:#000;stroke-width:1.4}
 .legend span{margin-right:14px;white-space:nowrap}
 .sw{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:4px;
     vertical-align:middle}
 .axis{stroke:#999;stroke-width:1}
 .tick{font-size:10px;fill:#666}
 .ttl{font-size:12px;fill:#333;font-weight:600}
 #strips{display:grid;grid-template-columns:repeat(4,1fr);gap:6px;max-width:1000px}
</style>
<h2>PCA of opportunity-normalized mutation spectra</h2>
<div class="legend" id="leg"></div>
<div id="pca"></div>
<h2>Per-virus mutation-type rates by Baltimore group</h2>
<div id="strips"></div>
<div id="tip"></div>
<script>
const D = __DATA__;
const tip = document.getElementById('tip');
const NS = 'http://www.w3.org/2000/svg';
const fmt = v => (v===null||v===undefined) ? '-' :
  (typeof v==='number' ? (Math.abs(v)<1 ? v.toFixed(5) : v.toLocaleString()) : v);

function show(e, r){
  let h = '<b>'+r.name+'</b>';
  for(const c of D.meta_cols){ if(c!=='segments') h += '<div>'+c.replace(/_/g,' ')+': '+fmt(r[c])+'</div>'; }
  if(r.segments) h += '<div style="margin-top:4px;color:#888">'+r.segments+'</div>';
  tip.innerHTML = h; tip.style.opacity = 1;
  const pad=14; let x=e.clientX+pad, y=e.clientY+pad;
  if(x+tip.offsetWidth>innerWidth) x=e.clientX-tip.offsetWidth-pad;
  if(y+tip.offsetHeight>innerHeight) y=e.clientY-tip.offsetHeight-pad;
  tip.style.left=x+'px'; tip.style.top=y+'px';
}
const hide = () => tip.style.opacity = 0;

// legend
document.getElementById('leg').innerHTML = D.groups.map(g =>
  '<span><i class="sw" style="background:'+D.colors[g]+'"></i>'+D.labels[g]+'</span>').join('');

function el(t,a){const n=document.createElementNS(NS,t);for(const k in a)n.setAttribute(k,a[k]);return n;}

// ---- PCA ----
(function(){
  const W=680,H=680,M=52;
  const xs=D.recs.map(r=>r.pc1), ys=D.recs.map(r=>r.pc2);
  const cx=(Math.min(...xs)+Math.max(...xs))/2, cy=(Math.min(...ys)+Math.max(...ys))/2;
  const half=0.54*Math.max(Math.max(...xs)-Math.min(...xs), Math.max(...ys)-Math.min(...ys));
  const X=v=>M+(v-(cx-half))/(2*half)*(W-2*M);
  const Y=v=>H-M-(v-(cy-half))/(2*half)*(H-2*M);
  const svg=el('svg',{width:W,height:H});
  svg.appendChild(el('line',{x1:M,y1:H-M,x2:W-M,y2:H-M,class:'axis'}));
  svg.appendChild(el('line',{x1:M,y1:M,x2:M,y2:H-M,class:'axis'}));
  const xl=el('text',{x:W/2,y:H-14,'text-anchor':'middle',class:'tick'});
  xl.textContent='PC1 ('+(D.evr[0]*100).toFixed(1)+'%)'; svg.appendChild(xl);
  const yl=el('text',{x:16,y:H/2,'text-anchor':'middle',class:'tick',
    transform:'rotate(-90 16 '+H/2+')'});
  yl.textContent='PC2 ('+(D.evr[1]*100).toFixed(1)+'%)'; svg.appendChild(yl);
  D.recs.forEach(r=>{
    const c=el('circle',{cx:X(r.pc1),cy:Y(r.pc2),r:5,fill:D.colors[r.g]||'#333',
      'fill-opacity':.85,stroke:'#fff','stroke-width':.6});
    c.addEventListener('mousemove',e=>show(e,r)); c.addEventListener('mouseleave',hide);
    svg.appendChild(c);
  });
  document.getElementById('pca').appendChild(svg);
})();

// ---- 4x4 strips ----
(function(){
  const B=['A','C','G','T'], host=document.getElementById('strips');
  const W=240,H=190,M={l:44,r:8,t:20,b:26};
  let seed=1; const rnd=()=>{seed=(seed*1103515245+12345)&0x7fffffff; return seed/0x7fffffff;};
  for(const fb of B) for(const tb of B){
    const cell=document.createElement('div');
    if(fb===tb){ host.appendChild(cell); continue; }
    const mt=fb+tb;
    const vals=D.recs.map(r=>r[mt]);
    const lo=Math.min(...vals), hi=Math.max(...vals), pad=(hi-lo)*0.08||1e-6;
    const Y=v=>H-M.b-(v-(lo-pad))/((hi+pad)-(lo-pad))*(H-M.t-M.b);
    const X=i=>M.l+(i+0.5)/D.groups.length*(W-M.l-M.r);
    const svg=el('svg',{width:W,height:H});
    const t=el('text',{x:W/2,y:14,'text-anchor':'middle',class:'ttl'});
    t.textContent=fb+'>'+tb; svg.appendChild(t);
    svg.appendChild(el('line',{x1:M.l,y1:H-M.b,x2:W-M.r,y2:H-M.b,class:'axis'}));
    svg.appendChild(el('line',{x1:M.l,y1:M.t,x2:M.l,y2:H-M.b,class:'axis'}));
    [lo,hi].forEach(v=>{const q=el('text',{x:M.l-5,y:Y(v)+3,'text-anchor':'end',class:'tick'});
      q.textContent=v.toPrecision(2); svg.appendChild(q);});
    const bw=(W-M.l-M.r)/D.groups.length;
    D.groups.forEach((g,i)=>{
      const lab=el('text',{x:X(i),y:H-M.b+12,'text-anchor':'middle',class:'tick'});
      lab.textContent=g; svg.appendChild(lab);
      const gs=D.recs.filter(r=>r.g===g);
      const ys=gs.map(r=>r[mt]).sort((a,b)=>a-b);
      if(ys.length){
        const med=ys[Math.floor(ys.length/2)];
        svg.appendChild(el('line',{x1:X(i)-bw*0.3,x2:X(i)+bw*0.3,y1:Y(med),y2:Y(med),
          stroke:'#222','stroke-width':1.3}));
      }
      gs.forEach(r=>{
        const c=el('circle',{cx:X(i)+(rnd()-0.5)*bw*0.6,cy:Y(r[mt]),r:2.8,
          fill:D.colors[g],'fill-opacity':.5});
        c.addEventListener('mousemove',e=>show(e,r)); c.addEventListener('mouseleave',hide);
        svg.appendChild(c);
      });
    });
    cell.appendChild(svg); host.appendChild(cell);
  }
})();
</script>"""
    with open(path, "w") as f:
        f.write(html.replace("__DATA__", json.dumps(payload)))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-s", "--spectra", required=True,
                    help="Opportunity-normalized spectrum TSV.")
    ap.add_argument("-b", "--baltimore", required=True,
                    help="Baltimore classification TSV (Sample, Baltimore_group[, genome_type]).")
    ap.add_argument("-o", "--output", default="spectra_pca.png",
                    help="Output figure path (default: spectra_pca.png).")
    ap.add_argument("--spectra-key", default="Sample")
    ap.add_argument("--balt-key", default="Sample")
    ap.add_argument("--balt-group-col", default="Baltimore_group")
    ap.add_argument("--weight-col", default="Total_Mutations",
                    help="Column used to weight segments when averaging (default: Total_Mutations).")
    ap.add_argument("--group-map", default=None,
                    help="Optional TSV of sample<TAB>group overriding built-in segment grouping.")
    ap.add_argument("--no-ellipse", action="store_true",
                    help="Do not draw a shaded hull polygon around each Baltimore group.")
    ap.add_argument("--standardize", action="store_true",
                    help="Z-score each feature before PCA (default: mean-center only).")
    ap.add_argument("--html", default=None,
                    help="Also write a self-contained interactive HTML with hover tooltips "
                         "(PCA + strips) showing virus metadata.")
    ap.add_argument("--dump-merged", default=None,
                    help="Optional path to write the merged per-virus spectra TSV.")
    args = ap.parse_args()

    spec = pd.read_csv(args.spectra, sep="\t")
    balt = pd.read_csv(args.baltimore, sep="\t")
    mut_cols = find_mutation_columns(spec.columns)
    if len(mut_cols) != 12:
        sys.stderr.write(f"Warning: expected 12 mutation columns, found {len(mut_cols)}\n")

    override = {}
    if args.group_map:
        gm = pd.read_csv(args.group_map, sep="\t", header=None)
        override = dict(zip(gm[0].astype(str), gm[1].astype(str)))

    spec["virus_group"] = spec[args.spectra_key].astype(str).map(lambda n: virus_group(n, override))

    # weight per segment
    if args.weight_col in spec.columns:
        w = pd.to_numeric(spec[args.weight_col], errors="coerce").fillna(0.0).to_numpy()
    else:
        sys.stderr.write(f"Weight column '{args.weight_col}' not found; using equal weights.\n")
        w = np.ones(len(spec))
    spec["_w"] = np.where(w > 0, w, 0.0)

    # weighted mean of the 12 features within each virus group
    def wavg(g):
        ww = g["_w"].to_numpy()
        tot = ww.sum()
        if tot <= 0:
            ww = np.ones(len(g)); tot = ww.sum()
        vals = (g[mut_cols].to_numpy() * ww[:, None]).sum(axis=0) / tot
        return pd.Series(vals, index=mut_cols)

    merged = spec.groupby("virus_group", sort=True).apply(wavg, include_groups=False)
    merged = merged.reset_index()
    # renormalize defensively so each row sums to 1
    rs = merged[mut_cols].sum(axis=1).replace(0, np.nan)
    merged[mut_cols] = merged[mut_cols].div(rs, axis=0)

    n_seg = spec.groupby("virus_group").size()
    multi = (n_seg > 1).sum()
    sys.stderr.write(f"{len(spec)} segments -> {len(merged)} viruses "
                     f"({multi} were multi-segment merges).\n")

    # attach Baltimore group via the original segment names
    keep = [args.balt_key, args.balt_group_col] + (
        ["genome_type"] if "genome_type" in balt.columns else [])
    seg2grp = spec[[args.spectra_key, "virus_group"]].rename(columns={args.spectra_key: args.balt_key})
    bmap = balt[keep].merge(seg2grp, on=args.balt_key, how="inner")
    # one Baltimore group per virus (warn on any disagreement)
    grp_balt = (bmap.groupby("virus_group")[args.balt_group_col]
                .agg(lambda s: s.mode().iat[0] if not s.mode().empty else "NA"))
    disagree = bmap.groupby("virus_group")[args.balt_group_col].nunique()
    for vg in disagree[disagree > 1].index:
        sys.stderr.write(f"Warning: mixed Baltimore groups within '{vg}'\n")
    gt = None
    if "genome_type" in bmap.columns:
        gt = bmap.groupby("virus_group")["genome_type"].agg(
            lambda s: s.mode().iat[0] if not s.mode().empty else "")

    merged = merged.merge(grp_balt.rename("Baltimore_group"), on="virus_group", how="left")
    merged["Baltimore_group"] = merged["Baltimore_group"].fillna("NA")
    if gt is not None:
        merged = merged.merge(gt.rename("genome_type"), on="virus_group", how="left")

    # per-virus metadata for tooltips
    meta = pd.DataFrame({"n_segments": spec.groupby("virus_group").size()})
    for col, out in [(args.weight_col, "total_mutations"), ("Number_Tips", "number_tips")]:
        if col in spec.columns:
            meta[out] = pd.to_numeric(spec[col], errors="coerce").groupby(
                spec["virus_group"]).sum()
    meta["segments"] = spec.groupby("virus_group")[args.spectra_key].agg(
        lambda s: ", ".join(map(str, s)))
    merged = merged.merge(meta.reset_index(), on="virus_group", how="left")

    if args.dump_merged:
        merged.to_csv(args.dump_merged, sep="\t", index=False, float_format="%.10g")

    # ---- PCA (numpy SVD, no sklearn dependency) ----
    X = merged[mut_cols].to_numpy(dtype=float)
    X = X - X.mean(axis=0)
    if args.standardize:
        sd = X.std(axis=0); sd[sd == 0] = 1.0
        X = X / sd
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    scores = U * S                       # projected coordinates
    ev = (S ** 2) / (len(X) - 1)
    evr = ev / ev.sum()

    # ---- plot ----
    k = len(evr)
    fig = plt.figure(figsize=(15, 22))
    gs = fig.add_gridspec(6, 4, height_ratios=[1, 1, 1, 1, 3, 3],
                          hspace=0.30, wspace=0.35)
    # 4x4 strip grid: rows = from-base, cols = to-base; diagonal blank
    strip_grid = {(r, c): fig.add_subplot(gs[r, c]) for r in range(4) for c in range(4)}
    ax = fig.add_subplot(gs[4:6, 0:3])    # PCA (left, square)
    axL = fig.add_subplot(gs[4, 3])       # loadings (top-right)
    axS = fig.add_subplot(gs[5, 3])       # scree (bottom-right)

    present = [g for g in BALT_ORDER if (merged["Baltimore_group"] == g).any()]
    labels = {}
    if "genome_type" in merged.columns:
        for g in present:
            gtv = merged.loc[merged["Baltimore_group"] == g, "genome_type"].dropna()
            labels[g] = f"{g} ({gtv.iloc[0]})" if len(gtv) else g
    else:
        labels = {g: g for g in present}

    # per-virus jitter strips, panel (from-base, to-base); diagonal blank
    _rng = np.random.default_rng(0)
    for r, fb in enumerate(BASES):
        for c, tb in enumerate(BASES):
            sax = strip_grid[(r, c)]
            if r == c:
                sax.axis("off")
                continue
            mt = fb + tb
            for i, g in enumerate(present):
                y = merged.loc[merged["Baltimore_group"] == g, mt].to_numpy(float)
                sax.scatter(i + _rng.uniform(-0.3, 0.3, len(y)), y, s=10, alpha=0.45,
                            color=BALT_COLORS.get(g, "#333"), edgecolor="none", zorder=2)
                sax.hlines(np.median(y), i - 0.35, i + 0.35, color="0.15", lw=1.2, zorder=3)
            sax.set_title(f"{fb}>{tb}", fontsize=10)
            sax.set_xticks(range(len(present)))
            sax.set_xticklabels(present, fontsize=6.5)
            sax.tick_params(axis="y", labelsize=7)
            sax.margins(x=0.05)
    # one centered label under the strip block
    pr = [strip_grid[(3, c)].get_position() for c in range(4)]
    xc = (pr[0].x0 + pr[-1].x1) / 2
    yb = min(p.y0 for p in pr)
    fig.text(xc, yb - 0.022, "Baltimore classification", ha="center", fontsize=11)

    for g in present:
        m = merged["Baltimore_group"] == g
        pts = scores[m.to_numpy(), :2]
        if not args.no_ellipse:
            hull = convex_hull(pts)
            if hull is not None:
                ax.fill(hull[:, 0], hull[:, 1], color=BALT_COLORS.get(g, "#333333"),
                        alpha=0.06, zorder=1, lw=0)
                ax.plot(np.append(hull[:, 0], hull[0, 0]),
                        np.append(hull[:, 1], hull[0, 1]),
                        color=BALT_COLORS.get(g, "#333333"), alpha=0.5, lw=1.0, zorder=1)
        ax.scatter(scores[m.to_numpy(), 0], scores[m.to_numpy(), 1],
                   s=42, alpha=0.85, edgecolor="white", linewidth=0.4,
                   color=BALT_COLORS.get(g, "#333333"), label=labels[g])
    ax.set_xlabel(f"PC1 ({evr[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({evr[1]*100:.1f}%)")
    _x0, _x1 = scores[:, 0].min(), scores[:, 0].max()
    _y0, _y1 = scores[:, 1].min(), scores[:, 1].max()
    _cx, _cy = (_x0 + _x1) / 2, (_y0 + _y1) / 2
    _half = 0.5 * max(_x1 - _x0, _y1 - _y0) * 1.08
    ax.set_xlim(_cx - _half, _cx + _half)
    ax.set_ylim(_cy - _half, _cy + _half)
    ax.set_aspect("equal", adjustable="box")
    ax.legend(title="Baltimore group", frameon=False, fontsize=9,
              loc="best", markerscale=1.1)

    # loadings: grouped bars, one group per mutation type, bars = PC1/PC2/PC3
    npc = min(3, Vt.shape[0])
    Lmat = Vt[:npc, :]                     # npc x 12
    xf = np.arange(len(mut_cols))
    width = 0.8 / npc
    pc_colors = ["#4c72b0", "#dd8452", "#55a868"]
    for j in range(npc):
        axL.bar(xf + (j - (npc - 1) / 2) * width, Lmat[j], width,
                label=f"PC{j+1}", color=pc_colors[j])
    axL.axhline(0, color="0.6", lw=0.8)
    axL.set_xticks(xf)
    axL.set_xticklabels([f"{m[0]}>{m[1]}" for m in mut_cols], rotation=90, fontsize=8)
    axL.set_ylabel("loading")
    axL.legend(frameon=False, fontsize=8, ncol=npc)

    # scree
    x = np.arange(1, k + 1)
    axS.bar(x, evr * 100, color="#4c72b0", alpha=0.9)
    axS.plot(x, np.cumsum(evr) * 100, "-o", color="#c44e52", ms=4, lw=1.2,
             label="cumulative")
    axS.set_xlabel("Principal component")
    axS.set_ylabel("Variance explained (%)")
    axS.set_xticks(x)
    axS.legend(frameon=False, fontsize=9)

    # align right-column panels to the PCA's drawn (square) vertical extent,
    # and shift them right so their y-labels clear the PCA
    fig.canvas.draw()
    p = ax.get_position()
    strip_right = max(strip_grid[(0, c)].get_position().x1 for c in range(4))
    rx0, rx1 = p.x1 + 0.05, strip_right
    gap = 0.03
    mid = (p.y0 + p.y1) / 2
    axL.set_position([rx0, mid + gap / 2, rx1 - rx0, p.y1 - (mid + gap / 2)])
    axS.set_position([rx0, p.y0, rx1 - rx0, (mid - gap / 2) - p.y0])

    fig.savefig(args.output, dpi=200)

    if args.html:
        write_html(args.html, merged, mut_cols, scores, evr, present,
                   BALT_COLORS, labels)
        sys.stderr.write(f"Wrote {args.html}\n")

    sys.stderr.write(f"PC1-3 explained variance: "
                     f"{evr[0]*100:.1f}%, {evr[1]*100:.1f}%, {evr[2]*100:.1f}%\n")
    sys.stderr.write(f"Wrote {args.output}\n")


if __name__ == "__main__":
    main()