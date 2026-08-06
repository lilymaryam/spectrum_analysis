#!/usr/bin/env python3
"""
Deduplicate a set of UShER/BTE phylogenetic trees by sample overlap.
For any pair of trees sharing samples, retain only the larger tree.
Samples present in all trees are assumed to be reference/root artifacts and ignored.

Usage:
    python dedup_trees.py -i *.pb *.pb.gz
    python dedup_trees.py -i *.pb.gz --min-tree-freq 0.9
"""

import bte
import os
import gzip
import shutil
import tempfile
import argparse
from collections import Counter
from itertools import combinations

from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn
from rich import box
from rich.text import Text

console = Console()


def strip_tree_ext(filename: str) -> str:
    for ext in (".pb.gz", ".pb", ".nwk.gz", ".nwk"):
        if filename.endswith(ext):
            return filename[: -len(ext)]
    return filename


def load_tree(fp: str) -> bte.MATree:
    if fp.endswith(".gz"):
        with tempfile.NamedTemporaryFile(suffix=".pb", delete=False) as tmp:
            with gzip.open(fp, "rb") as gz:
                shutil.copyfileobj(gz, tmp)
            tmp_path = tmp.name
        try:
            t = bte.MATree(tmp_path)
        finally:
            os.unlink(tmp_path)
    else:
        t = bte.MATree(fp)
    return t


def get_leaves(tree: bte.MATree) -> set:
    return {n.id for n in tree.depth_first_expansion() if n.is_leaf()}


def load_trees(filepaths: list) -> list:
    trees = []
    failed = []
    with Progress(
        SpinnerColumn(),
        TextColumn("[bold blue]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TextColumn("[dim]{task.fields[filename]}"),
        console=console,
    ) as progress:
        task = progress.add_task("Loading trees", total=len(filepaths), filename="")
        for fp in filepaths:
            progress.update(task, filename=os.path.basename(fp))
            try:
                t = load_tree(fp)
                leaves = get_leaves(t)
                trees.append((fp, t, leaves))
            except Exception as e:
                failed.append((fp, str(e)))
            progress.advance(task)

    if failed:
        console.print(f"\n[bold red]Failed to load {len(failed)} file(s):[/bold red]")
        for fp, err in failed:
            console.print(f"  [red]✗[/red] {os.path.basename(fp)}: [dim]{err}[/dim]")

    console.print(
        f"\n[green]✓[/green] Loaded [bold]{len(trees)}[/bold]/"
        f"[bold]{len(filepaths)}[/bold] trees successfully.\n"
    )
    return trees


def find_universal_samples(trees: list, min_tree_freq: float) -> set:
    counts = Counter()
    for _, _, leaves in trees:
        counts.update(leaves)
    threshold = min_tree_freq * len(trees)
    return {s for s, c in counts.items() if c >= threshold}


def print_tree_summary(trees: list, universal: set) -> None:
    table = Table(
        title="Tree Summary",
        box=box.ROUNDED,
        header_style="bold cyan",
        show_lines=False,
    )
    table.add_column("File", style="dim", no_wrap=True)
    table.add_column("Total samples", justify="right")
    table.add_column("Effective samples", justify="right")

    for fp, _, leaves in sorted(trees, key=lambda x: len(x[2]), reverse=True):
        effective = len(leaves - universal)
        table.add_row(
            os.path.basename(fp),
            str(len(leaves)),
            str(effective),
        )
    console.print(table)
    console.print()


def print_universal_samples(universal: set, min_tree_freq: float) -> None:
    if not universal:
        return
    panel_text = "\n".join(sorted(universal))
    console.print(Panel(
        f"[yellow]{panel_text}[/yellow]",
        title=f"[bold yellow]{len(universal)} sample(s) ignored "
              f"(≥{min_tree_freq*100:.0f}% of trees)[/bold yellow]",
        border_style="yellow",
        expand=False,
    ))
    console.print()


def print_overlap_report(trees: list, universal: set) -> None:
    pairs = [
        (fp1, fp2, (l1 - universal) & (l2 - universal), len(l1 - universal), len(l2 - universal))
        for (fp1, _, l1), (fp2, _, l2) in combinations(trees, 2)
    ]
    overlapping = [(fp1, fp2, ov, s1, s2) for fp1, fp2, ov, s1, s2 in pairs if ov]

    if not overlapping:
        console.print("[green]✓ No pairwise overlaps detected.[/green]\n")
        return

    table = Table(
        title="Pairwise Overlaps (effective samples only)",
        box=box.ROUNDED,
        header_style="bold cyan",
    )
    table.add_column("Tree A", style="dim")
    table.add_column("N(A)", justify="right")
    table.add_column("Tree B", style="dim")
    table.add_column("N(B)", justify="right")
    table.add_column("Shared", justify="right", style="bold red")

    for fp1, fp2, ov, s1, s2 in sorted(overlapping, key=lambda x: -len(x[2])):
        table.add_row(
            os.path.basename(fp1), str(s1),
            os.path.basename(fp2), str(s2),
            str(len(ov)),
        )
    console.print(table)
    console.print()


def select_non_overlapping(trees: list, universal: set) -> tuple:
    filtered = [(fp, tree, leaves - universal) for fp, tree, leaves in trees]
    sorted_trees = sorted(filtered, key=lambda x: len(x[2]), reverse=True)

    retained = []
    skipped = []
    retained_samples = set()

    for fp, tree, leaves in sorted_trees:
        overlap = leaves & retained_samples
        if overlap:
            skipped.append((fp, leaves, overlap))
        else:
            retained.append((fp, tree, leaves))
            retained_samples |= leaves

    return retained, skipped


def print_selection_report(retained: list, skipped: list) -> None:
    table = Table(
        title="Selection Results",
        box=box.ROUNDED,
        header_style="bold cyan",
        show_lines=False,
    )
    table.add_column("Decision", justify="center", width=8)
    table.add_column("File")
    table.add_column("Effective samples", justify="right")
    table.add_column("Overlap", justify="right")

    kept_names = {fp for fp, _, _ in retained}

    for fp, _, leaves in retained:
        table.add_row(
            "[bold green]KEEP[/bold green]",
            os.path.basename(fp),
            str(len(leaves)),
            "[dim]—[/dim]",
        )
    for fp, leaves, overlap in skipped:
        table.add_row(
            "[bold red]SKIP[/bold red]",
            Text(os.path.basename(fp), style="dim"),
            str(len(leaves)),
            str(len(overlap)),
        )

    console.print(table)
    console.print()


def print_final_list(retained: list) -> None:
    console.print(Panel(
        "\n".join(fp for fp, _, _ in retained),
        title=f"[bold green]Files to keep ({len(retained)})[/bold green]",
        border_style="green",
    ))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", nargs="+", required=True,
                        help="Input tree files (.pb, .pb.gz, or .nwk)")
    parser.add_argument("--min-tree-freq", type=float, default=1.0,
                        help="Exclude samples present in this fraction of trees "
                             "(default: 1.0 = all trees)")
    args = parser.parse_args()

    console.rule("[bold blue]UShER Tree Deduplicator[/bold blue]")
    console.print()

    trees = load_trees(args.input)

    universal = find_universal_samples(trees, args.min_tree_freq)
    print_universal_samples(universal, args.min_tree_freq)
    print_tree_summary(trees, universal)

    console.rule("[bold blue]Pairwise Overlap Analysis[/bold blue]")
    console.print()
    print_overlap_report(trees, universal)

    console.rule("[bold blue]Selection[/bold blue]")
    console.print()
    retained, skipped = select_non_overlapping(trees, universal)
    print_selection_report(retained, skipped)

    console.rule("[bold blue]Result[/bold blue]")
    console.print()
    print_final_list(retained)


if __name__ == "__main__":
    main()