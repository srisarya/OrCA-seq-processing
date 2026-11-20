#!/usr/bin/env python3
"""
make_filters.py — Generate a barbell filters.txt from `barbell inspect` output.

Usage:
    python make_filters.py inspect_output.txt [-o filters.txt]
    barbell inspect -i anno.tsv | python make_filters.py [-o filters.txt]
    python make_filters.py inspect_output.txt -v   # verbose: show accept/reject

A pattern is included in filters.txt only when it contains at least one complete
Ftag AND at least one complete Rtag that can form a valid amplicon pair:

    Forward amplicon:  Ftag[fw, …] … Rtag[rc, …]
    RC amplicon:       Rtag[fw, …] … Ftag[rc, …]

Cut markers (>> / <<) are inserted at amplicon boundaries so barbell trims the
adapters and retains only the insert.  Ligation / concat reads (N amplicon pairs
per read) receive numbered cut groups (>>1/<<1, >>2/<<2, …).

Fflank / Rflank elements are kept in the output line verbatim but never used as
pair anchors; orphan rc-orientation tags are likewise preserved but never paired.
"""

import re
import sys
import argparse
import fileinput
from dataclasses import dataclass


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Element:
    tag:         str   # Ftag | Rtag | Fflank | Rflank
    orientation: str   # fw | rc
    label:       str   # e.g. * or NB01 or ~group1
    position:    str   # e.g. @left(0..250)

    @property
    def is_full_tag(self) -> bool:
        return self.tag in ("Ftag", "Rtag")

    @property
    def is_opener(self) -> bool:
        """An opener starts an amplicon: Ftag[fw] or Rtag[fw]."""
        return self.is_full_tag and self.orientation == "fw"

    @property
    def closing_tag(self) -> tuple[str, str] | None:
        """The (tag, orientation) that closes an amplicon opened by this element."""
        if self.tag == "Ftag" and self.orientation == "fw":
            return ("Rtag", "rc")
        if self.tag == "Rtag" and self.orientation == "fw":
            return ("Ftag", "rc")
        return None

    def render(self, cut: str | None = None) -> str:
        fields = [self.orientation, self.label, self.position]
        if cut:
            fields.append(cut)
        return f"{self.tag}[{', '.join(fields)}]"


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

_HEADER_RE = re.compile(r"^\s*Pattern\s+\d+:\s+\d+\s+occurrences\s*$")
_BODY_RE   = re.compile(r"^\s{2,}(\S.+)$")
_ELEM_RE   = re.compile(
    r"(?P<tag>Ftag|Rtag|Fflank|Rflank)"
    r"\[(?P<orientation>fw|rc)"
    r",\s*(?P<label>[^,\]]+)"
    r",\s*(?P<position>@\w+\(\d+\.\.\d+\))"
    r"[^\]]*\]"          # ignore any existing cut specifiers / extra fields
)


def parse_inspect(lines: list[str]) -> list[str]:
    """Extract raw pattern strings from inspect output lines."""
    patterns = []
    expect_body = False
    for line in lines:
        if _HEADER_RE.match(line):
            expect_body = True
            continue
        if expect_body:
            m = _BODY_RE.match(line)
            if m:
                patterns.append(m.group(1).strip())
            expect_body = False
    return patterns


def parse_elements(pattern: str) -> list[Element]:
    """Parse a raw pattern string into a sequence of Element objects."""
    return [
        Element(
            tag=m.group("tag"),
            orientation=m.group("orientation"),
            label=m.group("label").strip(),
            position=m.group("position"),
        )
        for m in _ELEM_RE.finditer(pattern)
    ]


# ---------------------------------------------------------------------------
# Amplicon pair detection
# ---------------------------------------------------------------------------

def find_amplicon_pairs(elements: list[Element]) -> list[tuple[int, int]]:
    """
    Return (opener_idx, closer_idx) for each complete amplicon, left-to-right.

    Greedy: each opener claims the next available compatible closer.
    Fflank/Rflank and rc-orientation tags are skipped as openers.
    Each index is used at most once.
    """
    pairs: list[tuple[int, int]] = []
    used: set[int] = set()

    for i, el in enumerate(elements):
        if i in used or not el.is_opener:
            continue
        target = el.closing_tag          # (tag, orientation) we need
        for j in range(i + 1, len(elements)):
            if j in used:
                continue
            cel = elements[j]
            if (cel.tag, cel.orientation) == target:
                pairs.append((i, j))
                used.add(i)
                used.add(j)
                break

    return pairs


# ---------------------------------------------------------------------------
# Filter line construction
# ---------------------------------------------------------------------------

def make_filter_line(elements: list[Element], pairs: list[tuple[int, int]]) -> str:
    """
    Render elements as a barbell filter line with cut markers.

    Single pair  → >> / <<  (no group number needed)
    Multiple pairs → >>1/<<1, >>2/<<2, …  (cut groups for concat reads)
    """
    multi = len(pairs) > 1
    cut_map: dict[int, str] = {}
    
    for g, (opener, closer) in enumerate(pairs, start=1):
        cut_map[opener] = f">>{g}" if multi else ">>"
        cut_map[closer] = f"<<{g}" if multi else "<<"

    return "__".join(el.render(cut_map.get(i)) for i, el in enumerate(elements))


# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------

def process_patterns(
    patterns: list[str],
    verbose: bool = False,
) -> list[str]:
    """
    Convert raw inspect pattern strings to deduplicated filter lines.

    Acceptance criteria:
      1. Has at least one Ftag AND at least one Rtag  (fast pre-reject for the
         majority of single-tag patterns before the O(n²) pair search).
      2. At least one valid (opener, closer) amplicon pair can be formed.
    """
    seen: set[str] = set()
    result: list[str] = []

    for raw in patterns:
        elements = parse_elements(raw)

        has_ftag = any(e.tag == "Ftag" for e in elements)
        has_rtag = any(e.tag == "Rtag" for e in elements)
        if not (has_ftag and has_rtag):
            if verbose:
                print(f"  SKIP (no Ftag+Rtag): {raw}", file=sys.stderr)
            continue

        pairs = find_amplicon_pairs(elements)
        if not pairs:
            if verbose:
                print(f"  SKIP (no valid pair): {raw}", file=sys.stderr)
            continue

        line = make_filter_line(elements, pairs)
        if line not in seen:
            seen.add(line)
            result.append(line)
            if verbose:
                print(f"  KEEP ({len(pairs)} pair(s)): {line}", file=sys.stderr)

    return result


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate barbell filters.txt from `barbell inspect` output.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python make_filters.py inspect.txt\n"
            "  barbell inspect -i anno.tsv | python make_filters.py\n"
            "  python make_filters.py inspect.txt -o my_filters.txt -v\n"
        ),
    )
    parser.add_argument(
        "inspect_file",
        nargs="?",
        help="Inspect output file (omit or use - to read from stdin)",
    )
    parser.add_argument(
        "-o", "--output",
        default="filters.txt",
        help="Output file path (default: filters.txt)",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print accept/reject decisions to stderr",
    )
    args = parser.parse_args()

    # Read all input lines — fileinput handles file, stdin, and "-" uniformly
    src = [args.inspect_file] if args.inspect_file else []
    try:
        lines = list(fileinput.input(files=src or ("-",)))
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    patterns = parse_inspect(lines)
    if not patterns:
        print(
            "ERROR: No patterns found in input.\n"
            "Make sure you are passing the output of `barbell inspect`.",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"Parsed {len(patterns)} patterns.", file=sys.stderr)

    filter_lines = process_patterns(patterns, verbose=args.verbose)
    print(f"Generated {len(filter_lines)} filter lines.", file=sys.stderr)

    with open(args.output, "w") as fh:
        fh.writelines(line + "\n" for line in filter_lines)

    print(f"Written to: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
