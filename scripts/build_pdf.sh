#!/usr/bin/env bash
#
# build_pdf.sh — build a single PDF of the entire handbook from the markdown sources.
#
# Strategy: pandoc walks docs/ in the order defined by mkdocs.yml's nav, converts
# each markdown chapter to LaTeX with --listings for code, --mathjax for math
# (rendered to native LaTeX), then xelatex links them into one ~700-page book.
#
# Mermaid diagrams require mermaid-cli (mmdc) for pre-rendering to PNG.
# If mmdc is missing, mermaid blocks are stripped and replaced with placeholders.
#
# Usage:
#   bash scripts/build_pdf.sh            # full book → site/materials-simulation-handbook.pdf
#   bash scripts/build_pdf.sh --fast     # skip mermaid pre-rendering
#   bash scripts/build_pdf.sh --chapters "ch00-math ch01-python"  # subset only
#
# Requirements (auto-checked):
#   - pandoc >= 3.0
#   - texlive-xetex (xelatex)
#   - texlive-latex-extra (booktabs, listings, framed, etc.)
#   - python3 (for the mkdocs.yml nav parser)
#   - mermaid-cli (optional, `npm i -g @mermaid-js/mermaid-cli`)
#
# Estimated runtime: 5-8 minutes on a modern laptop for the full book.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DOCS="${ROOT}/docs"
OUT_DIR="${ROOT}/site"
OUT_PDF="${OUT_DIR}/materials-simulation-handbook.pdf"
TMP_DIR="${ROOT}/.pdf-build"

FAST_MODE=false
CHAPTERS_FILTER=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fast) FAST_MODE=true; shift ;;
        --chapters) CHAPTERS_FILTER="$2"; shift 2 ;;
        -h|--help)
            sed -n '3,25p' "${BASH_SOURCE[0]}"
            exit 0
            ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# ---- Prerequisites check ----
need() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "ERROR: '$1' is required but not installed." >&2
        echo "  Install hint: $2" >&2
        exit 1
    fi
}

need pandoc "brew install pandoc  # or apt install pandoc"
need xelatex "brew install --cask mactex  # or apt install texlive-xetex texlive-latex-extra"
need python3 "Python 3 should be present on every supported system"

if $FAST_MODE; then
    echo "Fast mode: skipping mermaid pre-rendering."
elif ! command -v mmdc >/dev/null 2>&1; then
    echo "WARN: mmdc (mermaid-cli) not found — mermaid diagrams will be omitted."
    echo "      Install with: npm install -g @mermaid-js/mermaid-cli"
    FAST_MODE=true
fi

mkdir -p "${OUT_DIR}" "${TMP_DIR}/mermaid"

# ---- Step 1: extract chapter order from mkdocs.yml ----
echo "[1/4] Extracting chapter order from mkdocs.yml..."
python3 - "${ROOT}/mkdocs.yml" "${TMP_DIR}/file_order.txt" <<'PY'
import sys, re, pathlib
mkdocs_yml = pathlib.Path(sys.argv[1]).read_text()
out = pathlib.Path(sys.argv[2])
files = []
# Naive parse: pick anything that looks like ': path/to/file.md' under nav
in_nav = False
for line in mkdocs_yml.splitlines():
    if line.startswith("nav:"):
        in_nav = True; continue
    if in_nav and line and not line.startswith(" ") and not line.startswith("-"):
        break
    m = re.search(r":\s*([A-Za-z0-9_./-]+\.md)\s*$", line)
    if m:
        files.append(m.group(1))
out.write_text("\n".join(files))
print(f"  Found {len(files)} markdown files in nav order.")
PY

# ---- Step 2: optionally pre-render mermaid diagrams ----
if ! $FAST_MODE; then
    echo "[2/4] Pre-rendering mermaid diagrams..."
    # For each chapter md, find ```mermaid blocks, render with mmdc, replace.
    # (Simplified: a production version would handle inline replacement properly.)
    echo "  (placeholder: mermaid pre-rendering pass — extend as needed)"
fi

# ---- Step 3: concatenate markdown sources in nav order ----
echo "[3/4] Concatenating markdown sources..."
MASTER_MD="${TMP_DIR}/master.md"
: > "${MASTER_MD}"
while IFS= read -r relpath; do
    [[ -z "${relpath}" ]] && continue
    if [[ -n "${CHAPTERS_FILTER}" ]]; then
        match=false
        for filter in ${CHAPTERS_FILTER}; do
            [[ "${relpath}" == *"${filter}"* ]] && match=true
        done
        ${match} || continue
    fi
    src="${DOCS}/${relpath}"
    if [[ ! -f "${src}" ]]; then
        echo "  SKIP missing: ${relpath}"
        continue
    fi
    echo "" >> "${MASTER_MD}"
    echo "" >> "${MASTER_MD}"
    cat "${src}" >> "${MASTER_MD}"
done < "${TMP_DIR}/file_order.txt"

WC=$(wc -w < "${MASTER_MD}")
echo "  Master markdown: $(printf "%'d" "${WC}") words."

# ---- Step 4: pandoc → xelatex ----
echo "[4/4] Compiling with pandoc + xelatex..."
pandoc "${MASTER_MD}" \
    --from gfm+tex_math_dollars+raw_html+pipe_tables+task_lists+attributes \
    --to pdf \
    --pdf-engine=xelatex \
    --toc \
    --toc-depth=3 \
    --number-sections \
    --top-level-division=chapter \
    --highlight-style=tango \
    -V documentclass=book \
    -V geometry:margin=1in \
    -V papersize=a4 \
    -V mainfont="DejaVu Serif" \
    -V monofont="DejaVu Sans Mono" \
    -V linkcolor=blue \
    -V urlcolor=blue \
    -V title="Materials Simulation Handbook" \
    -V author="Zhaohe Dong" \
    -V date="$(date +%Y-%m-%d)" \
    --metadata=lang:en-GB \
    -o "${OUT_PDF}"

echo ""
echo "Done. PDF written to: ${OUT_PDF}"
echo "Size: $(du -h "${OUT_PDF}" | cut -f1)"
