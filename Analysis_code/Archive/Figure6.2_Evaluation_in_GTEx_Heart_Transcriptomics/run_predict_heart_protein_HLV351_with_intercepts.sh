#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PAPER_ctOWAS_DIR:-/Users/admin/Library/CloudStorage/Dropbox/Paper_ctOWAS}"

export PAPER_ctOWAS_DIR="$PROJECT_DIR"

Rscript --vanilla "$SCRIPT_DIR/3_predict_celltype_heart_protein_HLV351_with_intercepts.R"
