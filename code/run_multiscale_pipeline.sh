#!/bin/bash
# ============================================================
# Multi-scale hazard model pipeline (51→52→53→54→55)
# 多尺度危险率模型流水线
#
# Usage:
#   R_PROFILE_USER="" bash code/run_multiscale_pipeline.sh
#   R_PROFILE_USER="" bash code/run_multiscale_pipeline.sh --skip-county
#
# All outputs go to outputs_multiscale/ — does NOT overwrite
# existing province-level results.
# ============================================================

set -euo pipefail

BASE_DIR="/home/dingchenchen/projects/bird-new-distribution-records/tasks/bird_hazard_model_effort_upgrade_v2"
OUT_DIR="outputs_multiscale"
SKIP_COUNTY=""

# Parse args
for arg in "$@"; do
  case $arg in
    --skip-county) SKIP_COUNTY="--skip-county" ;;
  esac
done

COMMON_ARGS="--base-dir $BASE_DIR --output-dir $OUT_DIR $SKIP_COUNTY"
LOG_DIR="$BASE_DIR/$OUT_DIR/logs"
mkdir -p "$LOG_DIR"

echo "[$(date '+%H:%M:%S')] === Pipeline START ==="
echo "  base_dir = $BASE_DIR"
echo "  skip_county = $SKIP_COUNTY"

# Step 1: glmmTMB M0-M4 at 4 scales
if [ -f "$BASE_DIR/$OUT_DIR/results/tables/table_multiscale_coefficients.csv" ]; then
  echo "[$(date '+%H:%M:%S')] 51 output exists, skipping 51"
else
  echo "[$(date '+%H:%M:%S')] === Step 1/5: 51_multiscale_full.R ==="
  R_PROFILE_USER="" Rscript "$BASE_DIR/code/51_multiscale_full.R" $COMMON_ARGS 2>&1 | tee -a "$LOG_DIR/pipeline.log"
  echo "[$(date '+%H:%M:%S')] 51 done (exit=$?)"
fi

# Step 2: Variance partitioning + RF importance
echo "[$(date '+%H:%M:%S')] === Step 2/5: 52_multiscale_varpart_rf.R ==="
R_PROFILE_USER="" Rscript "$BASE_DIR/code/52_multiscale_varpart_rf.R" $COMMON_ARGS 2>&1 | tee -a "$LOG_DIR/pipeline.log"
echo "[$(date '+%H:%M:%S')] 52 done (exit=$?)"

# Step 3: glmmTMB future prediction + maps
echo "[$(date '+%H:%M:%S')] === Step 3/5: 53_multiscale_future_prediction.R ==="
R_PROFILE_USER="" Rscript "$BASE_DIR/code/53_multiscale_future_prediction.R" $COMMON_ARGS 2>&1 | tee -a "$LOG_DIR/pipeline.log"
echo "[$(date '+%H:%M:%S')] 53 done (exit=$?)"

# Step 4: XGBoost + SHAP + future prediction
echo "[$(date '+%H:%M:%S')] === Step 4/5: 54_multiscale_xgboost.R ==="
R_PROFILE_USER="" Rscript "$BASE_DIR/code/54_multiscale_xgboost.R" $COMMON_ARGS 2>&1 | tee -a "$LOG_DIR/pipeline.log"
echo "[$(date '+%H:%M:%S')] 54 done (exit=$?)"

# Step 5: Publication choropleth maps
echo "[$(date '+%H:%M:%S')] === Step 5/5: 55_multiscale_choropleth_maps.R ==="
R_PROFILE_USER="" Rscript "$BASE_DIR/code/55_multiscale_choropleth_maps.R" $COMMON_ARGS 2>&1 | tee -a "$LOG_DIR/pipeline.log"
echo "[$(date '+%H:%M:%S')] 55 done (exit=$?)"

echo "[$(date '+%H:%M:%S')] === Pipeline COMPLETE ==="
echo "Outputs: $BASE_DIR/$OUT_DIR/"
