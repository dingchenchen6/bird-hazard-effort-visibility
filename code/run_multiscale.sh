#!/bin/bash
# run_multiscale.sh
# 服务器上执行多尺度风险集构建与建模的便捷脚本

set -e  # 遇到错误立即退出

# ============================================================
# 配置部分 - 请根据服务器实际情况修改
# ============================================================

# 基础工作目录
BASE_DIR="${BASE_DIR:-/home/dingchenchen/Documents/New records/bird-new-distribution-records/tasks/bird_hazard_model_effort_upgrade_v2}"

# 输出目录（相对于BASE_DIR）
OUTPUT_DIR="${OUTPUT_DIR:-outputs_multiscale}"

# R脚本路径
R_SCRIPT="${BASE_DIR}/code/50_multiscale_risk_build_and_model.R"

# 年份窗口
YEAR_MIN="${YEAR_MIN:-2002}"
YEAR_MAX="${YEAR_MAX:-2024}"

# 日志文件
LOG_DIR="${BASE_DIR}/${OUTPUT_DIR}/logs"
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/run_$(date +%Y%m%d_%H%M%S).log"

# ============================================================
# 函数定义
# ============================================================

print_header() {
    echo "=========================================="
    echo "多尺度风险集构建与建模"
    echo "执行脚本: $(basename ${R_SCRIPT})"
    echo "基础目录: ${BASE_DIR}"
    echo "输出目录: ${OUTPUT_DIR}"
    echo "年份窗口: ${YEAR_MIN}-${YEAR_MAX}"
    echo "日志文件: ${LOG_FILE}"
    echo "=========================================="
    echo ""
}

check_environment() {
    echo "[$(date +%H:%M:%S)] 检查环境..."

    # 检查基础目录
    if [ ! -d "${BASE_DIR}" ]; then
        echo "ERROR: 基础目录不存在: ${BASE_DIR}"
        exit 1
    fi

    # 检查R脚本
    if [ ! -f "${R_SCRIPT}" ]; then
        echo "ERROR: R脚本不存在: ${R_SCRIPT}"
        exit 1
    fi

    # 检查R是否可用
    if ! command -v Rscript &> /dev/null; then
        echo "ERROR: Rscript 命令未找到，请加载R模块或安装R"
        exit 1
    fi

    # 检查数据文件
    local data_dir="${BASE_DIR}/data/raw"
    if [ ! -d "${data_dir}" ]; then
        echo "WARNING: 数据目录不存在: ${data_dir}"
    else
        echo "[$(date +%H:%M:%S)] 检查关键数据文件..."
        for f in "events_100km_grid_assigned.csv" "climate_metrics_province_year.csv" "effort_panel_upgraded.csv"; do
            if [ ! -f "${data_dir}/${f}" ] && [ ! -L "${data_dir}/${f}" ]; then
                echo "  WARNING: ${f} 不存在"
            else
                echo "  OK: ${f}"
            fi
        done
    fi

    echo "[$(date +%H:%M:%S)] 环境检查完成"
    echo ""
}

run_analysis() {
    echo "[$(date +%H:%M:%S)] 开始执行分析..."
    echo "[$(date +%H:%M:%S)] 命令: Rscript ${R_SCRIPT} --base-dir ${BASE_DIR} --output-dir ${OUTPUT_DIR} --year-min ${YEAR_MIN} --year-max ${YEAR_MAX}"
    echo ""

    # 执行R脚本，同时输出到屏幕和日志
    Rscript "${R_SCRIPT}" \
        --base-dir "${BASE_DIR}" \
        --output-dir "${OUTPUT_DIR}" \
        --year-min "${YEAR_MIN}" \
        --year-max "${YEAR_MAX}" \
        2>&1 | tee -a "${LOG_FILE}"

    local exit_code=${PIPESTATUS[0]}

    if [ ${exit_code} -eq 0 ]; then
        echo ""
        echo "[$(date +%H:%M:%S)] 分析成功完成!"
    else
        echo ""
        echo "[$(date +%H:%M:%S)] ERROR: 分析失败，退出码: ${exit_code}"
        echo "请查看日志: ${LOG_FILE}"
        exit ${exit_code}
    fi
}

print_output_summary() {
    echo ""
    echo "=========================================="
    echo "输出文件摘要"
    echo "=========================================="

    local output_path="${BASE_DIR}/${OUTPUT_DIR}"

    if [ ! -d "${output_path}" ]; then
        echo "WARNING: 输出目录不存在"
        return
    fi

    echo "输出目录: ${output_path}"
    echo ""

    # 风险集
    echo "风险集数据:"
    for f in "${output_path}/data/derived/"*.csv "${output_path}/data/derived/"*.rds; do
        if [ -f "$f" ]; then
            local size=$(du -h "$f" | cut -f1)
            echo "  $f [${size}]"
        fi
    done

    echo ""
    echo "结果表格:"
    for f in "${output_path}/results/tables/"*.csv; do
        if [ -f "$f" ]; then
            local size=$(du -h "$f" | cut -f1)
            echo "  $f [${size}]"
        fi
    done

    echo ""
    echo "图表:"
    for f in "${output_path}/figures/main/"*; do
        if [ -f "$f" ]; then
            local size=$(du -h "$f" | cut -f1)
            echo "  $f [${size}]"
        fi
    done

    echo ""
    echo "日志:"
    ls -lh "${output_path}/logs/" 2>/dev/null || echo "  (无日志文件)"
}

# ============================================================
# 主程序
# ============================================================

main() {
    print_header

    check_environment

    # 创建输出目录
    mkdir -p "${BASE_DIR}/${OUTPUT_DIR}/"{data/derived,results/{tables,diagnostics},figures/{main,diagnostics},logs}

    # 记录开始时间
    local start_time=$(date +%s)

    # 运行分析
    run_analysis

    # 记录结束时间
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    local minutes=$((duration / 60))
    local seconds=$((duration % 60))

    echo ""
    echo "=========================================="
    echo "执行摘要"
    echo "=========================================="
    echo "开始时间: $(date -d @${start_time} '+%Y-%m-%d %H:%M:%S')"
    echo "结束时间: $(date -d @${end_time} '+%Y-%m-%d %H:%M:%S')"
    echo "总耗时: ${minutes}分${seconds}秒"
    echo ""

    print_output_summary

    echo ""
    echo "=========================================="
    echo "完成！"
    echo "查看详细日志: cat ${LOG_FILE}"
    echo "查看输出目录: ls -la ${BASE_DIR}/${OUTPUT_DIR}/"
    echo "=========================================="
}

# 执行主程序
main "$@"
