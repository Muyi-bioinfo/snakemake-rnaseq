#!/usr/bin/env bash
# ============================================================
# RNA-seq 流程容器构建脚本
# 支持 Docker 和 Apptainer 双引擎
# ============================================================
set -euo pipefail

IMAGE_NAME="${IMAGE_NAME:-rnaseq-pipeline}"
IMAGE_TAG="${IMAGE_TAG:-latest}"
FULL_NAME="${IMAGE_NAME}:${IMAGE_TAG}"
SIF_FILE="${IMAGE_NAME}.sif"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# ---- 颜色输出 ----
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

log_info()  { echo -e "${GREEN}[INFO]${NC}  $*"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC}  $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }
log_step()  { echo -e "${CYAN}[STEP]${NC} $*"; }

usage() {
    cat <<EOF
Usage: $0 <command> [options]

Commands:
  docker       Build Docker image
  singularity  Build Apptainer SIF from local .def file
  sif-from-docker   Convert local Docker image to SIF
  all          Build Docker, then convert to SIF
  pull         Pull pre-built image from Docker Hub → convert to SIF
  clean        Remove built images
  test         Test the built container (run tool version checks)

Options:
  IMAGE_NAME=name    Image name (default: rnaseq-pipeline)
  IMAGE_TAG=tag      Image tag (default: latest)
  REGISTRY=url       Docker registry for push/pull (default: none)

Examples:
  $0 docker                        # Build Docker image locally
  $0 singularity                   # Build SIF with apptainer (needs --fakeroot)
  $0 all                           # Docker + SIF
  $0 pull REGISTRY=docker.io/bio   # Pull & convert to SIF
  $0 test                          # Verify installed tools

Environment:
  IMAGE_NAME     Container image name (default: rnaseq-pipeline)
  IMAGE_TAG      Container tag (default: latest)
EOF
    exit 1
}

# ---- 参数解析 ----
CMD="${1:-}"
[[ -z "$CMD" ]] && usage
shift || true

# ---- 构建 Docker 镜像 ----
build_docker() {
    log_step "构建 Docker 镜像: ${FULL_NAME}"

    cd "$PROJECT_DIR"

    docker build \
        --tag "${FULL_NAME}" \
        --file containers/Dockerfile \
        .

    log_info "Docker 镜像构建完成: ${FULL_NAME}"

    # 显示镜像大小
    docker images "${IMAGE_NAME}:${IMAGE_TAG}"
}

# ---- 推送到 Registry ----
push_docker() {
    local registry="${REGISTRY:-}"
    if [[ -n "$registry" ]]; then
        local remote="${registry}/${FULL_NAME}"
        log_step "推送镜像到: ${remote}"
        docker tag "${FULL_NAME}" "${remote}"
        docker push "${remote}"
        log_info "推送完成: ${remote}"
    fi
}

# ---- 从 .def 文件构建 SIF ----
build_singularity() {
    local def_file="${SCRIPT_DIR}/apptainer.def"

    if [[ ! -f "$def_file" ]]; then
        log_error "Apptainer 定义文件不存在: ${def_file}"
        exit 1
    fi

    # 检测 apptainer 还是 singularity
    local cmd=""
    if command -v apptainer &>/dev/null; then
        cmd="apptainer"
    elif command -v singularity &>/dev/null; then
        cmd="singularity"
    else
        log_error "未找到 apptainer 或 singularity 命令"
        exit 1
    fi

    log_step "使用 ${cmd} 构建 SIF: ${SIF_FILE}"

    cd "$PROJECT_DIR"

    # 尝试 fakeroot，如果不行则提示需要 sudo
    if ${cmd} build --fakeroot "${SIF_FILE}" "${def_file}" 2>/dev/null; then
        log_info "SIF 构建完成: ${SIF_FILE}"
    else
        log_warn "--fakeroot 不可用，尝试使用 sudo..."
        sudo ${cmd} build "${SIF_FILE}" "${def_file}"
        log_info "SIF 构建完成: ${SIF_FILE}"
    fi
}

# ---- 从本地 Docker 镜像转换 SIF ----
docker_to_sif() {
    local cmd=""
    if command -v apptainer &>/dev/null; then
        cmd="apptainer"
    elif command -v singularity &>/dev/null; then
        cmd="singularity"
    else
        log_error "未找到 apptainer 或 singularity 命令"
        exit 1
    fi

    log_step "将 Docker 镜像转换为 SIF..."

    cd "$PROJECT_DIR"

    ${cmd} build "${SIF_FILE}" "docker-daemon://${FULL_NAME}"

    log_info "SIF 转换完成: ${SIF_FILE}"
}

# ---- 从 Docker Hub 拉取并转换 ----
pull_to_sif() {
    local registry="${REGISTRY:-docker.io}"
    local cmd=""
    if command -v apptainer &>/dev/null; then
        cmd="apptainer"
    elif command -v singularity &>/dev/null; then
        cmd="singularity"
    else
        log_error "未找到 apptainer 或 singularity 命令"
        exit 1
    fi

    local remote="docker://${registry}/${FULL_NAME}"
    log_step "从 ${remote} 拉取并转换为 SIF..."

    cd "$PROJECT_DIR"

    ${cmd} build "${SIF_FILE}" "${remote}"

    log_info "SIF 构建完成: ${SIF_FILE}"
}

# ---- 清理 ----
clean() {
    log_step "清理容器镜像..."
    docker rmi "${FULL_NAME}" 2>/dev/null || true
    rm -f "${SIF_FILE}"
    log_info "清理完成"
}

# ---- 测试 ----
test_container() {
    log_step "测试容器工具..."

    _run_tool() {
        local tool_cmd="$1"
        echo -n "  ${tool_cmd} ... "
        if docker run --rm "${FULL_NAME}" bash -c "${tool_cmd}" &>/dev/null; then
            echo -e "${GREEN}OK${NC}"
        else
            echo -e "${RED}FAIL${NC}"
        fi
    }

    _run_tool "fastqc --version"
    _run_tool "trimmomatic -version"
    _run_tool "hisat2 --version"
    _run_tool "samtools --version"
    _run_tool "featureCounts -v"
    _run_tool "multiqc --version"
    _run_tool "R --version"
    _run_tool 'Rscript -e "library(DESeq2);cat(packageVersion(\"DESeq2\"),\"\n\")"'
    _run_tool 'Rscript -e "library(clusterProfiler);cat(packageVersion(\"clusterProfiler\"),\"\n\")"'

    # 如果有 SIF 文件也测试一下
    if [[ -f "${SIF_FILE}" ]]; then
        log_step "测试 SIF 文件..."
        local cmd="apptainer"
        command -v apptainer &>/dev/null || cmd="singularity"

        echo -n "  fastqc --version ... "
        if ${cmd} exec "${SIF_FILE}" fastqc --version &>/dev/null; then
            echo -e "${GREEN}OK${NC}"
        else
            echo -e "${RED}FAIL${NC}"
        fi
    fi

    log_info "测试完成"
}

# ---- 主逻辑 ----
case "$CMD" in
    docker)
        build_docker
        ;;
    singularity)
        build_singularity
        ;;
    sif-from-docker)
        build_docker
        docker_to_sif
        ;;
    all)
        build_docker
        docker_to_sif
        ;;
    pull)
        pull_to_sif
        ;;
    clean)
        clean
        ;;
    test)
        test_container
        ;;
    *)
        usage
        ;;
esac
