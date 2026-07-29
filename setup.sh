#!/bin/bash
set -euo pipefail

log_info() {
    echo "$(date "+%Y-%m-%d %H:%M:%S,%3N") - INFO - $1"
}

log_error() {
    echo "$(date "+%Y-%m-%d %H:%M:%S,%3N") - ERROR - $1" >&2
}

LGVAR_INFO=$(pip show lgvar)
LOCATION=$(echo "$LGVAR_INFO" | grep '^Location:' | awk '{print $2}')
LGVAR_SRC="${LOCATION}/lgvar"
YML_PATH="${LGVAR_SRC}/environment.yml"

if [ ! -d "${LGVAR_SRC}" ]; then
    log_error "lgvar source folder not found! Check pip install"
    exit 1
fi

if [ ! -f "${YML_PATH}" ]; then
    log_error "environment.yml not exists under ${LGVAR_SRC}"
    exit 1
fi

if conda info --envs | awk '{print $1}' | grep -qx "LGVAR"; then
    log_info "Conda environment 'LGVAR' already exists. Skipping creation."
else
    log_info "Creating conda environment 'LGVAR' from ${YML_PATH}..."
    conda env create -f "${YML_PATH}" -q > /dev/null 2>&1
fi

eval "$(conda shell.bash hook)"
conda activate LGVAR

export PATH="${LGVAR_SRC}:$PATH"
log_info "Using 'conda activate LGVAR' to load corresponding dependencies"
log_info "Complete! You can run: LGVAR --help to test"
