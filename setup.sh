#!/bin/bash
set -euo pipefail

LGVAR_INFO=$(pip show lgvar)
LOCATION=$(echo "$LGVAR_INFO" | grep '^Location:' | awk '{print $2}')
LGVAR_SRC="${LOCATION}/lgvar"
YML_PATH="${LGVAR_SRC}/environment.yml"

if [ ! -d "${LGVAR_SRC}" ];then
    echo "ERROR: lgvar source folder not found! Check pip install"
    exit 1
fi
if [ ! -f "${YML_PATH}" ];then
    echo "ERROR: environment.yml not exists under ${LGVAR_SRC}"
    exit 1
fi

conda env create -f "${YML_PATH}"
conda activate LGVAR

export PATH="${LGVAR_SRC}:$PATH"
echo "Done! You can run: LGVAR --help to test"
