#!/usr/bin/env bash

set -euo pipefail
shopt -s dotglob
trap "exit" INT

# Directory where this script resides
THIS_DIR=$(
  cd "$(dirname "${BASH_SOURCE[0]}")"
  pwd
)

# Where the source code is
PROJECT_ROOT_DIR="$(realpath "${THIS_DIR}/..")"

source "${THIS_DIR}/lib/set_locales.sh"

source "${PROJECT_ROOT_DIR}/.env.example"
if [ -f "${PROJECT_ROOT_DIR}/.env" ]; then
  source "${PROJECT_ROOT_DIR}/.env"
fi

DOCKERHUB_ORG="nextstrain"
DOCKERHUB_PROJECT="nextclade_builder"
DOCKERHUB_REPO="${DOCKERHUB_ORG}/${DOCKERHUB_PROJECT}"

COMMIT_HASH=${CIRCLE_SHA1:=$(git rev-parse --short HEAD)}

docker build -f "${PROJECT_ROOT_DIR}/Dockerfile" \
  --build-arg UID="${UID:=$(id -u)}" \
  --build-arg GID="${GID:=$(id -g)}" \
  --tag "${DOCKERHUB_REPO}:latest" \
  --tag "${DOCKERHUB_REPO}:${COMMIT_HASH}" \
  .
