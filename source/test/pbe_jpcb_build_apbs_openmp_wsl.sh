#!/usr/bin/env bash
set -euo pipefail

version="${1:-3.4.1}"
jobs="${2:-4}"
root="${HOME}/.local/src/apbs-${version}"
build="${HOME}/.local/build/apbs-${version}-openmp"
prefix="${HOME}/.local/apbs-${version}-openmp"

if [[ ! -d "${root}/.git" ]]; then
  mkdir -p "$(dirname "${root}")"
  git clone --branch "v${version}" --recurse-submodules https://github.com/Electrostatics/apbs.git "${root}"
fi
git -C "${root}" submodule update --init --recursive
rm -rf "${build}"

cmake -S "${root}" -B "${build}" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${prefix}" \
  -DENABLE_OPENMP=ON \
  -DENABLE_PYTHON=OFF \
  -DENABLE_BEM=OFF \
  -DENABLE_GEOFLOW=OFF \
  -DENABLE_PBAM=OFF \
  -DENABLE_PBSAM=OFF \
  -DENABLE_TESTS=OFF \
  -DAPBS_STATIC_BUILD=OFF
cmake --build "${build}" --parallel "${jobs}"
cmake --install "${build}"

exe="${prefix}/bin/apbs"
[[ -x "${exe}" ]]
grep -q '^ENABLE_OPENMP:BOOL=ON$' "${build}/CMakeCache.txt"
cat > "${prefix}/apbs_openmp_build_manifest.txt" <<EOF
schema_version=1
apbs_version=${version}
source_commit=$(git -C "${root}" rev-parse HEAD)
cmake_enable_openmp=true
build_jobs=${jobs}
executable=${exe}
EOF
echo "APBS_OPENMP_EXECUTABLE=${exe}"
