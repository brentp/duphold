#!/bin/bash
set -euo pipefail

base=$(pwd)

git clone -b devel --depth 1000 git://github.com/nim-lang/nim nim
cd nim
# before unchecked was removed:
git checkout cc5b8c6
sh build_all.sh
export PATH=$(base)/nim/bin:$PATH

cd $base
git clone --depth 1 git://github.com/brentp/duphold.git
cd duphold
nimble install -y kexpr
nimble install -y genoiser

nimble install -y

nim c -d:release src/duphold.nim 
cp ./src/duphold /io
