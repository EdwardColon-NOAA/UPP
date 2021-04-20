#!/bin/bash
# Wen Meng 01/2020, Set up for cmake build.
#############################################

set -x

#Clean loaded modules
module purge

hostname
source ./detect_machine.sh
if [[ $(uname -s) == Darwin ]]; then
  readonly MYDIR=$(cd "$(dirname "$(greadlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
else
  readonly MYDIR=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
fi
PATHTR=${PATHTR:-$( cd ${MYDIR}/.. && pwd )}

#Load required modulefiles
module use $PATHTR/modulefiles
modulefile=${MACHINE_ID}
module load $modulefile
# module list

rm -rf build install
mkdir build && cd build
set +x
module list 2>&1 | tee modules_loaded.list
set -x
# cmake -DCMAKE_INSTALL_PREFIX=../install -DBUILD_WITH_WRFIO=ON ../..
cmake -DCMAKE_INSTALL_PREFIX=../install -DBUILD_WITH_WRFIO=ON ../.. 2>&1 | tee output.cmake
# make -j6 
make VERBOSE=1 -j 6 2>&1 | tee output.compile
make install

rm -rf $PATHTR/exec && mkdir $PATHTR/exec
cp $PATHTR/tests/install/bin/upp.x $PATHTR/exec/.
