#!/usr/bin/env bash
# (C) Copyright 2009-2016 ECMWF.
# (C) Copyright 2024 Meteorologisk Institutt.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

# Parameters
compare_script=$1
mpi=$2
omp=$3
test_name=$4
ref_name=$5

# Set number of OpenMP threads
export OMP_NUM_THREADS=${omp}

# Log and test outputs extensions
flog=${ref_name}.log.out
ftest=${ref_name}.test.out

if test -f "${ref_name}"; then
  # Run job, create log and test outputs
  mpirun -n ${mpi} ${test_name} | tee ${flog}
  if test "$?" = "0"; then
    grep -s 'Test     : ' ${flog} > ${ftest}

    test_size=`cat ${ftest} | wc -l`
    ref_size=`cat ${ref_name} | wc -l`
    if [ "${test_size}" != "0" ] || [ "${ref_size}" != "0" ]; then
      # Compare reference and log output
      yaml_name=${ref_name/"testref"/"testinput"}
      yaml_name=${yaml_name/".ref"/".yaml"}
      tol=`grep -s "float relative tolerance" ${yaml_name} | awk '{print $4}'`
      if test "${tol}" = ""; then
        tol=1.0e-12
      fi
      echo "Tolerance: "${tol}
      ${compare_script} ${flog} ${ref_name} $tol 0
    fi
  else
    exit $?
  fi
else
  # Run job
  mpirun -n ${mpi} ${test_name} | tee ${flog}
  exit $?
fi
