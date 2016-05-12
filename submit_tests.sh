#!/bin/bash

MATLAB="matlab"
BASE_DIR=$(pwd)
INSTANCES_DIR="${BASE_DIR}/Section_5_4_instances/unsolved_in_Reference_8"
RESULTS_DIR="${BASE_DIR}/Section_5_4_results/unsolved_in_Reference_8"

mkdir -p ${RESULTS_DIR}

cd ${INSTANCES_DIR}
for FILENAME in $(ls instance_*.mat); do
    FILENAME="${FILENAME%.*}"
    echo "${FILENAME}"
    sed -e "s/INSERT_FILE_NAME_HERE/${FILENAME}.mat/g" ${BASE_DIR}/submit_tests_template.m > tmp1.m
    sed -e "s#INSERT_BASE_DIR_HERE#${BASE_DIR}#g" tmp1.m > tmp2.m 
    sed -e "s#INSERT_RESULTS_DIR_HERE#${RESULTS_DIR}#g" tmp2.m > tmp3.m
    sed -e "s#INSERT_INSTANCE_DIR_HERE#${INSTANCES_DIR}#g" tmp3.m > ${FILENAME}.m
    rm -rf tmp?.m
    if [ ! -e "${RESULTS_DIR}/${FILENAME}.mat" ]
    then
        echo "${MATLAB} -r 'cd ${INSTANCES_DIR}; run ${FILENAME}.m; exit' > ${RESULTS_DIR}/${FILENAME}_out.txt" > instance.sh
        chmod +x instance.sh
        #qsub -q MANSCI instance.sh
        #qsub -q sandbox instance.sh
        qsub -q all.q instance.sh
        rm -rf instance.sh
    fi
done
cd ${BASE_DIR}

# Some sed lines above use pound sign (#) to avoid escaping forward
# slashes:
# http://superuser.com/questions/766595/properly-escaping-forward-slash-in-bash-script-for-usage-with-sed
