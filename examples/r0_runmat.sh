#!/bin/sh
# r0_runmat.sh
# ord_soumet r0_runmat.sh -v 1 -cpus 4x15x4 -t 7200 -cm 4G -notify complete
# S. Innocenti (sin007) - silvia.innocenti@ec.gc.ca
#
############################## 
#
echo "" 
echo " =========================================== "
echo "   run MBB matlab script `(date +'%Y/%m/%d %H:%M:%S')`"
echo " =========================================== " 
echo "" 


# # load the matlab package
. ssmuse-sh -x /fs/ssm/main/opt/matlab/matlab-R2020a

cd /home/sin007/boot_tide/examples/
matlab -nodisplay -batch ex1_generate_wl_resample
