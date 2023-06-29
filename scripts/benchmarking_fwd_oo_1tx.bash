#!/bin/bash
#
################################################################
# FWD
################################################################
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -f -m ../docs/vti_results/model_ISO_h10.mod -d ../docs/vti_results/Data_Z_1tx_1rx.dat -pd pred_mt_iso_oo.dat -es esol_mt_iso_oo.bin -cf ../docs/vti_results/oo_fwd_control_em1d
#
RESULT=$?
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > MT FWD ISO FAILS!"
    return $RESULT
    #
fi
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -f -m ../docs/vti_results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../docs/vti_results/Data_Z_1tx_1rx.dat -pd pred_mt_vti_oo.dat -es esol_mt_vti_oo.bin -cf ../docs/vti_results/oo_fwd_control_em1d
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > MT FWD VTI FAILS!"
    return $RESULT
    #
fi
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -f -m ../docs/vti_results/model_ISO_h10.mod -d ../docs/vti_results/Data_Ey_1tx_1rx.dat -pd pred_dip1d_oo.dat -es esol_dip1d_oo.bin -cf ../docs/vti_results/oo_fwd_control_dip1d
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > DIPOLE1D FWD FAILS!"
    return $RESULT
    #
fi
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -f -m ../docs/vti_results/model_ISO_h10.mod -d ../docs/vti_results/Data_Ey_1tx_1rx.dat -pd pred_em1d_iso_oo.dat -es esol_em1d_iso_oo.bin -cf ../docs/vti_results/oo_fwd_control_em1d
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > EM1D ISO FWD FAILS!"
    return $RESULT
    #
fi
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -f -m ../docs/vti_results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../docs/vti_results/Data_Ey_1tx_1rx.dat -pd pred_em1d_vti_oo.dat -es esol_em1d_vti_oo.bin -cf ../docs/vti_results/oo_fwd_control_em1d
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > EM1D VTI FWD FAILS!"
    return $RESULT
    #
fi
#
################################################################
# JMULT
################################################################
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -j -m ../docs/vti_results/model_ISO_h0.mod -pm ../docs/vti_results/model_ISO_h10.mod -d ../docs/vti_results/Data_Z_1tx_1rx.dat -jm jmhat_mt_iso_oo.dat -cf ../docs/vti_results/oo_fwd_control_em1d
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > MT JMULT ISO FAILS!"
    return $RESULT
    #
fi
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -j -m ../docs/vti_results/model_VTI_h0_v0.mod -pm ../docs/vti_results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../docs/vti_results/Data_Z_1tx_1rx.dat -jm jmhat_mt_vti_oo.dat -cf ../docs/vti_results/oo_fwd_control_em1d
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > MT JMULT VTI FAILS!"
    return $RESULT
    #
fi
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -j -m ../docs/vti_results/model_ISO_h0.mod -pm ../docs/vti_results/model_ISO_h10.mod -d ../docs/vti_results/Data_Ey_1tx_1rx.dat -jm jmhat_dip1d_oo.dat -cf ../docs/vti_results/oo_fwd_control_dip1d
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > DIPOLE1D JMULT FAILS!"
    return $RESULT
    #
fi
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -j -m ../docs/vti_results/model_ISO_h0.mod -pm ../docs/vti_results/model_ISO_h10.mod -d ../docs/vti_results/Data_Ey_1tx_1rx.dat -jm jmhat_em1d_iso_oo.dat -cf ../docs/vti_results/oo_fwd_control_em1d
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > EM1D ISO JMULT FAILS!"
    return $RESULT
    #
fi
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -j -m ../docs/vti_results/model_VTI_h0_v0.mod -pm ../docs/vti_results/VTI_H10_V10/model_VTI_h10_v10.mod -d ../docs/vti_results/Data_Ey_1tx_1rx.dat -jm jmhat_em1d_vti_oo.dat -cf ../docs/vti_results/oo_fwd_control_em1d
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > EM1D VTI JMULT FAILS!"
    return $RESULT
    #
fi
#
################################################################
# JMULT_T
################################################################
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -jt -m ../docs/vti_results/model_ISO_h0.mod -d pred_mt_iso_oo.dat -dm dsigma_mt_iso_oo.rho -cf ../docs/vti_results/oo_fwd_control_em1d
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > MT JMULT ISO FAILS!"
    return $RESULT
    #
fi
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -jt -m ../docs/vti_results/model_VTI_h0_v0.mod -d pred_mt_vti_oo.dat -dm dsigma_mt_vti_oo.rho -cf ../docs/vti_results/oo_fwd_control_em1d
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > MT JMULT VTI FAILS!"
    return $RESULT
    #
fi
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -jt -m ../docs/vti_results/model_ISO_h0.mod -d pred_dip1d_oo.dat -dm dsigma_dip1d_oo.rho -cf ../docs/vti_results/oo_fwd_control_dip1d
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > DIPOLE1D JMULT_T FAILS!"
    return $RESULT
    #
fi
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -jt -m ../docs/vti_results/model_ISO_h0.mod -d pred_em1d_iso_oo.dat -dm dsigma_em1d_iso_oo.rho -cf ../docs/vti_results/oo_fwd_control_em1d
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > EM1D ISO JMULT_T FAILS!"
    return $RESULT
    #
fi
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -jt -m ../docs/vti_results/model_VTI_h0_v0.mod -d pred_em1d_vti_oo.dat -dm dsigma_em1d_vti_oo.rho -cf ../docs/vti_results/oo_fwd_control_em1d
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > EM1D VTI JMULT_T FAILS!"
    return $RESULT
    #
fi
#
################################################################
# NLCG
################################################################
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -i -m ../docs/vti_results/model_ISO_h0.mod -d pred_mt_iso_oo.dat -ci ../docs/vti_results/oo_inv_control_nlcg -cf ../docs/vti_results/oo_fwd_control_em1d -o mt_iso_oo_nlcg
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > MT NLCG ISO FAILS!"
    return $RESULT
    #
fi
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -i -m ../docs/vti_results/model_VTI_h0_v0.mod -d pred_mt_vti_oo.dat -ci ../docs/vti_results/oo_inv_control_nlcg -cf ../docs/vti_results/oo_fwd_control_em1d -o mt_vti_oo_nlcg
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > MT NLCG VTI FAILS!"
    return $RESULT
    #
fi
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -i -m ../docs/vti_results/model_ISO_h0.mod -d pred_dip1d_oo.dat -ci ../docs/vti_results/oo_inv_control_nlcg -cf ../docs/vti_results/oo_fwd_control_dip1d -o dip1d_oo_nlcg
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > DIPOLE1D NLCG FAILS!"
    return $RESULT
    #
fi
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -i -m ../docs/vti_results/model_ISO_h0.mod -d pred_em1d_iso_oo.dat -ci ../docs/vti_results/oo_inv_control_nlcg -cf ../docs/vti_results/oo_fwd_control_em1d -o em1d_iso_oo_nlcg
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > EM1D ISO NLCG FAILS!"
    return $RESULT
    #
fi
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -i -m ../docs/vti_results/model_VTI_h0_v0.mod -d pred_em1d_vti_oo.dat -ci ../docs/vti_results/oo_inv_control_nlcg -cf ../docs/vti_results/oo_fwd_control_em1d -o em1d_vti_oo_nlcg
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > EM1D VTI NLCG FAILS!"
    return $RESULT
    #
fi
#
################################################################
# DCG
################################################################
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -i -m ../docs/vti_results/model_ISO_h0.mod -d pred_mt_iso_oo.dat -ci ../docs/vti_results/oo_inv_control_dcg -cf ../docs/vti_results/oo_fwd_control_em1d -o mt_iso_oo_dcg
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > MT DCG ISO FAILS!"
    return $RESULT
    #
fi
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -i -m ../docs/vti_results/model_VTI_h0_v0.mod -d pred_mt_vti_oo.dat -ci ../docs/vti_results/oo_inv_control_dcg -cf ../docs/vti_results/oo_fwd_control_em1d -o mt_vti_oo_dcg
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > MT DCG VTI FAILS!"
    return $RESULT
    #
fi
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -i -m ../docs/vti_results/model_ISO_h0.mod -d pred_dip1d_oo.dat -ci ../docs/vti_results/oo_inv_control_dcg -cf ../docs/vti_results/oo_fwd_control_dip1d -o dip1d_oo_dcg
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > DIPOLE1D DCG FAILS!"
    return $RESULT
    #
fi
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -i -m ../docs/vti_results/model_ISO_h0.mod -d pred_em1d_iso_oo.dat -ci ../docs/vti_results/oo_inv_control_dcg -cf ../docs/vti_results/oo_fwd_control_em1d -o em1d_iso_oo_dcg
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > EM1D ISO DCG FAILS!"
    return $RESULT
    #
fi
#
mpirun -np 2 --hostfile ../../hostfile_22.txt ./ModEM_MPI -i -m ../docs/vti_results/model_VTI_h0_v0.mod -d pred_em1d_vti_oo.dat -ci ../docs/vti_results/oo_inv_control_dcg -cf ../docs/vti_results/oo_fwd_control_em1d -o em1d_vti_oo_dcg
#
if [ "$RESULT" -ne "0" ]; then
    #
    echo " > EM1D VTI DCG FAILS!"
    return $RESULT
    #
fi
#