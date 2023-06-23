#!/bin/bash
#
# OO MPI
#
mpirun -np 22 ./ModEM_MPI -f -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey.dat -cf fwd_ctrl_dip1d.txt -pd pred_csem_dip1d_iso_mpi.dat
#
mpirun -np 22 ./ModEM_MPI -f -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey.dat -cf fwd_ctrl_em1d.txt -pd pred_csem_em1d_iso_mpi.dat
#
mpirun -np 22 ./ModEM_MPI -f -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_fix.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey.dat -cf fwd_ctrl_em1d.txt -pd pred_csem_em1d_ani_mpi.dat
#
# ON
#
cd ../../modem-on/src/
#
mpirun -np 22 ./Mod3DMT_STD -F ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod ../../modem_inputs/ANI_CSEM/Data_Ey.dat pred_csem_dip1d_iso_on.dat esoln control_dip1d
#
mpirun -np 22 ./Mod3DMT_STD -F ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod ../../modem_inputs/ANI_CSEM/Data_Ey.dat pred_csem_dip1d_iso_on.dat esoln control_em1d
#
mpirun -np 22 ./Mod3DMT_STD -F ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_fix.mod ../../modem_inputs/ANI_CSEM/Data_Ey.dat pred_csem_dip1d_iso_on.dat esoln control_em1d
#
# OO SERIAL
#
cd ../../modem-oo/src/
#
./ModEM_SERIAL -f -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey.dat -cf fwd_ctrl_dip1d.txt -pd pred_csem_dip1d_iso_serial.dat
#
./ModEM_SERIAL -f -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey.dat -cf fwd_ctrl_em1d.txt -pd pred_csem_em1d_iso_serial.dat
#
./ModEM_SERIAL -f -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_fix.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey.dat -cf fwd_ctrl_em1d.txt -pd pred_csem_em1d_ani_serial.dat
#
