!
!> Module with the sensitivity routines JMult, JMult_Tx, JMult_T, JMult_T_Tx 
!
module GlobalVariables
    !
    use Grid3D_SG
    !
    use ModelParameterCell_SG
    !
    use ModelOperator_MF
    !
    use ForwardSolverIT_DC
    !
    use ModelCovarianceRec
    !
    !> Global Variables
    class( Grid_t ), allocatable, target :: main_grid
    class( ModelParameter_t ), allocatable :: sigma0, pmodel
    class( ModelOperator_t ), allocatable :: model_operator
    !
    class( ForwardSolver_t ), allocatable, target :: forward_solver
    !
    class( ModelCovarianceRec_t ), allocatable :: model_cov
    !
    !> Program control variables
    character(:), allocatable :: control_file_name, model_file_name, pmodel_file_name, data_file_name, modem_job
    logical :: set_data_groups, has_control_file, has_model_file, has_pmodel_file, has_data_file, verbosis
    !
end module GlobalVariables
!