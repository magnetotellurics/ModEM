#!/bin/bash

OCCAM1DCSEM_ZIP_URL="https://marineemlab.ucsd.edu/Projects/Occam/1DCSEM/Occam1DCSEM_v3.13.zip"
export CSEM_MODEM_DIR="./3D_MT/CSEM_module"
export DIPOLE1D_DIR="Dipole1D"
ZIP_FNAME="occam1dcsem.zip"
OCCAM_DIR="Occam1DCSEM"

cmd_exists()
{
    command -v "$1" >/dev/null 2>&1
}

get_dl_cmd() {
    # Determine the download command, either Curl or Wget
    cmd=""
   
    if cmd_exists wget; then
        cmd="wget -O $ZIP_FNAME"
        echo $cmd
        return 0
    fi
 
    if cmd_exists curl; then
        cmd="curl -o $ZIP_FNAME"
        echo $cmd
        return 0
    fi

    return 1
}

check_for_dipole1d() {

    if [ -e $CSEM_MODEM_DIR/$DIPOLE1D_DIR/Dipole1D.f90 ]; then
        return 0
    fi
    return -1
}


download_dipole1d() {
    force=${1:-0}

    if [ $(check_for_dipole1d) -eq 0 ] && [ $force -eq 0 ]; then
        echo "Dipole1D already appears to be downloaded... skipping downloading"
        echo "To re-download it remove the files within: $CSEM_MODEM_DIR/$DIPOLE1D_DIR"
        echo "Or pass 1 to this function: 'download_dipole1d 1'"
        exit 0
    fi

    if [ ! -d $CSEM_MODEM_DIR ]; then
        echo "ERROR: We are not in the ./src directory of the CSEM branch (could not find $CSEM_MODEM_DIR/$DIPOLE1D_DIR)"
        echo "ERROR: Please ensure we are in 'ModEM-Model/src' directory"
        exit 1
    fi

    if [ ! -d "$CSEM_MODEM_DIR/$DIPOLE1D_DIR" ]; then
        echo "Making CSEM ModEM Directory: $CSEM_MODEM_DIR/$DIPOLE1D_DIR"
        mkdir "$CSEM_MODEM_DIR/$DIPOLE1D_DIR"
    fi

    pushd $CSEM_MODEM_DIR/$DIPOLE1D_DIR >/dev/null

    download_cmd=$(get_dl_cmd)
    if [ $? -ne 0 ]; then
        echo "ERROR: Neither 'wget' or 'curl' are available on this system"
        exit 1
    fi

    $download_cmd $OCCAM1DCSEM_ZIP_URL
    if [ $? -ne 0 ]; then
        popd >/dev/null
        echo "ERROR: there was an error trying to download from $OCCAM1DCSEM_ZIP_URL"
        exit 1
    fi

    unzip $ZIP_FNAME
    if [ $? -ne 0 ]; then
        popd >/dev/null
        echo "ERROR: There was an error unzipping "
        exit 1
    fi

    cp $OCCAM_DIR/Source/* .
    rm -rf $OCCAM_DIR

    popd >/dev/null
}

determine_if_dl_needed() {
    if [ "$csem_type" == "Dipole1D" ] || [ "$csem_type" == "dipole1d" ] \
        || [ "$csem_type" == "Dipole1D+EM1D" ] || [ "$csem_type" == "dipole1d+em1d" ]; then

        # Check to see if we should download Dipole1D
        check_for_dipole1d
        dipole_present=$?
        if [ $dipole_present != 0 ]; then
            echo "Dipole1D is not currently in $CSEM_MODEM_DIR/$DIPOLE1D_DIR"

            while true; do
                read -p "Would you like to have this script automatically download it now? [Yes/No]: " resp
                resp_lower=$(to_lower $resp)

                if [ "$resp_lower" == "yes" ] || [ "$resp_lower" == "y" ]; then
                    return 0
                elif [ "$resp_lower" == "no" ] || [ "$resp_lower" == "n" ]; then
                    echo "ERROR: Asked for Dipole1D CSEM but Dipole1D not in: $CSEM_MODEM_DIR/$DIPOLE1D_DIR"
                    echo "ERROR: Either place Dipole1D into the directory above manually, or use this tool to download it"
                    exit 1
                fi
            done
        fi
    fi
    return 1
}


