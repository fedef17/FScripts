#!/bin/bash

set -e

# Set paths and file names
VARIABLES_FILE="variables_longrunmip.txt"   # File containing list of variables (one per line)
OUTPUT_DIR="/home/fabiano/work_big/irods_move/LongRunMIP/"        # Directory to save merged files

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Check if variables file exists
if [[ ! -f "$VARIABLES_FILE" ]]; then
  echo "Error: Variables file '$VARIABLES_FILE' not found!"
  exit 1
fi

for mem in stabilization-hist-1990 stabilization-ssp585-2025 stabilization-ssp585-2050 stabilization-ssp585-2065 stabilization-ssp585-2080 stabilization-ssp585-2100; do
    # Read variables from the file and process each
    while IFS=" " read -r varname1 varname2 domain || [[ -n "$varname2" ]]; do
        echo "Processing variable: $varname2 for domain: $domain"
        DATA_DIR="/home/fabiano/work_big/irods_move/${mem}/r1i1p1f1/${domain}/${varname2}/"                 # Directory containing NetCDF files

        do_regrid=true
        if [ "$varname2" == "msftyz" ]; then do_regrid=false; fi

        # Output file for the merged data
        outdir="$OUTPUT_DIR/orig/${varname1}/"
        outdir_regrid="$OUTPUT_DIR/regrid/${varname1}/"
        mkdir -p $outdir
        mkdir -p $outdir_regrid

        if [ "$varname1" == "sic" ]; then
            output_file="${outdir}/${varname1}_mon_EC-Earth3_${mem}_1000.nc"
            output_file_regrid="${outdir_regrid}/${varname1}_mon_EC-Earth3_${mem}_1000_g025.nc"
        else
            output_file="${outdir}/${varname1}_ann_EC-Earth3_${mem}_1000.nc"
            output_file_regrid="${outdir_regrid}/${varname1}_ann_EC-Earth3_${mem}_1000_g025.nc"
        fi

        if [ -e "$output_file" ] && { [ -e "$output_file_regrid" ] || [ "$do_regrid" = false ]; }; then
            echo "Output file $output_file already exists. Skipping calculation."
        else
            mkdir -p $outdir/tmp/
            mkdir -p $outdir_regrid/tmp/

            
            first_file=$(ls "$DATA_DIR"/*.nc | head -n 1)
            if $do_regrid; then
                cdo genbil,grid_25.des ${first_file} ${outdir_regrid}/remapweights.nc
            fi

            for nc_file in "$DATA_DIR"/*.nc; do
                filename=$(basename "$nc_file")
                
                if [ "$varname1" == "sic" ]; then
                    ln -s $nc_file $outdir/tmp/$filename
                else
                    cdo yearmean $nc_file $outdir/tmp/$filename
                fi
                
                if $do_regrid; then
                    cdo remap,grid_25.des,${outdir_regrid}/remapweights.nc $outdir/tmp/$filename $outdir_regrid/tmp/$filename
                fi
            done

            # Use cdo to merge files and extract variable
            cdo cat $outdir/tmp/*nc "$output_file"
            if $do_regrid; then
                cdo cat $outdir_regrid/tmp/*nc "$output_file_regrid"
            fi

            # Check if the output file was created successfully
            if [[ -f "$output_file" ]]; then
                echo "Merged file created: $output_file"
            else
                echo "Error: Failed to create merged file for variable '$variable'."
                exit 1
            fi

            if do_regrid; then
                if [[ -f "$output_file_regrid" ]]; then
                    echo "Merged file created: $output_file_regrid"
                else
                    echo "Error: Failed to create merged file for variable '$variable'."
                    exit 1
                fi
            fi

            # Clean up temporary file
            rm -rf $outdir/tmp/
            if $do_regrid; then rm -rf $outdir_regrid/tmp/; fi
        fi


    done < "$VARIABLES_FILE"

    ## missing netTOA, surf, pr_ocn, rls
    echo "All variables processed successfully for mem ${mem}."
done
