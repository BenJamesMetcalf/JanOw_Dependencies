#!/bin/bash -l

#. /usr/share/Modules/init/bash
###This wrapper script validates the input arguments and creates the job-control.txt file which is needed to submit the qsub array job to the cluster.###

while getopts :i:f:s:o: option
do
    case $option in
        i) instr_dir=$OPTARG;;
	f) instr_file=$OPTARG;;
        #r) allDB_dir=$OPTARG;;
	s) species=$OPTARG;;
        o) output_dir=$OPTARG;;
    esac
done

###Check if batch directory and reference database directory arguments were given and if they exist###
type="blah"
if [[ -z "$instr_dir" && -z "$instr_file" ]]
then
    echo "Either the instrument directory (-i argument) or a list of flowcells in a text file (-f argument) needs to be given."
    exit 1
elif [[ -n "$instr_dir" && -n "$instr_file" ]]
then
    echo "Either the instrument directory (-i argument) or a list of flowcells in a text file (-f argument) needs to be given but not both."
    exit 1
else
    if [[ ! -z "$instr_dir" ]]
    then
	if [[ -d "$instr_dir" ]]
	then
            instr_dir=$(echo "$instr_dir" | sed 's/\/$//g')
            echo "The sequence directory is in the following location: $instr_dir"
	    type="dir"
	else
            echo "This sequence directory is not in the correct format or doesn't exist."
            echo "Make sure you provide the full directory path (/root/path/sequence_directory)."
            exit 1
	fi
    elif [[ ! -z "$instr_file" ]]
    then
	if [[ -s "$instr_file" ]]
	then
	    echo "The flowcell input file is in the following location: $instr_file"
	    type="file"
	else
	    echo "The flowcell input file is not in the correct format or doesn't exist."
	    echo "Make sure you provide the full directory path (/root/path/flowcell_file)."
	    exit 1
	fi
    fi
fi

#if [[ ! -z "$allDB_dir" ]]
#then
#    if [[ -d "$allDB_dir" ]]
#    then
#        allDB_dir=$(echo "$allDB_dir" | sed 's/\/$//g')
#        echo "The references directory is in the following location: $allDB_dir"
#    else
#        echo "This reference directory is not in the correct format or doesn't exist."
#        echo "Make sure you provide the full directory path (/root/path/reference_directory)."
#        exit 1
#    fi
#else
#    echo "No reference database directory path argument given."
#    exit 1
#fi

if [[ -z "$species" ]]
then
    echo "The species identifier has not been given. This argument must be 'SPN', 'GBS' or 'GAS'"
    exit 1
elif [[ "$species" != SPN && "$species" != GBS && "$species" != GAS ]]
then
    echo "The species input argument must be 'SPN', 'GBS' or 'GAS'"
    exit 1
fi

if [[ -z "$output_dir" ]]
then
    echo "The files will be output into the default directory '$species_Typing_Analysis'."
    if [[ ! -d ~/"$species"_Typing_Analysis ]]
    then
        mkdir ~/"$species"_Typing_Analysis
        out_dir="~/$species_Typing_Analysis"
        eval out_dir=$out_dir
        echo "The output directory has been created: $out_dir"
    else
        out_dir="~/$species_Typing_Analysis"
        eval out_dir=$out_dir
    fi
elif [[ ! -d "$output_dir" ]]
then
    output_dir=$(echo "$output_dir" | sed 's/\/$//g')
    mkdir "$output_dir"
    out_dir="$output_dir"
else
    output_dir=$(echo "$output_dir" | sed 's/\/$//g')
    out_dir="$output_dir"
fi

###Start Doing Stuff###
if [[ "$species" == "SPN" ]]
then
    cd "/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/External/SPN_Scripts_Reference/"
    comm="StrepLab-JanOw_SPN-wrapr.sh"
    allDB_dir="/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/External/SPN_Scripts_Reference/SPN_Reference_DB/"
elif [[ "$species" == "GAS" ]]
then
    cd "/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/External/GAS_Scripts_Reference"
    comm="StrepLab-JanOw_GAS-wrapr.sh"
    allDB_dir="/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/External/GAS_Scripts_Reference/GAS_Reference_DB/"
elif [[ "$species" == "GBS" ]]
then
    cd "/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/External/GBS_Scripts_Reference"
    comm="StrepLab-JanOw_GBS-wrapr.sh"
    allDB_dir="/scicomp/groups/OID/NCIRD/DBD/RDB/Strep_Lab/External/GBS_Scripts_Reference/GBS_Reference_DB/"
else
    echo "The input for the species '-s' arugment must be either 'SPN', 'GAS' or 'GBS'"
    exit 1
fi

if [[ "$type" == "dir" ]]
then
    instr_dir_star="${instr_dir}/*"
    for flowPath in $instr_dir_star
    do
	#echo "Flow Cell path is: $flowPath";
	flowName=$(basename "$flowPath" | sed 's/_new//g')
	flowID=$(echo "$flowName" | sed 's/.*-\([A-Za-z0-9]\+\).*/\1/g')
	instrID=$(echo "$flowName" | sed 's/.*_\(M032\(32\|20\)\)_.*/\1/g')
	#echo "Flow Cell: $flowName | Flowcell ID: $flowID | Instrument: $instrID"
	instrDir="/scicomp/instruments/18-B-429_Illumina-MiSeq-$instrID/$flowName/Data/Intensities/BaseCalls/"
	#ls "$instrDir"
	echo "screen -dmS $flowID bash $comm -s $instrDir -r $allDB_dir -o $out_dir/$flowName"
	screen -dmS "$flowID" bash "$comm" -s "$instrDir" -r "$allDB_dir" -o "$out_dir/$flowName"
    done
elif [[ "$type" == "file" ]]
then
    while read flowName
    do
	flowName=$(echo "$flowName" | sed 's/_new//g')
	flowID=$(echo "$flowName" | sed 's/.*-\([A-Za-z0-9]\+\).*/\1/g')
        instrID=$(echo "$flowName" | sed 's/.*_\(M032\(32\|20\)\)_.*/\1/g')
        #echo "Flow Cell: $flowName | Flowcell ID: $flowID | Instrument: $instrID"
        instrDir="/scicomp/instruments/18-B-429_Illumina-MiSeq-$instrID/$flowName/Data/Intensities/BaseCalls/"
        #ls "$instrDir"
        echo "screen -dmS $flowID bash $comm -s $instrDir -r $allDB_dir -o $out_dir/$flowName"
	screen -dmS "$flowID" bash "$comm" -s "$instrDir" -r "$allDB_dir" -o "$out_dir/$flowName"
    done < "$instr_file"
fi
