#!/bin/bash

ROOT=".root"

export add=0

makelist(){

    cat $1 | while read line 
    do

	#Thing to search for
	prev=${line%_*}
	prev=${prev%_*}
	prev=${prev%_*}
	prev=${prev##*/}
	prev=${prev#*_}

	add=0
	echo $prev

	while read secondline
	do

	    sline=${secondline%_*}
	    sline=${sline%_*}

	    if [[ *$sline* == *$prev* ]]; then
		add=1
		break;
	    fi
	  
	done <<<$(cat $2)


	if [[ $add == 0 ]]; then
	    echo "$line" >> $3
	fi

    done
	
}


rm -rf intrinsic_nue_icarus_reco_full_makeup.list;
rm -rf intrinsic_nue_uboone_reco_full_makeup.list;

#makelist test1.list test2.list makeup_test.list
makelist intrinsic_nue_icarus_reco_total.list intrinsic_nue_icarus_reco_full.list intrinsic_nue_icarus_reco_full_makeup.list

makelist intrinsic_nue_uboone_reco_total.list intrinsic_nue_uboone_reco_full.list intrinsic_nue_uboone_reco_full_makeup.list

#makelist intrinsic_nue_sbnd_reco_total.list intrinsic_nue_sbnd_reco_full.list intrinsic_nue_sbnd_reco_full_makeup.list


#makelist osc_nue_icarus_reco_total.list osc_nue_icarus_reco_full.list osc_nue_icarus_reco_full_makeup.list

#makelist osc_nue_uboone_reco_total.list osc_nue_uboone_reco_full.list osc_nue_uboone_reco_full_makeup.list

#makelist osc_nue_sbnd_reco_total.list osc_nue_sbnd_reco_full.list osc_nue_sbnd_reco_full_makeup.list

#makelist nu_icarus_reco_total.list nu_icarus_reco_full.list nu_icarus_reco_full_makeup.list

#makelist nu_uboone_reco_total.list nu_uboone_reco_full.list nu_uboone_reco_full_makeup.list

#makelist nu_sbnd_reco_total.list nu_sbnd_reco_full.list nu_sbnd_reco_full_makeup.list
