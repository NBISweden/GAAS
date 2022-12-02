#!/bin/bash

cd bin

if [[ ! ${PWD##*/} == "bin" ]];then
	echo "The script must be run in GAAS/bin folder"
	echo "currently here `pwd`"
	exit
fi


#While looking at all script within the repo we skip all Deprecated folder and what it contains
for i in $(find ../ -not \( -path */Deprecated -prune \) -not \( -path */terraform/* -prune \) -not \( -path */bin -prune \)  -not \( -path */blib -prune \) -name '*.pl' -o -name '*.sh' -o -name '*.py' -o -name '*.r' -o -name '*.R' -o -name '*.rb');do

        name=$(basename $i)

	# skip gaas_refresh_list.sh because must not be in the bin to avoid to be distributed
	if [[ $name == "gaas_refresh_bin.sh" || $name == "run_proxy.sh" || $name == "refresh_cert.sh" 
		|| $name == "get_cert_first.sh" || $name == "verify.sh"  ]] ; then
		continue
	fi

        # add new script
        if [[ ! -f ${name} ]];then

	    while true; do
	    	copyfile="no"
	    	read -p "The script ${name} does not exist in /bin. Do you want to add it? (yes/y or no/n)" yn
	    	case $yn in
	        	[Yy]* ) copyfile="yes"; break;;
	        	[Nn]* ) break;;
	        	* ) echo "Please answer yes or no.";;
	    	esac
	    done

	    if [[ $copyfile == "yes" ]];then
	    	cp $i .
	    fi

        else # upadte script
            echo -e "copy $i into /bin"
            cp $i .
        fi
done
