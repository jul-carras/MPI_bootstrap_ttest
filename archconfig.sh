#!/bin/bash

# archconfig.sh - generate an MPI hosts file, with slots computed from
#            physical or logical CPUs of available Arch servers
#
# usage: ./archconfig.sh <mode>
#       where <mode> = 'phy' OR 'log'

# set mode, physical or logical
mode=$1
case $mode in
	[Pp]hy)
	;&
	[Pp]hysical)
        mode="phy"
		echo -e "\nslots= physical CPUs\n"
	;;
	[Ll]og)
	;&
	[Ll]ogical)
        mode="log"
		echo -e "\nslots= logical CPUs\n"
	;;
	*)
		echo -e "\tusage: $0 <mode> (mode = 'phy' OR 'log')"
		exit 1
esac

# change hosts filename here
hostfile=my_hosts

# hosts to check
>$hostfile
local_hosts=(arch01 arch02 arch03 arch04 arch05 arch06 arch07 arch08 arch09 arch10)


for host in ${local_hosts[@]};
do
	echo `echo $host`

    # slots per machine
    tpc=`ssh $host "lscpu | grep Thread" | awk '{print $4}'`
    if [ "$tpc" \> 1 ] && [ "$mode" == "log" ]; then
        tpc=1
    fi

	# write to hostfile
	numprocs=`ssh $host "grep proc /proc/cpuinfo | wc -l"` && \
		echo $host slots=$((numprocs / tpc)) >> $hostfile
done
