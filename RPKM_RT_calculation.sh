function calcRPKM() {
    ##calculate rpkm, rpkm = 1e9*mappedReads/(totalMappedReads*genelen)
    ##usage: calcRPKM samfile transcriptName
    sam=$1
    transcriptID=$2
    mappdeReads=`cat ${sam} | awk -v tsID="$transcriptID" '{ if ( $3 == tsID ) count++ } END { print count}'`
    totalReads=`cat ${sam} | awk '{ if ( $0 !~ /^@/) count++ } END {print count}'`
    genelen=`cat $sam | awk -v pattern="^@SQ.*${transcriptID}" '{ if ( $0 ~ pattern )  print $3}'`
    genelen=${genelen##LN:}
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    echo "mappdeReads: $mappdeReads; totalReads: $totalReads; genelen: $genelen"
    echo "--------------------------------------------------------------------"
    rpkm=$(( 1000000000 * $mappdeReads / ( $totalReads * $genelen ) ))
    echo " ${transcriptID} rpkm in $sam is $rpkm"
    echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
}

function calcRT() {
    ##RT stop count, use sam directly..
    ##usage: calcRT samfile transcriptName
    sam=$1
    transcriptID=$2
    ##see http://www.math.utah.edu/docs/info/gawk_5.html for the ~ pattern usage, and pay attention to its disadvantage(expr with backslash should be escaped)
    genelen=`cat $sam | awk -v pattern="^@SQ.*${transcriptID}" '{ if ( $0 ~ pattern )  print $3}'`
    genelen=${genelen##LN:}
    declare -a RT_arr
    ##init
    for i in `seq 0 $(( $genelen - 1 ))`
    do
	RT_arr[$i]=0
    done
    while read line
    do
	if [[ $line == @* ]];then
	    continue
	fi
	f=( $line )
	if [[ ${f[2]} == $transcriptID ]];then
	    RT_arr[${f[3]}]=$(( ${RT_arr[${f[3]}]} + 1 ))
	fi
    done < $sam
    echo "${transcriptID} RT count: ${RT_arr[*]}"
}
