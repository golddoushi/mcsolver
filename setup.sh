#----------------------------------------------------------------------
# some simple codes for compiling c-lib for mcsolver
# written by Liang Liu@SZU
# May.11.2020 3:33 a.m.

CC=gcc  # c compiler

# set the path for python include and lib directories if you are sure
# otherwise, leave them empty. we would automatically search for them
INC=
LIB=
#----------------------------------------------------------------------
# No need to modify following 

if [ -z $INC ] || [ -z $LIB ]; then
echo automatically search relied python files
# first inquire if there is python3
PyC=$(which python3)
if [[ $PyC == *python3 ]] && [ -x $PyC ]; then
    ver1=$($PyC --version 2>&1|awk '{print $2}'|awk -F '.' '{print $1}')
    ver2=$($PyC --version 2>&1|awk '{print $2}'|awk -F '.' '{print $2}')
    ver3=$($PyC --version 2>&1|awk '{print $2}'|awk -F '.' '{print $3}')
    echo python3 found
else
    echo cannot find python3
    PyC=$(which python)
    if [[ $PyC == *python ]]; then
        echo maybe it is called python
        ver1=$($PyC --version 2>&1|awk '{print $2}'|awk -F '.' '{print $1}')
        ver2=$($PyC --version 2>&1|awk '{print $2}'|awk -F '.' '{print $2}')
        ver3=$($PyC --version 2>&1|awk '{print $2}'|awk -F '.' '{print $3}')
        echo python version is : $ver1.$ver2.$ver3
        if [ $ver1 -lt 3 ];then
            echo Error
            echo python3 is not installed or configured properly
            echo now exit
            exit
        fi
    else
        echo Error
        echo I cannot find python interpreter
        echo maybe you can build dynamic libs manually
        echo now i am sliding away, bye~
        exit
    fi
fi

# locate the parent dir for python
oldifs=$IFS
IFS='/'
path_dat=($PyC)
parentDir='/'
IFS=$oldifs
nnode=${#path_dat[@]}
let nnode=$nnode-2
ncount=1
while [ $ncount -lt $nnode ]; do
    parentDir=$parentDir${path_dat[$ncount]}/
    let ncount=$ncount+1
done
echo parent dir is: $parentDir

# locate include director
INC=$parentDir"include/"

postVer=$ver1.$ver2
if [ -f $INC"Python.h" ]; then
    echo Python.h found
elif [ -f $INC"python"$ver1.$ver2"/Python.h" ]; then
    INC=$INC"python"$ver1.$ver2
    echo Python.h found
elif [ -f $INC"python"$ver1"."$ver2"m/Python.h" ]; then
    INC=$INC"python"$ver1"."$ver2"m/"
    postVer=$postVer"m"
    echo Python.h found
else
    echo Error
    echo cannot find Python.h
    exit
fi
echo include dir is: $INC

# locate lib director
LIB=$parentDir"lib/"
if [ -f $parentDir"lib/libpython"$postVer".so" ]; then
    LIB=$parentDir"lib/"
    echo .so found in $LIB
elif [ -f $parentDir"lib64/libpython"$postVer".so" ]; then
    LIB=$parentDir"lib64/"
    echo .so found in $LIB
else
    echo Error
    echo cannot find libpython$postVer".so"
    exit
fi
else
    echo using custome setting for py inc. and lib.
fi

# start compile
echo start compile c-lib
$CC -std=c99 ./isingLib.c -shared -fPIC -o ./isinglib.so -I$INC -L$LIB -lpython3
$CC -std=c99 ./heisenbergLib.c -shared -fPIC -o ./heisenberglib.so -I$INC -L$LIB -lpython3
$CC -std=c99 ./xyLib.c -shared -fPIC -o ./xylib.so -I$INC -L$LIB -lpython3
if [ -f "./isinglib.so" ] && [ -f "./heisenberglib.so" ] && [ -f "./xylib.so" ]; then
    echo c-lib successfully built
else
    echo Error occured when compiling c-lib
    exit
fi

# try a small test
