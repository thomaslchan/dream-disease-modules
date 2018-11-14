#!/bin/bash

chmod +x Pascal
if [[ `uname` == 'Darwin' ]];then
    echo "not installing openBLAS; on OsX, there is no need to.." 
fi
if [[ `uname` != 'Darwin' ]];
then
echo "installing openBLAS"
echo "trying to unzip xianyi-OpenBLAS-v0.2.12-0-g7e4e195.zip : if this fails do the  unzipping manually and restart the script"
sleep 4
unzip xianyi-OpenBLAS-v0.2.12-0-g7e4e195.zip
cd xianyi-OpenBLAS-48f06dd/ 
make
make install PREFIX="../lib/openBLASlib"
cd ../
cd lib/openBLASlib/lib
ln -s libopenblas.so.0 libblas.so.3
ln -s libopenblas.so.0 liblapack.so.3
cd ../../../
fi
cd lib/fortranlibs
if [ `command -v gfortran | wc -l` -ne 0 ];then
    echo "compiling fortran libraries..."
    make clean
    make libmvtpack.so
    make libmvtpack.dylib
fi
if [ `command -v gfortran | wc -l` -eq 0 ];then
    echo "gfortran compiler not detected. "
    if [[ `uname` == 'Darwin' ]];then
        echo "using precompiled libraries.." 
    fi
    if [[ `uname` != 'Darwin' ]];
    then
        echo "<genescoring=max>-option might not work. Pascal might be slow. Consider installing gcc compiler family." 
    fi
fi
cd ../../
if [ `command -v java | wc -l` -eq 0 ];then
    echo "java VM not installed. Please install java VM (64-Bit) before proceeding." 
fi
if [ `command -v java | wc -l` -ne 0 ];then
    if [ `java -version 2>&1| grep "64\-Bit" | wc -l` -eq 0 ];then
        echo "java 32-Bit installed instead of 64-Bit. You might run into memory problems when running large data sets."
        echo "Consider installing the 64-Bit java VM."
    fi
fi

echo "done."
