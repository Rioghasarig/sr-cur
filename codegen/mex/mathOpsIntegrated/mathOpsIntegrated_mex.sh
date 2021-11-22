MATLAB="/global/software/sl-7.x86_64/modules/tools/matlab/r2017b"
Arch=glnxa64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/global/home/users/oekenta/.matlab/R2017b"
OPTSFILE_NAME="./setEnv.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for mathOpsIntegrated" > mathOpsIntegrated_mex.mki
echo "CC=$CC" >> mathOpsIntegrated_mex.mki
echo "CFLAGS=$CFLAGS" >> mathOpsIntegrated_mex.mki
echo "CLIBS=$CLIBS" >> mathOpsIntegrated_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> mathOpsIntegrated_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> mathOpsIntegrated_mex.mki
echo "CXX=$CXX" >> mathOpsIntegrated_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> mathOpsIntegrated_mex.mki
echo "CXXLIBS=$CXXLIBS" >> mathOpsIntegrated_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> mathOpsIntegrated_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> mathOpsIntegrated_mex.mki
echo "LDFLAGS=$LDFLAGS" >> mathOpsIntegrated_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> mathOpsIntegrated_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> mathOpsIntegrated_mex.mki
echo "Arch=$Arch" >> mathOpsIntegrated_mex.mki
echo "LD=$LD" >> mathOpsIntegrated_mex.mki
echo OMPFLAGS= >> mathOpsIntegrated_mex.mki
echo OMPLINKFLAGS= >> mathOpsIntegrated_mex.mki
echo "EMC_COMPILER=gcc" >> mathOpsIntegrated_mex.mki
echo "EMC_CONFIG=optim" >> mathOpsIntegrated_mex.mki
"/global/software/sl-7.x86_64/modules/tools/matlab/r2017b/bin/glnxa64/gmake" -B -f mathOpsIntegrated_mex.mk
