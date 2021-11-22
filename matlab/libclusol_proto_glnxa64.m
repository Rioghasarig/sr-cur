function [methodinfo,structs,enuminfo,ThunkLibName]=libclusol_proto
%LIBCLUSOL_PROTO Create structures to define interfaces found in 'clusol'.

%This function was generated by loadlibrary.m parser version  on Sun Nov 21 17:47:07 2021
%perl options:'clusol.i -outfile=libclusol_proto.m -thunkfile=libclusol_thunk_glnxa64.c -header=clusol.h'
ival={cell(1,0)}; % change 0 to the actual number of functions to preallocate the data.
structs=[];enuminfo=[];fcnNum=1;
fcns=struct('name',ival,'calltype',ival,'LHS',ival,'RHS',ival,'alias',ival,'thunkname', ival);
MfilePath=fileparts(mfilename('fullpath'));
ThunkLibName=fullfile(MfilePath,'libclusol_thunk_glnxa64');
% void clu1fac ( int64_t * m , int64_t * n , int64_t * nelem , int64_t * lena , int64_t * ap , int64_t * aq , int64_t * frank , int64_t * luparm , double * parmlu , double * a , int64_t * indc , int64_t * indr , int64_t * p , int64_t * q , int64_t * lenc , int64_t * lenr , int64_t * locc , int64_t * locr , int64_t * iploc , int64_t * iqloc , int64_t * ipinv , int64_t * iqinv , double * w , int64_t * inform ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu1fac'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'longPtr'};fcnNum=fcnNum+1;
% void clu6sol ( int64_t * mode , int64_t * m , int64_t * n , double * v , double * w , int64_t * lena , int64_t * luparm , double * parmlu , double * a , int64_t * indc , int64_t * indr , int64_t * p , int64_t * q , int64_t * lenc , int64_t * lenr , int64_t * locc , int64_t * locr , int64_t * inform ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu6sol'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr'};fcnNum=fcnNum+1;
% void clu8rpc ( int64_t * mode1 , int64_t * mode2 , int64_t * m , int64_t * n , int64_t * jrep , double * v , double * w , int64_t * lena , int64_t * luparm , double * parmlu , double * a , int64_t * indc , int64_t * indr , int64_t * p , int64_t * q , int64_t * lenc , int64_t * lenr , int64_t * locc , int64_t * locr , double * lv , int64_t * li , int64_t * lj , int64_t * lenlv , int64_t * inform , double * diag , double * vnorm ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu8rpc'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr'};fcnNum=fcnNum+1;
% void clu6mul ( int64_t * mode , int64_t * m , int64_t * n , double * v , double * w , int64_t * lena , int64_t * luparm , double * parmlu , double * a , int64_t * indc , int64_t * indr , int64_t * p , int64_t * q , int64_t * lenc , int64_t * lenr , int64_t * locc , int64_t * locr ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu6mul'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr'};fcnNum=fcnNum+1;
% void clu8adc ( int64_t * mode , int64_t * m , int64_t * n , double * v , double * w , int64_t * lena , int64_t * luparm , double * parmlu , double * a , int64_t * indc , int64_t * indr , int64_t * p , int64_t * q , int64_t * lenc , int64_t * lenr , int64_t * locc , int64_t * locr , int64_t * inform , double * diag , double * vnorm ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu8adc'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr'};fcnNum=fcnNum+1;
% void clu8adr ( int64_t * m , int64_t * n , double * w , int64_t * lena , int64_t * luparm , double * parmlu , double * a , int64_t * indc , int64_t * indr , int64_t * p , int64_t * q , int64_t * lenc , int64_t * lenr , int64_t * locc , int64_t * locr , int64_t * inform , double * diag ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu8adr'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'doublePtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr'};fcnNum=fcnNum+1;
% void clu8dlc ( int64_t * m , int64_t * n , int64_t * jdel , int64_t * lena , int64_t * luparm , double * parmlu , double * a , int64_t * indc , int64_t * indr , int64_t * p , int64_t * q , int64_t * lenc , int64_t * lenr , int64_t * locc , int64_t * locr , int64_t * inform ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu8dlc'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr'};fcnNum=fcnNum+1;
% void clu9diagu ( int64_t * m , int64_t * n , int64_t * lena , int64_t * luparm , double * diagU , double * a , int64_t * indc , int64_t * indr , int64_t * ip , int64_t * iq , int64_t * lenr , int64_t * locr , int64_t * inform ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu9diagu'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr'};fcnNum=fcnNum+1;
% void clu9mod ( int64_t * m , int64_t * n , double * beta , double * v , double * w , int64_t * lena , int64_t * luparm , double * parmlu , double * a , int64_t * indc , int64_t * indr , int64_t * p , int64_t * q , int64_t * lenc , int64_t * lenr , int64_t * locc , int64_t * locr , int64_t * inform ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu9mod'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr'};fcnNum=fcnNum+1;
% void clu8rpr ( int64_t * mode1 , int64_t * mode2 , int64_t * m , int64_t * n , int64_t * irep , double * v , double * w , double * wnew , int64_t * lena , int64_t * luparm , double * parmlu , double * a , int64_t * indc , int64_t * indr , int64_t * p , int64_t * q , int64_t * lenc , int64_t * lenr , int64_t * locc , int64_t * locr , int64_t * lenlv , int64_t * li , int64_t * lj , double * lv , int64_t * inform ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu8rpr'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'longPtr'};fcnNum=fcnNum+1;
% void clu9clr ( int64_t * m , int64_t * n , int64_t * lena , int64_t * luparm , double * parmlu , double * a , int64_t * indc , int64_t * indr , int64_t * p , int64_t * q , int64_t * lenr , int64_t * locc , int64_t * locr , double * c , double * lmax , int64_t * inform ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu9clr'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'doublePtr', 'longPtr'};fcnNum=fcnNum+1;
% void clu1or2 ( int64_t * n , int64_t * numa , int64_t * lena , double * a , int64_t * inum , int64_t * jnum , int64_t * lenc , int64_t * locc ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='clu1or2'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'longPtr', 'longPtr', 'longPtr', 'doublePtr', 'longPtr', 'longPtr', 'longPtr', 'longPtr'};fcnNum=fcnNum+1;
% void ccarrqr ( int32_t * m , int32_t * n , int32_t * blkMax , double * A , int32_t * lda , int32_t * jpvt , double * tau , int32_t * info ); 
fcns.thunkname{fcnNum}='voidvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrvoidPtrThunk';fcns.name{fcnNum}='ccarrqr'; fcns.calltype{fcnNum}='Thunk'; fcns.LHS{fcnNum}=[]; fcns.RHS{fcnNum}={'int32Ptr', 'int32Ptr', 'int32Ptr', 'doublePtr', 'int32Ptr', 'int32Ptr', 'doublePtr', 'int32Ptr'};fcnNum=fcnNum+1;
methodinfo=fcns;