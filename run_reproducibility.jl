include("./HiC_spector.jl");

r=20;

map_file1=ARGS[1];
map_file2=ARGS[2];

X1=readdlm(map_file1,Int64);
X2=readdlm(map_file2,Int64);

N=maximum([maximum(X1[:,1:2]),maximum(X2[:,1:2])]);

e=minimum([minimum(X1[:,1:2]),minimum(X2[:,1:2])]);

if e.>0
	M1=sparse(X1[:,1],X1[:,2],X1[:,3],N,N);
	M2=sparse(X2[:,1],X2[:,2],X2[:,3],N,N);
elseif e.==0;
	N=N+1;
	M1=sparse(X1[:,1]+1,X1[:,2]+1,X1[:,3],N,N);
	M2=sparse(X2[:,1]+1,X2[:,2]+1,X2[:,3],N,N);
end

if ~isequal(M1,M1');
	M1_tmp=M1-spdiagm(diag(M1));
	M1=M1+M1_tmp';
end

if ~isequal(M2,M2');
	M2_tmp=M2-spdiagm(diag(M2));
	M2=M2+M2_tmp';
end

Q,a1,b1=get_reproducibility(M1,M2,r);
	
println("size of maps:",size(M1,1),"\t","reproducibility score=",Q);

