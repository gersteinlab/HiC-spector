using HDF5; 
using JLD;
using MAT;
#using Graphs;
#using DataFrames;
#using CurveFit;
#using Distributions;


function get_Laplacian(M);
	K=vec(sum(M,1));
	i_nz=find(K.>0);
	D_nz=spdiagm(K[i_nz]);
	D_isq=spdiagm(1./sqrt(K[i_nz]));
	#L_nz=D_nz-M[i_nz,i_nz];
	#the smallest ev of L is 0.
	#in many networks, because of the existenc of singleton, we expect more than 1 zero ev..
	#if we do normalization, 0=lambda1<=lambda_1<=lambda_2.,,,<=2

	#to avoid matrix multiplication with 0, inf. actually
	#L_norm(i,j)=1 if i=j
	#L_norm(i,j)=-1/sqrt(deg(i)*deg(j)))
	#0 otherwise.
	Ln_nz=M[i_nz,i_nz]*D_isq;
	Ln_nz=I-D_isq*Ln_nz;
	n=size(M,1);
	#Ln=extend_mat(Ln_nz,i_nz,n);
	#check for Laplacian, K(v)-w(v,v) for diag, if 0, both 0..normalization goes to 1-0/0=0;
	#r=collect(1:n);
	#r=(r-1)*n+r;
	#Ln[r[Ln[r].==0]]=1;
	Ln_nz=(Ln_nz+Ln_nz')/2;
	return Ln_nz;
end

function get_reproducibility(M1,M2,num_evec);
	
	k1=sum(M1,1);
	k2=sum(M2,1);
	iz=find(k1+k2.>0);

	M1b=M1[iz,iz];
	M2b=M2[iz,iz];

	i_nz1=find(sum(M1b,2).>0);
	i_nz2=find(sum(M2b,2).>0);

	i_z1=find(sum(M1b,2).==0);
	i_z2=find(sum(M2b,2).==0);

	Ln1_nz1=get_Laplacian(M1b);
	Ln2_nz2=get_Laplacian(M2b);

	#(a1,b1)=eigs(Ln1_nz1,which=:SM,nev=num_evec);
	#(a2,b2)=eigs(Ln2_nz2,which=:SM,nev=num_evec);
	
	(a1,b1)=eig(full(Ln1_nz1));
	(a2,b2)=eig(full(Ln2_nz2));

	#ev from eigs diff from eig in the 4th decimal place
	#but more importantly, we found ev that are close to 0 in the 2th, 3rd place using
	#eigs. they are not the right one...
	#it seems ev diff. in the 4th decimal place..

	ord1=sortperm(a1)[1:num_evec];
	b1=b1[:,ord1];
	ord2=sortperm(a2)[1:num_evec];
	b2=b2[:,ord2];	

	b1_extend=zeros(size(M1b,1),num_evec);
	for i=1:num_evec
		b1_extend[i_nz1,i]=b1[:,i];
		x=b1_extend[:,i];
		x1=[x[2:end]' x[end]]';
		x2=[x[1] x[1:end-1]']';
		xx=(x+x1+x2)/3;
		b1_extend[i_z1,i]=xx[i_z1];		
	end

	b2_extend=zeros(size(M2b,1),num_evec);
	for i=1:num_evec
		b2_extend[i_nz2,i]=b2[:,i];
		x=b2_extend[:,i];
		x1=[x[2:end]' x[end]]';
		x2=[x[1] x[1:end-1]']';
		xx=(x+x1+x2)/3;
		b2_extend[i_z2,i]=xx[i_z2];
	end

	evd=zeros(num_evec);
	for i=1:num_evec;
		evd[i]=evec_distance(b1_extend[:,i],b2_extend[:,i]);
	end

	evs=abs(sqrt(2)-evd)/sqrt(2);

	return evd,evs,a1,a2;

end


function evec_distance(x,y);
	#as x and y are normalized in the first place, sqrt(d) makes sense, no need to scale with n
	d1=sum((x-y).^2);
	d2=sum((x+y).^2);
	if d1<d2
		d=d1;
	else 
		d=d2;
	end
	return sqrt(d);
end

function evec_similarity(x,y)
	d=evec_distance(x,y);
	max_d=sqrt(2);
	#this is verified by simulation up to certain accuracy..not proved yet
	s=abs(max_d-d)/max_d;

end
#it's very easy to transform evec_distance to evec_similarity,
#since early reprducib



