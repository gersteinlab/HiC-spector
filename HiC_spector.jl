using HDF5; 
using JLD;
using MAT;
#using Graphs;
using DataFrames;
using CurveFit;
using Interpolations;
#using Loess;
#using Distributions;:

#do power-iteration 
#for compartment, is ICED better in genera?
#can do we a few ipr analysis for leading ev..

function read_contact_map(input_file,chr_num,bin_size);
	#input is a 5 col. file with chr, pos, and contacts	
	chr_length=define_chr_size();
	X=readtable(input_file,separator='\t',header=false);
	chr2bins,bin2loc=generate_arbitrary_mapping_files(bin_size);
	ib=find(bin2loc[1,:].==chr_num-1);
	N=length(ib);
	chr_string=change_chr(chr_num);
	M=sparse(floor(Int64,X[iz,2]/bin_size)+1,floor(Int64,X[iz,4]/bin_size)+1,X[iz,5],N,N);
	M=M+0;

	return M;
end


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

function get_ipr(evec);
	ipr=1./sum(evec.^4,1);
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

	return evs,a1,a2;

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

	return s;

end
#it's very easy to transform evec_distance to evec_similarity

function knight_ruiz(M);
#adapted from the MATLAB code implemented in Knight and Ruiz, 
	M[isnan(M)]=0;
	L=size(M,1);
	iz=find(sum(M,2).>0);
	A=M[iz,iz];
	n=size(A,1); 
	e = ones(n,1); 
	res=[];
	delta = 0.1;
	x0 = e;
	tol = 1e-6;
	g=0.9; etamax = 0.1; # Parameters used in inner stopping criterion.
	
	eta = etamax;
	x = x0; rt = tol^2; v = x.*(A*x); rk = 1 - v;
	rho_km1=sum(rk.^2);
	rout = rho_km1; rold = rout;
	MVP = 0; # count matrix vector products.
	i = 0; # Outer iteration count.

	while rout > rt # Outer iteration
    	i = i + 1; k = 0; y = e;
    	innertol = maximum([eta^2*rout;rt]);
    	while rho_km1 > innertol #Inner iteration by CG
        	k = k + 1;
        	if k == 1
            	Z = rk./v; p=Z; rho_km1 = sum(rk.*Z);
        	else
            	beta=rho_km1/rho_km2;
            	p=Z + beta*p;
        	end
        	# Update search direction efficiently.
        	w = x.*(A*(x.*p)) + v.*p;
        	#w=squeeze(w,2);
        	alpha = rho_km1/sum(p.*w);
        	ap =squeeze(alpha*p,2);
        	# Test distance to boundary of cone.
        	ynew = y + ap;
        	if minimum(ynew) <= delta
            	if delta == 0
            		break
            	end
            	ind = find(ap .< 0);
            	gamma = minimum((delta - y[ind])./ap[ind]);
            	y = y + gamma*ap;
            	break
        	end
        	y = ynew;
        	rk = rk - alpha*w; rho_km2 = rho_km1; rho_km2=rho_km2[1];
        	Z = rk./v; rho_km1 = sum(rk.*Z);
    	end
    	x = x.*y; v = x.*(A*x);
    	rk = 1 - v; rho_km1 = sum(rk.*rk); rout = rho_km1;
    	MVP = MVP + k + 1;
    	# Update inner iteration stopping criterion.
    	rat = rout/rold; rold = rout; r_norm = sqrt(rout);
    	eta_o = eta; eta = g*rat;
    	if g*eta_o^2 > 0.1
        	eta = maximum([eta;g*eta_o^2]);
    	end
    	eta = maximum([minimum([eta;etamax]);0.5*tol/r_norm]);
    	#@sprintf("%3d %6d %.3e %.3e %.3e \n", i,k,r_norm,minimum(y),minimum(x));
        display(rout);
        #res=[res; r_norm];
	end
	#@printf("Matrix-vector products = %6d\n", MVP);
	x=squeeze(x,2);
	A2=A*diagm(x);
	A2=diagm(x)*A2;
	A_balance=extend_mat(A2,iz,L);
	A_balance=(A_balance+A_balance')/2;

	return x,A_balance;

end

function extend_mat(Z,iz,L);
    (u,v)=ind2sub(size(Z),find(Z.!=0));
    w=Z[find(Z)];
    #w=nonzeros(Z);
    u=iz[u];
    v=iz[v];
    Z_extend=sparse(u,v,w,L,L);
    Z_extend=full(Z_extend);
    return Z_extend;
end

#for simplicity, we provide the distance dependence for an individual chromsome, not the genome-wide version
function get_expect_vs_d_single_chr(W,chr2bins,bin_size);

	#unlike the previous version, we include all the pairs with 0 contact (but not dark) as data points..
	#as there are 0 we can't take log, smoothing ae done in linear scale.

	W=full(W);
	W[isnan(W)]=0;
	dark_bins=find(sum(W,1).==0);

	N=size(W,1);
	mmm=minimum(W[W.>0])/10;
	
	(u,v,w)=findnz(triu(W+mmm));
	dark_u=[x in dark_bins for x in u];
	dark_v=[x in dark_bins for x in v];

	iz=find((1-dark_u).*(1-dark_v).>0);

	u=u[iz];
	v=v[iz];
	w=w[iz]-mmm;

	d=float(v-u);
	d2=float(d);
	d2[d2.==0]=1/3;#this is the average distance for 2 points drawn from an uniform distribution between [0.1];
	d3=d2*bin_size;


	#x=log10(d3);
	#y=log10(w+mmm);

	xs,ys_smooth=local_smnoothing(d2,w);
	xs_aux=[round(i) for i in xs];
	xs_aux[1]=1/3;

	xs_all=collect(0:1.0:size(W,1)-1);xs_all[1]=1/3;
	
	ys_all=zeros(size(xs_all));
	for k=1:length(xs_all);
		ik=find(xs_aux.==xs_all[k]);
		if ~isempty(ik)
			ys_all[k]=ys_smooth[ik][1];
		end
	end	

	A_x=find(ys_all.>0);
	knots=(A_x,);
	itp=interpolate(knots,ys_smooth, Gridded(Linear()));

	A_nz=find(ys_all.==0);
	for i=1:length(A_nz);
		ys_all[A_nz[i]]=itp[A_nz[i]];
	end

	expect=ys_all;

	return xs_all*bin_size, expect;

end

function local_smnoothing(x,y);
	
	span=0.01;
	v=sortperm(x);
	x=x[v];
	y=y[v];
	ux=unique(x);
	uy_smooth=zeros(size(ux));
	n=Int(floor(length(x)*span/2));

	mm=zeros(size(x));
	L=2*n+1;
	i=n+1;
	st=1;
	ed=i+n;
	mm[i]=mean(y[st:ed]);
	for i=n+2:length(y)-n;
		#display(i);
    	ed=ed+1;
    	mm[i]=mm[i-1]+y[ed]/L-y[st]/L;
	    st=st+1;
	end
	for i=1:n
    	mm[i]=mean(y[1:n+i]);
	end
	for i=1:n;
    	mm[end-n+i]=mean(y[end-n+1-n+i:end]);
	end

	for i=1:length(ux);
    	iz=find(x.==ux[i]);
    	uy_smooth[i]=mean(mm[iz]);
	end
   
	return ux,uy_smooth;

end

# this is an early version. it doesn't include all zeros...
# but it seems that it's more conssitent with others do, like the compartment analysis in Stein2015
# fitting is also a bit better..
function get_expect_vs_d_single_chr_v0(W,chr2bins,bin_size);

	W=full(W);
	W[isnan(W)]=0;
	
	N=size(W,1);
	

	(u,v,w)=findnz(triu(W));
	d=float(v-u);
	d2=float(d);
	d2[d2.==0]=1/3;#this is the average distance for 2 points drawn from an uniform distribution between [0.1];
	d3=d2*bin_size;

	#model = loess(log10(d3),log10(w),span=0.01);
	#the loess fct is rather slow, and fail to work at some matrices (not sure why), we have replaced it by a simpler method

	x=log10(d3);
	y=log10(w);

	xs,ys_smooth=local_smnoothing(x,y);

	xs_all=collect(0:1.0:size(W,1)-1);xs_all[1]=1/3;
	xs_all=xs_all*bin_size;
	xs_all_aux=log10(xs_all);

	ys_all=zeros(size(xs_all));
	for k=1:length(xs_all_aux);
		ik=find(xs.==xs_all_aux[k]);
		if ~isempty(ik)
			ys_all[k]=ys_smooth[ik][1];
		end
	end	

	A_x=find(ys_all.>0);
	knots=(A_x,);
	itp=interpolate(knots,ys_smooth, Gridded(Linear()));

	A_nz=find(ys_all.==0);
	for i=1:length(A_nz);
		ys_all[A_nz[i]]=itp[A_nz[i]];
	end

	expect=10.^ys_all;

	return xs_all, expect;

end

function get_expect_vs_d_WG_v0(contact,chr2bins,bin_size);

	all_d2=Float64[];
	all_w=Float64[];
	Ltmp=zeros(23);
	for chr_num=1:23
	
		#display(chr_num);
		W=extract_chr(contact,chr2bins,chr_num);
		W=full(W);
		W[isnan(W)]=0;

		N=size(W,1);
		
		(u,v,w)=findnz(triu(W));
		
		d=float(v-u);
		d2=float(d);
		d2[d2.==0]=1/3;#this is the average distance for 2 points drawn from an uniform distribution between [0.1];
		
		all_d2=[all_d2;d2];
		all_w=[all_w;w];
		Ltmp[chr_num]=size(W,1);
	
	end

	all_d3=all_d2*bin_size;

	x=log10(all_d3);
	y=log10(all_w);

	xs,ys_smooth=local_smnoothing(x,y);

	xs_all=collect(0:1.0:maximum(Ltmp)-1);xs_all[1]=1/3;
	xs_all=xs_all*bin_size;
	xs_all_aux=log10(xs_all);

	ys_all=zeros(size(xs_all));
	for k=1:length(xs_all_aux);
		ik=find(xs.==xs_all_aux[k]);
		if ~isempty(ik)
			ys_all[k]=ys_smooth[ik][1];
		end
	end	

	A_x=find(ys_all.>0);
	knots=(A_x,);
	itp=interpolate(knots,ys_smooth, Gridded(Linear()));

	A_nz=find(ys_all.==0);
	for i=1:length(A_nz);
		ys_all[A_nz[i]]=itp[A_nz[i]];
	end

	expect=10.^ys_all;

	return xs_all, expect;

end




function get_expect_vs_d_WG(contact,chr2bins,bin_size);

	all_d2=Float64[];
	all_w=Float64[];
	Ltmp=zeros(23);
	for chr_num=1:23
	
		#display(chr_num);
		W=extract_chr(contact,chr2bins,chr_num);
		W=full(W);
		W[isnan(W)]=0;

		dark_bins=find(sum(W,1).==0);

		N=size(W,1);
		mmm=minimum(W[W.>0])/10;
	
		(u,v,w)=findnz(triu(W+mmm));
		dark_u=[x in dark_bins for x in u];
		dark_v=[x in dark_bins for x in v];

		iz=find((1-dark_u).*(1-dark_v).>0);

		u=u[iz];
		v=v[iz];
		w=w[iz]-mmm;

		d=float(v-u);
		d2=float(d);
		d2[d2.==0]=1/3;#this is the average distance for 2 points drawn from an uniform distribution between [0.1];
		d3=d2*bin_size;

		all_d2=[all_d2;d2];
		all_w=[all_w;w];
		Ltmp[chr_num]=size(W,1);
	
	end

	all_d3=all_d2*bin_size;

	xs,ys_smooth=local_smnoothing(all_d2,all_w);

	xs_aux=[round(i) for i in xs];
	xs_aux[1]=1/3;

	xs_all=collect(0:1.0:size(W,1)-1);xs_all[1]=1/3;
	
	ys_all=zeros(size(xs_all));
	for k=1:length(xs_all);
		ik=find(xs_aux.==xs_all[k]);
		if ~isempty(ik)
			ys_all[k]=ys_smooth[ik][1];
		end
	end	

	A_x=find(ys_all.>0);
	knots=(A_x,);
	itp=interpolate(knots,ys_smooth, Gridded(Linear()));

	A_nz=find(ys_all.==0);
	for i=1:length(A_nz);
		ys_all[A_nz[i]]=itp[A_nz[i]];
	end

	expect=ys_all;

	return xs_all*bin_size, expect;

end

#to find distance dependence, we should NOT iced the chr one by one.
#because in genome-wide scale dependance, we should keep the contacts in same base
#for individual genome dependence, may be ice again ot not both ok, but the difference is maginal, cor=0.99

function get_f_W(W,ys);

	N=size(W,1);
	W[isnan(W)]=0;
	dark_bins=find(sum(W,1).==0);
	num_dark=length(dark_bins);
	N_eff=N-num_dark;
	f_W=zeros(size(W));#what's f_W? it's a generation of ones(size(W));

	x=collect(1:N);

	for d=0:N-1
		f_W[1+d:N+1:end-d*N]=ys[d+1];
	end
	tmp=f_W-diagm(diag(f_W));
	f_W=f_W+tmp';
	#sum(f_W[1,:])=1 here..

	f_W[dark_bins,:]=0;
	f_W[:,dark_bins]=0;
	f_W=f_W/sum(f_W)*N_eff.^2;

	return f_W;

end


function get_compartment_A_B(W,f_W);

	iz=find(sum(W,2).>0);
	izz=find(sum(W,2).==0);
	Wn=W[iz,iz]./f_W[iz,iz];
	C=cor(Wn);
	(U,V)=eigs(C);
	i_max=indmax(U);
	ev=V[:,i_max:i_max+5];
	ev_whole=zeros(size(W,1),6);
	ev_whole[iz,:]=ev;
	ev_whole[izz,:]=NaN;

	(loc,span)=get_chunks_v2(sign(ev_whole[:,1]),1);#
	cpt=sign(ev_whole[loc,1]);

	return loc,span,ev_whole,cpt;

end


function generate_arbitrary_mapping_files(bin_size);

	num_of_chromosomes=25;
	chr2bins=zeros(2,num_of_chromosomes);
	chr_length=define_chr_size();
	chr_num_bins=round(Int64,floor(chr_length/bin_size))+1
	#chr_num_bins=int(floor(chr_length/bin_size))+1;
	chr2bins[2,:]=cumsum(chr_num_bins)'-1;
	chr2bins[1,1]=0;
	chr2bins[1,2:end]=chr2bins[2,1:end-1]+1;
	X=round(Int,chr2bins+1);
	bin2loc=zeros(3,X[2,end]);
	for c=1:25
		bin2loc[1,X[1,c]:X[2,c]]=c-1;
		bin2loc[2,X[1,c]:X[2,c]]=round(Int,collect(1:bin_size:chr_length[c]))';
		bin2loc[3,X[1,c]:X[2,c]]=[round(Int,collect(bin_size:bin_size:chr_length[c]))' chr_length[c]];
	end
	return round(Int64,chr2bins),round(Int64,bin2loc);
	
end


#id is the starting loc of a chunk, and d is the length it spans..
function get_chunks_v2(a,singleton=0);
	# adopt from a matlab code by Jiro Doke;
	 a                 = [NaN; a; NaN];
	 b                 = diff(a);
	 b1                = b;  # to be used in fullList (below)
	 ii                = trues(size(b));
	 ii[b.==0] = false;
	 b[ii]             = 1;
	 c                 = diff(b);
	 id                = find(c.==-1);

	 #Get single-element chunks also
	 if singleton.==1
	 	b1[id]          = 0;
	 	ii2             = find(b1[1:end-1]);
	 	d               = vcat(find(c.==1) - id + 1, ones(length(ii2)));
	 	id              = [id;ii2];
	 	v=sortperm(id);
	 	id=sort(id);
	 	#(id,tmp)        = sort(id);
	 	d               = d[v];
	 else 
	 	d               = find(c.==1) - id + 1;
	 end

	 return id,d;
end



function define_chr_size();
	chr_length=zeros(Int,25,1);
	chr_length[1]=249250621;
	chr_length[2]=243199373;
	chr_length[3]=198022430;
	chr_length[4]=191154276;
	chr_length[5]=180915260;
	chr_length[6]=171115067;
	chr_length[7]=159138663;
	chr_length[8]=146364022;
	chr_length[9]=141213431;
	chr_length[10]=135534747;
	chr_length[11]=135006516;
	chr_length[12]=133851895;
	chr_length[13]=115169878;
	chr_length[14]=107349540;
	chr_length[15]=102531392;
	chr_length[16]=90354753;
	chr_length[17]=81195210;
	chr_length[18]=78077248;
	chr_length[19]=59128983;
	chr_length[20]=63025520;
	chr_length[21]=48129895;
	chr_length[22]=51304566;
	chr_length[23]=155270560;#X
	chr_length[24]=59373566;#Y
	chr_length[25]=16571;
	return chr_length;
end

function change_chr(chr)
	if typeof(chr)==Float64||typeof(chr)==Int64;
			if chr==25;
				chr2="chrM";
			elseif chr==24;
				chr2="chrY";
			elseif chr==23;
				chr2="chrX";
			elseif chr==22;
				chr2="chr22";
			elseif chr==21;
				chr2="chr21";
			elseif chr==20;
				chr2="chr20";
			elseif chr==19;
				chr2="chr19";
			elseif chr==18;
				chr2="chr18";
			elseif chr==17;
				chr2="chr17";
			elseif chr==16;
				chr2="chr16";
			elseif chr==15;
				chr2="chr15";
			elseif chr==14;
				chr2="chr14";
			elseif chr==13;
				chr2="chr13";
			elseif chr==12;
				chr2="chr12";
			elseif chr==11;
				chr2="chr11";
			elseif chr==10;
				chr2="chr10";
			elseif chr==9;
				chr2="chr9";
			elseif chr==8;
				chr2="chr8";
			elseif chr==7;
				chr2="chr7";
			elseif chr==6;
				chr2="chr6";
			elseif chr==5;
				chr2="chr5";
			elseif chr==4;
				chr2="chr4";
			elseif chr==3;
				chr2="chr3";
			elseif chr==2;
				chr2="chr2";
			elseif chr==1;
				chr2="chr1";
			end

	elseif typeof(chr)==ASCIIString||typeof(chr)==SubString{ASCIIString}||typeof(chr)==UTF8String
			if chr=="chrM";
				chr2=25;
			elseif chr=="chrY";
				chr2=24;
			elseif chr=="chrX";
				chr2=23;
			elseif chr=="chr22";
				chr2=22;
			elseif chr=="chr21";
				chr2=21;
			elseif chr=="chr20";
				chr2=20;
			elseif chr=="chr19";
				chr2=19;
			elseif chr=="chr18";
				chr2=18;
			elseif chr=="chr17";
				chr2=17;
			elseif chr=="chr16";
				chr2=16;
			elseif chr=="chr15";
				chr2=15;
			elseif chr=="chr14";
				chr2=14;
			elseif chr=="chr13";
				chr2=13;
			elseif chr=="chr12";
				chr2=12;
			elseif chr=="chr11";
				chr2=11;
			elseif chr=="chr10";
				chr2=10;
			elseif chr=="chr9";
				chr2=9;
			elseif chr=="chr8";
				chr2=8;
			elseif chr=="chr7";
				chr2=7;
			elseif chr=="chr6";
				chr2=6;
			elseif chr=="chr5";
				chr2=5;
			elseif chr=="chr4";
				chr2=4;
			elseif chr=="chr3";
				chr2=3;
			elseif chr=="chr2";
				chr2=2;
			elseif chr=="chr1";
				chr2=1;
			end
	end
	return chr2;
end


