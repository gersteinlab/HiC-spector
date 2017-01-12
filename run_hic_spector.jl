include("/home/fas/gerstein/ky26/Github/HiC_spectra/HiC_spector.jl");
using PyPlot;

#####Calculating reproducibility scores using A549 data

hg19_info=define_hg19_genome();
bin_size=40000;

chr2bins,bin2loc=generate_arbitrary_mapping_files(hg19_info,bin_size);

map_file1="./data/A549/A549C-HindIII-R1_int.bed";
map_file2="./data/A549/A549D-HindIII-R2_int.bed";

num_evec=30;
r=20;

X1=readtable(map_file1,separator='\t',header=false);
X2=readtable(map_file2,separator='\t',header=false);

Q=zeros(23);

for chr_num=1:23;
	display(chr_num);

	ib=find(bin2loc[1,:].==chr_num-1);
	N=length(ib);

	chr_string=change_chr(hg19_info,chr_num);
		
	iz1=find(X1[:,1].==chr_string);
	M1=sparse(floor(Int64,X1[iz1,2]/bin_size)+1,floor(Int64,X1[iz1,4]/bin_size)+1,X1[iz1,5],N,N);
	iz2=find(X2[:,1].==chr_string);
	M2=sparse(floor(Int64,X2[iz2,2]/bin_size)+1,floor(Int64,X2[iz2,4]/bin_size)+1,X2[iz2,5],N,N);

	M1=full(M1);
	M2=full(M2);
	
	evs,a1,a2=get_reproducibility(M1,M2,num_evec);

	Q[chr_num]=mean(evs[1:r]);
end


#####Other analysis using MCF7 data

hg19_info=define_hg19_genome();
bin_size=250000;
chr2bins,bin2loc=generate_arbitrary_mapping_files(hg19_info,bin_size);

data_loc="./data/HiCStein-MCF7-WT.maps/";

all_gamma=zeros(23);
all_K=zeros(23);
all_w=[];
all_d2=[];

for chr_num=1:23
	display(chr_num);
	#chr_num=10;
	input_file=data_loc*"HiCStein-MCF7-WT__hg19__genome__C-250000-iced__"*change_chr(chr_num)*"__"*change_chr(chr_num)*"__cis.matrix";
	X=readtable(input_file,header=true,separator='\t');
	W=X[:,2:end];
	W=array(W);
	W[isnan(W)]=0;
	#d2,w, gamma, K=get_expect_vs_d_single_chr(W,chr2bins,bin_size);
	xs_all, expect=get_expect_vs_d_single_chr_v0(W,chr2bins,bin_size);

	all_d2=[all_d2;d2];
	all_w=[all_w;w];
	all_gamma[chr_num]=gamma;
	all_K[chr_num]=K;
end

all_d2=convert(Array{Float64,1},all_d2);
all_w=convert(Array{Float64,1},all_w);
tmp=linear_fit(log10(all_d2*250000),log10(all_w));
gamma=tmp[2];
K=10^tmp[1];

model = loess(log10(all_d2*250000),log10(all_w));
us=log10(all_d2*2500000);
vs=predict(model,us);



x=collect(4:.5:8.5);x=10.^x;
y=K*x.^gamma;
PyPlot.plot(log10(x),log10(y));

f_W=get_f_W(W,K,gamma);

compartment_file="/home/fas/gerstein/ky26/scratch/Hi-C_data/Stein_GB2015/Hi-C_MCF7_MCF10A_processed_HiCfiles/Compartments/HiCStein-MCF7-WT__hg19__genome__C-250000-iced__chr1__chr1__cis.matrix.compartments";

df=readtable(compartment_file,header=true,separator='\t');

gene_density=df[:geneDensity];
ev1=df[:eigen1];

loc,span,ev_whole,cpt=get_compartment_A_B(W,f_W);
my_ev=ev_whole[:,1];

y1=all_K[1]*x.^all_gamma[1];
f_W1=get_f_W(W,all_K[1],all_gamma[1]);

loc1,span1,ev_whole1,cpt1=get_compartment_A_B(W,f_W1);


