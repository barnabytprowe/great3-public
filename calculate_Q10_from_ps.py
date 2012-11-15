%    Code to generate Q values for power spectrum submissions to GREAT10
%    Copyright (C) 2011 Thomas Kitching
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.                                                                                                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Copyright T. D. Kitching 2011
%
%  Code to read in great10 variable shear catalogues and estimate Q
%  Adapted from analysis_shears (written August 2010 to December 2010) to read in power spectrum submissions 
%
%  * Written July 2011 to August 2011 
%  * Comments added (Sep 11)
%  * Bugs found when commenting :
%  - dll=l(j+1)-l(j)->log(l(j+1))-log(j) and log(dll)->dll (in sigma_sys)
%  - sigma_sys writting bug (line 299) values were factor 5 too small on leaderboard
%
%  * Notes :
%  - 2e-4 in tk talk on 26/9 used old plot should have been 5e-6
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LASTN = maxNumCompThreads('automatic'); %parallel matlab

lstr=int2str(1);
outfilename=strcat('/disk1/tdk/code/g10res/outs/q10.galaxy.tiers.',lstr,'.dat');
fout = fopen(outfilename,'a');
%___________________________ start _______________________________________%

%---------------------image dimensions -----------------------------------%

% harwired image numbers
g10postage=48;  %postage stamp size   
n      =100;    % number of galaxies on a side  
n_pix  =n.*n;   

% Survey parameters.
sky_size_deg=10.;              % effective sky size in degrees
sky_area=sky_size_deg.^2;      % Sky area.
nbin2  = 20;                   
logl_flag=0;                   % 1=log, 0=linear ***only linear implemented log implementation not complete

%Q factor numerator 
%see coordination wiki page http://great10coordination.pbworks.com/w/page/32500368/121110 for discussion and decision
QN=5e-6;      %from figure(99) in make_varshears.m
QN=QN.*1000.; %so good Q=1000

%define the l range  (indeces)
imin=2;
imax=9;

%--------------testing
topname='g10_gamma1';
pubname='g10_v1';
postfix='.dat';
nstartreal=1;
numnreal=200;
npsreal=1;
ncols=4;
ntiers=26;

%_________________________________________________________________________%

allQ=0.;

Qt=zeros(ntiers);

for tiernum=1:ntiers
    
sigmasys=0.;

%read in the file that maps the participant number to the g10 key 

%number in each tier 
ntot=npsreal*numnreal;
mapname=strcat('truth/',topname,'_t',int2str(tiernum),'.map');
fid = fopen(mapname,'r');
head = fgetl(fid); %get the header
[tr,count]=fscanf(fid,'%d %d %d %d',[4,ntot]);

avq=0.;

Q_10final=0.;

gtr=zeros(n,n);  
ein=zeros(n,n);  

ncount=0;

    % harwired image numbers
    g10postage=48;  %postage stamp size   
    n      =100;    % number of galaxies on a side  
    n_pix  =n.*n;   

    powernum=1;
nrealisation=1; %only one submissions per set (all shear power the same for each set)
    %read the true g10 data file ---------------------------------------------%
    truthnamein =strcat('truth/',topname,'_t',int2str(tiernum),'_',int2str(powernum),'_image_galaxy_',int2str(nrealisation),'.dat');

    fid = fopen(truthnamein,'r');

    head = fgetl(fid); %get the header
    [r,count]=fscanf(fid,'%29f',[29,10000]);

    g10postage=r(1,1);   
    n      =r(29,1)/g10postage;    
    n_pix  =n.*n;   
    ebulgetr=zeros(n,n);
    edisktr=zeros(n,n);
    psftr=zeros(n,n);
    x=r(27,:)+g10postage/2.;
    y=r(28,:)+g10postage/2.;
    for i=1:int32(n_pix)

        %note: -g2 below because of the great08 -g2 convention 
        %the make shear creates a +g2 power which is then written as -g2 for
        %the image code. This code then needs to re-convert to a -g2 for the
        %analysis (most shape codes will measure using the +g2 convention)
        gtr(int32(x(i)/g10postage+0.5),int32(y(i)/g10postage+0.5))=gtr(int32(x(i)/g10postage+0.5),int32(y(i)/g10postage+0.5))+complex(r(9,i),-r(10,i)); %shear field

    end
    fclose(fid);
     
    %------------------------------read the files in--------------------------%
    
    methodnamein=strcat('/disk2/tdk/code/g10res/galssub/1/dkirkby_set',int2str(tiernum),'.dat');
    
    gPowEE=0.;
    fid = fopen(methodnamein,'r');
    [r,count]=fscanf(fid,'%d %e %e',[3,8]);
    gPowEE=r(3,:);
    fclose(fid);
    %_________________________________________________________________________%

    %----------------------- set up the survey -------------------------------%

    deg2rad=pi./180.;               % Change degrees to radians.
    size_pix=sky_size_deg./n;       % angular Pixel size in degrees
    size_pix_rad=size_pix.*deg2rad; % angular Pixel size in radians

    % set up l ranges and resolutions
    L=sky_size_deg.*deg2rad;        % field size in rads
    Max_l_mode=2.*pi/size_pix_rad;  % find max l 
    Min_l_mode=2.*pi/L;             % find min l

    if (nbin2 >=n) 
        fprintf(' error nbin2 >= n : %f',nbin2);
        stop
    end 

    lbin=zeros(nbin2,1);
    if (logl_flag==1) 
    %log lbins: 
        dlogl=(log10(Max_l_mode/Min_l_mode)/(double(nbin2)-1.)); %log10 bin width
        lbin=Min_l_mode*10.0.^(dlogl.*(double(1:nbin2)-1))-1.+0.00001;
        lbin=lbin';
    end
    if (logl_flag==0) 
    %lin lbins
        dlogl=(Max_l_mode-Min_l_mode)/(double(nbin2)-1.); 
        lbin=Min_l_mode+(dlogl.*(double(1:nbin2)-1))-1.+0.00001;
        lbin=lbin';
    end

    nbin=int32(2.*n);    %convert 2* pixels to integer

    fsky=1.0;

    % Create a Complex wavevector 
    el1=2.*pi.*((1:n)-1.0-double(n)/2.)./L+1e-10; %need to hadd small constant to avoid ==0 causing NaNs in lrot
    lvec  =zeros(n,n);
    icoord=zeros(n,n);
    jcoord=zeros(n,n);
    for i1=1:n
        l1=el1(i1);
        for j1=1:n
            l2=el1(j1);     
            lvec(i1,j1)=l1+1i.*l2;
            icoord(i1,j1)=i1;
            jcoord(i1,j1)=j1;
        end
    end

    %need to remove any numerically small, but non-zero rounding errors in
    %the l-rotation matrix : (and the conjugate)
    lrot=lvec.*lvec./(lvec.*conj(lvec));
    llft=ifftshift(lrot);
    llft=ifft2(llft);
    lrot=(fft2(complex(real(llft)-mean(mean(real(llft))),imag(llft)-mean(mean(imag(llft)))))); 
   
    clrot=conj(lvec).*conj(lvec)./(lvec.*conj(lvec));
    cllft=ifftshift(clrot);
    cllft=ifft2(cllft);
    clrot=(fft2(complex(real(cllft)-mean(mean(real(cllft))),imag(cllft)-mean(mean(imag(cllft)))))); 
    
    % Estimate E and B modes assuming linear-KS.------------------------------%

    %estimated ellipticities
    ein=fftshift(ein);
    eft=fft2(ein);                    %FT of the ellipticity field
    
    kapi=clrot.*eft;
    kapi=ifftshift(kapi);
    kapi=ifft2(kapi);

    %the true shear values

    gtr=fftshift(gtr);
    tfieldft=fft2(gtr);                    %FT of the shear field
    tkapi=clrot.*tfieldft; %rotate the shear field to a kappa and beta
    tkapi=ifftshift(tkapi);                   %inverse FT back to real space 
    tkapi=ifft2(tkapi);

    %make power spectra

    tkapi=fftshift(tkapi);
    tkapft=fft2(real(tkapi));
    tbetft=fft2(imag(tkapi));
    
    tCEE_2  = real(tkapft).^2+imag(tkapft).^2;  %E mode power 
    tCBB_2  = real(tbetft).^2+imag(tbetft).^2;  %B mode power
    tCEB_2  = real(tkapft)*real(tbetft)-imag(tkapft)*imag(tbetft); %EB cross power

    % Angle average of power spectra in log10 bins.
    tPowEE     = zeros(nbin2,1); %allocate some measured power spectrum arrays
    tPowBB     = zeros(nbin2,1);
    tPowEB     = zeros(nbin2,1);
    ll         = zeros(nbin2,1);
    dll        = zeros(nbin2,1);                                                                                                                                                             
    %integrate over ell space
    for i1=1:n
        l1=el1(i1);
        for j1=1:n
            l2=el1(j1);

            l=sqrt(l1.^2+l2.^2);

            if (l<=Max_l_mode && l>=Min_l_mode) 
            if (logl_flag==1) 
                ibin=int32(log10((l+1.)./Min_l_mode)./dlogl)+1;
                tPowEE(ibin)  = tPowEE(ibin)+tCEE_2(i1,j1)/log(10.);
                tPowBB(ibin)  = tPowBB(ibin)+tCBB_2(i1,j1)/log(10.);
                tPowEB(ibin)  = tPowEB(ibin)+tCEB_2(i1,j1)/log(10.);
            end
            if (logl_flag==0) 
                ibin=int32((l+1-Min_l_mode)./dlogl)+1; 
                tPowEE(ibin)  = tPowEE(ibin)+tCEE_2(i1,j1)*l;
                tPowBB(ibin)  = tPowBB(ibin)+tCBB_2(i1,j1)*l;
                tPowEB(ibin)  = tPowEB(ibin)+tCEB_2(i1,j1)*l;     
            end
            ll(ibin)=l;
            
            if ibin > 1
                dll(ibin-1)=log(ll(ibin))-log(ll(ibin-1));
            end
            
            end
        end 
    end
    %normalise the integrals with the ML fft normalisation
    tPowEE  =  (tPowEE./(n.^4.*dlogl.*fsky));
    tPowBB  =  (tPowBB./(n.^4.*dlogl.*fsky));
    tPowEB  =  (tPowEB./(n.^4.*dlogl.*fsky));
    
    %_________________________________________________________________________%

    for j=imin:imax                                                                                                         
	    fprintf(' %f %f %f\n',ll(j),tPowEE(j,1),gPowEE(j-imin+1));                                                      
            Csys(1,1,j)=tPowEE(j)-gPowEE(j-imin+1);                                                                     
    end
    fprintf('\n\n');

    sigmasys=0.;
    for j=imin:imax
	sigmasys=sigmasys+(abs(Csys(1,1,j))).*dll(j); %has to be the same way that make_varshears calculates the QN
    end
    q10=QN./sigmasys;
    fprintf(' %d Q_10=%f %f\n',tiernum,q10,sigmasys); %Q of the mean

    Q_10final=q10;

    avq=avq+Q_10final;
    allQ=allQ+avq;
    %stop
    
    Qt(tiernum)=Q_10final;
        
end
fclose(fout);

allQ=allQ/ntiers;

fprintf(' Q=%f\n',allQ);

%print to a file
outfilename=strcat('q10.galaxy.power.lis');
fout = fopen(outfilename,'a');
%fprintf(fout,' %f %f\n',0.001/allQ,allQ);
fprintf(fout,' %f %f\n',QN./allQ,allQ);
fclose(fout);
%_________________________________________________________________________%
