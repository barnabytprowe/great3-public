"""Code to generate Q values for ellipticity catalogue submissions to GREAT10.

Based on (and ported from) analysis_shears_commented.m by Tom Kitching.
"""

#print Q per set files to a file
lstr = int2str(leadernum);
outfilename=strcat(['/Users/browe/g10_v1/qcode/results/' ...
                    'q10.galaxy.tiers.'],lstr,'.', methodtype, '.dat');
fout = fopen(outfilename, 'a');

%___________________________ start _______________________________________%

%---------------------image dimensions -----------------------------------%

% harwired image numbers
g10postage=48;     %postage stamp size (overwritten by read from truth)    
n         =100;    % number of galaxies on a side n_pix  =n.*n;   

% Survey parameters.
zs=1.0;                        % effective median depth of survey
sky_size_deg=10.;              % effective sky size in degrees
sky_area=sky_size_deg.^2;      % Sky area.
nbin2=20;                      % >= 20 is good for fine runs. number of ell bins in power
logl_flag=0;                   % 1=log, 0=linear ***only linear implemented log implementation not completed

%Q factor numerator 
%see coordination wiki page http://great10coordination.pbworks.com/w/page/32500368/121110 for discussion and decision
QN=5e-6;      %from figure(99) in make_varshears.m 
QN=QN.*1000.; %so good Q=1000

%define the l range  (indeces)
imin=2;
imax=9;

%-------------- set up method namesm number of tiers (sets) and number of images per set ----------%
topname='g10_gamma1';
pubname='g10_v1';
%methodtype='ksb'; %input from shell script (default is first entry = ksb)
postfix='.dat';
nstartreal=1;
numnreal=200; %number of iges per set 
npsreal=1; %dummy variable 
ncols=4; %number of expected columns in the file
ntiers=26; %number of tiers
fprintf(' leader=%d; version=%s; method=%s; number of images=%d \n',leadernum,topname,methodtype,numnreal);

%_________________________________________________________________________%

allQ=0.; %initialise total Q value
for tiernum=1:ntiers
    
sigmasys=0.; %initialise 

%read in the file that maps the participant number to the g10 key 

%number in each tier 
ntot=npsreal*numnreal;
mapname=strcat('truth/',topname,'_t',int2str(tiernum),'.map');
fid = fopen(mapname,'r');
head = fgetl(fid); %get the header
[tr,count]=fscanf(fid,'%d %d %d %d',[4,ntot]);
fclose(fid);

avq=0.; %average Q
for powernum=1:npsreal 
%redundant loop (npsreal=1) from older version of code, kept in in case flexibility required later

g2checkbest=0;
Q10_final=0.;

%checking catalogue formats (e1,e2 flips and x-y flips)
%g2checkmin=2; %1
%g2checkmax=2; %8

%checking catalogue formats (e1,e2 flips and x-y flips)
g2checkmin=1;                                                                
g2checkmax=8;


for g2check=g2checkmin:g2checkmax
    
    meantPowEE=zeros(nbin2,1); %initialise mean true E-mode power 
    meanPowEE=zeros(nbin2,1);  %initialise mean submitted E-mode power

    gtr=zeros(n,n);  %initlalise shear field
    ein=zeros(n,n);  %initialise submitted e field

    ncount=0;
for nrealisation=nstartreal:numnreal %number of images in a set (noise realisations)
%in this loop all the true shear and submitted ellipticities are read and avergaed over all images in a set 

    %find the truth key
    for t1=1:ntot
    if (tr(2,t1)==tiernum && tr(3,t1)==nrealisation && tr(4,t1)==powernum) 
        key=tr(1,t1);
    end
    end
    
    %read the true g10 data file ---------------------------------------------%
    truthnamein =strcat('truth/',topname,'_t',int2str(tiernum),'_',int2str(powernum),'_image_galaxy_',int2str(nrealisation),'.dat'); %construct the name of the truth table

    fid = fopen(truthnamein,'r'); %open truth table

    head = fgetl(fid); %get the header and ignore it 
    [r,count]=fscanf(fid,'%29e',[29,10000]); %read all columns 
    fclose(fid);

    g10postage=r(1,1); %size of postage stamps
    n         =r(29,1)/g10postage; %number of galaxies on a side of the image    
    n_pix     =n.*n;   
    ebulgetr=zeros(n,n); %bulge intrinsic
    edisktr=zeros(n,n);  %disk intrinsic
    x=r(27,:)+g10postage/2.; %positions of the galaxies note x(27 y(28 are bottom left corner of postage
    y=r(28,:)+g10postage/2.;
    for i=1:int32(n_pix)

        %note: -g2 below because of the great08 -g2 convention 
        %the make shear creates a +g2 power which is then written as -g2 for
        %the image code. This code then needs to re-convert to a -g2 for the
	%analysis (most shape codes will measure using the +g2 convention); redundant anyway because of g2checking loop
        gtr(int32(x(i)/g10postage+0.5),int32(y(i)/g10postage+0.5))=gtr(int32(x(i)/g10postage+0.5),int32(y(i)/g10postage+0.5))+complex(r(9,i),-r(10,i)); %shear field

	%bluge and disk (dont actually need this but here for completeness)
	%bulge
        e1btr=((1-r(17,i))./(1+r(17,i))).*cos(2.*r(18,i).*pi./180.);
        e2btr=((1-r(17,i))./(1+r(17,i))).*sin(2.*r(18,i).*pi./180.);
        ebulgetr(int32(x(i)/g10postage+0.5),int32(y(i)/g10postage+0.5))=complex(e1btr,e2btr);
        %disk
        e1dtr=((1-r(24,i))./(1+r(24,i))).*cos(2.*r(25,i).*pi./180.);
        e2dtr=((1-r(24,i))./(1+r(24,i))).*sin(2.*r(25,i).*pi./180.);
        edisktr(int32(x(i)/g10postage+0.5),int32(y(i)/g10postage+0.5))=complex(e1dtr,e2dtr);
    end
     
    %------------------------------read the submitted files --------------------------%
    
    methodnamein=strcat('./results/',pubname,'_objs_',int2str(tiernum),'_',int2str(key),'.',methodtype,postfix); %construct the name 

    fid = fopen(methodnamein,'r'); %open the file
    format=''; %generate the read format 
    for i=1:ncols 
        format=strcat(format,' %f');
    end
    [t,count]=fscanf(fid,format,[ncols,n_pix]); %read the submitted e values
    fclose(fid);

    if (count~=0) 
      %count==0 would mean a blank submission 

      %check for missing data
    if ((count/ncols)~=n_pix) 
        %fprintf('error not enough objects in the catalogue!! %d %d\n',count/4,n_pix);
        n_pix=count/ncols;
    end
        
    %x and y are the first two columns in units of pixels 
    x=t(1,:);
    y=t(2,:);
    %loop over all galaxies and make the e field (for each convention case) 
    for i=1:int32(n_pix)
	%x-y swaps
        if (g2check<=4) 
            jindex=int32(x(i)/g10postage+0.5);
            iindex=int32(y(i)/g10postage+0.5);
        else 
            jindex=int32(y(i)/g10postage+0.5);
            iindex=int32(x(i)/g10postage+0.5);
        end
        %check coordinates make sense
        if iindex>0 && iindex <=n && jindex>0 && jindex<=n 
	  %e1, e2 swaps
	    if (g2check==1 || g2check==5) 
                ein(iindex,jindex)=ein(iindex,jindex)+complex(t(3,i),t(4,i));
            end
            if (g2check==2 || g2check==6) 
                ein(iindex,jindex)=ein(iindex,jindex)+complex(t(3,i),-t(4,i));
            end
            if (g2check==3 || g2check==7) 
                ein(iindex,jindex)=ein(iindex,jindex)+complex(-t(3,i),-t(4,i));
            end
            if (g2check==4 || g2check==8) 
                ein(iindex,jindex)=ein(iindex,jindex)+complex(-t(3,i),t(4,i));
            end
        end        
    end %end loop over galaxies
    %fclose(fid)
    end %end if count~=0
      ncount=ncount+1; %count the successful images read per set (should be all images)
    end %end loop over images in set 

    %average over images in the set 
    ein=ein/ncount;
    gtr=gtr/ncount;

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

    % now follows power spectrum generation as described in GREAT10 Handbook appendix 
    % Create a Complex wavevector 
    el1=2.*pi.*((1:n)-1.0-double(n)/2.)./L+1e-10; %need to add small constant to avoid ==0 causing NaNs in lrot
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
   
    %conjugate 
    clrot=conj(lvec).*conj(lvec)./(lvec.*conj(lvec));
    cllft=ifftshift(clrot);
    cllft=ifft2(cllft);
    clrot=(fft2(complex(real(cllft)-mean(mean(real(cllft))),imag(cllft)-mean(mean(imag(cllft))))));
    
    %------------------------- Estimate E and B modes assuming linear-KS.------------------------------%

    %estimated ellipticities
    ein=fftshift(ein);
    eft=fft2(ein);                    %FT of the ellipticity field
    kapi=clrot.*eft;  %rotate the shear field to a kappa and beta
    kapi=ifftshift(kapi);  %inverse FT back to real space
    kapi=ifft2(kapi);

    %the true shear values
    gtr=fftshift(gtr);
    tfieldft=fft2(gtr);                    %FT of the shear field
    tkapi=clrot.*tfieldft; %rotate the shear field to a kappa and beta
    tkapi=ifftshift(tkapi);                   %inverse FT back to real space 
    tkapi=ifft2(tkapi);


    %now take kappa field and make power spectrum
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
	      %need l^2. One factror already taken into account by geometry (integration in azimuth with constant |l|)
	      %second factor multiplied below 
                ibin=int32((l+1-Min_l_mode)./dlogl)+1; 
                tPowEE(ibin)  = tPowEE(ibin)+tCEE_2(i1,j1)*l;
                tPowBB(ibin)  = tPowBB(ibin)+tCBB_2(i1,j1)*l;
                tPowEB(ibin)  = tPowEB(ibin)+tCEB_2(i1,j1)*l;     
            end
            ll(ibin)=l;
            end
        end 
    end
    %normalise the integrals with the Matlab fft normalisation
    tPowEE  =  (tPowEE./(n.^4.*dlogl));
    tPowBB  =  (tPowBB./(n.^4.*dlogl));
    tPowEB  =  (tPowEB./(n.^4.*dlogl));
    
    %the input ellipticities-------------------------------------------------%
    kapi=fftshift(kapi);
    gkapft=fft2(real(kapi));
    gbetft=fft2(imag(kapi));

    gCEE_2  = real(gkapft).^2+imag(gkapft).^2;  %E mode power 
    gCBB_2  = real(gbetft).^2+imag(gbetft).^2;  %B mode power
    gCEB_2  = real(gkapft)*real(gbetft)-imag(gkapft)*imag(gbetft); %EB cross power

    % Angle average of power spectra in log10 bins.
    gPowEE     = zeros(nbin2,1); %allocate some measured power spectrum arrays
    gPowBB     = zeros(nbin2,1);
    gPowEB     = zeros(nbin2,1);
    ll         = zeros(nbin2,1);

    %integrate over ell space
    for i1=1:n
        l1=el1(i1);
        for j1=1:n
            l2=el1(j1);

            l=sqrt(l1.^2+l2.^2);

            if (l<=Max_l_mode && l>=Min_l_mode) 
            if (logl_flag==1) 
                ibin=int32(log10((l+1.)./Min_l_mode)./dlogl)+1;
                gPowEE(ibin)  = gPowEE(ibin)+gCEE_2(i1,j1)/log(10.);
                gPowBB(ibin)  = gPowBB(ibin)+gCBB_2(i1,j1)/log(10.);
                gPowEB(ibin)  = gPowEB(ibin)+gCEB_2(i1,j1)/log(10.);
            end
            if (logl_flag==0) 
                ibin=int32((l+1-Min_l_mode)./dlogl)+1; 
                gPowEE(ibin)  = gPowEE(ibin)+gCEE_2(i1,j1)*l;
                gPowBB(ibin)  = gPowBB(ibin)+gCBB_2(i1,j1)*l;
                gPowEB(ibin)  = gPowEB(ibin)+gCEB_2(i1,j1)*l;     
            end
            ll(ibin)=l; 
            if ibin > 1
	    dll(ibin-1)=log(ll(ibin))-log(ll(ibin-1));
            end
            end
        end 
    end
    %normalise the integrals with the ML fft normalisation
    gPowEE  =  (gPowEE./(n.^4.*dlogl));
    gPowBB  =  (gPowBB./(n.^4.*dlogl));
    gPowEB  =  (gPowEB./(n.^4.*dlogl));

    %_________________________________________________________________________%

    %now make the C_sys
    meantPowEE=tPowEE; %rename (left in from old version of code)
    meanPowEE =gPowEE;

    %make C_sys at all here only restrict range at sigam_sys calculation 					 
    for j=1:nbin2
      Csys(1,1,j)=(meanPowEE(j,1)-meantPowEE(j,1));

      % -(0.1*0.1/40000.)); %will approx remove residual noise term
      % noise residual term removal changes Q values by 
      % approx e.g. 
      % Q=52.310702 to Q_10=52.070384 (set 1 leaderboard #13 entry)
    end   

    %now calculate sigma_sys
    sigmasys=0.;
    for j=imin:imax
     sigmasys=sigmasys+(abs(Csys(1,1,j))).*(dll(j)); 
     %has to be the same way that make_varshears calculates the QN (which has no 2pi) 
     %note that power is already multipled by l^2 
    end
    q10=QN./sigmasys;

    %checks for the g2check convention loop loop
    if (q10 >= Q10_final) 
        Q10_final=q10; %uses maximum Q from possible conventions
    end
    
    fprintf(' Q_10=%f\n',q10); %Q of the mean

end %loop over convention checks 

avq=avq+Q10_final; %add best Q from conventions to the running total

end %old loop 1:1 

fprintf(' %d set=%d Q=%f\n',g2checkbest,tiernum,avq);

fprintf(fout,' %d %f\n',tiernum,avq); %write Q for set to file

allQ=allQ+avq; 

end %loop over sets
fclose(fout);

allQ=allQ/ntiers; %average over sets 

fprintf(' Q=%f\n',allQ); %print Q

%print final Q and sigma_sys to a file
outfilename=strcat('q10.galaxy.',methodtype,'.lis');
fout = fopen(outfilename,'a');
%bug in writting sigma_sys values %fprintf(fout,' %f %f\n',0.001/allQ,allQ);
fprintf(fout,' %f %f\n',QN./allQ,allQ);
fclose(fout);
%_________________________________________________________________________%
