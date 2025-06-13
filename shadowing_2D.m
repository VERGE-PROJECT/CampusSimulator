function y=shadowing_2D(Nsize,Msize,pixel_size,sigma,dist_corr)
%Nsize, Msize   Map dimensions in pixels.
%pixel_size=1;  Measured in meters
%sigma=6;    Standard deviation of the shadowing in dB
%dist_corr:     Decorrelation distance in m

if sigma>0
%Initial image with randomly distributed normal variables
initial_noise=sigma*randn(Nsize,Msize);


%To generate the correlation function
dcorr=dist_corr/pixel_size;
size_filter=10*dcorr;

[X,Y]=meshgrid(-0.5*size_filter:1:0.5*size_filter);

R=2.^(-sqrt(X.*X+Y.*Y)/dcorr);

H=fft2(R);
filter=sqrt(abs(H));
fase=exp(1i*angle(H));
final_filter=filter.*fase;
final_filter=ifft2(final_filter);
suma=sum(sum(final_filter.*final_filter));
final_filter=final_filter./sqrt(suma);

result=conv2(initial_noise,final_filter);

y=result(0.5*size_filter:Nsize+0.5*size_filter-1,0.5*size_filter:Msize+0.5*size_filter-1);

else
    %Sigma is 0, so there is no shadowing.
    y=zeros(Nsize,Msize);
end
% figure;
% imshow(y,[-20,20]);colormap(jet);colorbar
end
