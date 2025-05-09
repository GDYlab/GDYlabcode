function cmp=stripedCMap(cmap,cb,cz,dc,strp)
%function cmp=stripedCMap(cmap,cb,cz,dc,strp)
% adds stripes and white zero;
% cmap(n,3) = some colormap
% cb = upper bound
% cz = upper lowest significant value (used to determine zero band)
% dc = stripe increment (above cz)
% strp = add stripes = true/false

        n=length(cmap(:,1));
% add white: cz corresponds to color index kz
% total number of colors:
        N=round(cb/(cb-cz)*n);
% zeros
        kz=round(cz/cb*N);
% stripe locations are
        ks=round((cz:dc:cb)/cb*N);
        cmp=ones(N,3);
       
% zero the map
        kM=round(0.02*n);
        vm=0.5;
        cmp(kz+1:N,:)=cmap;
% add stripes if requested
        if strp==0, return, end
        for k=1:length(ks)
                for k1=-kM:kM
                        kk=1+ks(k)+k1;
                        if kk>0 & kk<=N,
                                cmp(kk,:)= hsv2rgb(rgb2hsv(cmp(kk,:)).*[1,1,vm+sqrt(abs(k1)/kM)*(1-vm)]);
                        end            
                        kk=1-ks(k)+k1;
                        if kk>0 & kk<=2*N+1,
                                cmp(kk,:)= hsv2rgb(rgb2hsv(cmp(kk,:)).*[1,1,vm+sqrt(abs(k1)/kM)*(1-vm)]);
                        end            
                end    

        end


end

