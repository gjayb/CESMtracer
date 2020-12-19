function [weightedmean] = areaweightedmean(vartoaverage,areas,dim)
%AREAWEIGHTEDMEAN is built for CESM output. areas must have at least the
%first dimension the same length as vartoaverage. vartoaverage may be a 2-4
%dimensional matrix.
%   sample for mean vertical profile: 
%   plot(squeeze(nanmean(areaweightedmean(areaweightedmean(rho,taream,2),mean(taream,2),1),4)),-z)
ndimvar=length(size(vartoaverage));
ndimarea=length(size(areas));
areas2=areas;
if ndimarea<ndimvar
    switch ndimvar
        case 2
            [~,a]=size(vartoaverage);
            areas2=repmat(areas,1,a);
        case 3
            [~,a,b]=size(vartoaverage);
            if ndimarea==1
               areas2=repmat(areas,[1 a b]); 
            else
                areas2=repmat(areas,[1 1 b]);
            end
        case 4
            [~,a,b,c]=size(vartoaverage);
            switch ndimarea
                case 1
                    areas2=repmat(areas,[1 a b c]);
                case 2
                    areas2=repmat(areas,[1 1 b c]);
                case 3
                    areas2=repmat(areas,[1 1 1 c]);
            end
    end    
end

if size(vartoaverage)==size(areas2)
    areas2(isnan(vartoaverage))=NaN;
    weightedmean=nansum(vartoaverage.*areas2,dim)./nansum(areas2,dim);
else
    disp('error')
end
end

