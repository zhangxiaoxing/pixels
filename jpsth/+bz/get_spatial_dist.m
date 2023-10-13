function [avail,dist]=get_spatial_dist(reg1,reg2)
persistent reg_coord_map_
if isempty(reg_coord_map_)
    homedir=ephys.util.getHomedir('dtype','neupix');
    load(fullfile(replace(homedir,'SPKINFO','track_meta'),'reg2coord.mat'),'reg_coord_map');
    reg_coord_map_=reg_coord_map;
end
if reg_coord_map_.isKey(reg1) && reg_coord_map_.isKey(reg2)
    coord1=reg_coord_map_(char(reg1));
    coord2=reg_coord_map_(char(reg2));
    avail=true;
    dist=norm(coord1-coord2);
else
    avail=false;
    dist=intmax();
end
    
end