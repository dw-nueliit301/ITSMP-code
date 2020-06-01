%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ITSMP MATLAB code
% 20190628
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

diary('diary.txt');

inputfile='Inputdata.xlsx';
downstream_ID=1;
upstream_ID=2;
overnight_train_ID=1;
regular_train_ID=2;
inspection_train_ID=3;

objective_travel=1/3;
objective_planned_time=1/3;
objective_maintenance=1-objective_travel-objective_planned_time;
max_earliness=20*60;
max_tardiness=20*60;
extend_start_time=19*3600;
extend_end_time=9*3600;
horizon_start_time=20*3600;
horizon_end_time=8*3600;

sheetname='Station_info';
[station_info,station_name]=xlsread(inputfile,sheetname);
station_name(1,:)=[];
station_num=size(station_info,1);
station_track_info=zeros(sum(sum(station_info(:,2:3))),3);
vindex=0;
for station_idx=1:station_num
    for track_idx=1:station_info(station_idx,2)
        station_track_info(vindex+1,:)=[station_idx 2*track_idx-1 downstream_ID];
        vindex=vindex+1;
    end
    for track_idx=1:station_info(station_idx,3)
        station_track_info(vindex+1,:)=[station_idx 2*track_idx upstream_ID];
        vindex=vindex+1;
    end
end

sheetname='train_info';
[train_info,train_name]=xlsread(inputfile,sheetname);
train_num=size(train_info,1);
train_name(1,:)=[];

sheetname='down_arrive_time';
down_arrive_time=xlsread(inputfile,sheetname);
sheetname='down_depart_time';
down_depart_time=xlsread(inputfile,sheetname);
sheetname='up_arrive_time';
up_arrive_time=xlsread(inputfile,sheetname);
sheetname='up_depart_time';
up_depart_time=xlsread(inputfile,sheetname);
origin_arrive_time=[down_arrive_time;fliplr(up_arrive_time)];
origin_depart_time=[down_depart_time;fliplr(up_depart_time)];

sheetname='other_info';
other_info=xlsread(inputfile,sheetname);
accelerate_time=other_info(1);
decelerate_time=other_info(2);

train_section_running_info=zeros(train_num*station_num,5);
train_path=cell(train_num,1);
vindex=0;
for train_idx=1:train_num
    vindex1=find(~isnan(origin_arrive_time(train_idx,:)))';
    if train_info(train_idx,3)==downstream_ID
        train_path{train_idx}=vindex1;
    else
        train_path{train_idx}=flipud(vindex1);
    end
    for station_idx=1:length(train_path{train_idx})-1
        train_section_running_info(vindex+1,1)=train_idx;
        train_section_running_info(vindex+1,2)=train_path{train_idx}(station_idx);
        train_section_running_info(vindex+1,3)=train_path{train_idx}(station_idx+1);        
        if origin_depart_time(train_idx,train_path{train_idx}(station_idx))~=origin_arrive_time(train_idx,train_path{train_idx}(station_idx))
            stop_start=1;
        else
            stop_start=0;
        end
        if origin_depart_time(train_idx,train_path{train_idx}(station_idx+1))~=origin_arrive_time(train_idx,train_path{train_idx}(station_idx+1))
            stop_end=1;
        else
            stop_end=0;
        end
        running_time=origin_arrive_time(train_idx,train_path{train_idx}(station_idx+1))-origin_depart_time(train_idx,train_path{train_idx}(station_idx));
        if stop_start==1 && stop_end==1
            train_section_running_info(vindex+1,4)=running_time-accelerate_time-decelerate_time;
        elseif stop_start==1 && stop_end==0
            train_section_running_info(vindex+1,4)=running_time-accelerate_time;
        elseif stop_start==0 && stop_end==1
            train_section_running_info(vindex+1,4)=running_time-decelerate_time;
        else
            train_section_running_info(vindex+1,4)=running_time;
        end
        train_section_running_info(vindex+1,5)=train_section_running_info(vindex+1,4)+accelerate_time+decelerate_time;
        vindex=vindex+1;
    end
end
train_section_running_info(all(train_section_running_info==0,2),:)=[];
train_section_running=zeros(train_num,station_num,station_num,2);
for vindex=1:size(train_section_running_info,1)
    train_section_running(train_section_running_info(vindex,1),train_section_running_info(vindex,2),train_section_running_info(vindex,3),1)=train_section_running_info(vindex,4);
    train_section_running(train_section_running_info(vindex,1),train_section_running_info(vindex,2),train_section_running_info(vindex,3),2)=train_section_running_info(vindex,5);
end

vindex=origin_depart_time-origin_arrive_time;
min_station_dwell_time=zeros(station_num,1);
for station_idx=1:station_num
    vindex1=find(~isnan(vindex(:,station_idx)) & vindex(:,station_idx)~=0);
    if ~isempty(vindex1)==1
        min_station_dwell_time(station_idx)=min(vindex(vindex1,station_idx));
    end
end

sheetname='overnight_train_dwell_info';
[overnight_train_dwell_time,overnight_train_dwell_name]=xlsread(inputfile,sheetname);
overnight_train_dwell_name(1,:)=[];

sheetname='overnight_train_scale_info';
[overnight_train_scale_time,overnight_train_scale_name]=xlsread(inputfile,sheetname);
overnight_train_scale_name(1,:)=[];

sheetname='maintenance_duration';
maintenance_info=xlsread(inputfile,sheetname);
for segment_idx=1:size(maintenance_info,1)
    if maintenance_info(segment_idx,3)<=12*3600
        maintenance_info(segment_idx,3)=maintenance_info(segment_idx,3)+24*3600;
    end
    if maintenance_info(segment_idx,4)<=12*3600
        maintenance_info(segment_idx,4)=maintenance_info(segment_idx,4)+24*3600;
    end
end
maintenance_duration=[maintenance_info(:,1:2) maintenance_info(:,4)-maintenance_info(:,3)];    
maintenance_duration=full(sparse(maintenance_duration(:,1),maintenance_duration(:,2),maintenance_duration(:,3),station_num,station_num));    

train_station_dwell_info=zeros(train_num*station_num,7);
vindex=0;
for train_idx=1:train_num
    for station_idx=1:length(train_path{train_idx})
        train_station_dwell_info(vindex+1,1)=train_idx;
        train_station_dwell_info(vindex+1,2)=train_path{train_idx}(station_idx);
        if train_info(train_idx,2)==regular_train_ID
            train_station_dwell_info(vindex+1,3)=min_station_dwell_time(train_path{train_idx}(station_idx));
            train_station_dwell_info(vindex+1,4)=origin_arrive_time(train_idx,train_path{train_idx}(station_idx))-max_earliness;
            train_station_dwell_info(vindex+1,5)=origin_arrive_time(train_idx,train_path{train_idx}(station_idx))+max_tardiness;
            train_station_dwell_info(vindex+1,6)=origin_depart_time(train_idx,train_path{train_idx}(station_idx))-max_earliness;
            train_station_dwell_info(vindex+1,7)=origin_depart_time(train_idx,train_path{train_idx}(station_idx))+max_tardiness;
            train_station_dwell_info(vindex+1,8)=origin_arrive_time(train_idx,train_path{train_idx}(station_idx));
            train_station_dwell_info(vindex+1,9)=origin_depart_time(train_idx,train_path{train_idx}(station_idx));
        elseif train_info(train_idx,2)==overnight_train_ID
            for vindex1=1:length(overnight_train_dwell_time)
                if strcmp(overnight_train_dwell_name{vindex1,1},train_name{train_idx,1})==1 && strcmp(overnight_train_dwell_name{vindex1,2},station_name{train_path{train_idx}(station_idx),1})==1
                    train_station_dwell_info(vindex+1,3)=overnight_train_dwell_time(vindex1)*60;
                end
            end
            for vindex1=1:size(overnight_train_scale_time,1)
                if strcmp(overnight_train_scale_name{vindex1,1},train_name{train_idx,1})==1 && strcmp(overnight_train_scale_name{vindex1,2},station_name{train_path{train_idx}(1),1})==1
                    o_time1=overnight_train_scale_time(vindex1,1)*3600;
                    o_time2=overnight_train_scale_time(vindex1,2)*3600;
                    o_time3=overnight_train_scale_time(vindex1,3)*3600;
                    o_time4=overnight_train_scale_time(vindex1,4)*3600;
                end
                if strcmp(overnight_train_scale_name{vindex1,1},train_name{train_idx,1})==1 && strcmp(overnight_train_scale_name{vindex1,2},station_name{train_path{train_idx}(end),1})==1
                    d_time1=overnight_train_scale_time(vindex1,1)*3600;
                    d_time2=overnight_train_scale_time(vindex1,2)*3600;
                    d_time3=overnight_train_scale_time(vindex1,3)*3600;
                    d_time4=overnight_train_scale_time(vindex1,4)*3600;
                end
            end
            if station_idx==1
                train_station_dwell_info(vindex+1,4)=o_time1;
                train_station_dwell_info(vindex+1,5)=o_time2;
                train_station_dwell_info(vindex+1,6)=o_time3;
                train_station_dwell_info(vindex+1,7)=o_time4;
            elseif station_idx==length(train_path{train_idx})
                train_station_dwell_info(vindex+1,4)=d_time1;
                train_station_dwell_info(vindex+1,5)=d_time2;
                train_station_dwell_info(vindex+1,6)=d_time3;
                train_station_dwell_info(vindex+1,7)=d_time4;
            else
                for segment_idx=1:station_idx-1
                    o_time3=o_time3+train_section_running(train_idx,train_path{train_idx}(segment_idx),train_path{train_idx}(segment_idx+1),1);
                end
                for segment_idx=station_idx:length(train_path{train_idx})-1
                    d_time2=d_time2-train_section_running(train_idx,train_path{train_idx}(segment_idx),train_path{train_idx}(segment_idx+1),1);
                end
                train_station_dwell_info(vindex+1,4)=o_time3;
                train_station_dwell_info(vindex+1,5)=d_time2;
                train_station_dwell_info(vindex+1,6)=train_station_dwell_info(vindex+1,4);
                train_station_dwell_info(vindex+1,7)=train_station_dwell_info(vindex+1,5);
            end
        elseif train_info(train_idx,2)==inspection_train_ID
            train_station_dwell_info(vindex+1,3)=origin_depart_time(train_idx,train_path{train_idx}(station_idx))-origin_arrive_time(train_idx,train_path{train_idx}(station_idx));
            if station_idx~=length(train_path{train_idx})
                train_station_dwell_info(vindex+1,6)=other_info(10)+maintenance_duration(train_path{train_idx}(station_idx),train_path{train_idx}(station_idx+1))+other_info(9);
            else
                train_station_dwell_info(vindex+1,6)=other_info(10)+other_info(9)+other_info(7);
            end
            train_station_dwell_info(vindex+1,4)=train_station_dwell_info(vindex+1,6);
            train_station_dwell_info(vindex+1,5)=origin_arrive_time(train_idx,train_path{train_idx}(station_idx))+max_tardiness+60*60;
            train_station_dwell_info(vindex+1,7)=origin_depart_time(train_idx,train_path{train_idx}(station_idx))+max_tardiness+60*60;
            train_station_dwell_info(vindex+1,8)=origin_arrive_time(train_idx,train_path{train_idx}(station_idx));
            train_station_dwell_info(vindex+1,9)=origin_depart_time(train_idx,train_path{train_idx}(station_idx));  
        end
        vindex=vindex+1;
    end
end
train_station_dwell_info(all(train_station_dwell_info==0,2),:)=[];
train_station_dwell=zeros(train_num,station_num,7);
for vindex=1:size(train_station_dwell_info,1)
    for vindex1=1:7
        train_station_dwell(train_station_dwell_info(vindex,1),train_station_dwell_info(vindex,2),vindex1)=train_station_dwell_info(vindex,vindex1+2);
    end
end

train_passed_section=cell(train_num,1);
for train_idx=1:train_num
    train_passed_section{train_idx}=train_section_running_info(train_section_running_info(:,1)==train_idx,2:3);
end

sheetname='inspection_train';
[~,inspection_train]=xlsread(inputfile,sheetname);

sheetname='origin_track_occupation';
origin_track_occupation_info=xlsread(inputfile,sheetname);
origin_track_occupation=zeros(train_num,station_num,max(origin_track_occupation_info(:,3)));
for vindex=1:size(origin_track_occupation_info,1)
    origin_track_occupation(origin_track_occupation_info(vindex,1),origin_track_occupation_info(vindex,2),origin_track_occupation_info(vindex,3))=origin_track_occupation_info(vindex,4);
end

subscript_im=zeros(size(train_section_running_info,1)+train_num,2);
vindex=0;
for train_idx=1:train_num
    subscript_im(vindex+(1:length(train_path{train_idx})),:)=[train_idx*ones(length(train_path{train_idx}),1) train_path{train_idx}];
    vindex=vindex+length(train_path{train_idx});
end
subscript_im_prenum=full(sparse(subscript_im(:,1),subscript_im(:,2),(1:size(subscript_im,1))',train_num,station_num));

subscript_ijmn=zeros(train_num*train_num*station_num*station_num,4);
vindex=0;
for train_idx=1:train_num-1
    for train_idx1=train_idx+1:train_num
        overlap_station=intersect(train_path{train_idx},train_path{train_idx1});
        if train_info(train_idx,3)==downstream_ID && train_info(train_idx1,3)==downstream_ID
            for station_idx=1:length(overlap_station)-1
                if ~((train_station_dwell(train_idx,overlap_station(station_idx),5)+other_info(3)<=train_station_dwell(train_idx1,overlap_station(station_idx),4) &&...
                        train_station_dwell(train_idx,overlap_station(station_idx+1),3)+other_info(4)<=train_station_dwell(train_idx1,overlap_station(station_idx+1),2)) ||...
                        (train_station_dwell(train_idx1,overlap_station(station_idx),5)+other_info(3)<=train_station_dwell(train_idx,overlap_station(station_idx),4) &&...
                        train_station_dwell(train_idx1,overlap_station(station_idx+1),3)+other_info(4)<=train_station_dwell(train_idx,overlap_station(station_idx+1),2)))
                    subscript_ijmn(vindex+1,:)=[train_idx train_idx1 overlap_station(station_idx) overlap_station(station_idx+1)];
                    vindex=vindex+1;
                end
            end   
        elseif train_info(train_idx,3)==upstream_ID && train_info(train_idx1,3)==upstream_ID
            for station_idx=length(overlap_station):-1:2
                if ~((train_station_dwell(train_idx,overlap_station(station_idx),5)+other_info(3)<=train_station_dwell(train_idx1,overlap_station(station_idx),4) &&...
                        train_station_dwell(train_idx,overlap_station(station_idx-1),3)+other_info(4)<=train_station_dwell(train_idx1,overlap_station(station_idx-1),2)) ||...
                        (train_station_dwell(train_idx1,overlap_station(station_idx),5)+other_info(3)<=train_station_dwell(train_idx,overlap_station(station_idx),4) &&...
                        train_station_dwell(train_idx1,overlap_station(station_idx-1),3)+other_info(4)<=train_station_dwell(train_idx,overlap_station(station_idx-1),2)))
                    subscript_ijmn(vindex+1,:)=[train_idx train_idx1 overlap_station(station_idx) overlap_station(station_idx-1)];
                    vindex=vindex+1;
                end
            end
        end
    end
end
subscript_ijmn(all(subscript_ijmn==0,2),:)=[];    
subscript_ijmn_prenum=zeros(train_num,train_num,station_num,station_num);
for vindex=1:size(subscript_ijmn,1)
    subscript_ijmn_prenum(subscript_ijmn(vindex,1),subscript_ijmn(vindex,2),subscript_ijmn(vindex,3),subscript_ijmn(vindex,4))=vindex;
end

subscript_ijm=zeros(train_num*train_num*station_num,3);
vindex=0;
for train_idx=1:train_num
    for train_idx1=1:train_num
        if train_idx~=train_idx1
            overlap_station=intersect(train_path{train_idx},train_path{train_idx1});
            for station_idx=1:length(overlap_station)
                if ~(train_station_dwell(train_idx,overlap_station(station_idx),5)+other_info(5)<=train_station_dwell(train_idx1,overlap_station(station_idx),2) ||...
                        train_station_dwell(train_idx1,overlap_station(station_idx),5)+other_info(5)<=train_station_dwell(train_idx,overlap_station(station_idx),2))
                    subscript_ijm(vindex+1,:)=[train_idx train_idx1 overlap_station(station_idx)];
                    vindex=vindex+1;
                end   
            end
        end
    end
end
subscript_ijm(all(subscript_ijm==0,2),:)=[];
subscript_ijm_prenum=zeros(train_num,train_num,station_num);
for vindex=1:size(subscript_ijm,1)
    subscript_ijm_prenum(subscript_ijm(vindex,1),subscript_ijm(vindex,2),subscript_ijm(vindex,3))=vindex;
end

subscript_mn=zeros(2*(station_num-1),2);
vindex=0;
for station_idx=1:station_num-1
    subscript_mn(vindex+(1:2),:)=[station_idx station_idx+1;station_idx+1 station_idx];
    vindex=vindex+2;
end
subscript_mn_prenum=full(sparse(subscript_mn(:,1),subscript_mn(:,2),(1:size(subscript_mn,1))',station_num,station_num));

subscript_imn=zeros(train_num*station_num*station_num,3);
vindex=0;
for train_idx=1:train_num
    for section_idx=1:size(train_passed_section{train_idx},1)
        if ~(train_station_dwell(train_idx,train_passed_section{train_idx}(section_idx,2),3)+other_info(8)<=other_info(10) ||...
                other_info(11)+other_info(7)+other_info(9)<=train_station_dwell(train_idx,train_passed_section{train_idx}(section_idx,1),4))
            subscript_imn(vindex+1,:)=[train_idx train_passed_section{train_idx}(section_idx,:)];
            vindex=vindex+1;
        end
    end
end
subscript_imn(all(subscript_imn==0,2),:)=[];
subscript_imn_prenum=zeros(train_num,station_num,station_num);
for vindex=1:size(subscript_imn,1)
    subscript_imn_prenum(subscript_imn(vindex,1),subscript_imn(vindex,2),subscript_imn(vindex,3))=vindex;
end

subscript_mm=zeros(station_num*station_num,2);
vindex=0;
for station_idx=1:station_num-2
    subscript_mm(vindex+1,:)=[station_idx station_idx+1];
    vindex=vindex+1;
end
subscript_mm(all(subscript_mm==0,2),:)=[];
subscript_mm_prenum=full(sparse(subscript_mm(:,1),subscript_mm(:,2),(1:size(subscript_mm,1))',station_num,station_num));

subscript_imk=zeros(train_num*sum(sum(station_info(:,2:3))),3);
vindex=0;
for train_idx=1:train_num
    for station_idx=1:length(train_path{train_idx})
        track_set=station_track_info(station_track_info(:,1)==train_path{train_idx}(station_idx),2);
        for track_idx=1:length(track_set)
            subscript_imk(vindex+1,:)=[train_idx train_path{train_idx}(station_idx) track_set(track_idx)];
            vindex=vindex+1;
        end
    end
end
subscript_imk(all(subscript_imk==0,2),:)=[];
subscript_imk_prenum=zeros(train_num,station_num,max(station_track_info(:,2)));
for vindex=1:size(subscript_imk)
    subscript_imk_prenum(subscript_imk(vindex,1),subscript_imk(vindex,2),subscript_imk(vindex,3))=vindex;
end

vnum_x=size(subscript_im,1);
vnum_y=size(subscript_im,1);
vnum_w=size(subscript_im,1);
vnum_u=size(subscript_ijmn,1);
vnum_p=size(subscript_ijm,1);
vnum_z=size(subscript_mn,1);
vnum_e=size(subscript_imn,1);
vnum_v=size(subscript_imk,1);
vnum_alpha=size(subscript_mm,1);
vnum_delta=size(subscript_im,1);
vnum_epsilon=size(subscript_im,1);
vnum=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z+vnum_e+vnum_v+vnum_alpha+vnum_delta+vnum_epsilon;

Atemp=0;
Atem=0;

for train_idx=1:train_num
    for station_idx=1:length(train_path{train_idx})
        Atemp=Atemp+1;
        Atem=Atem+2; 
    end
end

for train_idx=1:train_num
    for section_idx=1:size(train_passed_section{train_idx},1)
        Atemp=Atemp+1;
        Atem=Atem+4;
    end
end
for train_idx=1:train_num
    for section_idx=1:size(train_passed_section{train_idx},1)
        Atemp=Atemp+1;
        Atem=Atem+4;
    end
end

for train_idx=1:train_num
    for station_idx=1:length(train_path{train_idx})
        Atemp=Atemp+1;
        Atem=Atem+3;
    end
end
for train_idx=1:train_num
    for station_idx=1:length(train_path{train_idx})
        Atemp=Atemp+1;
        Atem=Atem+3;
    end
end

for train_idx=1:train_num
    if train_info(train_idx,2)==regular_train_ID
        for station_idx=1:length(train_path{train_idx})
            Atemp=Atemp+1;
            Atem=Atem+2;
        end
    end
end
for train_idx=1:train_num
    if train_info(train_idx,2)==regular_train_ID
        for station_idx=1:length(train_path{train_idx})
            Atemp=Atemp+1;
            Atem=Atem+2;
        end
    end
end

for train_idx=1:train_num
    if train_info(train_idx,2)==regular_train_ID
        for station_idx=1:length(train_path{train_idx})
            Atemp=Atemp+1;
            Atem=Atem+2;
        end
    end
end
for train_idx=1:train_num
    if train_info(train_idx,2)==regular_train_ID
        for station_idx=1:length(train_path{train_idx})
            Atemp=Atemp+1;
            Atem=Atem+2;
        end
    end
end


for vindex=1:size(subscript_ijmn,1)
    Atemp=Atemp+1;
    Atem=Atem+3;
end

for vindex=1:size(subscript_ijmn,1)
    Atemp=Atemp+1;
    Atem=Atem+3;
end

for vindex=1:size(subscript_ijmn,1)
    Atemp=Atemp+1;
    Atem=Atem+3;
end

for vindex=1:size(subscript_ijmn,1)
    Atemp=Atemp+1;
    Atem=Atem+3;
end

for vindex=1:size(subscript_ijm,1)
    if train_info(subscript_ijm(vindex,1),3)==train_info(subscript_ijm(vindex,2),3)
        Atemp=Atemp+1;
        Atem=Atem+3;
    else
        Atemp=Atemp+1;
        if train_info(subscript_ijm(vindex,1),3)==downstream_ID
            location3=find(subscript_imk(:,1)==subscript_ijm(vindex,1) & subscript_imk(:,2)==subscript_ijm(vindex,3) & mod(subscript_imk(:,3),2)==0);
        else
            location3=find(subscript_imk(:,1)==subscript_ijm(vindex,1) & subscript_imk(:,2)==subscript_ijm(vindex,3) & mod(subscript_imk(:,3),2)==1);
        end
        Atem=Atem+3+length(location3);
    end
end

for vindex=1:size(subscript_ijm,1)
    if train_info(subscript_ijm(vindex,1),3)~=train_info(subscript_ijm(vindex,2),3)     
        Atemp=Atemp+1;
        if train_info(subscript_ijm(vindex,1),3)==downstream_ID
            location3=find(subscript_imk(:,1)==subscript_ijm(vindex,1) & subscript_imk(:,2)==subscript_ijm(vindex,3) & mod(subscript_imk(:,3),2)==0);
        else
            location3=find(subscript_imk(:,1)==subscript_ijm(vindex,1) & subscript_imk(:,2)==subscript_ijm(vindex,3) & mod(subscript_imk(:,3),2)==1);
        end
        Atem=Atem+3+length(location3);
    end
end

for train_idx=1:train_num
    for station_idx=1:length(train_path{train_idx})
        Atemp=Atemp+1;
        Atem=Atem+3;
    end
end

for vindex=1:size(subscript_ijm,1)
    if subscript_ijm(vindex,1)<subscript_ijm(vindex,2)
        track_set=station_track_info(station_track_info(:,1)==subscript_ijm(vindex,3),2);
        for track_idx=1:length(track_set)
            Atemp=Atemp+1;
            Atem=Atem+4;
        end
    end
end

for vindex=1:size(subscript_imn,1)
    Atemp=Atemp+1;
    Atem=Atem+3;
end

for vindex=1:size(subscript_imn,1)
    Atemp=Atemp+1;
    Atem=Atem+3;
end

for vindex=1:size(subscript_mm,1)
    Atemp=Atemp+1;
    Atem=Atem+3;
end

for vindex=1:size(subscript_mm,1)
    Atemp=Atemp+1;
    Atem=Atem+3;
end

for vindex=1:size(subscript_ijmn,1)
    location1=subscript_ijm_prenum(subscript_ijmn(vindex,1),subscript_ijmn(vindex,2),subscript_ijmn(vindex,3));
    if location1~=0
        Atemp=Atemp+1;
        Atem=Atem+2;
    end
end
for vindex=1:size(subscript_ijmn,1)
    location1=subscript_ijm_prenum(subscript_ijmn(vindex,2),subscript_ijmn(vindex,1),subscript_ijmn(vindex,3));
    if location1~=0
        Atemp=Atemp+1;
        Atem=Atem+2;
    end
end

for vindex=1:size(subscript_ijmn,1)
    if ismember(train_name{subscript_ijmn(vindex,1),1},inspection_train)==1 && strcmp(train_name{subscript_ijmn(vindex,2),4},'确认车')==0
        location1=subscript_imn_prenum(subscript_ijmn(vindex,2),subscript_ijmn(vindex,3),subscript_ijmn(vindex,4));%jmn
        if location1~=0
            Atemp=Atemp+1;
            Atem=Atem+2;
        end
    elseif strcmp(train_name{subscript_ijmn(vindex,1),4},'确认车')==0 && ismember(train_name{subscript_ijmn(vindex,2),1},inspection_train)==1
        location1=subscript_imn_prenum(subscript_ijmn(vindex,1),subscript_ijmn(vindex,3),subscript_ijmn(vindex,4));%imn
        if location1~=0
            Atemp=Atemp+1;
            Atem=Atem+2;
        end
    end
end

Aineq=zeros(Atem,3);
bineq=zeros(Atemp,1);
Atemp=0;
Atem=0;

for train_idx=1:train_num
    for station_idx=1:length(train_path{train_idx})
        Atemp=Atemp+1;
        location1=subscript_im_prenum(train_idx,train_path{train_idx}(station_idx));
        Arng=[location1;vnum_x+location1];
        Acoe=[1;-1];
        Aineq(Atem+(1:2),:)=[ones(2,1)*Atemp Arng Acoe];
        Atem=Atem+2;
        if train_info(train_idx,2)==regular_train_ID
            bineq(Atemp)=-1*min(train_station_dwell(train_idx,train_path{train_idx}(station_idx),1),...
                origin_depart_time(train_idx,train_path{train_idx}(station_idx))-origin_arrive_time(train_idx,train_path{train_idx}(station_idx)));
        elseif train_info(train_idx,2)==overnight_train_ID
            bineq(Atemp)=-1*train_station_dwell(train_idx,train_path{train_idx}(station_idx),1);
        end
    end
end

vindex=vnum_x+vnum_y;
for train_idx=1:train_num
    for section_idx=1:size(train_passed_section{train_idx},1)
        Atemp=Atemp+1;
        location1=subscript_im_prenum(train_idx,train_passed_section{train_idx}(section_idx,1));%下标为im
        location2=subscript_im_prenum(train_idx,train_passed_section{train_idx}(section_idx,2));%下标为in
        Arng=[location2;vnum_x+location1;vindex+location1;vindex+location2];
        Acoe=[-1;1;other_info(1);other_info(2)];
        Aineq(Atem+(1:4),:)=[ones(4,1)*Atemp Arng Acoe];
        Atem=Atem+4;
        bineq(Atemp)=-1*train_section_running(train_idx,train_passed_section{train_idx}(section_idx,1),train_passed_section{train_idx}(section_idx,2),1);
    end
end

for train_idx=1:train_num
    for section_idx=1:size(train_passed_section{train_idx},1)
        Atemp=Atemp+1;
        location1=subscript_im_prenum(train_idx,train_passed_section{train_idx}(section_idx,1));%下标为im
        location2=subscript_im_prenum(train_idx,train_passed_section{train_idx}(section_idx,2));%下标为in
        Arng=[location2;vnum_x+location1;vindex+location1;vindex+location2];
        Acoe=[1;-1;-1*other_info(1);-1*other_info(2)];
        Aineq(Atem+(1:4),:)=[ones(4,1)*Atemp Arng Acoe];
        Atem=Atem+4;
        bineq(Atemp)=train_section_running(train_idx,train_passed_section{train_idx}(section_idx,1),train_passed_section{train_idx}(section_idx,2),2);
    end
end

vindex=vnum_x+vnum_y;
for train_idx=1:train_num
    for station_idx=1:length(train_path{train_idx})
        Atemp=Atemp+1;
        location1=subscript_im_prenum(train_idx,train_path{train_idx}(station_idx));
        Arng=[location1;vnum_x+location1;vindex+location1];
        Acoe=[1;-1;1];
        Aineq(Atem+(1:3),:)=[ones(3,1)*Atemp Arng Acoe];
        Atem=Atem+3;
    end
end
for train_idx=1:train_num
    for station_idx=1:length(train_path{train_idx})
        Atemp=Atemp+1;
        location1=subscript_im_prenum(train_idx,train_path{train_idx}(station_idx));
        Arng=[location1;vnum_x+location1;vindex+location1];
        M1=min(other_info(7),train_station_dwell(train_idx,train_path{train_idx}(station_idx),5)-train_station_dwell(train_idx,train_path{train_idx}(station_idx),2));
        Acoe=[-1;1;-1*M1];
        Aineq(Atem+(1:3),:)=[ones(3,1)*Atemp Arng Acoe];
        Atem=Atem+3;
    end
end

vindex=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z+vnum_e+vnum_v+vnum_alpha;
for train_idx=1:train_num
    if train_info(train_idx,2)==regular_train_ID
        for station_idx=1:length(train_path{train_idx})
            Atemp=Atemp+1;
            location1=subscript_im_prenum(train_idx,train_path{train_idx}(station_idx));
            Arng=[location1;vindex+location1];
            Acoe=[1;-1];
            Aineq(Atem+(1:2),:)=[ones(2,1)*Atemp Arng Acoe];
            Atem=Atem+2;
            bineq(Atemp)=train_station_dwell(train_idx,train_path{train_idx}(station_idx),6);
        end
    end
end
for train_idx=1:train_num
    if train_info(train_idx,2)==regular_train_ID
        for station_idx=1:length(train_path{train_idx})
            Atemp=Atemp+1;
            location1=subscript_im_prenum(train_idx,train_path{train_idx}(station_idx));
            Arng=[location1;vindex+location1];
            Acoe=[-1;-1];
            Aineq(Atem+(1:2),:)=[ones(2,1)*Atemp Arng Acoe];
            Atem=Atem+2;
            bineq(Atemp)=-1*train_station_dwell(train_idx,train_path{train_idx}(station_idx),6);
        end
    end
end

vindex=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z+vnum_e+vnum_v+vnum_alpha+vnum_delta;
for train_idx=1:train_num
    if train_info(train_idx,2)==regular_train_ID
        for station_idx=1:length(train_path{train_idx})
            Atemp=Atemp+1;
            location1=subscript_im_prenum(train_idx,train_path{train_idx}(station_idx));
            Arng=[vnum_x+location1;vindex+location1];
            Acoe=[1;-1];
            Aineq(Atem+(1:2),:)=[ones(2,1)*Atemp Arng Acoe];
            Atem=Atem+2;
            bineq(Atemp)=train_station_dwell(train_idx,train_path{train_idx}(station_idx),7);
        end
    end
end
for train_idx=1:train_num
    if train_info(train_idx,2)==regular_train_ID
        for station_idx=1:length(train_path{train_idx})
            Atemp=Atemp+1;
            location1=subscript_im_prenum(train_idx,train_path{train_idx}(station_idx));
            Arng=[vnum_x+location1;vindex+location1];
            Acoe=[-1;-1];
            Aineq(Atem+(1:2),:)=[ones(2,1)*Atemp Arng Acoe];
            Atem=Atem+2;
            bineq(Atemp)=-1*train_station_dwell(train_idx,train_path{train_idx}(station_idx),7);
        end
    end
end

vindex1=vnum_x+vnum_y+vnum_w;
for vindex=1:size(subscript_ijmn,1)
    Atemp=Atemp+1;
    location1=subscript_im_prenum(subscript_ijmn(vindex,1),subscript_ijmn(vindex,3));
    location2=subscript_im_prenum(subscript_ijmn(vindex,2),subscript_ijmn(vindex,3));
    Arng=[vnum_x+location1;vnum_x+location2;vindex1+vindex];
    h_1=min(other_info(3),abs(origin_depart_time(subscript_ijmn(vindex,1),subscript_ijmn(vindex,3))-origin_depart_time(subscript_ijmn(vindex,2),subscript_ijmn(vindex,3))));
    M2=max(train_station_dwell(subscript_ijmn(vindex,1),subscript_ijmn(vindex,3),5),train_station_dwell(subscript_ijmn(vindex,2),subscript_ijmn(vindex,3),5))-...
        min(train_station_dwell(subscript_ijmn(vindex,1),subscript_ijmn(vindex,3),4),train_station_dwell(subscript_ijmn(vindex,2),subscript_ijmn(vindex,3),4))+h_1;
    Acoe=[1;-1;M2];
    Aineq(Atem+(1:3),:)=[ones(3,1)*Atemp Arng Acoe];
    Atem=Atem+3;
    bineq(Atemp)=M2-h_1;
end

for vindex=1:size(subscript_ijmn,1)
    Atemp=Atemp+1;
    location1=subscript_im_prenum(subscript_ijmn(vindex,1),subscript_ijmn(vindex,3));
    location2=subscript_im_prenum(subscript_ijmn(vindex,2),subscript_ijmn(vindex,3));
    Arng=[vnum_x+location1;vnum_x+location2;vindex1+vindex];

    h_1=min(other_info(3),abs(origin_depart_time(subscript_ijmn(vindex,1),subscript_ijmn(vindex,3))-origin_depart_time(subscript_ijmn(vindex,2),subscript_ijmn(vindex,3))));
    M2=max(train_station_dwell(subscript_ijmn(vindex,1),subscript_ijmn(vindex,3),5),train_station_dwell(subscript_ijmn(vindex,2),subscript_ijmn(vindex,3),5))-...
        min(train_station_dwell(subscript_ijmn(vindex,1),subscript_ijmn(vindex,3),4),train_station_dwell(subscript_ijmn(vindex,2),subscript_ijmn(vindex,3),4))+h_1;
    Acoe=[-1;1;-1*M2];
    Aineq(Atem+(1:3),:)=[ones(3,1)*Atemp Arng Acoe];
    Atem=Atem+3;
    bineq(Atemp)=-1*h_1;
end

vindex1=vnum_x+vnum_y+vnum_w;
for vindex=1:size(subscript_ijmn,1)
    Atemp=Atemp+1;
    location1=subscript_im_prenum(subscript_ijmn(vindex,1),subscript_ijmn(vindex,4));
    location2=subscript_im_prenum(subscript_ijmn(vindex,2),subscript_ijmn(vindex,4));
    Arng=[location1;location2;vindex1+vindex];
    h_2=min(other_info(4),abs(origin_depart_time(subscript_ijmn(vindex,1),subscript_ijmn(vindex,4))-origin_depart_time(subscript_ijmn(vindex,2),subscript_ijmn(vindex,4))));
    M3=max(train_station_dwell(subscript_ijmn(vindex,1),subscript_ijmn(vindex,4),3),train_station_dwell(subscript_ijmn(vindex,2),subscript_ijmn(vindex,4),3))-...
        min(train_station_dwell(subscript_ijmn(vindex,1),subscript_ijmn(vindex,4),2),train_station_dwell(subscript_ijmn(vindex,2),subscript_ijmn(vindex,4),2))+h_2;
    Acoe=[1;-1;M3];
    Aineq(Atem+(1:3),:)=[ones(3,1)*Atemp Arng Acoe];
    Atem=Atem+3;
    bineq(Atemp)=M3-h_2;
end

for vindex=1:size(subscript_ijmn,1)
    Atemp=Atemp+1;
    location1=subscript_im_prenum(subscript_ijmn(vindex,1),subscript_ijmn(vindex,4));%in
    location2=subscript_im_prenum(subscript_ijmn(vindex,2),subscript_ijmn(vindex,4));%jn
    Arng=[location1;location2;vindex1+vindex]; 
    h_2=min(other_info(4),abs(origin_depart_time(subscript_ijmn(vindex,1),subscript_ijmn(vindex,4))-origin_depart_time(subscript_ijmn(vindex,2),subscript_ijmn(vindex,4))));
    M3=max(train_station_dwell(subscript_ijmn(vindex,1),subscript_ijmn(vindex,4),3),train_station_dwell(subscript_ijmn(vindex,2),subscript_ijmn(vindex,4),3))-...
        min(train_station_dwell(subscript_ijmn(vindex,1),subscript_ijmn(vindex,4),2),train_station_dwell(subscript_ijmn(vindex,2),subscript_ijmn(vindex,4),2))+h_2;
    Acoe=[-1;1;-1*M3];
    Aineq(Atem+(1:3),:)=[ones(3,1)*Atemp Arng Acoe];
    Atem=Atem+3;
    bineq(Atemp)=-1*h_2;
end

vindex1=vnum_x+vnum_y+vnum_w+vnum_u;
vindex2=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z+vnum_e;
for vindex=1:size(subscript_ijm,1)
    M4=max(train_station_dwell(subscript_ijm(vindex,1),subscript_ijm(vindex,3),5),train_station_dwell(subscript_ijm(vindex,2),subscript_ijm(vindex,3),3))...
        -min(train_station_dwell(subscript_ijm(vindex,1),subscript_ijm(vindex,3),4),train_station_dwell(subscript_ijm(vindex,2),subscript_ijm(vindex,3),2))+other_info(5);
    if train_info(subscript_ijm(vindex,1),3)==train_info(subscript_ijm(vindex,2),3)
        Atemp=Atemp+1;
        location1=subscript_im_prenum(subscript_ijm(vindex,1),subscript_ijm(vindex,3));%im
        location2=subscript_im_prenum(subscript_ijm(vindex,2),subscript_ijm(vindex,3));%jm
        Arng=[location2;vnum_x+location1;vindex1+vindex];
        Acoe=[-1;1;M4];
        Aineq(Atem+(1:3),:)=[ones(3,1)*Atemp Arng Acoe];
        Atem=Atem+3;
        bineq(Atemp)=M4-other_info(5);
    else
        Atemp=Atemp+1;
        location1=subscript_im_prenum(subscript_ijm(vindex,1),subscript_ijm(vindex,3));%im
        location2=subscript_im_prenum(subscript_ijm(vindex,2),subscript_ijm(vindex,3));%jm
        if train_info(subscript_ijm(vindex,1),3)==downstream_ID
            location3=find(subscript_imk(:,1)==subscript_ijm(vindex,1) & subscript_imk(:,2)==subscript_ijm(vindex,3) & mod(subscript_imk(:,3),2)==0);
        else
            location3=find(subscript_imk(:,1)==subscript_ijm(vindex,1) & subscript_imk(:,2)==subscript_ijm(vindex,3) & mod(subscript_imk(:,3),2)==1);
        end
        Arng=[location2;vnum_x+location1;vindex1+vindex;vindex2+location3];
        Acoe=[-1;1;M4;M4*ones(length(location3),1)];
        Aineq(Atem+(1:3+length(location3)),:)=[ones(3+length(location3),1)*Atemp Arng Acoe];
        Atem=Atem+3+length(location3);
        bineq(Atemp)=2*M4-other_info(5);
    end
end

vindex1=vnum_x+vnum_y+vnum_w+vnum_u;
vindex2=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z+vnum_e;
for vindex=1:size(subscript_ijm,1)
    if train_info(subscript_ijm(vindex,1),3)~=train_info(subscript_ijm(vindex,2),3)
        M4=max(train_station_dwell(subscript_ijm(vindex,1),subscript_ijm(vindex,3),5),train_station_dwell(subscript_ijm(vindex,2),subscript_ijm(vindex,3),3))...
            -min(train_station_dwell(subscript_ijm(vindex,1),subscript_ijm(vindex,3),4),train_station_dwell(subscript_ijm(vindex,2),subscript_ijm(vindex,3),2))+other_info(6);        
        Atemp=Atemp+1;
        location1=subscript_im_prenum(subscript_ijm(vindex,1),subscript_ijm(vindex,3));
        location2=subscript_im_prenum(subscript_ijm(vindex,2),subscript_ijm(vindex,3));
        if train_info(subscript_ijm(vindex,1),3)==downstream_ID
            location3=find(subscript_imk(:,1)==subscript_ijm(vindex,1) & subscript_imk(:,2)==subscript_ijm(vindex,3) & mod(subscript_imk(:,3),2)==0);
        else
            location3=find(subscript_imk(:,1)==subscript_ijm(vindex,1) & subscript_imk(:,2)==subscript_ijm(vindex,3) & mod(subscript_imk(:,3),2)==1);
        end
        Arng=[location2;vnum_x+location1;vindex1+vindex;vindex2+location3];
        Acoe=[1;-1;-1*M4;M4*ones(length(location3),1)];
        Aineq(Atem+(1:3+length(location3)),:)=[ones(3+length(location3),1)*Atemp Arng Acoe];
        Atem=Atem+3+length(location3);
        bineq(Atemp)=M4-other_info(6);
    end
end

vindex=vnum_x+vnum_y;
vindex1=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z+vnum_e;
for train_idx=1:train_num
    for station_idx=1:length(train_path{train_idx})
        Atemp=Atemp+1;
        location1=subscript_im_prenum(train_idx,train_path{train_idx}(station_idx));%im
        location2=subscript_imk_prenum(train_idx,train_path{train_idx}(station_idx),downstream_ID);%imsigma_im
        location3=subscript_imk_prenum(train_idx,train_path{train_idx}(station_idx),upstream_ID);%imsigma_im
        Arng=[vindex+location1;vindex1+location2;vindex1+location3];
        Acoe=[1;1;1];
        Aineq(Atem+(1:3),:)=[ones(3,1)*Atemp Arng Acoe];
        Atem=Atem+3;
        bineq(Atemp)=1;
    end
end

vindex1=vnum_x+vnum_y+vnum_w+vnum_u;
vindex2=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z+vnum_e;
for vindex=1:size(subscript_ijm,1)
    if subscript_ijm(vindex,1)<subscript_ijm(vindex,2)
        track_set=station_track_info(station_track_info(:,1)==subscript_ijm(vindex,3),2);
        for track_idx=1:length(track_set)
            Atemp=Atemp+1;
            location1=subscript_ijm_prenum(subscript_ijm(vindex,2),subscript_ijm(vindex,1),subscript_ijm(vindex,3)); %jim
            location2=subscript_imk_prenum(subscript_ijm(vindex,1),subscript_ijm(vindex,3),track_set(track_idx));%imk
            location3=subscript_imk_prenum(subscript_ijm(vindex,2),subscript_ijm(vindex,3),track_set(track_idx));%jmk
            Arng=[vindex1+vindex;vindex1+location1;vindex2+location2;vindex2+location3];
            Acoe=[-1;-1;1;1];
            Aineq(Atem+(1:4),:)=[ones(4,1)*Atemp Arng Acoe];
            Atem=Atem+4;
            bineq(Atemp)=1;
        end
    end
end

vindex1=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p;
vindex2=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z;
for vindex=1:size(subscript_imn,1)
    Atemp=Atemp+1;
    location1=subscript_im_prenum(subscript_imn(vindex,1),subscript_imn(vindex,3));
    location2=subscript_mn_prenum(subscript_imn(vindex,2),subscript_imn(vindex,3));
    Arng=[location1;vindex1+location2;vindex2+vindex];
    M5=max(other_info(11),train_station_dwell(subscript_imn(vindex,1),subscript_imn(vindex,3),3))-...
        min(other_info(10),train_station_dwell(subscript_imn(vindex,1),subscript_imn(vindex,3),2))+other_info(8);
    Acoe=[1;-1;M5];
    Aineq(Atem+(1:3),:)=[ones(3,1)*Atemp Arng Acoe];
    Atem=Atem+3;
    bineq(Atemp)=M5-other_info(8);
end

for vindex=1:size(subscript_imn,1)
    Atemp=Atemp+1;
    location1=subscript_im_prenum(subscript_imn(vindex,1),subscript_imn(vindex,2));
    location2=subscript_mn_prenum(subscript_imn(vindex,2),subscript_imn(vindex,3));
    Arng=[vnum_x+location1;vindex1+location2;vindex2+vindex];
    M6=max(other_info(11),train_station_dwell(subscript_imn(vindex,1),subscript_imn(vindex,2),5))-...
        min(other_info(10),train_station_dwell(subscript_imn(vindex,1),subscript_imn(vindex,2),4))+other_info(9)+maintenance_duration(subscript_imn(vindex,2),subscript_imn(vindex,3));
    Acoe=[-1;1;-1*M6];
    Aineq(Atem+(1:3),:)=[ones(3,1)*Atemp Arng Acoe];
    Atem=Atem+3;
    bineq(Atemp)=-1*maintenance_duration(subscript_imn(vindex,2),subscript_imn(vindex,3))-other_info(9);
end

vindex1=vnum_x+vnum_y+vnum_w;
vindex2=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z;
for vindex=1:size(subscript_ijmn,1)
    if ismember(train_name{subscript_ijmn(vindex,1),1},inspection_train)==1 && strcmp(train_name{subscript_ijmn(vindex,2),4},'确认车')==0
        location1=subscript_imn_prenum(subscript_ijmn(vindex,2),subscript_ijmn(vindex,3),subscript_ijmn(vindex,4));%jmn
        if location1~=0
            Atemp=Atemp+1;
            Arng=[vindex1+vindex;vindex2+location1];
            Acoe=[-1;-1];
            Aineq(Atem+(1:2),:)=[ones(2,1)*Atemp Arng Acoe];
            Atem=Atem+2;
            bineq(Atemp)=-1;
        end
    elseif strcmp(train_name{subscript_ijmn(vindex,1),4},'确认车')==0 && ismember(train_name{subscript_ijmn(vindex,2),1},inspection_train)==1
        location1=subscript_imn_prenum(subscript_ijmn(vindex,1),subscript_ijmn(vindex,3),subscript_ijmn(vindex,4));%imn
        if location1~=0
            Atemp=Atemp+1;
            Arng=[vindex1+vindex;vindex2+location1];
            Acoe=[1;-1];
            Aineq(Atem+(1:2),:)=[ones(2,1)*Atemp Arng Acoe];
            Atem=Atem+2;
        end
    end
end

vindex1=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p;
vindex2=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z+vnum_e+vnum_v;
for vindex=1:size(subscript_mm,1)
    Atemp=Atemp+1;
    location1=subscript_mn_prenum(subscript_mm(vindex,1),subscript_mm(vindex,1)+1);%mn
    location2=subscript_mn_prenum(subscript_mm(vindex,2),subscript_mm(vindex,2)+1);%m'n'
    Arng=[vindex1+location1;vindex1+location2;vindex2+vindex];
    Acoe=[1;-1;-1];
    Aineq(Atem+(1:3),:)=[ones(3,1)*Atemp Arng Acoe];
    Atem=Atem+3;
end

for vindex=1:size(subscript_mm,1)
    Atemp=Atemp+1;
    location1=subscript_mn_prenum(subscript_mm(vindex,1),subscript_mm(vindex,1)+1);%mn
    location2=subscript_mn_prenum(subscript_mm(vindex,2),subscript_mm(vindex,2)+1);%m'n'
    Arng=[vindex1+location1;vindex1+location2;vindex2+vindex];
    Acoe=[-1;1;-1];
    Aineq(Atem+(1:3),:)=[ones(3,1)*Atemp Arng Acoe];
    Atem=Atem+3;
end

vindex1=vnum_x+vnum_y+vnum_w;
vindex2=vnum_x+vnum_y+vnum_w+vnum_u;
for vindex=1:size(subscript_ijmn,1)
    location1=subscript_ijm_prenum(subscript_ijmn(vindex,1),subscript_ijmn(vindex,2),subscript_ijmn(vindex,3));
    if location1~=0
        Atemp=Atemp+1;
        Arng=[vindex1+vindex;vindex2+location1];
        Acoe=[-1;1];
        Aineq(Atem+(1:2),:)=[ones(2,1)*Atemp Arng Acoe];
        Atem=Atem+2;
    end
end
for vindex=1:size(subscript_ijmn,1)
    location1=subscript_ijm_prenum(subscript_ijmn(vindex,2),subscript_ijmn(vindex,1),subscript_ijmn(vindex,3));
    if location1~=0
        Atemp=Atemp+1;
        Arng=[vindex1+vindex;vindex2+location1];
        Acoe=[1;1];
        Aineq(Atem+(1:2),:)=[ones(2,1)*Atemp Arng Acoe];
        Atem=Atem+2;
        bineq(Atemp)=1;
    end
end

Aineq=sparse(Aineq(:,1),Aineq(:,2),Aineq(:,3),Atemp,vnum);

Aeqtemp=0;
Aeqtem=0;

for train_idx=1:train_num
    for station_idx=1:length(train_path{train_idx})
        Aeqtemp=Aeqtemp+1;
        Aeqtem=Aeqtem+length(location1);
    end
end

for train_idx=1:train_num
    if train_info(train_idx,2)==regular_train_ID
        for station_idx=1:length(train_path{train_idx})
            Aeqtemp=Aeqtemp+1;
            if train_info(train_idx,3)==downstream_ID
                location1=find(subscript_imk(:,1)==train_idx & subscript_imk(:,2)==train_path{train_idx}(station_idx) & mod(subscript_imk(:,3),2)==0);
            else
                location1=find(subscript_imk(:,1)==train_idx & subscript_imk(:,2)==train_path{train_idx}(station_idx) & mod(subscript_imk(:,3),2)==1);
            end
            Aeqtem=Aeqtem+length(location1);
        end
    end
end

for station_idx=1:station_num-1
    Aeqtemp=Aeqtemp+1;
    Aeqtem=Aeqtem+2;
end

for train_idx=1:train_num
    if train_info(train_idx,2)==regular_train_ID
        for station_idx=1:length(train_path{train_idx})
            if train_station_dwell(train_idx,train_path{train_idx}(station_idx),6)<horizon_start_time...
                    || train_station_dwell(train_idx,train_path{train_idx}(station_idx),6)>horizon_end_time+24*3600
                Aeqtemp=Aeqtemp+1;
                Aeqtem=Aeqtem+1;
            end
        end
    end
end

for train_idx=1:train_num
    if train_info(train_idx,2)==regular_train_ID
        for station_idx=1:length(train_path{train_idx})
            if train_station_dwell(train_idx,train_path{train_idx}(station_idx),7)<horizon_start_time...
                    || train_station_dwell(train_idx,train_path{train_idx}(station_idx),7)>horizon_end_time+24*3600
                Aeqtemp=Aeqtemp+1;
                Aeqtem=Aeqtem+1;
            end
        end
    end
end

for vindex=1:size(subscript_imk,1)
    if train_info(subscript_imk(vindex,1),2)==regular_train_ID
        if origin_arrive_time(subscript_imk(vindex,1),subscript_imk(vindex,2))<horizon_start_time...
                || origin_arrive_time(subscript_imk(vindex,1),subscript_imk(vindex,2))>horizon_end_time+24*3600
            Aeqtemp=Aeqtemp+1;
            Aeqtem=Aeqtem+1;
        end
    end
end

for vindex=1:size(subscript_imn,1)
    if ismember(train_name{subscript_imn(vindex,1),1},inspection_train)
        Aeqtemp=Aeqtemp+1;
        Aeqtem=Aeqtem+1;
    end
end

Aeq=zeros(Aeqtem,3);
beq=zeros(Aeqtemp,1);
Aeqtemp=0;
Aeqtem=0;

vindex=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z+vnum_e;
for train_idx=1:train_num
    for station_idx=1:length(train_path{train_idx})
        Aeqtemp=Aeqtemp+1;
        location1=find(subscript_imk(:,1)==train_idx & subscript_imk(:,2)==train_path{train_idx}(station_idx));
        Arng=vindex+location1;
        Acoe=ones(length(location1),1);
        Aeq(Aeqtem+(1:length(location1)),:)=[ones(length(location1),1)*Aeqtemp Arng Acoe];
        Aeqtem=Aeqtem+length(location1);
        beq(Aeqtemp)=1;
    end
end

vindex1=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z;
for vindex=1:size(subscript_imn,1)
    if ismember(train_name{subscript_imn(vindex,1),1},inspection_train)
        Aeqtemp=Aeqtemp+1;
        Arng=vindex1+vindex;
        Acoe=1;
        Aeq(Aeqtem+(1:1),:)=[ones(1,1)*Aeqtemp Arng Acoe];
        Aeqtem=Aeqtem+1;
    end
end

vindex=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p;
for station_idx=1:station_num-1
    Aeqtemp=Aeqtemp+1;
    location1=subscript_mn_prenum(station_idx,station_idx+1);
    location2=subscript_mn_prenum(station_idx+1,station_idx);
    Arng=[vindex+location1;vindex+location2];
    Acoe=[1;-1];
    Aeq(Aeqtem+(1:2),:)=[ones(2,1)*Aeqtemp Arng Acoe];
    Aeqtem=Aeqtem+2;
end

vindex=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z+vnum_e;
for train_idx=1:train_num
    if train_info(train_idx,2)==regular_train_ID
        for station_idx=1:length(train_path{train_idx})
            Aeqtemp=Aeqtemp+1;
            if train_info(train_idx,3)==downstream_ID
                location1=find(subscript_imk(:,1)==train_idx & subscript_imk(:,2)==train_path{train_idx}(station_idx) & mod(subscript_imk(:,3),2)==0);
            else
                location1=find(subscript_imk(:,1)==train_idx & subscript_imk(:,2)==train_path{train_idx}(station_idx) & mod(subscript_imk(:,3),2)==1);
            end
            Arng=vindex+location1;
            Acoe=ones(length(location1),1);
            Aeq(Aeqtem+(1:length(location1)),:)=[ones(length(location1),1)*Aeqtemp Arng Acoe];
            Aeqtem=Aeqtem+length(location1);
        end
    end
end

for train_idx=1:train_num
    if train_info(train_idx,2)==regular_train_ID
        for station_idx=1:length(train_path{train_idx})
            if train_station_dwell(train_idx,train_path{train_idx}(station_idx),6)<horizon_start_time...
                    || train_station_dwell(train_idx,train_path{train_idx}(station_idx),6)>horizon_end_time+24*3600
                Aeqtemp=Aeqtemp+1;
                location1=subscript_im_prenum(train_idx,train_path{train_idx}(station_idx));
                Arng=location1;
                Acoe=1;
                Aeq(Aeqtem+(1:1),:)=[ones(1,1)*Aeqtemp Arng Acoe];
                Aeqtem=Aeqtem+1;
                beq(Aeqtemp)=train_station_dwell(train_idx,train_path{train_idx}(station_idx),6);
            end
        end
    end
end

for train_idx=1:train_num
    if train_info(train_idx,2)==regular_train_ID
        for station_idx=1:length(train_path{train_idx})
            if train_station_dwell(train_idx,train_path{train_idx}(station_idx),7)<horizon_start_time...
                    || train_station_dwell(train_idx,train_path{train_idx}(station_idx),7)>horizon_end_time+24*3600
                Aeqtemp=Aeqtemp+1;
                location1=subscript_im_prenum(train_idx,train_path{train_idx}(station_idx));
                Arng=vnum_x+location1;
                Acoe=1;
                Aeq(Aeqtem+(1:1),:)=[ones(1,1)*Aeqtemp Arng Acoe];
                Aeqtem=Aeqtem+1;
                beq(Aeqtemp)=train_station_dwell(train_idx,train_path{train_idx}(station_idx),7);
            end
        end
    end
end
vindex1=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z+vnum_e;

for vindex=1:size(subscript_imk,1)
    if train_info(subscript_imk(vindex,1),2)==regular_train_ID
        if origin_arrive_time(subscript_imk(vindex,1),subscript_imk(vindex,2))<horizon_start_time...
                || origin_arrive_time(subscript_imk(vindex,1),subscript_imk(vindex,2))>horizon_end_time+24*3600
            Aeqtemp=Aeqtemp+1;
            Arng=vindex1+vindex;
            Acoe=1;
            Aeq(Aeqtem+(1:1),:)=[ones(1,1)*Aeqtemp Arng Acoe];
            Aeqtem=Aeqtem+1;
            beq(Aeqtemp)=origin_track_occupation(subscript_imk(vindex,1),subscript_imk(vindex,2),subscript_imk(vindex,3));
        end
    end
end

Aeq=sparse(Aeq(:,1),Aeq(:,2),Aeq(:,3),Aeqtemp,vnum);

lb=zeros(vnum,1);
ub=ones(vnum,1);
for vindex=1:size(subscript_im,1)
    lb(vindex)=train_station_dwell(subscript_im(vindex,1),subscript_im(vindex,2),2);
    ub(vindex)=train_station_dwell(subscript_im(vindex,1),subscript_im(vindex,2),3);
end

for vindex=1:size(subscript_im,1)
    lb(vnum_x+vindex)=train_station_dwell(subscript_im(vindex,1),subscript_im(vindex,2),4);
    ub(vnum_x+vindex)=train_station_dwell(subscript_im(vindex,1),subscript_im(vindex,2),5);
end

vindex1=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p;
for vindex=1:size(subscript_mn,1)
    lb(vindex1+vindex)=other_info(10);
    ub(vindex1+vindex)=other_info(11);
end

vindex=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z+vnum_e+vnum_v;
ub(vindex+(1:vnum_alpha))=(other_info(11)-other_info(10))*ones(vnum_alpha,1);

vindex=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z+vnum_e+vnum_v+vnum_alpha;
ub(vindex+(1:vnum_delta+vnum_epsilon))=max(max_earliness,max_tardiness)*ones(vnum_delta+vnum_epsilon,1);

c=zeros(vnum,1);
for train_idx=1:train_num
    if train_info(train_idx,2)==overnight_train_ID
        location1=subscript_im_prenum(train_idx,train_path{train_idx}(1));
        c(location1)=-1*objective_travel;
        location1=subscript_im_prenum(train_idx,train_path{train_idx}(end));
        c(vnum_x+location1)=objective_travel;
    end
end
vindex=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z+vnum_e+vnum_v+vnum_alpha;
for train_idx=1:train_num
    if train_info(train_idx,2)==regular_train_ID
        for station_idx=1:length(train_path{train_idx})
            location1=subscript_im_prenum(train_idx,train_path{train_idx}(station_idx));
            c(vindex+location1)=objective_planned_time;
            c(vindex+vnum_delta+location1)=objective_planned_time;
        end
    end
end
vindex=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z+vnum_e+vnum_v;
c(vindex+(1:vnum_alpha))=objective_maintenance*ones(vnum_alpha,1);
ctype=[repmat('C',1,vnum_x+vnum_y) repmat('I',1,vnum_w+vnum_u+vnum_p) repmat('C',1,vnum_z) repmat('B',1,vnum_e+vnum_v) repmat('C',1,vnum_alpha+vnum_delta+vnum_epsilon)];

cplex=Cplex();
cplex.Model.sense='minimize';
cplex.Model.obj=c;
cplex.Model.lb=lb;
cplex.Model.ub=ub;
cplex.Model.ctype=ctype;
cplex.Model.A=[Aineq;Aeq];
cplex.Model.lhs=[-Inf*ones(Atemp,1);beq];
cplex.Model.rhs=[bineq;beq];
cplex.Param.output.clonelog.Cur=0;
cplex.Param.mip.tolerances.mipgap.Cur=0.000001;
cplex.Param.timelimit.Cur=3600*4;
cplex.solve();

actual_arrive_time=[subscript_im cplex.Solution.x(1:vnum_x)];
actual_arrive_time=round(full(sparse(actual_arrive_time(:,1),actual_arrive_time(:,2),actual_arrive_time(:,3),train_num,station_num)));
actual_depart_time=[subscript_im cplex.Solution.x(vnum_x+(1:vnum_y))];
actual_depart_time=round(full(sparse(actual_depart_time(:,1),actual_depart_time(:,2),actual_depart_time(:,3),train_num,station_num)));
vindex=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p;
actual_maintenance_start=[subscript_mn cplex.Solution.x(vindex+(1:vnum_z))];
vindex=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z+vnum_e+vnum_v+vnum_alpha;
actual_arrive_time_deviation=[subscript_im cplex.Solution.x(vindex+(1:vnum_delta))];
actual_depart_time_deviation=[subscript_im cplex.Solution.x(vindex+vnum_delta+(1:vnum_epsilon))];
vindex=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z+vnum_e;
actual_track_occupation=[subscript_imk cplex.Solution.x(vindex+(1:vnum_v))];
objective_function_value=cplex.Solution.objval;
best_objective_function_value=cplex.Solution.bestobjval;

total_deviation=0;
for train_idx=1:train_num
    if train_info(train_idx,2)==regular_train_ID
        for station_idx=1:length(train_path{train_idx})
            total_deviation=total_deviation+abs(actual_arrive_time(train_idx,train_path{train_idx}(station_idx))-origin_arrive_time(train_idx,train_path{train_idx}(station_idx)));
            total_deviation=total_deviation+abs(actual_depart_time(train_idx,train_path{train_idx}(station_idx))-origin_depart_time(train_idx,train_path{train_idx}(station_idx)));
        end
    end
end

total_travelling_time=0;
for train_idx=1:train_num
    if train_info(train_idx,2)==overnight_train_ID
        total_travelling_time=total_travelling_time+actual_depart_time(train_idx,train_path{train_idx}(end))-actual_arrive_time(train_idx,train_path{train_idx}(1));
    end
end

vindex=vnum_x+vnum_y+vnum_w+vnum_u+vnum_p+vnum_z+vnum_e+vnum_v;
total_maintenance_deviation=sum(cplex.Solution.x(vindex+(1:vnum_alpha)));
optimized_solution=[objective_function_value;best_objective_function_value;total_deviation;total_travelling_time;total_maintenance_deviation];

outfile='output.xlsx';
sheetname='optimized_solution';
xlswrite(outfile,optimized_solution,sheetname,'A1');
sheetname='actual_arrive_time';
xlswrite(outfile,actual_arrive_time,sheetname,'A1');
sheetname='actual_depart_time';
xlswrite(outfile,actual_depart_time,sheetname,'A1');
sheetname='actual_maintenance_start';
xlswrite(outfile,actual_maintenance_start,sheetname,'A1');
sheetname='actual_arrive_time_deviation';
xlswrite(outfile,actual_arrive_time_deviation,sheetname,'A1');
sheetname='actual_depart_time_deviation';
xlswrite(outfile,actual_depart_time_deviation,sheetname,'A1');
sheetname='actual_track_occupation';
xlswrite(outfile,actual_track_occupation,sheetname,'A1');

set(gcf,'unit','centimeters','position',[0,0,16,8])
for vindex=1:size(actual_maintenance_start,1)
    if actual_maintenance_start(vindex,1)<actual_maintenance_start(vindex,2)
        rectangle('Position',[actual_maintenance_start(vindex,3),station_num-actual_maintenance_start(vindex,1),maintenance_duration(actual_maintenance_start(vindex,1),actual_maintenance_start(vindex,2)),1],'LineWidth',0.25,'LineStyle','none','FaceColor',[0.85,0.85,0.85]);
        hold on;
    end
end

for station_idx=1:size(station_info,1)
    x_alxe=[extend_start_time extend_end_time+24*3600];
    y_alxe=[station_idx station_idx];
    if station_idx==1 || station_idx==length(station_info)
        plot(x_alxe,y_alxe,'g','linewidth',0.5)
        hold on;
    else
        plot(x_alxe,y_alxe,'g','linewidth',0.25)
        hold on;
    end
    hold on;
end
for time_idx=extend_start_time:10*60:extend_end_time+24*3600
    x_alxe=[time_idx time_idx];
    y_alxe=[1 size(station_info,1)];
    if time_idx==extend_start_time || time_idx==extend_start_time
        plot(x_alxe,y_alxe,'g','linewidth',0.5)
        hold on;
    else if mod(time_idx,60*60)==0
            plot(x_alxe,y_alxe,'g','linewidth',0.5)
            hold on;
        else if mod(time_idx,30*60)==0
                plot(x_alxe,y_alxe,'--g','linewidth',0.25)
                hold on;
            else
                plot(x_alxe,y_alxe,'g','linewidth',0.1)
                hold on;
            end
        end
    end
end
for train_idx=1:train_num
    if strcmp(train_name{train_idx,4},'动车组')==1
        color_type='r';
    else if strcmp(train_name{train_idx,4},'确认车')==1
            color_type='m';
        else
            color_type='b';
        end
    end
    for section_idx=1:size(train_passed_section{train_idx},1)
        time1=actual_depart_time(train_idx,train_passed_section{train_idx}(section_idx,1));
        time2=actual_arrive_time(train_idx,train_passed_section{train_idx}(section_idx,2));
        position1=station_num+1-train_passed_section{train_idx}(section_idx,1);
        position2=station_num+1-train_passed_section{train_idx}(section_idx,2);
        if time1<=extend_start_time && time2>=extend_start_time
            if train_info(train_idx,3)==downstream_ID
                x_alxe=[extend_start_time time2];
                y_alxe=[min(position1,position2)+(time2-extend_start_time)/(time2-time1) position2];
            else
                x_alxe=[extend_start_time time2];
                y_alxe=[min(position1,position2)+(extend_start_time-time1)/(time2-time1) position2];
            end
            plot(x_alxe,y_alxe,color_type,'linewidth',0.1)
            hold on;
        elseif time1<=extend_end_time+24*3600 && time2>=extend_end_time+24*3600
            if train_info(train_idx,3)==downstream_ID
                x_alxe=[time1 extend_end_time+24*3600];
                y_alxe=[position1 min(position1,position2)+(time2-extend_end_time-24*3600)/(time2-time1)];
            else
                x_alxe=[time1 extend_end_time+24*3600];
                y_alxe=[position1 min(position1,position2)+(extend_end_time+24*3600-time1)/(time2-time1)];
            end
            plot(x_alxe,y_alxe,color_type,'linewidth',0.1)
            hold on;
            
        elseif ~(time1>extend_end_time+24*3600 || time2<extend_start_time)
            x_alxe=[time1 time2];
            y_alxe=[position1 position2];
            plot(x_alxe,y_alxe,color_type,'linewidth',0.1)
            hold on;
        end     
    end
end

for train_idx=1:train_num
    if strcmp(train_name{train_idx,4},'动车组')==1
        color_type='r';
    else if strcmp(train_name{train_idx,4},'确认车')==1
            color_type='m';
        else
            color_type='b';
        end
    end
    for station_idx=1:length(train_path{train_idx})
        time1=actual_arrive_time(train_idx,train_path{train_idx}(station_idx));
        time2=actual_depart_time(train_idx,train_path{train_idx}(station_idx));
        if time1<time2
            if time1<=extend_start_time && time2>=extend_start_time
                x_alxe=[extend_start_time time2];
                y_alxe=[station_num+1-train_path{train_idx}(station_idx) station_num+1-train_path{train_idx}(station_idx)];
                plot(x_alxe,y_alxe,color_type,'linewidth',0.1)
                hold on;
            elseif time1<=extend_end_time+24*3600 && time2>=extend_end_time+24*3600
                x_alxe=[time1 extend_end_time+24*3600];
                y_alxe=[station_num+1-train_path{train_idx}(station_idx) station_num+1-train_path{train_idx}(station_idx)];
                plot(x_alxe,y_alxe,color_type,'linewidth',0.1)
                hold on;
            elseif ~(time1>extend_end_time+24*3600 || time2<extend_start_time)
                x_alxe=[time1 time2];
                y_alxe=[station_num+1-train_path{train_idx}(station_idx) station_num+1-train_path{train_idx}(station_idx)];
                plot(x_alxe,y_alxe,color_type,'linewidth',0.1)
                hold on;
            end
        end
    end
end

for train_idx=1:train_num
    if strcmp(train_name{train_idx,4},'动车组')==1
        color_type='r';
    else if strcmp(train_name{train_idx,4},'确认车')==1
            color_type='m';
        else
            color_type='b';
        end
    end
    for section_idx=1:size(train_passed_section{train_idx},1)
        time1=actual_depart_time(train_idx,train_passed_section{train_idx}(section_idx,1));
        time2=actual_arrive_time(train_idx,train_passed_section{train_idx}(section_idx,2));
        position1=station_num+1-train_passed_section{train_idx}(section_idx,1);
        position2=station_num+1-train_passed_section{train_idx}(section_idx,2);
        if time1<=extend_start_time && time2>=extend_start_time
            if train_info(train_idx,3)==downstream_ID
                x_alxe=[extend_start_time time2];
                y_alxe=[min(position1,position2)+(time2-extend_start_time)/(time2-time1) position2];
                x=sum(x_alxe)/2;
                y=sum(y_alxe)/2;
                t1=(y_alxe(end)-y_alxe(1))/(station_num-1);
                t2=(x_alxe(end)-x_alxe(1))/(extend_end_time+24*3600-extend_start_time);
                angle_pi=atan2(t1,t2)*180/pi;
                text(x,y,num2str(2*train_idx-1),'Rotation',angle_pi,'FontSize',0.4,'FontName','Times New Roman','HorizontalAlignment','left','Color','k');
            else
                x_alxe=[extend_start_time time2];
                y_alxe=[min(position1,position2)+(extend_start_time-time1)/(time2-time1) position2];
                x=sum(x_alxe)/2;
                y=sum(y_alxe)/2;
                t1=(y_alxe(end)-y_alxe(1))/(station_num-1);
                t2=(x_alxe(end)-x_alxe(1))/(extend_end_time+24*3600-extend_start_time);
                angle_pi=atan2(t1,t2)*180/pi;
                text(x,y,num2str(2*(train_idx-size(down_arrive_time,1))),'Rotation',angle_pi,'FontSize',0.4,'FontName','Times New Roman','HorizontalAlignment','left','Color','c');
            end
        elseif time1<=extend_end_time+24*3600 && time2>=extend_end_time+24*3600
            if train_info(train_idx,3)==downstream_ID
                x_alxe=[time1 extend_end_time+24*3600];
                y_alxe=[position1 min(position1,position2)+(time2-extend_end_time-24*3600)/(time2-time1)];
                x=sum(x_alxe)/2;
                y=sum(y_alxe)/2;
                t1=(y_alxe(end)-y_alxe(1))/(station_num-1);
                t2=(x_alxe(end)-x_alxe(1))/(extend_end_time+24*3600-extend_start_time);
                angle_pi=atan2(t1,t2)*180/pi;
                text(x,y,num2str(2*train_idx-1),'Rotation',angle_pi,'FontSize',0.4,'FontName','Times New Roman','HorizontalAlignment','left','Color','k');
            else
                x_alxe=[time1 extend_end_time+24*3600];
                y_alxe=[position1 min(position1,position2)+(extend_end_time+24*3600-time1)/(time2-time1)];
                x=sum(x_alxe)/2;
                y=sum(y_alxe)/2;
                t1=(y_alxe(end)-y_alxe(1))/(station_num-1);
                t2=(x_alxe(end)-x_alxe(1))/(extend_end_time+24*3600-extend_start_time);
                angle_pi=atan2(t1,t2)*180/pi;
                text(x,y,num2str(2*(train_idx-size(down_arrive_time,1))),'Rotation',angle_pi,'FontSize',0.4,'FontName','Times New Roman','HorizontalAlignment','left','Color','c');
            end
        elseif ~(time1>extend_end_time+24*3600 || time2<extend_start_time)
            x_alxe=[time1 time2];
            y_alxe=[position1 position2];
            x=sum(x_alxe)/2;
            y=sum(y_alxe)/2;
            t1=(y_alxe(end)-y_alxe(1))/(station_num-1);
            t2=(x_alxe(end)-x_alxe(1))/(extend_end_time+24*3600-extend_start_time);
            angle_pi=atan2(t1,t2)*180/pi;
            if train_info(train_idx,3)==downstream_ID
                text(x,y,num2str(2*train_idx-1),'Rotation',angle_pi,'FontSize',0.4,'FontName','Times New Roman','HorizontalAlignment','left','Color','k');
            else
                text(x,y,num2str(2*(train_idx-size(down_arrive_time,1))),'Rotation',angle_pi,'FontSize',0.4,'FontName','Times New Roman','HorizontalAlignment','left','Color','c');
            end
        end     
    end
end

for station_idx=1:length(station_name)
    x=extend_end_time+24*60*60+0.2*60*60;
    y=length(station_name)-station_idx+1;
    text(x,y,num2str(station_idx),'Rotation',0,'FontSize',4,'FontName','Times New Roman','HorizontalAlignment','left');
    
    x=extend_start_time-0.2*60*60;
    y=length(station_name)-station_idx+1;
    text(x,y,num2str(station_idx),'Rotation',0,'FontSize',4,'FontName','Times New Roman','HorizontalAlignment','right');
end
for time_idx=extend_start_time/3600:1:extend_end_time/3600+24
    x=time_idx*3600;
    y=length(station_info)+1;
    if time_idx<24
        text(x,y,[num2str(time_idx) ':00'],'Rotation',0,'FontSize',4,'FontName','Times New Roman','HorizontalAlignment','center');
    else
        text(x,y,['0' num2str(time_idx-24) ':00'],'Rotation',0,'FontSize',4,'FontName','Times New Roman','HorizontalAlignment','center');
    end
    x=time_idx*3600;
    y=0;
    if time_idx<24
        text(x,y,[num2str(time_idx) ':00'],'Rotation',0,'FontSize',4,'FontName','Times New Roman','HorizontalAlignment','center');
    else
        text(x,y,['0' num2str(time_idx-24) ':00'],'Rotation',0,'FontSize',4,'FontName','Times New Roman','HorizontalAlignment','center');
    end
end

diary off;
