%%%%%%%%%

% N Z A % 

% grid cells: 228
% pairs: 749


%%%%%%%%%


function find_pairs_of_cells_mosimol


params.dir_load =...
    'C:\Noam\muscimol\DB_musc_MEC\';
    %'C:\Users\noam\Desktop\proj\muscimol\DB_musc_MEC\';

params.dir_save =....
    'C:\Noam\muscimol\noam_output\';
    %'C:\Users\noam\Desktop\proj\muscimol\noam_output';
%
tic
 files = dir(strcat(params.dir_load,'DB*.mat'));
 for i=1:length(files) 
    data{i} = load(strcat(params.dir_load,files(i).name));
    fprintf('%.1f \n',(i*100/length(files)));
 end
 toc
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data = load(strcat(params.dir_load,'ALL_CELLS.mat')); 
%data = data.data;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

grid_cells = 0;
for i=1:length(data) %a = load(strcat(params.dir_load,files(i).name));
    a = data(i); a = a{1}.db;
    if (strcmp(a.cell_type_subj, 'GRID') || strcmp(a.cell_type_subj, 'CONJ')) % '_GRID' ?????      
        for j=i+1:length(data)
            b = data(j); b = b{1}.db;
            if ( strcmp(a.rat, b.rat) && strcmp(a.date, b.date)...
                    && (strcmp(a.cell_type_subj, 'GRID') || strcmp(a.cell_type_subj, 'CONJ'))  ...
                    && strcmp(a.area, b.area)  ...  %area????
                    && ( a.tetrode ~= b.tetrode || a.cell ~= b.cell) ...
                )
                %fprintf('same rat:%s date:%s a:%d,%d b:%d,%d \n', a.rat, ...
                    %a.date, a.tetrode, a.cell, b.tetrode, b.cell);
                
                grid_cells = grid_cells + 1;
                pairs_by_file_index(grid_cells,1) = i; 
                pairs_by_file_index(grid_cells,2) = j; 
            end
        end
    else
        disp (a.cell_type_subj);
    end
end
disp(grid_cells);
disp('fin!');

%spike plot
%%
a = data{1}.db.B(2) % B session is combination of all recordings. has 3 structs: before, during, after muscomol injection
plot( mean([a.pos_data.x1, a.pos_data.x2] ,2) , mean([a.pos_data.y1, a.pos_data.y2] ,2) );
hold on;
plot( mean([a.spike_data.x, a.spike_data.x2] ,2) , mean([a.spike_data.y, a.spike_data.y2] ,2), 'r.' );

c1.x =  data{1}.db.B(2).mean([a.pos_data.x1, a.pos_data.x2], 2); 

%%
for i=1:length(pairs_by_file_index)
    
end


            %{
            if strcmp(rat{i},rat{j}) && strcmp(date{i},date{j}) && strcmp(session{i},session{j})...
                    && (tetrode(i)~=tetrode(j) || cell(i)~=cell(j)) && ...
                    (p_gridness(i)<0.05  && p_gridness(j)<0.05)
            %} 

end

            