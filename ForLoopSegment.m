load 'Experimental Files\AllSurfaceVolumes.mat'
TotalNumberOfSurfaces = size(AllSurfaceVolumeData,1);
tSurface = transpose(0:TotalNumberOfSurfaces);
tNodes = zeros(TotalNumberOfSurfaces,1);
tEndPoints = zeros(TotalNumberOfSurfaces,1);
parfor idx = 1:TotalNumberOfSurfaces
    
%     vDataItem = vSurface.GetSingleMask(idx - 1,0,0,0,600,600,600,300,300,300); %% IMDIALTE IMAGE??
% %     SurfaceData(idx,1) = vDataItem;
% %     vDataItem = SurfaceData{idx};
%     FloatData = vDataItem.GetDataVolumeFloats(0,0);
    logicalFloats = AllSurfaceVolumeData{idx};
    % isosurface(allMaskData,1/2)
% axis([0 600 0 600 0 64])
    
    %% Dilation of data
    se = strel('cuboid',[2 1 1]);
    dilatedImage = imdilate(logicalFloats,se);
    erodedImage = imerode(dilatedImage,se);
    
    %%
    %Insert error catching here?
    try 
        skel = Skeleton3D(erodedImage);

        w = size(skel,1);
        l = size(skel,2);
        h = size(skel,3);

        % initial step: condense, convert to voxels and back, detect cells
        [~,node,link] = Skel2Graph3D(skel,0);

        % total length of network
        wl = sum(cellfun('length',{node.links}));
        
        % converts the network graph back into a cleaned-up voxel skeleton image
        skel2 = Graph2Skel3D(node,link,w,l,h);
        
        if isempty(find(skel2,1))
            
        end

        [~,node2,link2] = Skel2Graph3D(skel2,0); %WILL CAUSE ERROR IF NO SKELETON IS AVAILABLE DUE TO LOW BRANCHING COMPLEXITY

        % calculate new total length of network
        wl_new = sum(cellfun('length',{node2.links}));

        % iterate the same steps until network length changed by less than 0.5%
        while(wl_new~=wl)

            wl = wl_new;   

             skel2 = Graph2Skel3D(node2,link2,w,l,h);
             [A2,node2,link2] = Skel2Graph3D(skel2,0);

             wl_new = sum(cellfun('length',{node2.links}));
        end;
        
        % view cell
        showCellSkeleton(erodedImage,link2,node2)
%         fileName = sprintf('%s_%d','Images\Surface',idx);
%         saveas(gcf,fileName);
        vNode2ep = [node2.ep];
        Nodes = vNode2ep == 0;
        EndPoints = vNode2ep == 1;
        tNodes(idx,1) = sum(Nodes(:));
        tEndPoints(idx,1) = sum(EndPoints(:));
        
    catch ME
        
        if isempty(find(skel2,1))
            
            try
                
                showCellSkeleton(erodedImage,link,node)
                vNodeep = [node.ep];
                Nodes = vNodeep == 0;
                EndPoints = vNodeep == 1;
                tNodes(idx,1) = sum(Nodes(:));
                tEndPoints(idx,1) = sum(EndPoints(:));
                
            catch ME
               
                vNodeep = [node.ep];
                Nodes = vNodeep == 0;
                EndPoints = vNodeep == 1;
                tNodes(idx,1) = sum(Nodes(:));
                tEndPoints(idx,1) = sum(EndPoints(:));
                
            end
            
        elseif (strcmp('MATLAB:badsubscript',ME.identifier)) && ~isempty(find(skel2,1))%SET THE ERROR CATCH WITHIN THE SHOWCELLSKELETON.M
            
            vNode2ep = [node2.ep];
            Nodes = vNode2ep == 0;
            EndPoints = vNode2ep == 1;
            tNodes(idx,1) = sum(Nodes(:));
            tEndPoints(idx,1) = sum(EndPoints(:));
            
        else
            
            vNode2ep = [node2.ep];
            Nodes = vNode2ep == 0;
            EndPoints = vNode2ep == 1;
            tNodes(idx,1) = sum(Nodes(:));
            tEndPoints(idx,1) = sum(EndPoints(:));
            errordlg(ME.message);
            uiwait; 
            rethrow(ME);
            
        end
        
   
    
    end
set(gcf,'NumberTitle','off') %don't show the figure number
titleName = sprintf('%s_%d','Surface ',idx - 1);
set(gcf,'Name',titleName);
fileName = sprintf('%s_%d','Images\Surface',idx-1);
saveas(gcf,fileName);  

end
%Annotate each fig so we know which surface it is from
%Save each figure to table.
dataTable = table(tSurface,tNodes,tEndPoints)
writetable(dataTable,'Images\BranchAnalysis.xlsx');
