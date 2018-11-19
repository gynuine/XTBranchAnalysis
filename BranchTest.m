
%% Extract Vantage Spot Counts
%
%    <CustomTools>
%      <Menu>
%        <Submenu name="Custom Scriptss">
%        <Item name="Extract Vantage Spot Count" icon="Matlab" tooltip="Extracts spot count from Vantage.xls file.">
%          <Command>MatlabXT::XTExtractData(%i)</Command>
%        </Item>
%        </Submenu>
%      </Menu>
%    </CustomTools>
%
%% Imaris Bridge

function BranchTest(aImarisApplicationID)

function anObjectID = GetObjectID

        vServer = vImarisLib.GetServer;
        vNumberOfObjects = vServer.GetNumberOfObjects;
        anObjectID = vServer.GetObjectID(vNumberOfObjects - 1);

end

if ~isa(aImarisApplicationID, 'Imaris.IApplicationPrxHelper')
    
    javaaddpath ImarisLib.jar;
    vImarisLib = ImarisLib;
    vImarisApplication = vImarisLib.GetApplication(GetObjectID);
    
else
    
    vImarisApplication = aImarisApplicationID;
    
end

vCurrentSurpassSelection = vImarisApplication.GetSurpassSelection;
vSurface = vImarisApplication.GetFactory.ToSurfaces(vCurrentSurpassSelection);
vSurpassScene = vImarisApplication.GetSurpassScene;

if isequal(vSurpassScene, [])
  errordlg('Please load an image','Error!','modal');
  return;
end

TotalNumberOfSurfaces = vSurface.GetNumberOfSurfaces();
SurfaceData = cell(TotalNumberOfSurfaces,1);

tSurface = transpose(1:TotalNumberOfSurfaces);
tNodes = zeros(TotalNumberOfSurfaces,1);
tEndPoints = zeros(TotalNumberOfSurfaces,1);

%% Which voxel value is best? Optimised for speed and accuracy
% vDataItem1 = vSurface.GetSingleMask(1,0,0,0,512,512,512,250,250,100);
% FloatData1 = vDataItem1.GetDataVolumeFloats(0,0);
% vDataItem2 = vSurface.GetSingleMask(1,0,0,0,512,512,512,250,250,150);%%Best option
% FloatData2 = vDataItem2.GetDataVolumeFloats(0,0);
% vDataItem3 = vSurface.GetSingleMask(1,0,0,0,512,512,512,300,300,150);
% FloatData3 = vDataItem3.GetDataVolumeFloats(0,0);
% vDataItem4 = vSurface.GetSingleMask(1,0,0,0,512,512,512,300,300,300);
% FloatData4 = vDataItem4.GetDataVolumeFloats(0,0);
% subplot(2,2,1); isosurface(FloatData1,1/2);
% subplot(2,2,2); isosurface(FloatData2,1/2);
% subplot(2,2,3); isosurface(FloatData3,1/2);
% subplot(2,2,4); isosurface(FloatData4,1/2);

%% Extract Single surface volume from imaris and append into cell array

AllSurfaceVolumeData = cell(TotalNumberOfSurfaces,1);
for surfaceIndex = 1:TotalNumberOfSurfaces
    vDataItem = vSurface.GetSingleMask(surfaceIndex - 1 ,0,0,0,512,512,512,250,250,150);
    FloatData = vDataItem.GetDataVolumeFloats(0,0);
    AllSurfaceVolumeData{surfaceIndex,1} = FloatData;
end

%% Show surfacedata
% maskData = vSurface.GetMask(0,0,0,512,512,512,600,600,64,0);
% allMaskData = maskData.GetDataVolumeFloats(0,0);
% isosurface(allMaskData,1/2)
% axis([0 600 0 600 0 64])

%% Image skeletonization & clean up.
for idx = 1:TotalNumberOfSurfaces
    
%     vDataItem = vSurface.GetSingleMask(idx - 1,0,0,0,600,600,600,300,300,300); %% IMDIALTE IMAGE??
% %     SurfaceData(idx,1) = vDataItem;
% %     vDataItem = SurfaceData{idx};
%     FloatData = vDataItem.GetDataVolumeFloats(0,0);
    logicalFloats = AllSurfaceVolumeData{idx};
    
    
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
        fileName = sprintf('%s_%d','Images\Surface',idx);
        saveas(gcf,fileName);
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
end
%Annotate each fig so we know which surface it is from
%Save each figure to table.
dataTable = table(tSurface,tNodes,tEndPoints)

end