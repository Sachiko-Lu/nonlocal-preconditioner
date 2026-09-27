function out = geometry_triangulate_polytope4d(vertices,opts)

if nargin < 2 || isempty(opts)
    opts = struct();
end
validateattributes(vertices, {'numeric'}, {'2d','real','finite','ncols',4}, ...
    mfilename, 'vertices');
rankTol = local_option(opts,'rankTolerance',2e-12);
volumeTol = local_option(opts,'volumeTolerance',2e-8);

out = local_empty_output(vertices);
if size(vertices,1) < 5
    return
end
affineRank = rank(vertices-mean(vertices,1), ...
    rankTol*max(1,norm(vertices,'fro')));
if affineRank < 4
    return
end

try
    [facets,hullVolume] = convhulln(vertices,{'Qt','Qx'});
catch firstFailure
    try
        [facets,hullVolume] = convhulln(vertices,{'Qt','QJ'});
    catch
        rethrow(firstFailure)
    end
end
facets = unique(sort(facets,2),'rows','stable');
center = mean(vertices,1);

if size(vertices,1) == 5
    simplexVertices = reshape(vertices,1,5,4);
    simplexVolumes = local_simplex_volume(vertices);
else
    nFacets = size(facets,1);
    simplexVertices = zeros(nFacets,5,4);
    simplexVolumes = zeros(nFacets,1);
    for k = 1:nFacets
        R = [center; vertices(facets(k,:),:)];
        simplexVertices(k,:,:) = R;
        simplexVolumes(k) = local_simplex_volume(R);
    end
    keepThreshold = volumeTol*max(hullVolume,realmin)/max(1,nFacets);
    keep = simplexVolumes > keepThreshold;
    simplexVertices = simplexVertices(keep,:,:);
    simplexVolumes = simplexVolumes(keep);
end

sumVolume = sum(simplexVolumes);
relativeError = abs(sumVolume-hullVolume)/max(hullVolume,realmin);
if relativeError > volumeTol
    error('geometry:TriangulationVolumeMismatch', ...
        ['Centroid-cone simplex volume differs from the convex-hull ', ...
         'volume by %.3e (allowed %.3e).'],relativeError,volumeTol);
end

out.vertices = vertices;
out.facets = facets;
out.center = center;
out.simplexVertices = simplexVertices;
out.simplexVolumes = simplexVolumes;
out.volume = sumVolume;
out.hullVolume = hullVolume;
out.volumeRelativeError = relativeError;
out.isFullDimensional = true;
end

function volume = local_simplex_volume(R)
volume = abs(det(R(2:5,:)-R(1,:)))/factorial(4);
end

function out = local_empty_output(vertices)
out = struct('vertices',vertices,'facets',zeros(0,4), ...
    'center',zeros(0,4),'simplexVertices',zeros(0,5,4), ...
    'simplexVolumes',zeros(0,1),'volume',0,'hullVolume',0, ...
    'volumeRelativeError',0,'isFullDimensional',false);
end

function value = local_option(opts,name,defaultValue)
if isfield(opts,name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = defaultValue;
end
validateattributes(value, {'numeric'}, ...
    {'scalar','real','finite','positive'}, mfilename, name);
end
