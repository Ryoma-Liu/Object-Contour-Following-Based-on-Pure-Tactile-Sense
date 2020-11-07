function clearMeshes(manip)
%clearMeshes Function to reduce the Sawyer mesh size to ensure faster visualization

% Copyright 2018 The MathWorks, Inc.

% hideMeshes = {'controller_box', 'pedestal_feet', 'pedestal'};
hideMeshes = {'pedestal'};
% hideMeshes = {'controller_box', 'pedestal_feet', 'pedestal', 'right_arm_base_link'};

for i = 1:length(hideMeshes)

    body = getBody(manip, hideMeshes{i});

    if ~isempty(body.Visuals)
        clearVisual(body);
    end

end

end