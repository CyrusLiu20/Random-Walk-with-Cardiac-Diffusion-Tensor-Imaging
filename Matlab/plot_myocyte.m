clear;
loadfile;
clf(figure(2));
figure(2);
plotGeometry(myocytes);
% Axis
xlabel("X-axis");
ylabel("Y-axis");
zlabel("Z-axis");
for i_particle = 1:size(histories, 2)
%     if valid(i_particle)
        plot3(histories(i_particle).position(:,1), histories(i_particle).position(:,2), histories(i_particle).position(:,3), linewidth = 0.5);
        hold on;
%     end
end

%Plot Geometry
function plotGeometry(geometry)
    hold on
    for i=1:length(geometry)
        vertices = double(geometry(i).Vertices);
        faces = double(geometry(i).Faces);
        poly_tri = triangulation(faces,vertices);
        tri_1 = trisurf(poly_tri, 'FaceAlpha', 0.5, 'FaceColor',"#FF8080", 'EdgeColor', 'none');
        tri_1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end