%Plot Geometry
plotGeometry(myocytes)

function plotGeometry(geometry)
hold on
% plotting every myocytes
for i=1:length(geometry)
    vertices = double(geometry(i).Vertices);
    faces = double(geometry(i).Faces);
    poly_tri = triangulation(faces,vertices);
    tri_1 = trisurf(poly_tri, 'FaceColor',"#FF8080", 'EdgeColor', 'none');
    tri_1.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
end