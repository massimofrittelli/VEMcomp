clearvars
close all

[P, h, boundarynode, EGamma, Elements] = plot_mesh_step_3_new(33);

file = fopen('mesh4.voro','w');

boundaryFaceCounter = 0;
internalFaceCounter = 0;
for i=1:length(Elements)
    Element = Elements(i);
    Pind = Element.Pind;
    npoints = length(Element.P);
    nfaces = length(Element.Faces);
    fprintf(file, 'VoronoiCell %d\n%d\n', i-1, npoints);
    for j=1:npoints
        fprintf(file, '%.16f %.16f %.16f\n', Element.P(j,1), Element.P(j,2), Element.P(j,3));
    end
    fprintf(file,'\n%d\n', nfaces);
    for j=1:nfaces
        [~,p1,p2] = intersect(Pind, Element.Faces(j).Pind,'stable');
        Pind_loc = p1(p2) - 1;
        for k=1:length(Pind_loc)
            % ho cambiato l'orientamento delle facce
           fprintf(file, '%d ', Pind_loc(length(Pind_loc)-k+1)); 
        end
        fprintf(file, '#\n');
    end
    fprintf(file, '\n');
    for j=1:nfaces
       boundaryFaceCounter = boundaryFaceCounter - Element.Faces(j).boundary;
       internalFaceCounter = internalFaceCounter + (1- Element.Faces(j).boundary);
       marker = Element.Faces(j).boundary*boundaryFaceCounter + (1- Element.Faces(j).boundary)*internalFaceCounter;
       fprintf(file, '%d ', marker); 
    end
    fprintf(file,'\n\n');
end

fclose(file);