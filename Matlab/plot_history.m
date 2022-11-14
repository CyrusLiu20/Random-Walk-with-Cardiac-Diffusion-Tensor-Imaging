clear
loadfile;
clf(figure(1));
figure(1);

% Which particle we want to pot
i_particle = 27;
% How many particles at a given time frame
lasting = 30;

% Domain
xmin = min(histories(i_particle).position(:, 1));
xmax = max(histories(i_particle).position(:, 1));
ymin = min(histories(i_particle).position(:, 2));
ymax = max(histories(i_particle).position(:, 2));
zmin = min(histories(i_particle).position(:, 3));
zmax = max(histories(i_particle).position(:, 3));

% Plot every time step
for step = 1:15:size(histories(i_particle).position,1)-lasting
    hold off;
    scatter3(histories(i_particle).position(step:step+lasting-1,1), histories(i_particle).position(step:step+lasting-1,2), histories(i_particle).position(step:step+lasting-1,3), 10, "filled", "blue");
    hold on;
    scatter3(histories(i_particle).position(step+lasting,1), histories(i_particle).position(step+lasting,2), histories(i_particle).position(step+lasting,3), 10, "filled", "red");
    hold on;
    legend(["Time step : ", num2str(step+lasting)], Location="best")

    % Axis
    xlabel("X-axis");
    ylabel("Y-axis");
    zlabel("Z-axis");

    % Limits
    xlim([xmin, xmax]);
    ylim([ymin, ymax]);
    zlim([zmin, zmax]);

    hold on;
    pause(0.0001);

end
hold off;
plot3(histories(i_particle).position(:,1), histories(i_particle).position(:,2), histories(i_particle).position(:,3), linewidth = 0.5);

% Axis
xlabel("X-axis");
ylabel("Y-axis");
zlabel("Z-axis");

Limits
xlim([xmin, xmax]);
ylim([ymin, ymax]);
zlim([zmin, zmax]);