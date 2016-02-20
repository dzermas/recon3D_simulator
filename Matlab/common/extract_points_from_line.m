function xy = extract_points_from_line(x_point_edges, line_params, npoints)

x = linspace(x_point_edges(1), x_point_edges(2), npoints);
y = line_params(1)*x + line_params(2);

xy = [x' y'];