#include <utilities.h>

void powerSpacedGrid(double low, double high, double curv, grid_type& grid) {
	int n = grid.size();
	linspace(0.0, 1.0, n, grid);

	for (int i=0; i<n; ++i) {
		grid[i] = low + (high - low) * pow(grid[i], 1.0 / curv);
	}
}

void adjustPowerSpacedGrid(grid_type& grid) {
	if (grid.size() >= 10)
		for (int i=0; i<9; ++i)
			grid[i] = (i - 1) * grid[9] / (10.0 - 1.0);
}
