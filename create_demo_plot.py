import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Create some demo data with multiple degrees and braid lengths
demo_data = [
    {'Braid_Length': 8, 'Degree': 3, 'Original_Time': 0.5, 'Optimized_Time': 0.4},
    {'Braid_Length': 8, 'Degree': 4, 'Original_Time': 0.8, 'Optimized_Time': 0.6},
    {'Braid_Length': 8, 'Degree': 5, 'Original_Time': 0.96, 'Optimized_Time': 2.08},
    {'Braid_Length': 9, 'Degree': 3, 'Original_Time': 0.7, 'Optimized_Time': 0.5},
    {'Braid_Length': 9, 'Degree': 4, 'Original_Time': 1.1, 'Optimized_Time': 0.9},
    {'Braid_Length': 9, 'Degree': 5, 'Original_Time': 1.44, 'Optimized_Time': 1.35},
    {'Braid_Length': 10, 'Degree': 3, 'Original_Time': 0.9, 'Optimized_Time': 0.7},
    {'Braid_Length': 10, 'Degree': 4, 'Original_Time': 1.3, 'Optimized_Time': 1.1},
    {'Braid_Length': 10, 'Degree': 5, 'Original_Time': 1.34, 'Optimized_Time': 1.88},
    {'Braid_Length': 11, 'Degree': 3, 'Original_Time': 1.1, 'Optimized_Time': 0.9},
    {'Braid_Length': 11, 'Degree': 4, 'Original_Time': 1.5, 'Optimized_Time': 1.3},
    {'Braid_Length': 11, 'Degree': 5, 'Original_Time': 1.37, 'Optimized_Time': 1.90},
    {'Braid_Length': 12, 'Degree': 3, 'Original_Time': 1.3, 'Optimized_Time': 1.1},
    {'Braid_Length': 12, 'Degree': 4, 'Original_Time': 1.7, 'Optimized_Time': 1.5},
    {'Braid_Length': 12, 'Degree': 5, 'Original_Time': 1.92, 'Optimized_Time': 1.90},
]

def create_surface_plot(results, output_filename="demo_computation_time_difference_plot.png"):
    """Create a surface plot showing average computation time differences"""
    print("Creating surface plot...")
    
    # Group results by braid length and degree
    grouped_data = {}
    for result in results:
        key = (result['Braid_Length'], result['Degree'])
        if key not in grouped_data:
            grouped_data[key] = []
        time_diff = result['Original_Time'] - result['Optimized_Time']
        grouped_data[key].append(time_diff)
    
    # Calculate average differences
    avg_data = {}
    for key, values in grouped_data.items():
        avg_data[key] = sum(values) / len(values)
    
    # Extract unique braid lengths and degrees
    braid_lengths = sorted(set(key[0] for key in avg_data.keys()))
    degrees = sorted(set(key[1] for key in avg_data.keys()))
    
    # Create meshgrid
    X, Y = np.meshgrid(degrees, braid_lengths)
    Z = np.zeros_like(X, dtype=float)
    
    # Fill Z matrix
    for i, bl in enumerate(braid_lengths):
        for j, deg in enumerate(degrees):
            if (bl, deg) in avg_data:
                Z[i, j] = avg_data[(bl, deg)]
            else:
                Z[i, j] = 0
    
    # Create 3D surface plot
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    surface = ax.plot_surface(X, Y, Z, cmap='RdYlBu', alpha=0.8)
    
    ax.set_xlabel('Degree')
    ax.set_ylabel('Braid Length')
    ax.set_zlabel('Time Difference (Original - Optimized) [seconds]')
    ax.set_title('Performance Comparison: Original vs Optimized Implementation\n(Positive values = optimization is faster)')
    
    # Add colorbar
    fig.colorbar(surface, shrink=0.5, aspect=5)
    
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Surface plot saved to {output_filename}")
    plt.close()

if __name__ == "__main__":
    create_surface_plot(demo_data)
    print("Demo surface plot created successfully!")