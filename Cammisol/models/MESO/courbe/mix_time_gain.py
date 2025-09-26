import matplotlib.pyplot as plt
import numpy as np
import os

script_dir = os.path.dirname(os.path.abspath(__file__))

# Filtered Data
processors = np.array([1, 2, 4, 8, 16, 32, 64])
average_time = np.array([345.4, 216.3, 106.9, 74.9, 65.7, 42.4, 39])

# Compute performance gain
performance_gain = (average_time[0] - average_time) / average_time[0] * 100

# --- Chart 1: Average Time and Performance Gain ---
fig1, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(processors, average_time, marker='o', linestyle='--', color='b', markersize=8, linewidth=2, label="Average Time (s)")
ax1.set_xscale('log', base=2)
ax1.set_xticks(processors)
ax1.set_xticklabels(processors, rotation=45, ha="right")
ax1.set_xlabel("Number of Processors")
ax1.set_ylabel("Average Time (s)", color='b')
ax1.set_title("Performance Analysis: Average Time and Performance Gain")
ax1.grid(True, linestyle=":", linewidth=0.7, alpha=0.7)

# Create a secondary y-axis for performance gain
ax2 = ax1.twinx()
ax2.plot(processors, performance_gain, marker='s', linestyle='-', color='r', markersize=6, linewidth=1.5, label="Performance Gain (%)")
ax2.set_ylabel("Performance Gain (%)", color='r')

# Combine legends from both axes
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc='center right', bbox_to_anchor=(1.0, 0.5))

# Save and show Chart 1
plt.savefig(os.path.join(script_dir, "performance_time_gain.png"), dpi=300, bbox_inches='tight')
plt.show()