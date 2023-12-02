import matplotlib.pyplot as plt

# Create a figure
fig, ax = plt.subplots()

# Define the y-axis positions for each row of data
y_pos = range(len(df.index))

# Loop through each row of data
for i, row in df.iterrows():
    # Calculate the x-axis positions of the start and end points
    x_start = row['start']
    x_end = row['end']
    
    # Calculate the height of the block
    block_height = 1
    
    # Draw the block
    ax.barh(y_pos[i], x_end - x_start, left=x_start, height=block_height)
    
# Set the y-axis tick labels to be the gene names
ax.set_yticks(y_pos)
ax.set_yticklabels(df['Gene'])

# Set the x-axis limits
ax.set_xlim(0, df['end'].max())

# Invert the y-axis so that the top row is at the top of the plot
ax.invert_yaxis()

# Show the plot
plt.show()

