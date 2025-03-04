import numpy as np
import matplotlib.pyplot as plt

# Load array
image_data = np.load("ST/images/raw_array_slice1_BN9966.npy")

# Display image
plt.imshow(image_data)
plt.axis('off')
plt.show()
