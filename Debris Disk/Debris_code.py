import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Set up the figure and axis - 800x800 AU grid in km --> 1.2*10^11 x 1.2*10^11 km
fig, ax = plt.subplots()
ax.set_xlim(-6*(10**10), 6*(10**10))
ax.set_ylim(-6*(10**10), 6*(10**10))
ax.set_aspect('equal')

# Number of dust particles
num_particles = 100000
particles, = ax.plot([], [], 'bo', markersize=0.5)  # Blue dust particles

# Common starting parameters for all particles
a = 9*(10**9)  # Semi-major axis (same for all particles) 60 AU
e = np.zeros(num_particles)  # Eccentricity 0 initially


# Generate random true anomalies for each particle
f = np.random.rand(num_particles) * 2 * np.pi   # True anomaly (uniform distribution between 0 and 2π)

# Function to calculate the max value of beta for each particle
def calculate_beta_max(e, f):
    return ((1 - e**2) / (2 * (1 + e * np.cos(f)))) * 0.997

# Calculate beta_max for each particle
beta_max = calculate_beta_max(e, f)

# Generate random beta values within the specified range for each particle using a uniform distribution
beta_min = 0.001
beta = np.random.uniform(beta_min, beta_max, num_particles)

# Function to calculate the new eccentricity for each particle
def calculate_e_prime(e, beta, f):
    return np.sqrt(e**2 + 2 * beta * e * np.cos(f) + beta**2) / (1 - beta)

# Function to calculate the new radius for each particle
def calculate_radius(a, e, f):
    return ((a*(1 - e**2 )*(1 - beta)) / (1 - e**2 - 2*beta*(1 + e*np.cos(f))))


# Generate random velocities for each particle with a mean of 4 and variance of 0.1 (km/s)
velocities = np.random.normal(4, 0.1, num_particles)

# Function to update the animation
def update(frame):
    theta = frame * velocities * 0.01  # Update angles for orbiting with different velocities
    e_prime = calculate_e_prime(e, beta, f)  # Calculate modified eccentricity
    r = calculate_radius(a, e, f) * (1 - e_prime**2) / (1 + e_prime * np.cos(theta))
    y = -1 * r * np.cos(theta)
    x = r * np.sin(theta)
    print(x, y)
    particles.set_data(x, y)
    return particles,

# Create the animation
ani = FuncAnimation(fig, update, frames=5000, interval=50, blit=True)

plt.show()

# Model uses random beta values based on uniform normal distribution as opposed to Dohnanyi size distribution
