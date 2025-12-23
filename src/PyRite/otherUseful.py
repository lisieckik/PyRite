import numpy as np

def is_pixel_in_ellipse(image_size, center, a, b, theta, scale = 1, ap = 0):
    # Generate the grid of coordinates
    a = a*scale +ap
    b = b*scale +ap
    y, x = np.indices(image_size)
    theta = np.deg2rad(theta)

    # Shift the coordinates to the ellipse center
    x_shifted = x - center[0]
    y_shifted = y - center[1]

    # Apply the rotation matrix to the shifted coordinates
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    x_rot = x_shifted * cos_theta + y_shifted * sin_theta
    y_rot = -x_shifted * sin_theta + y_shifted * cos_theta

    # Ellipse equation (x_rot / a)^2 + (y_rot / b)^2 <= 1
    ellipse_mask = ((x_rot / a) ** 2 + (y_rot / b) ** 2) <= 1
    return ellipse_mask
