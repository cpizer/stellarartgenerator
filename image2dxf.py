import numpy as np
import cv2 as cv
import ezdxf

#User inputs down below
covariance_radius_px = 3
points_radius_mm = 0.05
#End of user inputs

im_nparray = cv.imread("m101.png")
im_gray_nparray = cv.cvtColor(im_nparray, cv.COLOR_BGR2GRAY)

cartesian_points = []
rand_mean = [0, 0]
cov = [[covariance_radius_px, 0], [0, covariance_radius_px]]
for x in range(im_gray_nparray.shape[0]):
    for y in range(im_gray_nparray.shape[1]):
        points_x, points_y = np.random.multivariate_normal(rand_mean, cov, 1000).T
        for i in range(len(points_x)):
            cartesian_points.append([points_x[i], points_y[i]])

#Write the data to a dxf-File
doc = ezdxf.new("R2010")
msp = doc.modelspace()
for point in cartesian_points:
    msp.add_circle((point[0], point[1]), points_radius_mm)
doc.saveas("image_generated_pointcloud.dxf")