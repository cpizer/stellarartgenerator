import numpy as np
import astropy.coordinates as coord
import gaia.tap
import numpy as np
from astropy.table import Table
import ezdxf
from astropy import units as u
from astropy.coordinates import Angle
from astroquery.simbad import Simbad

#User inputs!!!
result_table = Simbad.query_object("gam Cas") #Name des Objekts (Existenz evtl. auf https://simbad.u-strasbg.fr/simbad/sim-fid kontrollieren)
no_of_sources = 500 #Anzahl der Sterne in der dxf-Datei
roi_diameter_deg = 10 #Sternbild ca. 10-20 Grad, Kugelsternhaufen 0.1-0.25 Grad
no_layers = 5 #Die Anzahl der Layer, welchen die Sterne je nach Helligkeit zugeordnet werden
create_circles = True #Auf False setzen, falls anstatt Kreisen nur Punkte gezeichnet werden sollen
max_radius_dxf = 1 #Der Radius des hellsten Sterns/der hellsten Quelle. Die Querschnittsflaeche aller weiteren Sterne wird auf die Querschnittsflaeche dieses Sterns normiert
canvas_height_mm = 240  #Die Hoehe der Leinwand
canvas_width_mm = 300   #Die Breite der Leinwand
#End of user inputs!!!

roi_center_ra = Angle(result_table['RA'][0] + 'h').to(u.degree).to_value()
roi_center_dec = Angle(result_table['DEC'][0] + 'd').to_value()

#spherical_date -> [ra, dec]
#central_date -> [ra, dec]
def gnomonic_spherical_to_cartesian(spherical_date, central_date):
    _phi = spherical_date[1] * np.pi / 180  #latitude - dec
    _lambda = spherical_date[0] * np.pi / 180 #longitude - ra
    _phi_1 = central_date[1] * np.pi / 180  #center_latitude - dec
    _lambda_0 = central_date[0] * np.pi / 180 #center_longitude - ra
    #The actual projection
    _cos_c = (np.sin(_phi_1)*np.sin(_phi))+(np.cos(_phi_1)*np.cos(_phi)*np.cos(_lambda-_lambda_0))
    _x = (np.cos(_phi)*np.sin(_lambda-_lambda_0))/_cos_c
    _y = (np.cos(_phi_1)*np.sin(_phi)-np.sin(_phi_1)*np.cos(_phi)*np.cos(_lambda-_lambda_0))/_cos_c
    return [_x, _y]

def euclidean_norm(cartesian_date):
    return np.sqrt(cartesian_date[0]**2 + cartesian_date[1]**2)

params = ['source_id', 'ra', 'dec', 'phot_g_mean_mag']
catalog = "gaiadr3.gaia_source_lite"
parameters = ''
for p in params:
    parameters = parameters + p + ", "  # string out of parameters
parameters = parameters[:-2]  # cut away last comma
projection_center = [roi_center_ra, roi_center_dec]

search_string = ("SELECT TOP {} {} FROM {} "
        "WHERE 1=CONTAINS(POINT('ICRS', {}.ra, {}.dec), CIRCLE('ICRS', {}, {}, {})) "
        "ORDER BY phot_g_mean_mag ASC".format(no_of_sources, parameters, catalog, catalog, catalog, roi_center_ra, roi_center_dec, roi_diameter_deg))

sources_table = gaia.tap.query(search_string)

#Remove all masked values from the results
if sources_table.has_masked_values:
    sources_table.remove_rows(np.where([c.data for c in sources_table.mask.itercols()])[-1])

#Ensure, that the number of sources is a multiple of the layer count
if not (len(sources_table) % no_layers == 0):
    sources_table = sources_table[0:(-(len(sources_table) % no_layers))]

magnitude_column = sources_table['phot_g_mean_mag']
magnitudes_list = []
if create_circles:
    radii_list = []
    max_magnitude_area = max_radius_dxf**2
    max_magnitude = magnitude_column[0]
for tmp_source in sources_table:
    tmp_magnitude = tmp_source['phot_g_mean_mag']
    magnitudes_list.append(tmp_magnitude)
    if create_circles:
        magnitude_delta = max_magnitude - tmp_magnitude
        tmp_magnitude_area = max_magnitude_area * (2.512**magnitude_delta)
        tmp_radius_dxf = np.sqrt(tmp_magnitude_area)
        radii_list.append(tmp_radius_dxf)

hist, bin_edges = np.histogram(magnitudes_list, no_layers)

layer_min_ind = 0
layer_max_ind = 0
max_euclidean_distance = 0
#project all points using a gnomonic projection
cartesian_points_layerwise = []
for bin_no in range(no_layers):
    layer_points = []
    if bin_no > 0:
        layer_min_ind = layer_min_ind + hist[bin_no-1]
    layer_max_ind = layer_max_ind + hist[bin_no]
    tmp_subtable = sources_table[layer_min_ind:layer_max_ind]
    for tmp_row in tmp_subtable:
        spherical_date = [tmp_row['ra'], tmp_row['dec']]
        cartesian_date = gnomonic_spherical_to_cartesian(spherical_date, projection_center)
        if euclidean_norm(cartesian_date) > max_euclidean_distance:
            max_euclidean_distance = euclidean_norm(cartesian_date)
        #print("{} -> {} (|x|={})".format(spherical_date, cartesian_date, euclidean_norm(cartesian_date)))
        layer_points.append(cartesian_date)
    cartesian_points_layerwise.append(layer_points)

#norm the maximum distance of all points to the desired maximum
desired_max_euclidean_distance = np.sqrt(canvas_height_mm**2 + canvas_width_mm**2) / 2
for layer_ind in range(len(cartesian_points_layerwise)):
    for i in range(len(cartesian_points_layerwise[layer_ind])):
        cartesian_points_layerwise[layer_ind][i] = [-1 * desired_max_euclidean_distance * cartesian_points_layerwise[layer_ind][i][0] / max_euclidean_distance, desired_max_euclidean_distance * cartesian_points_layerwise[layer_ind][i][1] / max_euclidean_distance]

#Write the data to a dxf-File
doc = ezdxf.new("R2010")
msp = doc.modelspace()
point_ind = 0
for layer_ind in range(len(cartesian_points_layerwise)):
    tmp_layer_name = "Layer{}".format(layer_ind+1)
    doc.layers.add(name=tmp_layer_name)
    for point in cartesian_points_layerwise[layer_ind]:
        #check if the point lies outside the canvas
        if abs(point[0]) > abs(canvas_width_mm/2) or abs(point[1]) > abs(canvas_height_mm/2):
            continue
        if create_circles:
            msp.add_circle((point[0], point[1]), radii_list[point_ind], dxfattribs={"layer": tmp_layer_name})
            point_ind = point_ind + 1
        else:
            msp.add_point((point[0], point[1]), dxfattribs={"layer": tmp_layer_name})
msp.add_line((-canvas_width_mm/2, canvas_height_mm/2), (canvas_width_mm/2, canvas_height_mm/2), dxfattribs={"layer": "canvas"})
msp.add_line((-canvas_width_mm/2, -canvas_height_mm/2), (canvas_width_mm/2, -canvas_height_mm/2), dxfattribs={"layer": "canvas"})
msp.add_line((canvas_width_mm/2, canvas_height_mm/2), (canvas_width_mm/2, -canvas_height_mm/2), dxfattribs={"layer": "canvas"})
msp.add_line((-canvas_width_mm/2, canvas_height_mm/2), (-canvas_width_mm/2, -canvas_height_mm/2), dxfattribs={"layer": "canvas"})
doc.saveas("gaia_based_image.dxf")