// Authors: Michiel Pezij & Harm-Jan Benninga
// Date: 14 June 2018


// --- Instructions ---
//  - For the value of a pixel, go to Inspector (tab on right bar) and click on the map viewer
//  - Turn off/on maps in the viewer (below) at 'Layers' (right)
//  - To make a soil moisture map for another area:
//      1. Draw a new polygon in the viewer (below) at 'Geometry Imports' (left)
//      2. Define 'Area_small' (line 16) as this new polygon

// shared partners

// --- Imports ---
var SAT_5cm = ee.Image("users/hjbenninga/VGE_SAT_5cm"),
    WP_5cm = ee.Image("users/hjbenninga/VGE_WP_5cm"),
    study_area = /* color: #d63000 */ee.Geometry.Polygon(
        [[[5.353, 52.469],
          [5.353, 53.063],
          [6.25, 53.063],
          [6.25, 52.469]]]),
    Netherlands = /* color: #ffc82d */ee.Geometry.Polygon(
        [[[3.3, 53.7],
          [3.3, 50.7],
          [7.3, 50.7],
          [7.3, 53.7]]]),
    Raam = /* color: #98ff00 */ee.Geometry.Polygon(
        [[[5.62, 51.775],
          [5.62, 51.5],
          [5.96, 51.5],
          [5.96, 51.775]]]),
    Noordoostpolder_polder = /* color: #d63000 */ee.Geometry.Polygon(
        [[[5.8069610595703125, 52.66097536464507],
          [5.82275390625, 52.644521658618494],
          [5.876655578613281, 52.66597273265816],
          [5.869791361455782, 52.672426747254946],
          [5.860862731933594, 52.682002097753795]]]);



// --- Settings ---
var EXPORT_image = true;    // Make 'true' to export an image. Then go to 'Tasks' (tab on right bar) to start the export to Google Drive.

var Area = Netherlands;
var Area_small = Noordoostpolder_polder;            // Change this to plot for another area
var Date_image = ee.Date('2018-06-04T00:00:00');    // Change this date to plot and export image for another date


// --- Load Sentinel-1 C-band SAR Ground Range collection (log scaling, VV co-polar) ---
var collection_S1_STATS = ee.ImageCollection('COPERNICUS/S1_GRD').filterBounds(Area)
.filter(ee.Filter.listContains('transmitterReceiverPolarisation','VV'))
.filter(ee.Filter.eq('instrumentMode','IW'))
.filterDate('2016-03-01', '2018-03-01');    // Statistics are calculated over two complete hydrological years

var collection_S1_TOTAL = ee.ImageCollection('COPERNICUS/S1_GRD').filterBounds(Area)
.filter(ee.Filter.listContains('transmitterReceiverPolarisation','VV'))
.filter(ee.Filter.eq('instrumentMode','IW'))
.filterDate('2014-10-01', Date.now());

// --- Incidence angle correction ---
var n = 1;            // normalization coefficient
var angle_ref = 30;   // reference angle

// define incidence angle correction function
var incidence_angle_correction_function = function(image) {
  var image_m2m2 = image.expression(
    '10**(image/10)', {
      'image': image.select('VV')
    });
  
  var image_m2m2_cor = image_m2m2.expression(
    'sigma0*((cos(pi/180*angle_ref)**n)/(cos(pi/180*angle)**n))', {
      'n': n,
      'angle_ref': angle_ref,
      'sigma0': image_m2m2,
      'angle': image.select('angle'),
      'pi': Math.PI
    });
    
  var output_image = image_m2m2_cor.expression(
    'log10(sigma0_cor)*10', {
      'sigma0_cor': image_m2m2_cor
    });  
    
  return output_image.set('system:time_start', image.get('system:time_start'));
};

// apply incidence angle correction function
var collection_S1_STATS_IC_cor = collection_S1_STATS.map(incidence_angle_correction_function);
var collection_S1_TOTAL_IC_cor = collection_S1_TOTAL.map(incidence_angle_correction_function);

// --- Mask Sentinel-1 observations ---
var min_value = -20;        // Minimum value that is valid
var max_value = -0.001;     // Maximum value that is valid

// apply threshold values to image collections
var collection_S1_STATS_mask = collection_S1_STATS_IC_cor.map(function(img) {
                   var mask = img.gte(min_value);
                   var new_img = img.updateMask(mask);
                   
                   var mask2 = new_img.lte(max_value);
                   return new_img.updateMask(mask2);
                   });

var collection_S1_TOTAL_mask = collection_S1_TOTAL_IC_cor.map(function(img) {
                   var mask = img.gte(min_value);
                   var new_img = img.updateMask(mask);
                   
                   var mask2 = new_img.lte(max_value);
                   return new_img.updateMask(mask2);
                   });

// --- Get statistics as input to Change Detection ---
//var max_collection = collection_S1_STATS_mask.reduce(ee.Reducer.max())    // Maximum in each pixel
//var min_collection = collection_S1_STATS_mask.reduce(ee.Reducer.min())    // Minimum in each pixel

var max_collection = collection_S1_STATS_mask.reduce(ee.Reducer.percentile([97.5]));   // 97.5% percentile to exclude outliers
var min_collection = collection_S1_STATS_mask.reduce(ee.Reducer.percentile([2.5]));    // 2.5% percentile to exclude outliers

var count_collection = collection_S1_STATS_mask.reduce(ee.Reducer.count());

// --- Define change detection function ---

var change_detection_function = function(image) {
  var output_image = image.expression(
    '(s1_im - s1_min) / (s1_max - s1_min)', {
      's1_im': image,
      's1_min': min_collection,
      's1_max': max_collection
    });
  return output_image.set('system:time_start', image.get('system:time_start'));
};

// apply change detection
var cd_s1 = collection_S1_TOTAL_mask.map(change_detection_function);

// --- Scale between Wilting Point (WP) and Saturation (SAT) ---
// Wilting point and Saturation soil moisture content are adopted from BOFEK2012.

var WP_SAT_scale_function = function(image) {
  var output_image = image.expression(
    '(MAX - MIN) * index + MIN', {
      'index': image,
      'MIN': WP_5cm,
      'MAX': SAT_5cm
    });
  return output_image.set('system:time_start', image.get('system:time_start'));
};

// apply scaling between WP and SAT
var cd_s1_volumetric = cd_s1.map(WP_SAT_scale_function);
print(cd_s1_volumetric);

// --- Map display ---

Map.centerObject(Area, 6);

Map.addLayer(SAT_5cm, {min: 0, max: 0.7, opacity:1, palette: ['ff1c05', 'fff705','4dff03','07ffe8','0501ff']},
  'Saturation soil moisture [m^3/m^3]');          // Add map with saturation soil moisture (from BOFEK2012)
Map.addLayer(WP_5cm, {min: 0, max: 0.7, opacity:1, palette: ['ff1c05', 'fff705','4dff03','07ffe8','0501ff']}, 
  'Wilting point soil moisture [m^3/m^3]');       // Add map with wilting point soil moisture (from BOFEK2012)

var count_collection_Study_area = count_collection.clip(Area);
Map.addLayer(count_collection_Study_area, {min: 0, max: 500, opacity:1, palette: ['LightBlue','blue']},
  'Number of images');                            // Add map with number of Sentinel-1 images over the Netherlands

// map soil moisture image 
var cd_s1_volumetric_date_image = ee.Image(cd_s1_volumetric.filterDate(Date_image,Date.now())
  .filterBounds(Area_small).first());             // Select image that covers area of interest
  
var cd_s1_volumetric_date_image_Study_area = cd_s1_volumetric_date_image.clip(Area_small);   // Clip to image to area of interest
print(cd_s1_volumetric_date_image_Study_area);

Map.addLayer(cd_s1_volumetric_date_image_Study_area, {min: 0, max: 0.7, opacity:1,
  palette: ['ff1c05', 'fff705','4dff03','07ffe8','0501ff']}, 'Volumetric soil moisture [m^3/m^3]');


// --- Plot figure soil moisture in time ---

var SoilMoisture_TimeSeries = ui.Chart.image.series(cd_s1_volumetric, Area_small, ee.Reducer.mean())
  .setOptions({
    hAxis:{title:'Date'},
    vAxis:{title:'Volumetric soil moisture'}
  });

print(SoilMoisture_TimeSeries);


// --- Export the image ---

var date_object = ee.Date(cd_s1_volumetric_date_image_Study_area.get('system:time_start'));
var date_string = date_object.format("YYYYMMdd_HHmm");

print('Timestamp image: ', date_string);

if(EXPORT_image === true) {
  Export.image.toDrive({
    image: cd_s1_volumetric_date_image_Study_area,
    description: 'S1_Vol_SoilMoisture_' + date_string.getInfo(),
    scale: 10,          // In meter
  region: Area_small
  });
} else {
  print('No export');
}
