# Script to automatically run african elephant habitat monitoring
import ee
import main
import datetime
ee.Initialize()

# region of interest
# region needs to be a geometry.
# region = ee.FeatureCollection('users/path/to/feature').geometry().bounds()
region = ee.Geometry.Polygon([[[22.37055462536107, -19.69234130304949],
                                   [23.161822166438526, -19.675148989974225],
                                   [23.519800725057106, -18.180985057390387],
                                   [21.87293615648901, -17.80809895124535],
                                   [21.43371056179063, -19.463056253246073]]])

# parameters for fire as burn
# fire history length, list [integer, string],
# integer = unit of time to to collected fire data
# string = time unit of 'month','day','year'
fireHistoryLength = [-6, 'month']

# parameters for vegetation products
# list [integer, string],
# duration to include in anomalie mapping
# integer = unit of time to include
# string = time unit of 'year'
vegetationHistoryLength = [5,'year']

# parameters for water products
# list [integer, string],
# duration to include in water availability mapping
# integer = unit of time to include
# string = time unit of 'month','day','year'
waterHistoryLength = [-2,'month']

# Export information
bucket = 'wwf-gee-bucket'
# export to asset base path
assetbase = "projects/wwf-de/African_Elephants/"
# Image collection names that are within the base path
fireIC = 'fire'
waterIC = 'water'
vegetationIC = 'veg'

# # export to cloud
fireFolder = 'AfricanElephants/fireAssets'
waterFolder = 'AfricanElephants/waterAssets'
vegetationFolder = 'AfricanElephants/vegetationAssets'


# Date time parameters
py_date = datetime.datetime.utcnow()

# move back one day to the last date of previous month
ed = ee.Date(py_date).advance(-16,'day')
# ed = ee.Date('2020-09-30')

# Fire products
def make_fire():
    fire = main.Fire().burnOut(ed, -2, 'month').toInt16()

    historyFireStartDate = ed.advance(fireHistoryLength[0],fireHistoryLength[1])
    historyFire = main.Fire().historyFire(historyFireStartDate,ed)

    fire = fire.addBands(historyFire).select(['confidence', 'allFires', 'past', 'distance10k', 'historyFire'])

    return fire
# Vegetation products
def nd(image):
    out = image.normalizedDifference(['nir', 'red']).copyProperties(image, ["system:time_start"])
    return out

def make_vegetation():
    ic = main.landsat().preprocess(startDate=ee.Date('2000-01-01'),endDate=ed,studyArea=region).select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2']).map(nd)

    m = ed.get('month').getInfo()
    y = ed.get('year').getInfo()
    vegetation = main.Vegetation().byRegion(m, y, ic, region, vegetationHistoryLength)
    return vegetation

def make_water():
    # Water products
    sd = ee.Date(ed).advance(waterHistoryLength[0],waterHistoryLength[1])
    sen2 = main.sentinel2().preprocess(startDate=sd,endDate=ed,studyArea=region).select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])
    img = main.Water().wlc(sen2, region, startDate=sd,endDate=ed)
    water = main.Water().wlcexpression(img, region)
    
    return water
# Exporting
datestr = py_date.strftime("%Y%m%d")
waterDesc = '%s%s' % ('water_',datestr)
vegDesc = '%s%s' % ('vegetation_',datestr)
fireDesc = '%s%s' % ('fire_',datestr)

def export_imgs():
    main.base().exportMapToCloud(water,waterDesc,region,bucket,prefix=waterFolder,scale=10)
    main.base().exportMapToCloud(vegetation,vegDesc,region,bucket,prefix=vegetationFolder,scale=30)
    main.base().exportMapToCloud(fire,fireDesc,region,bucket,prefix=fireFolder,scale=500)

    main.base().exportMapToAsset(water,waterDesc,region,assetbase+waterIC,scale=10)
    main.base().exportMapToAsset(vegetation,vegDesc,region,assetbase+vegetationIC,scale=30)
    main.base().exportMapToAsset(fire,fireDesc,region,assetbase+fireIC,scale=500)

if __name__ == "__main__":
   
    fire = make_fire()
    vegetation = make_vegetation()
    water = make_water()

    # export_imgs()