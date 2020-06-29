# Script to automatically run african elephant habitat monitoring
import ee
import main
import datetime
ee.Initialize()

# Date time parameters
# parameters for ndvi anomallies
py_date = datetime.datetime.utcnow()

sd = ee.Date(py_date).advance(-2,'month')
ed = ee.Date(py_date)

# region of interest
region = ee.Geometry.Polygon([[[22.37055462536107, -19.69234130304949],
                                   [23.161822166438526, -19.675148989974225],
                                   [23.519800725057106, -18.180985057390387],
                                   [21.87293615648901, -17.80809895124535],
                                   [21.43371056179063, -19.463056253246073]]])

# parameters for fire as burn
# fire history length, list [integer, string],
# integer = unit of time to to collected fire data
# string = time unit of 'months','days','years'
fireHistoryLength = [-6, 'months']

# parameters for vegetation products
# list [integer, string],
# duration to include in anomalie mapping
# integer = unit of time to include
# string = time unit of 'years'
vegetationHistoryLength = [5,'year']

# Export information
bucket = 'wwf-gee-bucket'
fireAsset = ''
# # export to cloud
fireFolder = 'AfricanElephants/fireAssets'
waterFolder = 'AfricanElephants/waterAssets'
vegetationFolder = 'AfricanElephants/vegetationAssets'


# Fire products
fire = main.Fire().burnOut(ed.advance(-1,'month'), -2, 'month').toInt16()
historyFireStartDate = ed.advance(fireHistoryLength[0],fireHistoryLength[1])
historyFire = main.Fire().historyFire(historyFireStartDate,ed)

fire = fire.addBands(historyFire).select(['confidence', 'allFires', 'past', 'distance10k', 'historyFire'])

print(fire.bandNames().getInfo(),historyFire.bandNames().getInfo())

# Vegetation products
def nd(image):
    out = image.normalizedDifference(['nir', 'red']).copyProperties(image, ["system:time_start"])
    return out

ic = main.landsat().preprocess(startDate=ee.Date('2000-01-01'),endDate=ed,studyArea=region).select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2']).map(nd)

m = sd.get('month').getInfo()
y = sd.get('year').getInfo()
vegetation = main.Vegetation().byRegion(m, y, ic, region, vegetationHistoryLength)

# Water products
sen2 = main.sentinel2().preprocess(startDate=sd,endDate=ed,studyArea=region).select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])
img = main.Water().wlc(sen2)
water = main.Water().wlcexpression(img, region)


# main.base().exportMapToCloud(water,'watertest',region,bucket,prefix=waterFolder)
# main.base().exportMapToCloud(vegetation,'vegetationtest',region,bucket,prefix=vegetationFolder,scale=30)
# main.base().exportMapToCloud(fire,'firetest',region,bucket,prefix=fireFolder,scale=500)

# export to asset
# assetbase = "users/TEST"
# main.base().exportMapToAsset(water,'watertestass',region,assetbase)
# main.base().exportMapToAsset(vegetation,'vegetationtest2',region,assetbase,scale=30)
# main.base().exportMapToAsset(fire,'firetestass',region,assetbase,scale=500)