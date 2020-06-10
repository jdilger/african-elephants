# Script to automatically run african elephant habitat monitoring
import ee
import main
import datetime
ee.Initialize()

# TODO change sd/ed to pull current data
py_date = datetime.datetime.utcnow()

sd = ee.Date(py_date).advance(-2,'month')
ed = ee.Date(py_date)

# Export information
bucket = 'goldminehack'

region = ee.Geometry.Polygon([[[22.37055462536107, -19.69234130304949],
                                   [23.161822166438526, -19.675148989974225],
                                   [23.519800725057106, -18.180985057390387],
                                   [21.87293615648901, -17.80809895124535],
                                   [21.43371056179063, -19.463056253246073]]])

# Fire products
# todo: when check added for 0 fires (current/previous) update ed
fire = main.Fire().burnOut(ed.advance(-1,'month'), -2, 'month').toInt16()

# Vegetation products
def nd(image):
    out = image.normalizedDifference(['nir', 'red']).copyProperties(image, ["system:time_start"])
    return out

ic = main.landsat().preprocess(startDate=ee.Date('2000-01-01'),endDate=ed,studyArea=region).select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2']).map(nd)
print(ic.size().getInfo())
m = sd.get('month').getInfo()
y = sd.get('year').getInfo()
vegetation = main.Vegetation().byRegion(m, y, ic, region)

# Water products
sen2 = main.sentinel2().preprocess(startDate=sd,endDate=ed,studyArea=region).select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])
img = main.Water().wlc(sen2)
water = main.Water().wlcexpression(img, region)

# # export to cloud
# main.base().exportMapToCloud(water,'watertest',region,'goldminehack',prefix='wwfWater')
# main.base().exportMapToCloud(vegetation,'vegetationtest',region,'goldminehack',prefix='wwfVegetation',scale=30)
# main.base().exportMapToCloud(fire,'firetest',region,'goldminehack',prefix='wwfFire',scale=500)
# export to asset
assetbase = "users/TEST"
# main.base().exportMapToAsset(water,'watertestass',region,assetbase)
# main.base().exportMapToAsset(vegetation,'vegetationtestass',region,assetbase,scale=30)
main.base().exportMapToAsset(fire,'firetestass',region,assetbase,scale=500)