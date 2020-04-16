#init
import ee
import time

ee.Initialize()
# functions
def lsTOA(img):
    return ee.Algorithms.Landsat.TOA(img)
def rescale(img, thresholds):
    return img.subtract(thresholds[0]).divide(thresholds[1] - thresholds[0])


def s2CloudMask(img):
    score = ee.Image(1.0)
    qa = img.select("QA60").unmask()

    score = score.min(rescale(img.select(['B2']), [0.1, 0.5]))
    score = score.min(rescale(img.select(['B1']), [0.1, 0.3]))
    score = score.min(rescale(img.select(['B1']).add(img.select(['B10'])), [0.15, 0.2]))
    score = score.min(rescale(img.select(['B4']).add(img.select(['B3'])).add(img.select('B2')), [0.2, 0.8]))

    #Clouds are moist
    ndmi = img.normalizedDifference(['B8A','B11'])
    score=score.min(rescale(ndmi, [-0.1, 0.1]))

    # However, clouds are not snow.
    ndsi = img.normalizedDifference(['B3', 'B11'])
    score=score.min(rescale(ndsi, [0.8, 0.6]))
    score = score.multiply(100).byte().lte(cloudThresh)
    mask = score.And(qa.lt(1024)).rename(['cloudMask'])\
                .clip(img.geometry())
    img = img.addBands(mask)
    return img.divide(10000).set("system:time_start", img.get("system:time_start"),
                                 'CLOUD_COVER',img.get('CLOUD_COVERAGE_ASSESSMENT'),
                                 'SUN_AZIMUTH',img.get('MEAN_SOLAR_AZIMUTH_ANGLE'),
                                 'SUN_ZENITH',img.get('MEAN_SOLAR_ZENITH_ANGLE'),
                                 'scale',ee.Number(10))


def lsCloudMask(img):
    blank = ee.Image(0)
    scored = ee.Algorithms.Landsat.simpleCloudScore(img)
    clouds = blank.where(scored.select(['cloud']).lte(cloudThresh), 1)\
                .clip(img.geometry())
    img = img.addBands(clouds.rename(['cloudMask']))
    return img.updateMask(clouds).set("system:time_start", img.get("system:time_start"),
                   "SUN_ZENITH",ee.Number(90).subtract(img.get('SUN_ELEVATION')),
                   "scale",ee.Number(30))
def pondClassifier(shape):
    waterList = waterCollection.filterBounds(shape.geometry())\
                .sort('system:time_start',False)
    latest = ee.Image(waterList.first())
    # dates = ee.Array(waterCollection.aggregate_array('system:time_start'))

    # true = ee.Image(ee.Algorithms.If(latest.geometry().contains(shape.geometry),latest,
    #         waterList.get(1)))

    avg = latest.reduceRegion(
        reducer=ee.Reducer.mean(),
        scale=10,
        geometry=shape.geometry(),
        bestEffort=True
    )

    try:
        val = ee.Number(avg.get('water'))
        mVal = ee.Number(avg.get('cloudShadowMask'))

        test = ee.Number(ee.Algorithms.If(mVal.lt(0.5),ee.Number(3),ee.Number(0)))
        cls = test.add(val.gt(0.25))

        cls = cls.add(val.gt(0.75))
    except:
        val = random.choice(range(2))
        cls = ee.Number(val)

    return ee.Feature(shape).set({'pondCls': cls.int8()})
def mergeCollections(l8, s2, studyArea, t1, t2):
    lc8rename = l8.filterBounds(studyArea).filterDate(t1, t2).map(lsTOA).filter(
        ee.Filter.lt('CLOUD_COVER', 75)).map(lsCloudMask).select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7','cloudMask'],
                                                                 ['blue', 'green', 'red', 'nir', 'swir1', 'swir2','cloudMask'])

    st2rename = s2.filterBounds(studyArea).filterDate(t1, t2).filter(
        ee.Filter.lt('CLOUD_COVERAGE_ASSESSMENT', 75)).map(s2CloudMask).select(
        ['B2', 'B3', 'B4', 'B8', 'B11', 'B12','cloudMask'],
        ['blue', 'green', 'red', 'nir', 'swir1', 'swir2','cloudMask'])#.map(bandPassAdjustment)

    return ee.ImageCollection(lc8rename.merge(st2rename))
def simpleTDOM2(collection, zScoreThresh, shadowSumThresh, dilatePixels):
    def darkMask(img):
        zScore = img.select(shadowSumBands).subtract(irMean).divide(irStdDev)
        irSum = img.select(shadowSumBands).reduce(ee.Reducer.sum())
        TDOMMask = zScore.lt(zScoreThresh).reduce(ee.Reducer.sum()).eq(2).And(irSum.lt(shadowSumThresh)).Not()
        TDOMMask = TDOMMask.focal_min(dilatePixels)
        img = img.addBands(TDOMMask.rename(['TDOMMask']))
        # Combine the cloud and shadow masks
        combinedMask = img.select('cloudMask').mask().And(img.select('TDOMMask'))\
            .rename('cloudShadowMask');
        return img.addBands(combinedMask).updateMask(combinedMask)

    shadowSumBands = ['nir','swir1','swir2']
    irStdDev = collection.select(shadowSumBands).reduce(ee.Reducer.stdDev())
    irMean = collection.select(shadowSumBands).mean()

    collection = collection.map(darkMask)

    return collection
def calcWaterIndex(img):
    mndwi = img.expression('0.1511*B1 + 0.1973*B2 + 0.3283*B3 + 0.3407*B4 + -0.7117*B5 + -0.4559*B7',{
        'B1': img.select('blue'),
        'B2': img.select('green'),
        'B3': img.select('red'),
        'B4': img.select('nir'),
        'B5': img.select('swir1'),
        'B7': img.select('swir2'),
    }).rename('mndwi')
    return mndwi.addBands(img.select('cloudShadowMask')).copyProperties(img, ["system:time_start", "CLOUD_COVER"])

def waterClassifier(img):
    THRESHOLD = ee.Number(-0.1304)

    water = ee.Image(0).where(img.select('cloudShadowMask'), img.select('mndwi').gt(THRESHOLD))

    result = img.addBands(water.rename(['water']))
    return result.copyProperties(img, ["system:time_start","CLOUD_COVER"])
def addArea(feature):
    return feature.set('area',feature.area());
#lines  430 -~460
#https://github.com/SERVIR/WaterWatch/blob/master/tethysapp/waterwatch/utilities.py
studyArea = ee.Geometry.Rectangle([-15.866, 14.193, -12.990, 16.490])
lc8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_RT')
st2 = ee.ImageCollection('COPERNICUS/S2')
ponds = ee.FeatureCollection('projects/servir-wa/services/ephemeral_water_ferlo/ferlo_ponds')\
                .map(addArea).filter(ee.Filter.gt("area",10000))
today = time.strftime("%Y-%m-%d")

iniTime = ee.Date('2015-01-01')
endTime = ee.Date(today)

dilatePixels = 2;
cloudHeights = ee.List.sequence(200,5000,500);
zScoreThresh = -0.8;
shadowSumThresh = 0.35;
cloudThresh = 10

mergedCollection = mergeCollections(lc8, st2, studyArea, iniTime, endTime).sort('system:time_start', False)

mergedCollection = simpleTDOM2(mergedCollection, zScoreThresh, shadowSumThresh, dilatePixels)#.map(cloudProject)

mndwiCollection = mergedCollection.map(calcWaterIndex)

waterCollection = mndwiCollection.map(waterClassifier)

ponds_cls = ponds.map(pondClassifier)

visParams = {'min': 0, 'max': 3, 'palette': 'red,yellow,green,gray'}

pondsImg = ponds_cls.reduceToImage(properties=['pondCls'],
                                   reducer=ee.Reducer.first())

pondsImgID = pondsImg.getMapId(visParams)

img = mergedCollection.median().clip(studyArea)

mndwiImg = mndwiCollection.median().clip(studyArea)

gfs = ee.ImageCollection('NOAA/GFS0P25')
cfs = ee.ImageCollection('NOAA/CFSV2/FOR6H').select(['Precipitation_rate_surface_6_Hour_Average'],['precip'])
elv = ee.Image('USGS/SRTMGL1_003')

print(mndwiImg.bandNames().getInfo())