import ee
from utils import *
import sun_angles
import view_angles


class base(object):
    def __init__(self):
        self.date = 123
        self.studyArea = ee.FeatureCollection('somehwrere').geometry()

        # Vegetation
        self.ecoregions = ee.FeatureCollection('RESOLVE/ECOREGIONS/2017')

        # Fire
        self.modisFireAqua = ee.ImageCollection('MODIS/006/MYD14A2').select([0])
        self.modisFireTerra = ee.ImageCollection('MODIS/006/MOD14A2').select([0])

        # self.startDate = ee.Date('2019-01-01')
        self.startDate = ee.Date('2019-04-01')
        self.endDate = ee.Date('2019-06-01')
        self.studyArea = ee.Geometry.Polygon([[[22.58124445939474, -18.13023785466269],
                                               [22.58124445939474, -18.203308698548458],
                                               [22.68012141251974, -18.203308698548458],
                                               [22.68012141251974, -18.13023785466269]]])
        self.metadataCloudCoverMax = 80
        self.maskSR = True
        self.cloudMask = True
        self.brdfCorrect = False
        self.shadowMasking = False
        self.cloudScoreThresh = 20

        self.toaOrSR = False

    def exportMapToAsset(self, img, desc, region, assetbase, **kwargs):
        scale = kwargs.get('scale', 20)

        task_ordered = ee.batch.Export.image.toAsset(image=img,
                                                     description=desc,
                                                     assetId='%s/%s' % (assetbase,desc),
                                                     region=region.getInfo()['coordinates'],
                                                     maxPixels=1e13,
                                                     # crs=self.epsg,
                                                     scale=scale)

        task_ordered.start()

    def exportMapToCloud(self, img, desc, region,bucket,**kwargs):

        scale = kwargs.get('scale', 20)
        prefix = kwargs.get('prefix',None)

        if prefix:
            fileNamePrefix = "%s/%s" % (prefix,desc)
        else:
            fileNamePrefix = desc

        task_ordered = ee.batch.Export.image.toCloudStorage(image=img,
                                                     description=desc,
                                                     bucket=bucket,
                                                     fileNamePrefix=fileNamePrefix,
                                                     region=region.getInfo()['coordinates'],
                                                     maxPixels=1e13,
                                                     # crs=self.epsg,
                                                     scale=scale)

        task_ordered.start()


class Fire(base):
    # todo: add check for no fires in current/previous months
    def __init__(self):
        super(Fire, self).__init__()

    def reclassify(self, img):
        remapped = img.remap([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], [0, 0, 0, 1, 1, 1, 1, 2, 3, 4]).rename(['confidence'])
        d = ee.Date(img.get('system:time_start'))
        y = ee.Image(d.get('year')).int16().rename(['year'])
        m = ee.Image(d.get('month')).int16().rename(['month']);
        day = ee.Image(d.get('day')).int16().rename(['day']);
        binary = remapped.select('confidence').gte(2).rename('binary')
        out = remapped.addBands(y).addBands(m).addBands(day).addBands(binary);
        out = out.updateMask(remapped.gte(2));

        return out

    def burnOut(self, sd, step, unit):
        """takes in a startdate as ee.Date() -maybe change to string later step as an interger
        unit of time as string (e.g. 'day','month','year') """
        cM = sd.get('month')
        cY = sd.get('year')
        currentFires = self.getFire(cY, cM)

        pD = sd.advance(step, unit)
        pM = pD.get('month')
        pY = pD.get('year')
        pastFires = self.getFire(pY, pM)

        mask = pastFires.select('binary')
        newFires = currentFires.where(mask.eq(1), 0).selfMask()
        print(currentFires.bandNames().getInfo(),'mask fire')
        allCurrentFires = currentFires.select('binary').rename('allFires')

        kernel = ee.Kernel.euclidean(100000, 'meters')
        distance = allCurrentFires.select('allFires').distance(kernel, False).rename('distance10k')

        return newFires.addBands(allCurrentFires).addBands(pastFires.select(['binary'], ['past'])).addBands(distance)

    def getFire(self, targetYear, targetMonth):
        # Bring in MYD14/MOD14
        modisFire = self.modisFireTerra.merge(self.modisFireAqua)
        singleMonth = modisFire.filter(ee.Filter.calendarRange(targetYear, targetYear, 'year')).filter(
            ee.Filter.calendarRange(targetMonth, targetMonth, 'month'))
        print(singleMonth.size().getInfo())
        # Recode it, and find the year, month, and day- then add it to the map
        singleMonth = singleMonth.map(self.reclassify);
        sum_denisty = singleMonth.select('binary').sum().rename('denisty')

        return singleMonth.mosaic().addBands(sum_denisty)


class Vegetation(base):
    def __init__(self):
        super(Vegetation, self).__init__()

    def monthlyNDVI(self, m, y, ic, geometry):
        """
        Monthly NDVI(NDVIi, m, y) for each Monitoring Unit (MU) i
        in month m and year y is obtained by averaging the 3 dekadal values in each month
        full text: https://www.frontiersin.org/articles/10.3389/fenvs.2019.00187/full
        @param self:
        @param m: month as integer
        @param y: year as integer
        @param ic: NDVI image collection
        @return: Four band image with 2 month NDVI mean, anomaly, standard anomaly, and Vegetative Control Index
        """

        ic = ic.filter(ee.Filter.calendarRange(m, m, 'month'))

        month_i_ic = ic.filter(ee.Filter.calendarRange(y, y, 'year'))

        month_i_mean = month_i_ic.mean().rename('NDVI_mean')

        month_i_mean_std = month_i_mean.reduceRegion(
            **{'reducer': ee.Reducer.stdDev(), 'geometry': geometry, 'scale': 30, 'bestEffort': True,
               'maxPixels': 1e13}).get('NDVI_mean')

        baseline_ic = ic.filter(ee.Filter.calendarRange(y, y - 10, 'year'))
        baseline_mean = baseline_ic.mean()

        aandvi = month_i_mean.subtract(baseline_mean).float().rename('AANDVI')

        sandvi = aandvi.divide(ee.Image.constant(month_i_mean_std).float()).rename('SANDVI')

        vci_min = baseline_mean.reduceRegion(
            **{'reducer': ee.Reducer.min(), 'geometry': geometry, 'scale': 30, 'bestEffort': True, 'maxPixels': 1e13})
        vci_max = baseline_mean.reduceRegion(
            **{'reducer': ee.Reducer.max(), 'geometry': geometry, 'scale': 30, 'bestEffort': True, 'maxPixels': 1e13})
        vci_min = ee.Image.constant(vci_min.get('nd')).float()
        vci_max = ee.Image.constant(vci_max.get('nd')).float()
        vci = month_i_mean.subtract(vci_min).divide(vci_max.subtract(vci_min)).rename('VCI')

        return ee.Image.cat([month_i_mean, aandvi, sandvi, vci])

    def byRegion(self, m, y, ic, region):
        eco = self.ecoregions.filterBounds(region)
        biomes = ee.List(eco.aggregate_array('BIOME_NUM')).distinct()

        def monthByRegion(b):
            a = eco.filter(ee.Filter.eq('BIOME_NUM', ee.Number(b)))
            c = ee.Feature(region).difference(ee.Feature(a.union(1).first()))
            return self.monthlyNDVI(m, y, ic, c.geometry())

        ndviByBiome = biomes.map(monthByRegion)

        # mosaic together
        out = ee.ImageCollection(ndviByBiome).mosaic()
        return out


class Water(base):
    def __init__(self):
        super(Water, self).__init__()

    def wlc(self, collection, **kwargs):
        valid_kwagrs = ['startDate', 'endDate']
        if len(kwargs.keys()) > 0 and len([i for i in valid_kwagrs if i in kwargs]) != len(kwargs.keys()):
            return print('Only valid arguments are startDate and endDate check that all key word are correct:{}'.format(
                kwargs.keys()))
        sd = kwargs.get('startDate', self.startDate)
        ed = kwargs.get('endDate', self.endDate)

        def forecastscleaning(img):
            chour = ee.Date(img.get('creation_time')).get('hour')
            fhour = ee.Date(img.get('forecast_time')).get('hour')
            fday = ee.Date(img.get('forecast_time')).get('day')
            cday = ee.Date(img.get('creation_time')).get('day')
            eq = cday.eq(fday)
            return img.set('fhour', fhour, 'chour', chour, 'fday', fday, 'cday', cday, 'test', eq)

        gfs = ee.ImageCollection("NOAA/GFS0P25").select(['temperature_2m_above_ground']).filterDate(sd, ed) \
            .filter(ee.Filter.calendarRange(12, 12, 'hour')).map(forecastscleaning) \
            .filter(ee.Filter.eq('fhour', 12)) \
            .filter(ee.Filter.eq('chour', 12)) \
            .filter(ee.Filter.eq('test', 1)).median()

        chirps = ee.ImageCollection("UCSB-CHG/CHIRPS/PENTAD").filterBounds(self.studyArea).filterDate(ee.Date(sd),
                                                                                                      ee.Date(ed)).sum()
        cfs = ee.ImageCollection('NOAA/CFSV2/FOR6H').select(['Precipitation_rate_surface_6_Hour_Average'],
                                                            ['precip']).filterDate(sd, ed).filterBounds(
            self.studyArea).sum()
        smap = ee.ImageCollection("NASA_USDA/HSL/SMAP_soil_moisture").select('ssm').filterDate(sd, ed).sum()
        chirps_spi = self.spi(ee.ImageCollection("UCSB-CHG/CHIRPS/PENTAD").filterBounds(self.studyArea), sd, ed).rename(
            'chirps_spi')

        img = collection.map(self.waterindicies).median()
        img = ee.Image.cat([img, gfs, chirps, cfs, smap, chirps_spi]).set('system:time_start', sd, 'sd', sd, 'ed', ed)
        return img

    def wlcexpression(self, img, region):
        img = img.select(['tcw', 'chirps_spi', 'ssm', 'mndwi', 'nwi', 'ndmi', 'temperature_2m_above_ground', 'ndwi'])
        img = self.normalizeBands(img, region)
        exout = img.expression("(b1*f1) +(b2 * f2) + (b3 *f3) + (b4 * f4)+ (b5 * f5)+ (b6 * f6)+ (b7 * f7)+ (b8 * f8)",
                               {'b1': img.select(0), 'b2': img.select(1),
                                'b3': img.select(2), 'b4': img.select(3),
                                'b5': img.select(4), 'b6': img.select(5),
                                'b7': ee.Image(1).subtract(img.select(6)),
                                # temp inverse  relation with water avalibility
                                'b8': img.select(7),
                                'f1': 7, 'f2': 6, 'f3': 3,
                                'f4': 4, 'f5': 5, 'f6': 5,
                                'f7': 4, 'f8': 7})
        return exout

    def waterindicies(self, image):
        ndwi = image.normalizedDifference(['green', 'nir']).rename('ndwi')
        ndmi = image.normalizedDifference(['nir', 'swir1']).rename('ndmi')
        mndwi = image.normalizedDifference(['green', 'swir1']).rename('mndwi')
        nwi = image.expression('((b-(n+s+w))/(b+(n+s+w))*100)', {
            'b': image.select('blue'),
            'n': image.select('nir'),
            's': image.select('swir1'),
            'w': image.select('swir2')}).rename('nwi')
        # add tesselcap wetness
        tcw = image.expression('0.1511*B1 + 0.1973*B2 + 0.3283*B3 + 0.3407*B4 + -0.7117*B5 + -0.4559*B7', {
            'B1': image.select('blue'),
            'B2': image.select('green'),
            'B3': image.select('red'),
            'B4': image.select('nir'),
            'B5': image.select('swir1'),
            'B7': image.select('swir2'),
        }).rename('tcw')
        # var factors = ee.ImageCollection.fromImages([waterlss2.select('mndwi'),chirps_spi,smap,mndwi,nwi,ndmi,gfs,ndwi]).map(function(f){ return n(f,aoi)})
        imglist = [ndwi, ndmi, mndwi, nwi, tcw]

        return ee.Image.cat(imglist)

    def spi(self, col, sd, ed):
        nsd = sd.advance(-6, 'month')

        spi6 = col.filterDate(nsd, sd)

        mean = spi6.mean()

        std = spi6.reduce(ee.Reducer.stdDev())

        img = col.filterDate(sd, ed).median()
        img = img.subtract(mean).divide(std).set('system:time_start', sd)
        return ee.Image(img)

    def normalizeBands(self, img,region):
        def mm(plist):
            plist = ee.List(plist)
            # // get bandname
            fbn = ee.String(plist.get(0)).slice(0, -4)
            imgMin = ee.Image.constant(mmvalues.get(plist.get(0)))
            imgMax = ee.Image.constant(mmvalues.get(plist.get(1)))
            normImg = img.select(fbn).subtract(imgMin).divide(imgMax.subtract(imgMin));
            return normImg.toFloat()

        def rename(i):
            imin = ee.String(i).cat(preMin)
            imax = ee.String(i).cat(preMax)
            return ee.List([imin, imax])

        mmvalues = img.reduceRegion(reducer=ee.Reducer.minMax(), scale=300, maxPixels=1e13, bestEffort=True, \
                                    geometry=region, tileScale=6)
        bnog = img.bandNames()
        preMax = '_max'
        preMin = '_min'
        bn = bnog.map(rename)

        return ee.ImageCollection(bn.map(mm)).toBands().rename(bnog)


class landsat(base):
    def __init__(self):
        super(landsat, self).__init__()
        self.divideBands = ee.List(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])
        self.bandNamesLandsat = ee.List(
            ['blue', 'green', 'red', 'nir', 'swir1', 'thermal', 'swir2', 'sr_atmos_opacity', 'pixel_qa'])
        self.sensorBandDictLandsatSR = ee.Dictionary({'L8': ee.List([1, 2, 3, 4, 5, 7, 6, 9, 10]), \
                                                      'L7': ee.List([0, 1, 2, 3, 4, 5, 6, 7, 9]), \
                                                      'L5': ee.List([0, 1, 2, 3, 4, 5, 6, 7, 9]), \
                                                      'L4': ee.List([0, 1, 2, 3, 4, 5, 6, 7, 9])})

    def loadls(self, startDate, endDate, studyArea):
        landsat8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR').filterDate(startDate,
                                                                            endDate).filterBounds(studyArea)
        landsat8 = landsat8.filterMetadata('CLOUD_COVER', 'less_than', self.metadataCloudCoverMax)
        landsat8 = landsat8.select(self.sensorBandDictLandsatSR.get('L8'), self.bandNamesLandsat)

        landsat5 = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR').filterDate(startDate,
                                                                            endDate).filterBounds(studyArea)
        landsat5 = landsat5.filterMetadata('CLOUD_COVER', 'less_than', self.metadataCloudCoverMax)
        landsat5 = landsat5.select(self.sensorBandDictLandsatSR.get('L5'), self.bandNamesLandsat).map(
            self.defringe)

        landsat7 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR').filterDate(startDate,
                                                                            endDate).filterBounds(studyArea)
        landsat7 = landsat7.filterMetadata('CLOUD_COVER', 'less_than', self.metadataCloudCoverMax)
        landsat7 = landsat7.select(self.sensorBandDictLandsatSR.get('L7'), self.bandNamesLandsat)

        return landsat5.merge(landsat7).merge(landsat8)

    def preprocess(self, **kwargs):
        valid_kwagrs = ['startDate', 'endDate','studyArea']
        if len(kwargs.keys()) > 0 and len([i for i in valid_kwagrs if i in kwargs]) != len(kwargs.keys()):
            return print('Only valid arguments are startDate, studyArea, and endDate check that all key word are correct:{}'.format(
                kwargs.keys()))

        startDate = kwargs.get('startDate', self.startDate)
        endDate = kwargs.get('endDate', self.endDate)
        studyArea = kwargs.get('studyArea',self.studyArea)

        landsat = self.loadls(startDate=startDate, endDate=endDate, studyArea=studyArea)
        if landsat.size().getInfo() > 0:
            if self.maskSR == True:
                print("removing clouds")
                landsat = landsat.map(self.CloudMaskSRL8)

            landsat = landsat.map(self.scaleLandsat)

            if self.cloudMask == True:
                print("removing some more clouds")
                landsat = landsat.map(self.maskClouds)
            if self.brdfCorrect == True:
                landsat = landsat.map(self.brdf)

            return landsat

    def defringe(self, img):
        # threshold for defringing landsat5 and 7
        fringeCountThreshold = 279

        k = ee.Kernel.fixed(41, 41,
                            [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

        m = ee.Image(img).mask().reduce(ee.Reducer.min())
        sum = m.reduceNeighborhood(ee.Reducer.sum(), k, 'kernel')
        mask = sum.gte(fringeCountThreshold)

        return img.updateMask(mask)

    def CloudMaskSRL8(self, img):
        """apply cf-mask Landsat"""
        QA = img.select("pixel_qa")

        shadow = QA.bitwiseAnd(8).neq(0);
        cloud = QA.bitwiseAnd(32).neq(0);
        return img.updateMask(shadow.Not()).updateMask(cloud.Not()).copyProperties(img)

    def maskClouds(self, img):
        """
        Computes spectral indices of cloudyness and take the minimum of them.

        Each spectral index is fairly lenient because the group minimum
        is a somewhat stringent comparison policy. side note -> this seems like a job for machine learning :)
        originally written by Matt Hancher for Landsat imageryadapted to Sentinel by Chris Hewig and Ian Housman
        """

        score = ee.Image(1.0);
        # Clouds are reasonably bright in the blue band.
        blue_rescale = img.select('blue').subtract(ee.Number(0.1)).divide(ee.Number(0.3).subtract(ee.Number(0.1)))
        score = score.min(blue_rescale);

        # Clouds are reasonably bright in all visible bands.
        visible = img.select('red').add(img.select('green')).add(img.select('blue'))
        visible_rescale = visible.subtract(ee.Number(0.2)).divide(ee.Number(0.8).subtract(ee.Number(0.2)))
        score = score.min(visible_rescale);

        # Clouds are reasonably bright in all infrared bands.
        infrared = img.select('nir').add(img.select('swir1')).add(img.select('swir2'))
        infrared_rescale = infrared.subtract(ee.Number(0.3)).divide(ee.Number(0.8).subtract(ee.Number(0.3)))
        score = score.min(infrared_rescale);

        # Clouds are reasonably cool in temperature.
        temp_rescale = img.select('thermal').subtract(ee.Number(300)).divide(ee.Number(290).subtract(ee.Number(300)))
        score = score.min(temp_rescale);

        # However, clouds are not snow.
        ndsi = img.normalizedDifference(['green', 'swir1']);
        ndsi_rescale = ndsi.subtract(ee.Number(0.8)).divide(ee.Number(0.6).subtract(ee.Number(0.8)))
        score = score.min(ndsi_rescale).multiply(100).byte().rename('cloudScore');
        mask = score.lt(self.cloudScoreThresh).rename(['cloudMask']);
        img = img.updateMask(mask).addBands([mask]).addBands(score);

        return img;

    def brdf(self, img):

        # import sun_angles
        # import view_angles

        def _apply(image, kvol, kvol0):
            blue = _correct_band(image, 'blue', kvol, kvol0, f_iso=0.0774, f_geo=0.0079, f_vol=0.0372)
            green = _correct_band(image, 'green', kvol, kvol0, f_iso=0.1306, f_geo=0.0178, f_vol=0.0580)
            red = _correct_band(image, 'red', kvol, kvol0, f_iso=0.1690, f_geo=0.0227, f_vol=0.0574)
            nir = _correct_band(image, 'nir', kvol, kvol0, f_iso=0.3093, f_geo=0.0330, f_vol=0.1535)
            swir1 = _correct_band(image, 'swir1', kvol, kvol0, f_iso=0.3430, f_geo=0.0453, f_vol=0.1154)
            swir2 = _correct_band(image, 'swir2', kvol, kvol0, f_iso=0.2658, f_geo=0.0387, f_vol=0.0639)
            return replace_bands(image, [blue, green, red, nir, swir1, swir2])

        def _correct_band(image, band_name, kvol, kvol0, f_iso, f_geo, f_vol):
            """fiso + fvol * kvol + fgeo * kgeo"""
            iso = ee.Image(f_iso)
            geo = ee.Image(f_geo)
            vol = ee.Image(f_vol)
            pred = vol.multiply(kvol).add(geo.multiply(kvol)).add(iso).rename(['pred'])
            pred0 = vol.multiply(kvol0).add(geo.multiply(kvol0)).add(iso).rename(['pred0'])
            cfac = pred0.divide(pred).rename(['cfac'])
            corr = image.select(band_name).multiply(cfac).rename([band_name])
            return corr

        def _kvol(sunAz, sunZen, viewAz, viewZen):
            """Calculate kvol kernel.
            From Lucht et al. 2000
            Phase angle = cos(solar zenith) cos(view zenith) + sin(solar zenith) sin(view zenith) cos(relative azimuth)"""

            relative_azimuth = sunAz.subtract(viewAz).rename(['relAz'])
            pa1 = viewZen.cos() \
                .multiply(sunZen.cos())
            pa2 = viewZen.sin() \
                .multiply(sunZen.sin()) \
                .multiply(relative_azimuth.cos())
            phase_angle1 = pa1.add(pa2)
            phase_angle = phase_angle1.acos()
            p1 = ee.Image(PI().divide(2)).subtract(phase_angle)
            p2 = p1.multiply(phase_angle1)
            p3 = p2.add(phase_angle.sin())
            p4 = sunZen.cos().add(viewZen.cos())
            p5 = ee.Image(PI().divide(4))

            kvol = p3.divide(p4).subtract(p5).rename(['kvol'])

            viewZen0 = ee.Image(0)
            pa10 = viewZen0.cos() \
                .multiply(sunZen.cos())
            pa20 = viewZen0.sin() \
                .multiply(sunZen.sin()) \
                .multiply(relative_azimuth.cos())
            phase_angle10 = pa10.add(pa20)
            phase_angle0 = phase_angle10.acos()
            p10 = ee.Image(PI().divide(2)).subtract(phase_angle0)
            p20 = p10.multiply(phase_angle10)
            p30 = p20.add(phase_angle0.sin())
            p40 = sunZen.cos().add(viewZen0.cos())
            p50 = ee.Image(PI().divide(4))

            kvol0 = p30.divide(p40).subtract(p50).rename(['kvol0'])

            return (kvol, kvol0)

        date = img.date()
        footprint = determine_footprint(img)
        print(footprint.getInfo()), 'footprint'
        (sunAz, sunZen) = sun_angles.create(date, footprint)
        (viewAz, viewZen) = view_angles.create(footprint)
        (kvol, kvol0) = _kvol(sunAz, sunZen, viewAz, viewZen)
        return _apply(img, kvol.multiply(PI()), kvol0.multiply(PI()))

    def scaleLandsat(self, img):
        """Landast is scaled by factor 0.0001 """
        thermal = img.select(ee.List(['thermal'])).multiply(0.1)
        scaled = ee.Image(img).select(self.divideBands).multiply(ee.Number(0.0001))

        return img.select(['pixel_qa']).addBands(scaled).addBands(thermal)


class sentinel2(base):
    def __init__(self):
        super(sentinel2, self).__init__()
        self.s2BandsIn = ee.List(
            ['QA60', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12'])
        self.s2BandsOut = ee.List(
            ['QA60', 'cb', 'blue', 'green', 'red', 're1', 're2', 're3', 'nir', 're4', 'waterVapor', 'cirrus', 'swir1',
             'swir2'])
        self.divideBands = ee.List(
            ['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 're4', 'cb', 'cirrus', 'swir1', 'swir2', 'waterVapor'])

        # contractPixels: The radius of the number of pixels to contract (negative buffer) clouds and cloud shadows by. Intended to eliminate smaller cloud
        #    patches that are likely errors (1.5 results in a -1 pixel buffer)(0.5 results in a -0 pixel buffer)
        # (1.5 or 2.5 generally is sufficient)
        self.contractPixels = 1.5

        # dilatePixels: The radius of the number of pixels to dilate (buffer) clouds
        # and cloud shadows by. Intended to include edges of clouds/cloud shadows
        # that are often missed (1.5 results in a 1 pixel buffer)(0.5 results in a 0 pixel buffer)
        # (2.5 or 3.5 generally is sufficient)
        self.dilatePixels = 1.5
        self.zScoreThresh = -0.8
        self.shadowSumThresh = 30

    def loads2(self, start, end, studyArea):
        collection = 'COPERNICUS/S2'
        if self.toaOrSR:  collection = 'COPERNICUS/S2_SR'
        s2s = ee.ImageCollection(collection).filterDate(start, end) \
            .filterBounds(studyArea) \
            .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', self.metadataCloudCoverMax)) \
            # .filter(ee.Filter.lt('CLOUD_COVERAGE_ASSESSMENT', self.metadataCloudCoverMax)) \

        return s2s

    def preprocess(self, **kwargs):
        valid_kwagrs = ['startDate', 'endDate', 'studyArea']
        if len(kwargs.keys()) > 0 and len([i for i in valid_kwagrs if i in kwargs]) != len(kwargs.keys()):
            return print('Only valid arguments are studyArea, startDate and endDate check that all key word are correct:{}'.format(
                kwargs.keys()))

        startDate = kwargs.get('startDate', self.startDate)
        endDate = kwargs.get('endDate', self.endDate)
        studyArea = kwargs.get('studyArea', self.studyArea)

        s2 = self.loads2(startDate, endDate, studyArea)

        if self.toaOrSR:
            self.s2BandsIn = self.s2BandsIn.remove('B10')
            self.s2BandsOut = self.s2BandsOut.remove('cirrus')
            self.divideBands = self.s2BandsOut.remove('cirrus')
            self.cloudMask = False

        if s2.size().getInfo() > 0:
            s2 = s2.map(self.scaleS2)

            if self.shadowMasking == True:
                s2 = self.maskShadows(s2, studyArea)
            s2 = s2.select(self.s2BandsIn, self.s2BandsOut)

            if self.maskSR == True:
                print("use QA band for cloud Masking")
                s2 = s2.map(self.QAMaskCloud)

            if self.cloudMask == True:
                print("sentinel cloud score...")
                s2 = s2.map(self.sentinelCloudScore)
                s2 = self.cloudMasking(s2)

            if self.brdfCorrect == True:
                print("apply brdf correction..")
                s2 = s2.map(self.brdf)

            return s2

    def scaleS2(self, img):
        divideBands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B11', 'B12']
        bandNames = img.bandNames()
        otherBands = bandNames.removeAll(divideBands)

        others = img.select(otherBands)
        out = img.select(divideBands).divide(10000)
        return out.addBands(others).copyProperties(img,
                                                   ['system:time_start', 'system:footprint', 'MEAN_SOLAR_ZENITH_ANGLE',
                                                    'MEAN_SOLAR_AZIMUTH_ANGLE']).set("centroid",
                                                                                     img.geometry().centroid())

    def maskShadows(self, collection, studyArea):

        def TDOM(image):
            zScore = image.select(shadowSumBands).subtract(irMean).divide(irStdDev)
            irSum = image.select(shadowSumBands).reduce(ee.Reducer.sum())
            TDOMMask = zScore.lt(self.zScoreThresh).reduce(ee.Reducer.sum()).eq(2) \
                .And(irSum.lt(self.shadowSumThresh)).Not()
            TDOMMask = TDOMMask.focal_min(self.dilatePixels)
            return image.addBands(TDOMMask.rename(['TDOMMask']))

        def mask(image):
            outimg = image.updateMask(image.select(['TDOMMask']))
            return outimg

        shadowSumBands = ['B8', 'B11']

        allCollection = ee.ImageCollection('COPERNICUS/S2').filterBounds(studyArea).filter(
            ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 30))

        # Get some pixel-wise stats for the time series
        irStdDev = allCollection.select(shadowSumBands).reduce(ee.Reducer.stdDev())
        irMean = allCollection.select(shadowSumBands).reduce(ee.Reducer.mean())

        # Mask out dark dark outliers
        collection_tdom = collection.map(TDOM)

        return collection_tdom.map(mask)

    def sentinelCloudScore(self, img):
        """
        Computes spectral indices of cloudyness and take the minimum of them.

        Each spectral index is fairly lenient because the group minimum
        is a somewhat stringent comparison policy. side note -> this seems like a job for machine learning :)

        originally written by Matt Hancher for Landsat imagery
        adapted to Sentinel by Chris Hewig and Ian Housman
        """

        def rescale(img, thresholds):
            """
            Linear stretch of image between two threshold values.
            """
            return img.subtract(thresholds[0]).divide(thresholds[1] - thresholds[0])

        # cloud until proven otherwise
        score = ee.Image(1)
        blueCirrusScore = ee.Image(0)

        # clouds are reasonably bright
        blueCirrusScore = blueCirrusScore.max(rescale(img.select(['blue']), [0.1, 0.5]))
        blueCirrusScore = blueCirrusScore.max(rescale(img.select(['cb']), [0.1, 0.5]))
        blueCirrusScore = blueCirrusScore.max(rescale(img.select(['cirrus']), [0.1, 0.3]))
        score = score.min(blueCirrusScore)

        score = score.min(rescale(img.select(['red']).add(img.select(['green'])).add(img.select('blue')), [0.2, 0.8]))
        score = score.min(rescale(img.select(['nir']).add(img.select(['swir1'])).add(img.select('swir2')), [0.3, 0.8]))

        # clouds are moist
        ndsi = img.normalizedDifference(['green', 'swir1'])
        score = score.min(rescale(ndsi, [0.8, 0.6]))
        score = score.multiply(100).byte();
        score = score.clamp(0, 100);

        return img.addBands(score.rename(['cloudScore']))

    def QAMaskCloud(self, img):
        bandNames = img.bandNames()
        otherBands = bandNames.removeAll(self.divideBands)
        others = img.select(otherBands)

        qa = img.select('QA60').int16();

        img = img.select(self.divideBands)

        # Bits 10 and 11 are clouds and cirrus, respectively.
        cloudBitMask = int(math.pow(2, 10));
        cirrusBitMask = int(math.pow(2, 11));

        # Both flags should be set to zero, indicating clear conditions.
        mask = qa.bitwiseAnd(cloudBitMask).eq(0).And(qa.bitwiseAnd(cirrusBitMask).eq(0));

        img = img.updateMask(mask).addBands(others)

        # Return the masked and scaled data.
        return img

    def cloudMasking(self, collection):

        def maskClouds(img):
            cloudMask = img.select(['cloudScore']).lt(self.cloudScoreThresh) \
                .focal_min(self.dilatePixels) \
                .focal_max(self.contractPixels) \
                .rename(['cloudMask'])

            bandNames = img.bandNames()
            otherBands = bandNames.removeAll(self.divideBands)
            others = img.select(otherBands)

            img = img.select(self.divideBands).updateMask(cloudMask)

            return img.addBands(cloudMask).addBands(others);

        # Find low cloud score pctl for each pixel to avoid comission errors
        # minCloudScore = collection.select(['cloudScore']).reduce(ee.Reducer.percentile([self.cloudScorePctl]));

        collection = collection.map(maskClouds)

        return collection

    def brdf(self, img):

        def _apply(image, kvol, kvol0):
            blue = _correct_band(image, 'blue', kvol, kvol0, f_iso=0.0774, f_geo=0.0079, f_vol=0.0372)
            green = _correct_band(image, 'green', kvol, kvol0, f_iso=0.1306, f_geo=0.0178, f_vol=0.0580)
            red = _correct_band(image, 'red', kvol, kvol0, f_iso=0.1690, f_geo=0.0227, f_vol=0.0574)
            re1 = _correct_band(image, 're1', kvol, kvol0, f_iso=0.2085, f_geo=0.0256, f_vol=0.0845)
            re2 = _correct_band(image, 're2', kvol, kvol0, f_iso=0.2316, f_geo=0.0273, f_vol=0.1003)
            re3 = _correct_band(image, 're3', kvol, kvol0, f_iso=0.2599, f_geo=0.0294, f_vol=0.1197)
            nir = _correct_band(image, 'nir', kvol, kvol0, f_iso=0.3093, f_geo=0.0330, f_vol=0.1535)
            re4 = _correct_band(image, 're4', kvol, kvol0, f_iso=0.2907, f_geo=0.0410, f_vol=0.1611)
            swir1 = _correct_band(image, 'swir1', kvol, kvol0, f_iso=0.3430, f_geo=0.0453, f_vol=0.1154)
            swir2 = _correct_band(image, 'swir2', kvol, kvol0, f_iso=0.2658, f_geo=0.0387, f_vol=0.0639)
            return replace_bands(image, [blue, green, red, re1, re2, re3, nir, re4, swir1, swir2])

        def _correct_band(image, band_name, kvol, kvol0, f_iso, f_geo, f_vol):
            """fiso + fvol * kvol + fgeo * kgeo"""
            iso = ee.Image(f_iso)
            geo = ee.Image(f_geo)
            vol = ee.Image(f_vol)
            pred = vol.multiply(kvol).add(geo.multiply(kvol)).add(iso).rename(['pred'])
            pred0 = vol.multiply(kvol0).add(geo.multiply(kvol0)).add(iso).rename(['pred0'])
            cfac = pred0.divide(pred).rename(['cfac'])
            corr = image.select(band_name).multiply(cfac).rename([band_name])
            return corr

        def _kvol(sunAz, sunZen, viewAz, viewZen):
            """Calculate kvol kernel.
            From Lucht et al. 2000
            Phase angle = cos(solar zenith) cos(view zenith) + sin(solar zenith) sin(view zenith) cos(relative azimuth)"""

            relative_azimuth = sunAz.subtract(viewAz).rename(['relAz'])
            pa1 = viewZen.cos() \
                .multiply(sunZen.cos())
            pa2 = viewZen.sin() \
                .multiply(sunZen.sin()) \
                .multiply(relative_azimuth.cos())
            phase_angle1 = pa1.add(pa2)
            phase_angle = phase_angle1.acos()
            p1 = ee.Image(PI().divide(2)).subtract(phase_angle)
            p2 = p1.multiply(phase_angle1)
            p3 = p2.add(phase_angle.sin())
            p4 = sunZen.cos().add(viewZen.cos())
            p5 = ee.Image(PI().divide(4))

            kvol = p3.divide(p4).subtract(p5).rename(['kvol'])

            viewZen0 = ee.Image(0)
            pa10 = viewZen0.cos() \
                .multiply(sunZen.cos())
            pa20 = viewZen0.sin() \
                .multiply(sunZen.sin()) \
                .multiply(relative_azimuth.cos())
            phase_angle10 = pa10.add(pa20)
            phase_angle0 = phase_angle10.acos()
            p10 = ee.Image(PI().divide(2)).subtract(phase_angle0)
            p20 = p10.multiply(phase_angle10)
            p30 = p20.add(phase_angle0.sin())
            p40 = sunZen.cos().add(viewZen0.cos())
            p50 = ee.Image(PI().divide(4))

            kvol0 = p30.divide(p40).subtract(p50).rename(['kvol0'])

            return (kvol, kvol0)

        date = img.date()
        footprint = ee.List(img.geometry().bounds().bounds().coordinates().get(0));
        (sunAz, sunZen) = sun_angles.create(date, footprint)
        (viewAz, viewZen) = view_angles.create(footprint)
        (kvol, kvol0) = _kvol(sunAz, sunZen, viewAz, viewZen)

        bandNames = img.bandNames()
        otherBands = bandNames.removeAll(self.divideBands)
        others = img.select(otherBands)

        img = ee.Image(_apply(img, kvol.multiply(PI()), kvol0.multiply(PI())))

        return img

    def sentinelCloudScore(self, img):
        """
        Computes spectral indices of cloudyness and take the minimum of them.

        Each spectral index is fairly lenient because the group minimum
        is a somewhat stringent comparison policy. side note -> this seems like a job for machine learning :)

        originally written by Matt Hancher for Landsat imagery
        adapted to Sentinel by Chris Hewig and Ian Housman
        """

        def rescale(img, thresholds):
            """
            Linear stretch of image between two threshold values.
            """
            return img.subtract(thresholds[0]).divide(thresholds[1] - thresholds[0])

        # cloud until proven otherwise
        score = ee.Image(1)
        blueCirrusScore = ee.Image(0)

        # clouds are reasonably bright
        blueCirrusScore = blueCirrusScore.max(rescale(img.select(['blue']), [0.1, 0.5]))
        blueCirrusScore = blueCirrusScore.max(rescale(img.select(['cb']), [0.1, 0.5]))
        blueCirrusScore = blueCirrusScore.max(rescale(img.select(['cirrus']), [0.1, 0.3]))
        score = score.min(blueCirrusScore)

        score = score.min(rescale(img.select(['red']).add(img.select(['green'])).add(img.select('blue')), [0.2, 0.8]))
        score = score.min(rescale(img.select(['nir']).add(img.select(['swir1'])).add(img.select('swir2')), [0.3, 0.8]))

        # clouds are moist
        ndsi = img.normalizedDifference(['green', 'swir1'])
        score = score.min(rescale(ndsi, [0.8, 0.6]))
        score = score.multiply(100).byte();
        score = score.clamp(0, 100);

        return img.addBands(score.rename(['cloudScore']))


if __name__ == "__main__":

    ee.Initialize()
    # tests

    ndvi_tests = True
    fire_test = False
    collection_test = False
    water_tests = False

    region = ee.Geometry.Polygon([[[22.37055462536107, -19.69234130304949],
                                   [23.161822166438526, -19.675148989974225],
                                   [23.519800725057106, -18.180985057390387],
                                   [21.87293615648901, -17.80809895124535],
                                   [21.43371056179063, -19.463056253246073]]])

    if water_tests:
        t = sentinel2().preprocess().select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])
        # g = Water().wlc(t)
        # wtf = Water().linearScale(g.select([0]), base().studyArea)
        # print('ok')
        # b = Water().wlcexpression(g, base().studyArea)
        try:
            img = Water().wlc(t)

            b = Water().wlcexpression(img, region)
            print(b.bandNames().getInfo())
        except Exception as e:
            print(e.args)
            print('fml')
        try:
            str == Water().wlc(t, supDate=111)
            print('passes wrong kwargs')
        except:
            print('fails wrong kwargs')
        try:
            print(Water().wlc(t).bandNames().getInfo())

            print('passes wlc with defaults')
        except:
            print('wlc fails with defaults')
        try:
            sd = ee.Date('2018-01-01')
            ed = ee.Date('2018-03-01')
            img = Water().wlc(t, startDate=sd, endDate=ed)
            sd.getInfo()['value'] == img.get('sd').getInfo()['value']
            print('passes useing custom date')
        except:
            print('failed to use correct date')
        try:
    #         test scaling bands
            b = t.first()
            scale = Water().normalizeBands(b,region)
            if isinstance(scale,ee.Image) == 0:
                raise BaseException()
            print('scaling all bands returns image')
        except Exception as e:
            print('failed scalling bands test')
            print(e.args,e.message)
        try:
            img = Water().wlc(t)
            b = Water().wlcexpression(img, region)
            base().exportMap(b,'wlctest20_20170401_20170601_toa',region)
        except Exception as e:
            print('failed exporting')
            print(e.args)


    if collection_test:
        startDate = ee.Date('2017-04-01')
        endDate = ee.Date('2017-05-01')
        t = sentinel2().loads2(start=startDate, end=endDate, studyArea=region).select(
            ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])
        # brdf = sentinel2().brdf(t.first())
        # print(brdf)
        try:
            t = sentinel2().preprocess().select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])
            t2 = landsat().preprocess().select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])
            t2.first().bandNames().getInfo() == t.first().bandNames().getInfo()
            print('passed loading collections')
        except:
            print('failed loading collections')
        try:
            t = sentinel2().preprocess().select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])
            t2 = landsat().preprocess().select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])
            merged_img = t.merge(t2).min()
            print('passed merging collections')
        except:
            print('failed merging collections')
        try:

            sentinel2().maskShadows(t, region)
            print(t.size().getInfo())
            print(t.first().bandNames().getInfo())
        except:
            print('no')

    if fire_test:
        sd = ee.Date('2019-02-01')

        t = Fire().burnOut(sd, -2, 'month')
        print(t.bandNames().getInfo())

        pass
    if ndvi_tests:
        ic = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR")
        m = 2
        y = 2019


        def maskL8sr(image):
            cloudShadowBitMask = 1 << 3;
            cloudsBitMask = 1 << 5;
            qa = image.select('pixel_qa');
            mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0).And(qa.bitwiseAnd(cloudsBitMask).eq(0))
            img = image.updateMask(mask).divide(10000).select("B[0-9]*")
            out = img.normalizedDifference(['B5', 'B4']).copyProperties(image, ["system:time_start"])
            return out


        ic = ic.map(maskL8sr).filterBounds(region)
        print(ic.size().getInfo())
        a = Vegetation().byRegion(m, y, ic, region)
        print(a.bandNames().getInfo())
        try:
            t2 = landsat().preprocess().select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])
            v = Vegetation().byRegion(m,y,t2,region)
            print(v.bandNames().getInfo())
        except Exception as e:
            print('failed using ls collection')
            print(e.args)
