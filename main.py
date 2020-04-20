class base(object):
    def __init__(self):
        self.date = 123
        self.studyArea = ee.FeatureCollection('somehwrere').geometry()
        self.ecoregions = ee.FeatureCollection('RESOLVE/ECOREGIONS/2017')


class Fire(base):
    def __inti__(self):
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
        currentFires = getFire(cY, cM)

        pD = sd.advance(step, unit)
        pM = pD.get('month')
        pY = pD.get('year')
        pastFires = getFire(pY, pM)

        mask = pastFires.select('binary')
        newFires = currentFires.where(mask.eq(1), 0).selfMask()
        allCurrentFires = currentFires.select('binary').rename('allFires')

        return newFires.addBands(allCurrentFires)

    def getFire(self, targetYear, targetMonth):
        # Bring in MYD14/MOD14
        modisFireAqua = ee.ImageCollection('MODIS/006/MYD14A2').select([0]);
        modisFireTerra = ee.ImageCollection('MODIS/006/MOD14A2').select([0]);
        modisFire = modisFireTerra.merge(modisFireAqua);
        singleMonth = modisFire.filter(ee.Filter.calendarRange(targetYear, targetYear, 'year')).filter(
            ee.Filter.calendarRange(targetMonth, targetMonth, 'month'));

        # Recode it, and find the year, month, and day- then add it to the map
        singleMonth = singleMonth.map(self.reclassify);
        sum_denisty = singleMonth.select('binary').sum().rename('denisty')

        return singleMonth.mosaic().addBands(sum_denisty)


class Vegetation(base):
    def __inti__(self):
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
            {'reducer': ee.Reducer.stdDev(), 'geometry': geometry, 'scale': 30, 'bestEffort': true, 'maxPixels': 1e13}).get(
            'NDVI_mean')

        baseline_ic = ic.filter(ee.Filter.calendarRange(y, y - 10, 'year'))
        baseline_mean = baseline_ic.mean()

        aandvi = month_i_mean.subtract(baseline_mean).float().rename('AANDVI')
        sandvi = aandvi.divide(ee.Image.constant(month_i_mean_std).float()).rename('SANDVI')

        vci_min = baseline_mean.reduceRegion(
            {'reducer': ee.Reducer.min(), 'geometry': geometry, 'scale': 30, 'bestEffort': true, 'maxPixels': 1e13})
        vci_max = baseline_mean.reduceRegion(
            {'reducer': ee.Reducer.max(), 'geometry': geometry, 'scale': 30, 'bestEffort': true, 'maxPixels': 1e13})
        vci_min = ee.Image.constant(vci_min.get('nd')).float()
        vci_max = ee.Image.constant(vci_max.get('nd')).float()
        vci = month_i_mean.subtract(vci_min).divide(vci_max.subtract(vci_min)).rename('VCI')

        return ee.Image.cat([month_i_mean, aandvi, sandvi, vci])

    def byRegion(self, geometry):
        biomes = ee.List(self.ecoregions.filterBounds(geometry).aggregate_array('BIOME_NUM')).distinct()

        def monthByRegion(self, b):
            a = self.ecoregions.filterBounds(geometry).filter(ee.Filter.eq('BIOME_NUM', ee.Number(b)))
            c = ee.Feature(geometry).difference(ee.Feature(a.union(1).first()))
            return monthlyNDVI(m, y, ic, c.geometry())

        ndviByBiome = biomes.map(monthByRegion)
        return ndviByBiome


class Water(base):
    def __inti__(self):
        super(Water, self).__init__()

    def waterindicies(image):
        # T NDMI
        # todo: do we want to keep band names or remap in import?
        tndmi = image.normalizedDifference(['B8', 'B11']).add(0.5).sqrt().rename('tndmi')
        tndvi = image.normalizedDifference(['NIR', 'RED']).add(0.5).sqrt().rename('tndvi')
        # NMDI - normalized multiband drought
        # https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2007GL031021
        # B8A - (B11 - B12) / B8A + (B11 + B12)
        nmdi = image.divide(1000)
        nmdi = image.expression("(B8A - (B11 - B12)) / (B8A + (B11 + B12))", {'B8A': image.select('B8A'),
                                                                              'B11': image.select('B11'),
                                                                              'B12': image.select('B12')
                                                                              }).rename('nmdi')
        return ee.Image.cat([image, tndmi, tndvi, nmdi])


if __name__ == "__main__":
    import ee

    ee.Initialize()
    # tests
    ic = ee.ImageCollection('ee.ImageCollection("LANDSAT/LC08/C01/T1_SR')
    m = 2
    y = 2019


    def maskL8sr(image):
        cloudShadowBitMask = 1 << 3;
        cloudsBitMask = 1 << 5;
        qa = image.select('pixel_qa');
        mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0).And(qa.bitwiseAnd(cloudsBitMask).eq(0))
        return image.updateMask(mask).divide(10000).select("B[0-9]*").copyProperties(image, ["system:time_start"])


    ic = ic.map(maskL8sr)
    region = ee.Geometry.Polygon([[[22.37055462536107, -19.69234130304949],
                                   [23.161822166438526, -19.675148989974225],
                                   [23.519800725057106, -18.180985057390387],
                                   [21.87293615648901, -17.80809895124535],
                                   [21.43371056179063, -19.463056253246073]]])

    a = Vegetation().monthlyNDVI(m, y, ic, region)
    print(a.bandNames().getInfo())
