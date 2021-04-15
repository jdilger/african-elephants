from run import *
import unittest
from dateutil.relativedelta import *

class runTests(unittest.TestCase):

    def test_check_valid_region(self):
        self.assertNotIsInstance(region, ee.FeatureCollection, msg="Input region needs to be a geometry not feature collection. Try casting to geometry with .geometry()")
        self.assertIsInstance(region, ee.Geometry, msg="Input region needs to be a geometry.")

    def test_region_is_accesible(self):
        # region = ee.FeatureCollection('users/abitathafrog/abitathafrog_11_27_2019') # not shared
        # region = ee.FeatureCollection('users/abitathafrog/stratRandomSamples-201911') #shared
        try:
            ee.data.computeValue(region)
        except:
            self.fail('Failed to read region asset. Check path is correct and that it is readable from current account.')

    def test_fire_params(self):
        self.assertEqual(len(fireHistoryLength),2,msg='fire history length requires 2 inputs for durration and time unit')
        self.assertIsInstance(fireHistoryLength[0],int, msg='fire history lengths first value should be an integer')
        self.assertIsInstance(fireHistoryLength[1],str,msg='fire history lengths second value need to be month, day, or year,')
        self.assertIn(fireHistoryLength[1],['month','day','year'],msg='fire history lengths second value need to be month, day, or year,')
    
    def test_vegetation_params(self):
        self.assertEqual(len(vegetationHistoryLength),2,msg='vegetation history length requires 2 inputs for durration and time unit')
        self.assertIsInstance(vegetationHistoryLength[0],int, msg='vegetation history lengths first value should be an integer')
        self.assertIsInstance(vegetationHistoryLength[1],str,msg='vegetation history lengths second value need to be month, day, or year,')
        self.assertIn(vegetationHistoryLength[1],['month','day','year'],msg='vegetation history lengths second value need to be month, day, or year,')

    def test_water_params(self):
        self.assertEqual(len(waterHistoryLength),2,msg='water history length requires 2 inputs for durration and time unit')
        self.assertIsInstance(waterHistoryLength[0],int, msg='water history lengths first value should be an integer')
        self.assertIsInstance(waterHistoryLength[1],str,msg='water history lengths second value need to be month, day, or year,')
        self.assertIn(waterHistoryLength[1],['month','day','year'],msg='water history lengths second value need to be month, day, or year,')

    def test_check_run_date(self):
        enddate = py_date+relativedelta(days=-12)
        self.assertIn(enddate.day, [30,31,1])

    def test_image_outputs(self):
        self.assertIsInstance(make_fire(),ee.Image,msg='failed building fire output')
        self.assertIsInstance(make_vegetation(),ee.Image,msg='failed building vegetation output')
        self.assertIsInstance(make_water(),ee.Image,msg='failed building water output')

    def test_export_asset(self):
        # assetbase
        try:
            ee.data.getAssetAcl(assetbase)
        except:
            self.fail("Gee asset base does not exist or doesn't allow this operation. Check path is correct and is writeable by account.")

    def test_export_cloud(self):
        # bucket
        # bucket = bucket
        self.assertIsInstance(bucket,str,msg='bucket should be a string of your buckets name. e.g. "BUCKET" for the bucket "gs://BUCKET/" ')
        self.assertNotIn("gs://",bucket,msg="bucket does not need the full path only the bucket name. remove gs://")

    def test_make_vegeation(self):
        def nd(image):
            out = image.normalizedDifference(['nir', 'red']).copyProperties(image, ["system:time_start"])
            return out

        ic = main.landsat().preprocess(startDate=ee.Date('2000-01-01'),endDate=ed,studyArea=region).select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2']).map(nd)
        self.assertIsNot(ic.size().getInfo(),0, msg="The input imagecollection has 0 images.")
        self.assertEqual(len(ic.first().bandNames().getInfo()),1,msg="More than one band is in the input image collection.")
        
        m = ed.get('month').getInfo()
        y = ed.get('year').getInfo()
        
        vegetation = main.Vegetation().byRegion(m, y, ic, region, vegetationHistoryLength)
        print(vegetation.bandNames().getInfo())

STARTDATE = '2020-01-01'
ENDDATE = '2020-12-01'
EE_ENDDATE = ee.Date(ENDDATE)

class Fire(unittest.TestCase):
    def setUp(self):
        self.fire = main.Fire()
        self.MODIS_IMG = ee.Image("MODIS/006/MYD14A2/2016_01_01").select([0])
        self.testGeo = ee.Geometry.Polygon(
        [[[-79.10665790098028, -1.6728296592564953],
          [-79.10665790098028, -2.5182035157244504],
          [-77.80367931830006, -2.5182035157244504],
          [-77.80367931830006, -1.6728296592564953]]], None, False)

    def test_fire_burnout(self):
        burn_out = self.fire.burnOut(EE_ENDDATE, -2, 'month').toInt16()
        
        self.assertEqual(burn_out.bandNames().getInfo(),['confidence', 'year', 'month', 'day', 'binary', 'unix', 'denisty', 'allFires', 'past', 'distance10k'])
    
    def test_get_fire(self):
        cM = EE_ENDDATE.get('month')
        cY = EE_ENDDATE.get('year')
        get_fire = self.fire.getFire(cY,cM)

        self.assertEqual(get_fire.bandNames().getInfo(),['confidence', 'year', 'month', 'day', 'binary', 'unix', 'denisty'])

    def test_fire_reclassify_hist(self):
        reclassified = self.fire.reclassify(self.MODIS_IMG).select('confidence')
        hist = reclassified.reduceRegion(**{'reducer':ee.Reducer.frequencyHistogram().unweighted(), 'geometry':self.testGeo, 'scale':1000,'maxPixels':1e13}).getInfo()

        self.assertEqual(hist,{'confidence': {'2': 2, '3': 10, '4': 6, 'null': 13593}})

    def test_fire_reclassify(self):
        reclassified = self.fire.reclassify(self.MODIS_IMG).bandNames().getInfo()
        self.assertEqual(reclassified,['confidence', 'year', 'month', 'day', 'binary', 'unix'])
    

# class Vegetation(unittest.TestCase):

if __name__ == '__main__':
    unittest.main()