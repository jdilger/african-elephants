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

class Fire(unittest.TestCase):
    def setUp(self):
        self.fire = main.Fire()
    def test_fire(self):
        self.assertTrue(False)
if __name__ == '__main__':
    unittest.main()