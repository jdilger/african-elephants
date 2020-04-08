class Base(object):
    def __init__(self):
        self.date = 
        self.forestMask 

class Fire(base):
	def __inti__(self):
		super(Fire, self).__init__()
		
	def reclassify(self,img):
		remapped = img.remap([0,1,2,3,4,5,6,7,8,9],[0,0,0,1,1,1,1,2,3,4]).rename(['confidence'])
		d = ee.Date(img.get('system:time_start'))
		y = ee.Image(d.get('year')).int16().rename(['year'])
		m = ee.Image(d.get('month')).int16().rename(['month']);
		day = ee.Image(d.get('day')).int16().rename(['day']);
		binary = remapped.select('confidence').gte(2).rename('binary')
		out = remapped.addBands(y).addBands(m).addBands(day).addBands(binary);
		out = out.updateMask(remapped.gte(2));

		return out


	def getFire(self,targetYear,targetMonth):
		#Bring in MYD14/MOD14
		modisFireAqua = ee.ImageCollection('MODIS/006/MYD14A2').select([0]);
		modisFireTerra = ee.ImageCollection('MODIS/006/MOD14A2').select([0]);
		modisFire = modisFireTerra.merge(modisFireAqua);
		singleMonth = modisFire.filter(ee.Filter.calendarRange(targetYear,targetYear,'year'))/
		                        .filter(ee.Filter.calendarRange(targetMonth,targetMonth,'month'));

		#Recode it, and find the year, month, and day- then add it to the map
		singleMonth = singleMonth.map(self.reclassify);
		sum_denisty = singleMonth.select('binary').sum().rename('denisty')

		return singleMonth.mosaic().addBands(sum_denisty)

	def burnOut(self,sd,step,unit):
		"""takes in a startdate as ee.Date() -maybe change to string later step as an interger
		unit of time as string (e.g. 'day','month','year') """
		cM = sd.get('month')
		cY = sd.get('year')
		currentFires = getFire(cY,cM)

		pD = sd.advance(step, unit)
		pM = pD.get('month')
		pY = pD.get('year')
		pastFires = getFire(pY,pM)

		mask = pastFires.select('binary')
		newFires = currentFires.where(mask.eq(1),0).selfMask()
		allCurrentFires = currentFires.select('binary').rename('allFires')


		return newFires.addBands(allCurrentFires)
  
  

class Vegetation(base):
	def __inti__(self):
		super(Vegetation, self).__init__()

	