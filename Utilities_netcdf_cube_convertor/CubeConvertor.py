"""
This script is a library to support converting an existing netcdf file to 
Esri Space Time Cube, which contains 3 classes: 
    CubeInfoRetriever: Retrieve Space Time Cube required info from netcdf file.
    CubeBuilder      : Construct a Space Time Cube based on the retrieved info.
    CubeSpatialRef   : Restructure the Spatial Reference for Space Time Cube.

Created by Leo Chin, 1/10/2015
"""
################################## Imports ####################################
import arcpy as ARCPY
import datetime as DT
import netCDF4 as NET
import numpy as NUM
import SSUtilities as UTILS
import time as TIME
from os import path as PATH

################################# Global Dictionary ###########################
TimeUnitConverter = {"SECONDS":1, "MINUTES":60, "HOURS": 3600, "DAYS": 86400,
                     "WEEKS": 604800}

#### Enable Overwrite Output ####
ARCPY.env.overwriteOutput = True

################################# Classes #####################################
class CubeInfoRetriever(object):
    """
    INPUT:
    
    OUTPUT:
    
    ATTRIBUTE:
    """
    def __init__(self, rawFile):
        #### Runtime Check for Input NetCDF File Existence ####
        if PATH.exists(rawFile):
            extension = PATH.splitext(rawFile)[1]
        else:
            msg = 'Input NetCDF File is not existed'
            ARCPY.AddError(msg)
            raise SystemExit
        
        #### Initialize Global Variable ####
        if extension != ".nc":
            raise SystemExit()
        else:
            self.variables = []
            self.extentInDegree = ARCPY.Extent(0,0,0,0)
            self.extent = ARCPY.Extent(0,0,0,0)
            self.numTime = 0
            self.numRows = 0
            self.numCols = 0
            self.cellSizeInfo = "1 Meter"
            self.timeSizeInfo = "1 Day"
            self.cellSizeInfoLon = "0.000001 Degree"
            self.cellSizeInfoLat = "0.000001 Degree"
            self.timeData = []
            self.startEndTimeStr = ""
            if ARCPY.env.outputCoordinateSystem:
                self.spatialRef = ARCPY.env.outputCoordinateSystem
                if self.spatialRef.projectionCode == 0:
                    msg = "Environment Output Coordination System is Not Projected. "
                    msg += "Set Spatial Reference to Default WGS 1984 World Mercator"
                    ARCPY.AddWarning(msg)
                    self.spatialRef = ARCPY.SpatialReference(3395)
            else:
                #### Default Projection: WGS 1984 World Mercator ####
                self.spatialRef = ARCPY.SpatialReference(3395)
                msg = "Default Spatial Reference is WGS 1984 World Mercator"
                ARCPY.AddWarning(msg)
    
        #### Initialize Local Variable ####
        rawDataset = NET.Dataset(rawFile)
        [minLon, minLat, maxLon, maxLat] = [0, 0, 0, 0]
        [flagTime, flagLon, flagLat] = [0, 0, 0]
        
        #### Checking Input NetCDF File Version ####
        if "NETCDF4" not in rawDataset.file_format:
            msg = "Your File is "+ rawDataset.file_format +". "
            msg += "Space Time Cube Is in NetCDF4 Format, An Conversion Will Be Made"
            ARCPY.AddWarning(msg)

        #### Retrieving Dimension Info ####
        [tDim, xDim, yDim] = [None, None, None]
        for dim in rawDataset.dimensions:
            dInfo = rawDataset.dimensions[dim]
            if 'TIME' in dim.upper():
                self.numTime = int(dInfo.__len__())
                tDim = dim
            elif 'LON' in dim.upper():
                self.numCols = int(dInfo.__len__())
                xDim = dim
            elif 'LAT' in dim.upper():
                self.numRows = int(dInfo.__len__())
                yDim = dim
            else:
                pass
        #### Runtime Check for Dimensions ####
        if tDim == None or xDim == None or yDim == None:
            msg = "Unable to retrieve valid dimensions in your file"
            ARCPY.AddError(msg)
            rawDataset.close()
            raise SystemExit()

        #### Retrieving Variable Info ####
        for var in rawDataset.variables:
            vInfo = rawDataset.variables[var]
            dims = str(vInfo.dimensions).upper()
            if hasattr(vInfo, "missing_value"):
                missingValue = vInfo.missing_value
            else:
                missingValue = -9999
            correctDim = ("TIME" in dims) and ("LON" in dims) and ("LAT" in dims)
            if (vInfo.ndim == 3) and correctDim:
                data = vInfo[:]
                data[data==missingValue] = -9999
                self.variables.append([var, vInfo.long_name, vInfo.dtype, data])
            #### Retrieve Time Units ####
            elif ("TIME" in dims) and (vInfo.ndim == 1):
                flagTime = 1
                self.numTime = len(vInfo[:])
                if len(vInfo[:]) >= 10 :
                    self.timeData = vInfo[:]
                    if hasattr(vInfo, "units"):
                        timeUnit = vInfo.units
                    else:
                        timeUnit = "seconds since 1800-1-1 00:00:00"
                else:
                    msg = 'Less than 10 time slices. Failed to convert to cube'
                    ARCPY.AddError(msg)
                    rawDataset.close()
                    raise SystemExit()
            
            #### Retrieve Extent in Degree and in Projection ####
            elif ("LON" in dims) and (vInfo.ndim == 1):
                flagLon = 1
                self.numCols = len(vInfo[:])
                temp = [vInfo[0], vInfo[-1]]
                cellSizeLonLat = abs(vInfo[0] - vInfo[1])
                halfCellSizeLonLat = cellSizeLonLat * .5
                minLon = min(temp) - halfCellSizeLonLat 
                maxLon = max(temp) + halfCellSizeLonLat
                self.cellSizeInfoLon = str(vInfo[1] - vInfo[0]) + " Degrees"            
            elif ("LAT" in dims) and (vInfo.ndim == 1):
                flagLat = 1
                self.numRows = len(vInfo[:])
                temp = [vInfo[0], vInfo[-1]]
                cellSizeLonLat = abs(vInfo[0] - vInfo[1])
                halfCellSizeLonLat = cellSizeLonLat * .5
                minLat = min(temp) - halfCellSizeLonLat 
                maxLat = max(temp) + halfCellSizeLonLat
                self.cellSizeInfoLat = str(vInfo[0] - vInfo[1]) + " Degrees"
            else:
                pass
            
        #### Ensure More than One Variable ####
        numVars = len(self.variables)
        if numVars == 0:
            msg = "No available variables can be converted."
            ARCPY.AddError(msg)
            rawDataset.close()
            raise SystemExit()
        else:
            msg = 'There are total {0} variable(s) in your file converted to cube'.format(str(numVars))
            ARCPY.AddWarning(msg)
           
        #### Project Extent to PCS ####        
        self.extentInDegree = ARCPY.Extent(minLon, minLat, maxLon, maxLat)
        self.extentInDegree.spatialReference = self.spatialRef.GCS
        self.extent = self.extentInDegree.projectAs(self.spatialRef)
        
        #### Retrieve Cell Size Info String ####
        cellSize = (self.extent.XMax - self.extent.XMin) / (self.numCols - 1)
        cellUnit = UTILS.distanceUnitInfo[self.spatialRef.linearUnitName.upper()][0]
        self.cellSizeInfo = str(cellSize) + " " + cellUnit
        cellSizeLon = (self.extentInDegree.XMax - self.extentInDegree.XMin) / (self.numCols - 1)
        cellSizeLat = (self.extentInDegree.YMax - self.extentInDegree.YMin) / (self.numRows - 1)
        self.cellSizeInfoLon = str(cellSizeLon) + " Degrees" 
        self.cellSizeInfoLat = str(cellSizeLat) + " Degrees"
        
        #### Retrieve Time Size Info ####
        try:
            [unit, since, date, time] = timeUnit.split()
        except:
            [unit, since, date] = timeUnit.split()
            time = '00:00:00'
        dtString = date + " " + time
        
        ## Generate Time Size Info String ##
        tmpSize = float(self.timeData[1] - self.timeData[0])
        self.timeSizeInfo = str(tmpSize) + " " + unit.title()

        try:        
            #### Create Date Time Object ####
            startTimeObj = DT.datetime.strptime(dtString, '%Y-%m-%d %H:%M:%S')
        
            self.timeData *= TimeUnitConverter[unit.upper()]
    
            ## Start Start Time ##
            tmpSize = float(self.timeData[0])
            tmpDTObject = startTimeObj + DT.timedelta(seconds = tmpSize)
            startStartTime = tmpDTObject.isoformat().replace("T", " ")
    
            ## Start End Time ##
            tmpSize = float(self.timeData[1])
            tmpDTObject = startTimeObj + DT.timedelta(seconds = tmpSize)
            startEndTime = tmpDTObject.isoformat().replace("T", " ")
            ## End Start Time ##
            tmpSize = float(self.timeData[-2])
            tmpDTObject = startTimeObj + DT.timedelta(seconds = tmpSize)
            endStartTime = tmpDTObject.isoformat().replace("T", " ")
            ## End End Time ##
            tmpSize = float(self.timeData[-1])
            tmpDTObject = startTimeObj + DT.timedelta(seconds = tmpSize)
            endEndTime = tmpDTObject.isoformat().replace("T", " ")
            ## Start End Time String ##
            self.startEndTimeInfo = startStartTime + ";" + startEndTime + ";" 
            self.startEndTimeInfo += endStartTime + ";" + endEndTime
            ## Adjust self.timeData[] to Start from 0 ##
            self.timeData -= self.timeData[0]
        except OverflowError:
            msg = "The date value is out of range. "
            msg += "Check the time variable in your netcdf."
            ARCPY.AddError(msg)
            rawDataset.close()
            raise SystemExit()
        
        #### Verifying Info Retrieved from File Is Sufficient for Cube ####
        if flagTime == 0 or ((flagLon * flagLat) < 1):
            rawDataset.close()
            msg = "Can Not Retrieve Enough Dimension Information, Conversion Failed"
            ARCPY.AddError(msg)
            raise SystemExit()

        rawDataset.close()

class CubeBuilder(object):
    """
    INPUT:
    
    OUTPUT:
    
    ATTRIBUTE:
    """
    
    def __init__(self, location, numRows, numCols, numTime, extent = None, 
                 extentInDegree = None, cellSizeInfo = "", userCellSizeInfo = "",
                 cellSizeInfoLon = "", cellSizeInfoLat = "", timeSizeInfo = "", 
                 startEndTimeInfo = "", variables = [], timeData = [], 
                 spatialRef = None, mask = None):
        #### Set Input Parameters to Global ####
        UTILS.assignClassAttr(self, locals())
        
        #### Initialize Output Cube File ####
        try:
            self.cube = NET.Dataset(location, 'w')
        except:
            ARCPY.AddIDMessage("ERROR", 110014)
            raise SystemExit()
        
        #### Restructure Spatial Reference by CubeSpatialRef Class ####
        self.spatialRef = CubeSpatialRef(spatialRef, extentInDegree)
        
        #### Cell Size and Cell Unit ####
        self.cellSize, self.cellUnit = cellSizeInfo.split()
        self.timeSize, self.timeUnit = timeSizeInfo.split()
        if self.userCellSizeInfo:
            self.userCellSize, self.userCellUnit = self.userCellSizeInfo.split()
        else:
            self.userCellSize = self.cellSize
            self.userCellUnit = self.cellUnit.title()
            
        #### Cell Size and Unit in Longitude and Latitude ####
        self.cellSizeLon = float(self.cellSizeInfoLon.split()[0])
        self.cellSizeLat = float(self.cellSizeInfoLat.split()[0])
        self.cellSize = float(self.cellSize)
        self.userCellSize = float(self.userCellSize)
        
        #### Start-End Time Info ####
        self.startStartTime, self.startEndTime, self.endStartTime, self.endEndTime = self.startEndTimeInfo.split(";")
        
        #### Set Global Attribute for Output Cube File ####
        self.setGlobalAttr()
        
        #### Create Dimension and Its Related Variables for Output Cube File ####
        self.createDimension()
        
        #### Process Time Variable in Output Cube File ####
        self.defineTemporalRef()
        
        #### Process Spatial Related Variables in Output Cube File ####
        self.defineSpatialRef()
        
        for var in self.variables:
            self.create3DVariable(self.cube, var[0], var[1], var[2], var[3])
        
        #### Create Mask Variable for First Variable ####
        self.createMaskVariable(self.cube, self.variables[0][0])
        
        #### Create Supplementary Variable for Cube ####
        self.createSupplementVariable(self.cube)
        
        self.cube.close()
               
    def createDimension(self):        
        self.cube.createDimension('time', self.numTime)
        self.cube.createDimension('lon', self.numCols)
        self.cube.createDimension('lat', self.numRows)
        self.cube.createDimension('x', self.numCols)
        self.cube.createDimension('y', self.numRows)
        
        self.time = self.cube.createVariable('time', 'f8', ('time'))
        self.projection = self.cube.createVariable('projection', 'i4')
        self.x = self.cube.createVariable('x', 'f8', ('x'))
        self.y = self.cube.createVariable('y', 'f8', ('y'))
        self.lat = self.cube.createVariable('lat', 'f8', ('lat'))
        self.lon = self.cube.createVariable('lon', 'f8', ('lon'))
        
    def setGlobalAttr(self, ssdo = None):          
        self.cube.description = 'Space-Time Pattern Mining Cube'
        self.cube.history = 'Created by ' + TIME.ctime(TIME.time())
        self.cube.source = 'Space-Time Pattern Mining Tool, Esri ArcGIS Pro'
        self.cube.cell_size = self.cellSize
        self.cube.cell_unit = self.cellUnit
        self.cube.user_cell_size = self.userCellSize
        self.cube.user_cell_unit = self.userCellUnit
        self.cube.time_size = self.timeSize
        self.cube.time_unit = self.timeUnit
        self.cube.time_step_label = self.timeSizeInfo
        self.cube.start_time = self.startStartTime
        self.cube.end_time = self.startEndTime
        self.cube.end_start_time = self.endStartTime
        self.cube.end_end_time = self.endEndTime
        self.cube.xy_length = UTILS.prettyUnits(self.userCellSize, self.userCellUnit)
        extentStr = str(self.extent.XMin) + " " + str(self.extent.YMin) + " "
        extentStr += str(self.extent.XMax) + " " + str(self.extent.YMax)
        self.cube.extent = extentStr
        self.cube.projection_authority_code = self.spatialRef.PCSCode
        self.cube.esri_pe_string = self.spatialRef.peString
        if ssdo:
            self.cube.raw_pe_string = ssdo.spatialRef.exportToString()
        else:
            self.cube.raw_pe_string = self.cube.esri_pe_string
    
    def defineTemporalRef(self):        
        timeSize = int(self.timeData[1])
       
        self.time.long_name = 'time'
        self.time.standard_name = 'time'
        self.time.units =  'seconds since '+ self.startStartTime
        self.time.calendar = 'gregorian'
        self.time._CoordinateAxisType = 'Time'
        self.time._ChunkSize = timeSize
        self.time.type = 'dimension'
        self.time[:] = self.timeData
        
    def defineSpatialRef(self):
        cubeSpatialRef = self.spatialRef        
        cellSize = self.cellSize
        cellSizeLon = self.cellSizeLon
        cellSizeLat = self.cellSizeLat

        
        mappingType = cubeSpatialRef.gridMapping
        self.projection.grid_mapping_name = mappingType
        if mappingType == "lambert_azimuthal_equal_area":
            self.projection.longitude_of_projection_origin = cubeSpatialRef.lonOfProjOrigin
            self.projection.latitude_of_projection_origin = cubeSpatialRef.latOfProjOrigin
            self.projection.false_easting = cubeSpatialRef.falseEasting
            self.projection.false_northing = cubeSpatialRef.falseNorthing

        elif mappingType == "lambert_conformal":
            self.projection.standard_parallel = cubeSpatialRef.standardParallel
            self.projection.longitude_of_central_meridian = cubeSpatialRef.lonOfCentralMeridian
            self.projection.latitude_of_projection_origin = cubeSpatialRef.latOfProjOrigin
            self.projection.false_easting = cubeSpatialRef.falseEasting
            self.projection.false_northing = cubeSpatialRef.falseNorthing

        elif mappingType == "lambert_cylindrical_equal_area":
            self.projection.longitude_of_central_meridian = cubeSpatialRef.lonOfCentralMeridian
            self.projection.scale_factor_at_projection_origin = cubeSpatialRef.scaleFactorAtProjOrigin
            self.projection.false_easting = cubeSpatialRef.falseEasting
            self.projection.false_northing = cubeSpatialRef.falseNorthing

        elif mappingType == "mercator":
            self.projection.longitude_of_projection_origin = cubeSpatialRef.latOfProjOrigin
            self.projection.scale_factor_at_projection_origin = cubeSpatialRef.scaleFactorAtProjOrigin
            self.projection.false_easting = cubeSpatialRef.falseEasting
            self.projection.false_northing = cubeSpatialRef.falseNorthing

        elif mappingType == "orthographic":
            self.projection.longitude_of_projection_origin = cubeSpatialRef.lonOfProjOrigin
            self.projection.latitude_of_projection_origin = cubeSpatialRef.latOfProjOrigin
            self.projection.false_easting = cubeSpatialRef.falseEasting
            self.projection.false_northing = cubeSpatialRef.falseNorthing

        elif mappingType == "polar_stereographic":
            self.projection.straight_vertical_longitude_from_pole = cubeSpatialRef.lonOfProjOrigin
            self.projection.latitude_of_projection_origin = cubeSpatialRef.latOfProjOrigin
            self.projection.scale_factor_at_projection_origin = cubeSpatialRef.scaleFactorAtProjOrigin
            self.projection.false_easting = cubeSpatialRef.falseEasting
            self.projection.false_northing = cubeSpatialRef.falseNorthing

        elif mappingType == "rotated_pole":
            self.projection.grid_north_pole_longitude = cubeSpatialRef.lonOfProjOrigin
            self.projection.grid_north_pole_latitude = cubeSpatialRef.latOfProjOrigin
            self.x.units = cubeSpatialRef.linearUnitName.split("_")[0]
            self.x.long_name = 'grid_longitude'
            self.x.standard_name = 'grid_longitude'
            self.y.units = self.x.units
            self.y.long_name = 'grid_latitude'
            self.y.standard_name = 'grid_latitude'

        elif mappingType == "stereographic":
            self.projection.longitude_of_projection_origin = cubeSpatialRef.lonOfProjOrigin
            self.projection.latitude_of_projection_origin = cubeSpatialRef.latOfProjOrigin
            self.projection.scale_factor_at_projection_origin = cubeSpatialRef.scaleFactorAtProjOrigin
            self.projection.false_easting = cubeSpatialRef.falseEasting
            self.projection.false_northing = cubeSpatialRef.falseNorthing

        elif mappingType == "transverse_mercator":
            self.projection.scale_factor_at_central_meridian = cubeSpatialRef.scaleFactorAtCentralMeridian
            self.projection.longitude_of_central_meridian = cubeSpatialRef.lonOfCentralMeridian
            self.projection.latitude_of_projection_origin = cubeSpatialRef.latOfProjOrigin
            self.projection.false_easting = cubeSpatialRef.falseEasting
            self.projection.false_northing = cubeSpatialRef.falseNorthing

        elif mappingType == "vertical_perspective":
            self.projection.latitude_of_projection_origin = cubeSpatialRef.latOfProjOrigin
            self.projection.longitude_of_projection_origin = cubeSpatialRef.lonOfProjOrigin
            self.projection.perspective_point_height = ""
            self.projection.false_easting = cubeSpatialRef.falseEasting
            self.projection.false_northing = cubeSpatialRef.falseNorthing

        self.projection.esri_pe_string = cubeSpatialRef.peString

        #### Specify Parameters for x ####
        self.x.units = cubeSpatialRef.linearUnitName.split("_")[0]
        self.x.long_name = 'x coordinate of projection'
        self.x.standard_name = 'projection_x_coordinate'
        self.x.type = 'dimension'
        xRange = NUM.arange(self.extent.XMin, self.extent.XMax + cellSize, cellSize)
        self.x[:] = xRange[0:self.numCols] + cellSize * .5

        #### Specify Parameters for y ####
        self.y.units = self.x.units
        self.y.long_name = 'y coordinate of projection'
        self.y.standard_name = 'projection_y_coordinate'
        self.y.type = 'dimension'
        yRange = NUM.arange(self.extent.YMax, self.extent.YMin - cellSize, -cellSize)
        self.y[:] = yRange[0:self.numRows] - cellSize * .5
        
        #### Define Lat/Lon Info ####
        self.lat.units = 'degrees_north'
        self.lat.long_name = 'latitude coordinate'
        self.lat.standard_name = 'latitude'
        self.lat.type = 'dimension'
        latRange = NUM.arange(self.extentInDegree.YMax, self.extentInDegree.YMin - cellSizeLat, -cellSizeLat)
        self.lat[:] = latRange[0:self.numRows] - cellSizeLat * .5

        self.lon.units = 'degrees_east'
        self.lon.long_name = 'longitude coordinate'
        self.lon.standard_name = 'longitude'
        self.lon.type = 'dimension'
        lonRange = NUM.arange(self.extentInDegree.XMin, self.extentInDegree.XMax + cellSizeLon, cellSizeLon)
        self.lon[:] = lonRange[0:self.numCols] - cellSizeLon * .5

    def create3DVariable(self, dataset, varName, longName = "",varType = NUM.float64,
                         varData = []):
        if NUM.array(varData).shape != (self.numTime, self.numRows, self.numCols): 
            msg = 'Dimension did not match. Failed to create variable: {0}'.format(varName)
            ARCPY.AddWarning(msg)
        else:
            if varName in dataset.variables:
                msg = 'Variable {0} is existed in the cube. Overwrite the existing values'.format(varName)
                ARCPY.AddWarning(msg)
                var = dataset.variables[varName]
            else:
                if longName == "":
                    longName = varName
                if 'f' in str(varType):
                    missingVal = -9999.
                else:
                    missingVal = -9999
    
                var = dataset.createVariable(varName, varType, ('time', 'y', 'x'))
                var.coordinates = 'time x y'
                var.long_name = longName
                var.standard_name = varName
                var.grid_mapping = 'projection'
                var.esri_pe_string = dataset.esri_pe_string
                var.missing_value = missingVal
                var.type = 'variable'
            var[:] = varData
     
    def create2DVariable(self, dataset, varName, longName = "",varType = NUM.float64,
                         varData = []):
        if NUM.array(varData).shape != (self.numRows, self.numCols): 
            msg = 'Dimension did not match. Failed to create variable: {0}'.format(varName)
            ARCPY.AddWarning(msg)
        else:
            if varName in dataset.variables:
                msg = 'Variable {0} is existed in the cube. Overwrite the existing values'.format(varName)
                ARCPY.AddWarning(msg)
                var = dataset.variables[varName]
            else:
                if longName == "":
                    longName = varName
                if 'f' in str(varType):
                    missingVal = -9999.
                else:
                    missingVal = -9999
    
                var = dataset.createVariable(varName, varType, ('y', 'x'))
                var.coordinates = 'x y'
                var.long_name = longName
                var.standard_name = varName
                var.grid_mapping = 'projection'
                var.esri_pe_string = dataset.esri_pe_string
                var.missing_value = missingVal
                var.type = 'variable'
            var[:] = varData       

    def createMaskVariable(self, dataset, varName):
        maskName = 'MASK_' + varName.upper()
        if varName not in dataset.variables:
            msg = 'Can not find the variable. Failed to create mask'
            ARCPY.AddWarning(msg)
        else:
            varData = dataset.variables[varName][:]
            maskArray = varData.sum(axis = 0)
            maskArray[maskArray > 0] = 1
            maskArray[maskArray <= 0] = 0
            if maskArray.shape != (self.numRows, self.numCols):
                msg = 'Variable dimension did not match. Failed to create mask'
                ARCPY.AddWarning(msg)
            else:
                longName = 'PORCESSION_BINARY_MASK_FOR_' + varName.upper()                
                self.create2DVariable(dataset, maskName, longName, 'i', maskArray)
     
    def createSupplementVariable(self, dataset):
        sliceSize = self.numRows * self.numCols
        
        #### Time Step ID ####
        timeIDData = NUM.repeat(NUM.arange(0, self.numTime), sliceSize)
        timeIDData.shape = (self.numTime, self.numRows, self.numCols)
        self.create3DVariable(dataset, 'time_step_ID', 'TIME_STEO_ID', 'i4', timeIDData)
                
        #### Location ID ####
        locationData = NUM.tile(NUM.arange(0, sliceSize), self.numTime)
        locationData.shape = (self.numTime, self.numRows, self.numCols)
        self.create3DVariable(dataset, 'location_ID', 'LOCATION_ID', 'i4', locationData)
        
class CubeSpatialRef(object):
    """
    INPUT:
    
    OUTPUT:
    
    ATTRIBUTE:
    """
    def __init__(self, spatialRef, ssdo = None, extent = None, 
                 extentInDegree = None):
        self.useLatLon = False
        if spatialRef.projectionName != '':
            #### Specify the Spatial Reference Parameters ####
            self.gridMapping = spatialRef.projectionName.lower()
            self.lonOfProjOrigin = spatialRef.GCS.longitude
            self.PCSCode = spatialRef.PCSCode
            self.latOfProjOrigin = self.latOfOrgin(self.PCSCode)
            self.falseEasting = spatialRef.falseEasting
            self.falseNorthing = spatialRef.falseNorthing
            self.standardParallel = spatialRef.GCS.standardParallel1
            self.lonOfCentralMeridian = spatialRef.GCS.centralMeridian
            self.scaleFactorAtProjOrigin = spatialRef.scaleFactor
            self.scaleFactorAtCentralMeridian = self.scaleFactorAtProjOrigin
            if extentInDegree == None:
                self.extentInDegree = self.addExtentInDegree(ssdo, extent)
            else:
                self.extentInDegree = extentInDegree
            self.meterPerUnit = spatialRef.metersPerUnit
            self.linearUnitName = spatialRef.linearUnitName
            
    def latOfOrgin(self, PCScode):
        """
        INPUT:
            PCScode(str): Projection code

        OUTPUT
            latOfOrigin(float): Latitude of the projection origin(1)

        NOTE:
            (1)Used when arcpy can not find origin property.  Value is retrieved
               from projection definition.
        """

        import re as RE
        spatialRef = ARCPY.SpatialReference(PCScode)
        self.peString = spatialRef.exportToString()
        self.peString = self.peString.split(";")[0]
        lat_of_origin_regex = "Latitude_Of_Origin',(\d+\.?\d*)]"
        match = RE.search(lat_of_origin_regex, self.peString)
        self.peString = self.peString.replace("'", "\\\"")
        if match:
            return float(match.groups()[0])
        else:
            #### When Can't Find Latitude of Origin, Return 0 ####
            return 0.

    def addExtentInDegree(self, ssdo, extent):
        """
        INPUT:
            ssdo (obj): SSDataObject
            extent (NUM arr): Extent Object

        OUTPUT:
            extentInDegree (obj): extent object(1)

        NOTE:
            (1) Used PCS extent to creat a temp layer and project it to GCS
                and then obtain GCS extent (in degrees)
        """
        gcsProjection = ssdo.spatialRef.GCS
        if extent is None:
            extent = ssdo.extent

        gcs = extent.projectAs(gcsProjection)
        extentDegree = ARCPY.Extent(gcs.XMin, gcs.YMin, gcs.XMax, gcs.YMax,
                                    gcs.ZMin, gcs.ZMax, gcs.MMin, gcs.MMax)

        return extentDegree
    