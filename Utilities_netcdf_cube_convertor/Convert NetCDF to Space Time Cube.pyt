import arcpy as ARCPY

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Convert NetCDF to Space Time Cube"
        self.alias = ""
        self.helpContext = 0

        #### List of tool classes associated with this toolbox ####
        self.tools = [Convert2Cube]

class Convert2Cube(object):
    """
    Create Space-Time Cube. Aggregate point counts or summarize variables to a
    space-time cube.
    METHOD:
        __init__(): Define tool name and class info
        getParameterInfo(): Define parameter definitions in tool
        isLicensed(): Set whether tool is licensed to execute
        updateParameters():Modify the values and properties of parameters
                           before internal validation is performed
        updateMessages(): Modify the messages created by internal validation
                          for each tool parameter.
        execute(): Runtime script for the tool
    """
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Convert NetCDF to Space Time Cube"
        self.description = "Restructure NetCDF4 to Space Time Cube"
        self.canRunInBackground = False
        self.helpContext = 0

    def getParameterInfo(self):
        """Define parameter definitions"""
        #### Local Imports ####
        import os as OS
        import sys as SYS
        #### Define Parameters ####
        param0 = ARCPY.Parameter(displayName="Input NetCDF File",
                                 name="in_netcdf",
                                 datatype="DEFile",
                                 parameterType="Required",
                                 direction="Input")
        param0.filter.list = ['nc']

        param1 = ARCPY.Parameter(displayName="Output Space Time Cube",
                                 name="output_cube",
                                 datatype="DEFile",
                                 parameterType="Required",
                                 direction="Output")
        param1.filter.list = ['nc']

        params = [param0, param1]

        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        from os import path as PATH
        if parameters[0].value or parameters[0].altered:
            inputNC = parameters[0].value.value
            if PATH.exists(inputNC):
                filePathNoExt, ext = PATH.splitext(inputNC)
                filePathNoExt += "_Cube"
                parameters[1].value = PATH.join(filePathNoExt + ext)
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        from os import path as PATH
        if parameters[0].value or parameters[0].altered:
            inputNC = parameters[0].value.value
            if PATH.exists(inputNC) != True:
                parameters[0].setIDMessage("ERROR", 110003, inputNC)
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        import CubeConvertor as CC
        
        #### Input Parameters ####
        inputNC = parameters[0].valueAsText
        outputCube = parameters[1].valueAsText
        
        #### Retrieve Cube Structure Info from NetCDF File ####
        cube = CC.CubeInfoRetriever(inputNC)
        
        #### Build Cube Based on Info Retrieved from NetCDF File ####
        cb = CC.CubeBuilder(outputCube, cube.numRows, cube.numCols, 
                            cube.numTime, cube.extent, cube.extentInDegree, 
                            cube.cellSizeInfo, None, cube.cellSizeInfoLon, 
                            cube.cellSizeInfoLat, cube.timeSizeInfo, cube.startEndTimeInfo, 
                            cube.variables, cube.timeData, cube.spatialRef, mask = None)

        return

